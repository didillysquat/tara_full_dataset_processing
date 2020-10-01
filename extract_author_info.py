#!/usr/bin/env python3

"""
Purpose of this script is to produce a plain text list of authors and associated afiiliations based on a
TARA-PACIFIC_authors-lists template.
It will take the path to the .xlsx formatted TARA-PACIFIC_authors-lists, and the name of the sheet that the extraction
should happen from.
"""

import pandas as pd
import os
import sys
import difflib
from collections import defaultdict
import requests
import json

# For getting the superscript and subscript numbers
# https://stackoverflow.com/questions/8651361/how-do-you-print-superscript-in-python


SUBSCRIPT_MAP = {
    "0": "₀", "1": "₁", "2": "₂", "3": "₃", "4": "₄", "5": "₅", "6": "₆",
    "7": "₇", "8": "₈", "9": "₉"}

class AuthorInfoExtraction:
    def __init__(self, fixes):
        """
        If fixes is true then we will make our best attempt at fixing the affiliation mess by implementing
        clustering of affiliations and checking for string similarity
        """
        self.fixes = fixes
        self.excel_path = sys.argv[1]
        self.target_sheet_name = sys.argv[2]
        self.output_dir = sys.argv[3]
        self.author_categories = ['First author(s)', 'Contributing authors list #1',
                                  'Contributing authors list #2', 'Consortium Coordinators', 'Scientific Directors',
                                  'Contributing authors list #3']
        # Make a df where index is lastname with first letter of initial appended, and has cols
        # ['last name', 'first name', 'first name initial(s)']
        self.author_info_df = self._make_author_info_df()
        # Generate the author order where each author is represented by the index version of their name used in the
        # self.name_df created above.
        self.author_order = self._make_author_order()

        # Now create the author to affiliation number dictionary
        # And the affiliation number to affiliation string dictionary
        (
            self.affiliation_list,
            self.affil_str_to_affil_num_dict,
            self.affil_num_to_affil_str_dict,
            self.author_to_affil_num_list_dict
        ) = self._make_affiliations_dicts()

        self.creator_array = self._make_creator_array()

        # personal access token is required for making uploads and publishing on Zenodo using their API
        with open(os.path.join(os.path.dirname(sys.argv[0]), 'personal_access_token.txt'), 'r') as f:
            self.access_token = f.readline().rstrip()

    def _make_creator_array(self):
        """
        Make the creator array that will be passed to the zenodo meta data item.
        This should be a list, that contains one dictionary per person with the keys
        name: Name of creator in the format Family name, Given names
        affiliation: Affiliation of creator (optional).
        orcid: ORCID identifier of creator (optional).
        """
        creator_array = []
        for author in self.author_order:
            author_dict = {}
            author_dict['name'] = f'{self.author_info_df.at[author, "last name"]}, {self.author_info_df.at[author, "first name"]}'
            if self.author_info_df.at[author, 'affiliation'] != 'not-provided':
                author_dict['affiliation'] = self.author_info_df.at[author, 'affiliation']
            if self.author_info_df.at[author, 'ORCID'] != 'not-provided':
                author_dict['orcid'] = self.author_info_df.at[author, 'ORCID']
            creator_array.append(author_dict)
        return creator_array

    def do_zenodo_submission(self):
        """
        Use the class objects to create the Zenodo submission using their API documented here:
        https://developers.zenodo.org/#quickstart-upload
        """
        filename = 'TARA_PACIFIC_METAB_ITS2_v1.zip'
        path = '/Users/humebc/Google_Drive/projects/tara/its2_for_zenodo/TARA_PACIFIC_METAB_ITS2_v1.zip'
        params = {'access_token': self.access_token}
        headers = {"Content-Type": "application/json"}

        # Create a blank deposition
        r = requests.post('https://zenodo.org/api/deposit/depositions', params = {'access_token': self.access_token}, json={})
        status_code = r.status_code
        r_json = r.json()
        bucket_url = r.json()["links"]["bucket"]
        deposition_id = r.json()['id']

        # Add a file to the deposition
        # We pass the file object (fp) directly to the request as the 'data' to be uploaded.
        # The target URL is a combination of the buckets link with the desired filename seperated by a slash.
        with open(path, "rb") as fp:
            r = requests.put(
                "%s/%s" % (bucket_url, filename),
                data=fp,
                # No headers included in the request, since it's a raw byte request
                params=params,
            )
        r_json = r.json()

        # Now add datda to the deposition
        data = {
            'metadata': {
                'title': 'TARA Pacific ITS2 Symbiodiniaceae data release version 1',
                'upload_type': 'dataset',
                'description': 'This data is the result of the primary analysis of the Symbiodiniaceae ITS2 rDNA gene '
                               'sequencing data collected from all islands as part of the Tara Pacific expedition.'
                               'The analysis was conducted using SymPortal.',
                'access_right': 'embargoed',
                'license': 'CC-BY-4.0',
                'embargo_date' : '2025-10-01',
                'communities' : [{'identifier': 'tarapacific'}],
                'version' : '1',
                'language': 'eng',
                'creators': self.creator_array
            }
        }

        r = requests.put(f'https://zenodo.org/api/deposit/depositions/{deposition_id}', params = {'access_token': self.access_token}, data = json.dumps(data),
        headers = headers)

        # now its ready to be published but we will not do this as part of the script
        # so that the user is forced to check over the submission and verfiy that everything is correct.

        foo = 'bar'



    def output_author_info(self):
        # First output the two variations of the author lists
        author_string_list_w_o_affiliation_numbers = []
        for author in self.author_order:
            author_string_list_w_o_affiliation_numbers.append(f'{self.author_info_df.at[author, "last name"]}, {self.author_info_df.at[author, "first name initial(s)"]}')
        author_string_w_o_affiliation_numbers = '; '.join(author_string_list_w_o_affiliation_numbers)

        author_string_list_w_affiliation_numbers = []
        for author in self.author_order:
            #first get the affiliation string
            affil_super_script_list = []
            for affil_num in self.author_to_affil_num_list_dict[author]:
                affil_super_script_list.append(self._superscript(affil_num))
            sup_affil_string = '˒'.join(affil_super_script_list)
            author_string_list_w_affiliation_numbers.append(
                f'{self.author_info_df.at[author, "last name"]}, {self.author_info_df.at[author, "first name initial(s)"]}{sup_affil_string}')
        author_string_w_affiliation_numbers = '; '.join(author_string_list_w_affiliation_numbers)

        # author_string_w_o_affiliation_numbers = '; '.join([f'{self.author_info_df.at[author, "last name"]}, {self.author_info_df.at[author, "first name initial(s)"]}' for author in self.author_order])
        # author_string_w_affiliation_numbers = '; '.join([f'{self.author_info_df.at[author, "last name"]}, {self.author_info_df.at[author, "first name initial(s)"]}{','.join([self._superscript(self.author_to_affil_num_dict[author]])}' for author in self.author_order])

        # Then output the affiliations
        affiliations_one_line = '; '.join([f'{affil_num + 1}-{self.affil_num_to_affil_str_dict[affil_num + 1]}' for affil_num in range(len(self.affiliation_list))])
        affiliations_new_lines = ';\n'.join([f'{affil_num + 1}-{self.affil_num_to_affil_str_dict[affil_num + 1]}' for affil_num in range(len(self.affiliation_list))])

        if self.fixes:
            with open(os.path.join(self.output_dir, f'author_string_w_o_affiliation_numbers_w_fixes'), 'w') as f:
                for line in author_string_w_o_affiliation_numbers:
                    f.write(line)

            with open(os.path.join(self.output_dir, f'author_string_w_affiliation_numbers_w_fixes'), 'w') as f:
                for line in author_string_w_affiliation_numbers:
                    f.write(line)

            with open(os.path.join(self.output_dir, f'affiliations_one_line_w_fixes'), 'w') as f:
                for line in affiliations_one_line:
                    f.write(line)

            with open(os.path.join(self.output_dir, f'affiliations_new_lines_w_fixes'), 'w') as f:
                for line in affiliations_new_lines:
                    f.write(line)
        else:
            with open(os.path.join(self.output_dir, f'author_string_w_o_affiliation_numbers_w_o_fixes'), 'w') as f:
                for line in author_string_w_o_affiliation_numbers:
                    f.write(line)

            with open(os.path.join(self.output_dir, f'author_string_w_affiliation_numbers_w_o_fixes'), 'w') as f:
                for line in author_string_w_affiliation_numbers:
                    f.write(line)

            with open(os.path.join(self.output_dir, f'affiliations_one_line_w_o_fixes'), 'w') as f:
                for line in affiliations_one_line:
                    f.write(line)

            with open(os.path.join(self.output_dir, f'affiliations_new_lines_w_o_fixes'), 'w') as f:
                for line in affiliations_new_lines:
                    f.write(line)

        foo = 'bar'

    def _cluster_affiliation(self, affil_string):
        if 'Konstanz' in affil_string:
            return 'Department of Biology, University of Konstanz, 78457 Konstanz, Germany'
        if 'Genoscope' in affil_string:
            return 'Génomique Métabolique, Genoscope, Institut François Jacob, CEA, CNRS, Univ Evry, ' \
                   'Université Paris-Saclay, 91057 Evry, France'
        if 'School of Marine Sciences, University of Maine' in affil_string:
            return 'School of Marine Sciences, University of Maine, Orono, 04469, Maine, USA'
        if 'Department of Biology, Institute of Microbiology and Swiss Institute of Bioinformatics' in affil_string:
            return 'Department of Biology, Institute of Microbiology and Swiss Institute of Bioinformatics, ' \
                   'Vladimir-Prelog-Weg 4, ETH Zürich, CH-8093 Zürich, Switzerland'
        if 'PSL Research University' in affil_string:
            return 'PSL Research University: EPHE-UPVD-CNRS, USR 3278 CRIOBE, Laboratoire d’Excellence CORAIL, ' \
                   'Université de Perpignan, 52 Avenue Paul Alduy, 66860 Perpignan Cedex, France'
        if 'Oregon State University, Department of Microbiology' in affil_string:
            return 'Oregon State University, Department of Microbiology, 220 Nash Hall, 97331 Corvallis OR USA'
        if 'Sorbonne Université, CNRS, Station Biologique de Roscoff, AD2M, UMR 7144, ECOMAP' in affil_string:
            return 'Sorbonne Université, CNRS, Station Biologique de Roscoff, AD2M, UMR 7144, ' \
                   'ECOMAP 29680 Roscoff, France'
        return affil_string

    def _make_affiliations_dicts(self):
        """
        Produce the affiliation list in order and three dicts
        We want to have each affiliation listed only once and numbered in order of the authors
        The affiliation numbers should start at 1
        """
        affiliation_list = []
        affil_str_to_affil_num_dict = {}
        affil_num_to_affil_str_dict = {}
        author_to_affil_num_list_dict = defaultdict(list)

        for author in self.author_order:

            if self.fixes:
                # Get the affiliation of the author
                # Four authors have two affiliations
                if author == 'de VargasC':
                    author_affil_list = [self._cluster_affiliation("Sorbonne Université, CNRS, Station Biologique de Roscoff, AD2M, UMR 7144, ECOMAP 29680 Roscoff, France"), self._cluster_affiliation("Research Federation for the study of Global Ocean Systems Ecology and Evolution, FR2022/ Tara Oceans-GOSEE, 3 rue Michel-Ange, 75016 Paris, France")]
                elif author in ['ForcioliD', 'FurlaP']:
                    author_affil_list = [self._cluster_affiliation("Université Côte d’Azur, CNRS, Inserm, IRCAN, France"), self._cluster_affiliation("Department of Medical Genetics, CHU of Nice, France")]
                elif author == 'CassarN':
                    author_affil_list = [self._cluster_affiliation("Division of Earth and Ocean Sciences, Duke University, Durham, USA"), self._cluster_affiliation("Laboratoire des Sciences de l’Environnement Marin (LEMAR), UMR 6539 UBO/CNRS/IRD/IFREMER, Institut Universitaire Européen de la Mer (IUEM), Brest, France")]
                else:
                    try:
                        author_affil_list = [self._cluster_affiliation(self.author_info_df.at[author, 'affiliation'])]
                    except KeyError as e:
                        raise RuntimeWarning(f'An affiliation could not be found for {author}\n'
                                         f'No affiliation will be associated.')

                for affil_string in author_affil_list:
                    if affil_string not in affiliation_list:
                        # Check that there is not a high similarity to an affiliation that is already in the affil list
                        for existing_affil in affiliation_list:
                            if difflib.SequenceMatcher(None, affil_string, existing_affil).ratio() > 0.6:
                                raise RuntimeError(f"{affil_string} and {existing_affil} seem to be very similar.")
                        # No similarity Error raised
                        affiliation_list.append(affil_string)
                        affil_str_to_affil_num_dict[affil_string] = len(affiliation_list)
                        author_to_affil_num_list_dict[author].append(len(affiliation_list))
                        affil_num_to_affil_str_dict[len(affiliation_list)] = affil_string
                    else:
                        author_to_affil_num_list_dict[author].append(affil_str_to_affil_num_dict[affil_string])
            else:
                # Get the affiliation of the author
                # Four authors have two affiliations

                try:
                    author_affil_list = [self.author_info_df.at[author, 'affiliation']]
                except KeyError as e:
                    raise RuntimeWarning(f'An affiliation could not be found for {author}\n'
                                         f'No affiliation will be associated.')

                for affil_string in author_affil_list:
                    if affil_string not in affiliation_list:
                        affiliation_list.append(affil_string)
                        affil_str_to_affil_num_dict[affil_string] = len(affiliation_list)
                        author_to_affil_num_list_dict[author].append(len(affiliation_list))
                        affil_num_to_affil_str_dict[len(affiliation_list)] = affil_string
                    else:
                        author_to_affil_num_list_dict[author].append(affil_str_to_affil_num_dict[affil_string])
        return affiliation_list, affil_str_to_affil_num_dict, affil_num_to_affil_str_dict, author_to_affil_num_list_dict

    def _make_author_order(self):
        df = pd.read_excel(io=self.excel_path, sheet_name=self.target_sheet_name, header=0)
        # Drop any rows that contain only nan, these may be the 'filtering' rows that excel inserts
        df.dropna(axis='index', how='all', inplace=True)
        # Drop any rows that have a score of 0
        df = df[df['sum'] > 0]
        df = df.fillna(value=0)
        # Keep only the useful cols
        df = df[
            ['last name', 'first name', 'First author(s)', 'Contributing authors list #1',
             'Contributing authors list #2', 'Consortium Coordinators', 'Scientific Directors',
             'Contributing authors list #3', 'sum']
        ]
        # Create the same index as for the name, orchid and affiliation df
        name_index = [last_name + first_name[0] for last_name, first_name in zip(df['last name'], df['first name'])]
        # Check that the last names are unique
        assert (len(name_index) == len(set(name_index)))
        df.index = name_index

        # The author order is then generated by going column by column
        # within each column going in order of top to bottom.
        # An author will only be added once obviously.
        # Where authors appear in two authorship categories, they will be placed into the first category they
        # appear in the order of
        # ['First author(s)', 'Contributing authors list #1',
        # 'Contributing authors list #2', 'Consortium Coordinators']
        # UNLESS the author appears in the Scientific Directors category or the Contributing authors list #3 in which
        # case they will be placed in this category.

        author_order_list = []
        for author_category in self.author_categories:

            if author_category not in  ['Scientific Directors', 'Contributing authors list #3']:
                # Then we are working with one of the first 4 author categories
                for author in df.index:
                    if author not in author_order_list and df.at[author, author_category] > 0 and df.at[author, 'Scientific Directors'] == 0 and df.at[author, 'Contributing authors list #3'] == 0:
                        # Then the author has not yet been placed into the author list
                        # The author is listed in the given author_category
                        # and the author is not listed in the scientific directors or contributing #3 categories
                        author_order_list.append(author)
            else:
                for author in df.index:
                    # Then we are in the scientific director or contributing #3 categories
                    if author not in author_order_list and df.at[author, author_category] > 0:
                        author_order_list.append(author)

        assert(len(author_order_list) == len(df.index))

        return author_order_list

    def _make_author_info_df(self):
        df = pd.read_excel(io=self.excel_path, sheet_name='Template', header=0)
        # Drop any rows that contain only nan, these may be the 'filtering' rows that excel inserts
        df.dropna(axis='index', how='all', inplace=True)
        df = df[['last name', 'first name', 'first name initial(s)', 'affiliation', 'ORCID']]
        # last names are not unique so create a key from the last name and fist initial
        name_index = [last_name + first_name[0] for last_name, first_name in zip(df['last name'], df['first name'])]
        # Check that the last names are unique
        assert(len(name_index) == len(set(name_index)))
        df.index = name_index
        if self.fixes:
            df.loc['ClayssenQ'] = ['Clayssen', 'Quentin', 'C.', 'Génomique Métabolique, Genoscope, Institut François Jacob, CEA, CNRS, Univ Evry, Université Paris-Saclay, 91057 Evry, France', 'not-provided']
            df.at['McMindsR', 'affiliation'] = 'Center of Modeling, Simulation and Interactions, Université Côte d′Azur, Nice, France'
        else:
            df.at['PogoreutzC', 'affiliation'] = 'Department of Biology, University of Konstanz, 78457 Konstanz, Germany'
            df.loc['ClayssenQ'] = ['Clayssen', 'Quentin', 'C.', 'not-provided', 'not-provided']
        return df

    def _superscript(self, number_to_convert):
        sup_map = {
            "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶",
            "7": "⁷", "8": "⁸", "9": "⁹"}
        num = str(number_to_convert)
        sup_num = ''
        for c in num:
            sup_num += sup_map[c]
        return sup_num


aie = AuthorInfoExtraction(fixes=False)
aie.output_author_info()
# At this point we have all of the objects we need to create the Zenodo submission object.
aie.do_zenodo_submission()
