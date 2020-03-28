"""
Script to cover all of the processing for the TARA ITS2 samples from the corals and from the surface waters
It will use a large part of the code that was used when we had only the data from the first three islands
available. This previous script was call tara_processing_script.py

To start with we will want to create an info df for the samples. This in turn can be used to generate a
datasheet that we will then be able to use to run the samples through SymPortal.

This is probaly a good stage to get to for a first attempt
"""
import os
import pickle
import sys
import pandas as pd
import subprocess
from collections import defaultdict
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import json
import numpy as np
import time
from bs4 import BeautifulSoup
from multiprocessing import Pool, current_process
import requests
import datetime
from pathlib import Path

class ITS2Processing:
    # This whole process became a little complicated when we realised that there were sometimes multiple sets
    # of fastq files per barcode_id/marker pairing. There were three categories of why this happened
    # 1 - mutiple lanes of sequencing on same sequencing run of same library - green
    # 2 - multiple sequencing runs of same library - yellow
    # 2 - different extraction of PCR methods (i.e. different libraries) for a single barcode_id/marker pair
    # We have decided for sake of simplicity and documentation that we will just use one pair of fastq files
    # per barcode_id/marker pairing. If we did decide to merge pairs of fastq files, we would need to create a new
    # file name and this could cause errors in documentaiton as we would have created a file that was not tracked
    # by genoscope of by the wider Tara Pacific meta data.

    # In this script, we join the code from two scripts that we have been using. One to parse the genoscope website
    # to see which fastq files are available for each sample, and one to create the datasheet for loading the data
    # into symportal.
    # At the end of this script, we want to have a datasheet for loading the coral samples, a datasheet for loading
    # the non-coral samples. We also want to have a dataframe for all of the seq files that we will be using and
    # similarly a dataframe for all of the seq files that we will not be using, i.e. that need removing from the
    # SymPortal web page. We need to have a careful think about what info will be useful in each of the dataframes.
    # I think the best way to do this is to have a single dataframe that has a will or will not use column.
    # Columns should be:
    # barcode_id, read_set, fastq_fwd_file_name, fastq_fwd_full_path, fastq_rev_file_name, fastq_rev_full_path,
    # will_use, genoscope_directory, has_replicates, access_date_time, replicate_class, replicate_color_code,
    # file_size_compressed

    # To make this dataframe we will need to use the info tables that Julie provided. We will also need to have a
    # list of the fastq files that Julie says it is OK to use for the red_cases.
    # For all read cases except 2, Julie has given us a specific set of fastq files to use. So for all red cases
    # except 2, we can check to see that we have been given the ones to keep. We can explicity check the name of the
    # two that we do not have info for. For these two we will keep the largest pair. For all yellow and green cases
    # we will just keep largest of the fastq files.

    #TODO we need to update so that we are using the new tara_sample_provenance table.
    # We also need to make sure that we are collecting the data from the provenance table
    # rather than from the directory structure as the directory structure may be wrong.
    # We have made a decision to keep only one of the sequencing files if there are multiple
    # pairs of fastq.gz files per barcode_id and marker gene. There were several different cases
    # in which this situation arises. In the majority of cases, we have sequencing data from two
    # different lanes but from the same library. In theory these fastq files could be safely joined
    # together but then we would be creating a new sequencing file that would have a name that is not
    # tracked by genoscope or stephanes meta data. So we will avoid that situation by just using the largest
    # pair of fastq files. The second situation was where two sequencing runs had been performed on the same
    # barcode id. The third situation was where a different methodology had been used (either extraction or PCR)
    # for the various sequencing pairs. In all cases we will just use one fastq pair. For the final cases I mentioned
    # here, there are some specific sequencing files that we should not be using in some cases. However, in other cases
    # Julie has said that we can use either of the sequencing sets. In these cases, we will use the largest pair again
    # if there are still mutiple fastq pairs after not considering the pairs that Julie said not to use.
    def __init__(self):
        self.base_directory_of_sequencing_data = "/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20200116_tara_pacific_its2/"
        self.input_dir = os.path.abspath(os.path.join('.', 'input'))
        self.output_dir = os.path.abspath(os.path.join('.', 'output'))
        os.makedirs(self.input_dir, exist_ok=True)
        os.makedirs(self.output_dir, exist_ok=True)
        self.sample_provenance_path = os.path.join(self.input_dir, 'sample_provenance_20200201.csv')
        self.sample_provenance_df = self._make_sample_provenance_df()
        self.readset_info_dir = "/home/humebc/projects/tara/replication_testing/readset_csvs"
        self.readset_df = self._make_readset_info_dir()
        self.cache_dir = os.path.abspath(os.path.join('.', 'cache'))
        # Two dictionaries that will hold the information for creating the dataframes form that will become the
        # symportal datasheets for doing the loading
        self.coral_sp_datasheet_df_dict = {}
        self.non_coral_sp_datasheet_df_dict = {}
        # The dict that will hold the info for the information df to document which files we used and which we did not
        self.output_information_df = {}
        # We can download the files that we are going to keep while we're at it
        # We should save them to a single directory
        self.seq_file_download_directory = "/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20200326_tara_its2_data"
        # Here we need to start our walk of the remote directories and
        self.remote_base_dir = "https://www.genoscope.cns.fr/sadc/tarapacific/METABARCODING/ITS2/ITS2_SYM_VAR_5.8S2_SYM_VAR_REV/"
        self.exe_path = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.authorisation_tup = self._make_auth_tup()
        self.headers = {'User-Agent': 'Benjamin Hume', 'From': 'benjamin.hume@kaust.edu.sa'}

    def _make_readset_info_dir(self):
        # read in the three sepearate csv files
        coral_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "coral_readset_info.csv"), skiprows=[0],
                                       names=['readset', 'primers', 'barcode_id', 'pcr_sample_name',
                                              'dna_sample_name'], dtype={'readset': str, 'primers': str, 'barcode_id': str, 'pcr_sample_name': str, 'dna_sample_name': str})
        sed_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "ssed_readset_info.csv"), skiprows=[0],
                                     names=['readset', 'primers', 'barcode_id', 'pcr_sample_name',
                                            'dna_sample_name'], dtype={'readset': str, 'primers': str, 'barcode_id': str, 'pcr_sample_name': str, 'dna_sample_name': str})
        fish_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "fish_readset_info.csv"), names=['readset', 'primers', 'barcode_id', 'pcr_sample_name', 'pcr_fl_sample_name',
                                   'dna_sample_name'], dtype={'readset': str, 'primers': str, 'barcode_id': str, 'pcr_sample_name': str, 'dna_sample_name': str, 'pcr_fl_sample_name': str})
        # fish_readset_df.drop(columns='PCR FL sample name', inplace=True)
        # fish_readset_df.columns = ['readset', 'primers', 'barcode_id', 'pcr_sample_name', 'pcr_fl_sample_name',
        #                            'dna_sample_name']
        plankton_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "plankton_readset_info.csv"), names=['readset', 'primers', 'barcode_id', 'pcr_sample_name', 'pcr_fl_sample_name',
                                   'dna_sample_name'], dtype={'readset': str, 'primers': str, 'barcode_id': str, 'pcr_sample_name': str, 'dna_sample_name': str, 'pcr_fl_sample_name': str})
        # # plankton_readset_df.drop(columns='PCR FL sample name', inplace=True)
        # plankton_readset_df.columns = ['readset', 'primers', 'barcode_id', 'pcr_sample_name', 'pcr_fl_sample_name',
        #                                'dna_sample_name']
        df = pd.concat([coral_readset_df, sed_readset_df, fish_readset_df, plankton_readset_df])
        return df.set_index('readset', drop=True)

    def _make_auth_tup(self):
        auth_path = os.path.join(self.exe_path, 'auth.txt')
        with open(auth_path, 'r') as f:
            auth_lines = [line.rstrip() for line in f]
        return (auth_lines[0], auth_lines[1])

    def start_walking(self):
        soup = BeautifulSoup(requests.get(self.remote_base_dir, auth=self.authorisation_tup, headers=self.headers).text,
                             features="html.parser", )
        worker_base_dirs = [link.string for link in soup.find_all('a') if ((link.string not in ['Name', 'Last modified',
                                                                                                'Size', 'Description',
                                                                                                'Parent Directory',
                                                                                                'NEGATIVE_CONTROLS/',
                                                                                                'NEGATIVE_CONTROL/']) and (
                                                                                       '/'.join([
                                                                                                    self.remote_base_dir.strip(
                                                                                                        '/'),
                                                                                                    link.string]) not in self.remote_base_dir))]
        worker_base_dirs = [os.path.join(self.remote_base_dir, _) for _ in worker_base_dirs]

        # NB We were originally mapping the rep_walker_list directly to the ReplicationWalkerWork class and running
        # its _walk function from within the __init__. However this was causing problems when running
        # for some of the markers and giving us an issue about errors returning the error and recursion (Pool error).
        # Instead we now pass in instances of the ReplicationWalkerWorker class and then run its _walk method, returning
        # only the list.
        # Create a ReplicationWalker for every worker_base_dir
        rep_walker_list = []
        for w_dir in worker_base_dirs:
            rep_walker_list.append(ReplicationWalkerWorker(w_dir, prov_df=self.sample_provenance_df, readset_df=self.readset_df))
        with Pool(20) as p:
            self.error_df_list_of_lists = p.map(self._run_walk_on_rep_walker_item, rep_walker_list)

        # At this point we will have the info required to make the sp_data sheets and output the
        # information dataframe
        # Maybe let's run it using just one thread to start with to do an initial debug.
        foo = 'bar'

    @staticmethod
    def _run_walk_on_rep_walker_item(rep_walker_class_instance):
        return rep_walker_class_instance._walk()

    def _make_sample_provenance_df(self):
        # The BARCODE ID is in the sequencing file name so we will be able to use this to get the latitute
        df = pd.read_csv(self.sample_provenance_path, skiprows=[1])
        df.set_index(keys='BARCODE ID', drop=True, inplace=True)
        return df

class ReplicationWalkerWorker:
    def __init__(self, remote_base_dir, prov_df=None, readset_df=None):
        self.remote_base_dir = remote_base_dir
        self.readset_info_dir = "/home/humebc/projects/tara/replication_testing/readset_csvs"
        if not readset_df is None:
            self.readset_df = readset_df
        else:
            self.readset_df = self._make_readset_info_dir()
        self.input_dir = os.path.abspath(os.path.join('.', 'input'))
        self.sample_provenance_path = os.path.join(self.input_dir, 'sample_provenance_20200201.csv')
        if not prov_df is None:
            self.sample_provenance_df = prov_df
        else:
            self.sample_provenance_df = self._make_sample_provenance_df()
        self.current_remote_dir = self.remote_base_dir
        self.done_list = set()
        self.done_and_empty_list = set()
        self.headers = {'User-Agent': 'Benjamin Hume', 'From': 'benjamin.hume@kaust.edu.sa'}
        self.exe_path = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.authorisation_tup = self._make_auth_tup()
        self.fastq_gz_list_current = []
        self.links_to_visit_current = []
        self.last_fork_location = None
        self.home_dir_reached = None
        self.s = requests.Session()
        self.s.auth = self.authorisation_tup
        self.s.headers = self.headers
        # Two dictionaries that will hold the information for creating the dataframes form that will become the
        # symportal datasheets for doing the loading
        self.coral_sp_datasheet_df_dict = {}
        self.non_coral_sp_datasheet_df_dict = {}
        # The list that will hold the info for the information df to document which files we used and which we did not
        self.output_information_list = []
        self.seq_file_download_directory = "/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20200326_tara_its2_data"
        # For the seq files that were replicated due to different methodologies being used, i.e. red cases
        # Julie gave us a list of the files that we should be using (one per barcode id). These files are listed here
        # (fwd files only). These are the strings from the readsets that we need to keep. They are not exact
        # UID of readsets, and there may be serveral sequencing replicates. We will have to search for the biggest to
        # keep.
        self.fwd_readset_strings_to_keep = [
            'HKNVMBCX2.12BA157',
            'HGY2FBCX2.12BA013',
            'HKNVMBCX2.12BA133',
            'BG8KK.12BA056',
            'BG8KK.12BA013',
            'BG8KK.12BA293',
            'BG8KK.12BA115'
        ]
        # She gave us the reads to keep for all but two of the barcodes. The two barcodes that she didn't give us
        # keep reads for are:
        self.no_keep_info_red_barcodes_list = ['CO-0002385', 'CO-0004425']

    def _walk(self):
        # This is the core unit of logic processing. Here we are visiting directories one by one and gathering
        # the information from them and downloading if necessary
        while True:
            print(f'{current_process()}: Current directory: {self.current_remote_dir}')
            soup = BeautifulSoup(self.s.get(self.current_remote_dir).text, features="html.parser")
            self.links_to_visit_current = [link.string for link in soup.find_all('a') if
                                           ((link.string not in ['Name', 'Last modified', 'Size', 'Description',
                                                                 'Parent Directory', 'NEGATIVE_CONTROLS/']) and
                                            ('/'.join([self.current_remote_dir.strip('/'),
                                                       link.string]) not in self.done_list))]
            self.fastq_gz_list_current = [link.string for link in self.links_to_visit_current if
                                          'fastq.gz' in link.string]
            self.links_to_visit_current = [link for link in self.links_to_visit_current if
                                           link not in self.fastq_gz_list_current]
            if len(self.links_to_visit_current) > 1:
                self.last_fork_location = self.current_remote_dir
            else:
                self.done_and_empty_list.add(self.current_remote_dir)
            if self.current_remote_dir == self.remote_base_dir and not self.links_to_visit_current:
                break

            if self.fastq_gz_list_current:
                try:
                    self._document_fastq_files()
                except IndexError as e:
                    print(e)
                    foo = 'bar'

            self.done_list.add(self.current_remote_dir)
            if self.links_to_visit_current:
                self.current_remote_dir = os.path.join(self.current_remote_dir, self.links_to_visit_current[0])
            else:
                if self.last_fork_location:
                    self.current_remote_dir = self.last_fork_location
                    self.last_fork_location = None
                else:
                    while True:
                        self.current_remote_dir = os.path.dirname(self.current_remote_dir)
                        if self.current_remote_dir == self.remote_base_dir or self.current_remote_dir + '/' == self.remote_base_dir:
                            if self.current_remote_dir not in self.done_and_empty_list and self.current_remote_dir + '/' not in self.done_and_empty_list:
                                self.walking_complete = False
                                break
                            else:
                                self.walking_complete = True
                                break
                        if self.current_remote_dir not in self.done_and_empty_list and self.current_remote_dir + '/' not in self.done_and_empty_list:
                            self.walking_complete = False
                            break
                    if self.walking_complete:
                        break

        return (self.coral_sp_datasheet_df_dict, self.non_coral_sp_datasheet_df_dict, self.output_information_list)

    def _document_fastq_files(self):
        # Then we can count how many there are, add current dir to done
        # and continue walking the directories
        if len(self.fastq_gz_list_current) > 2:
            self._handle_multiple_fastq_files()
        else:
            if not len(self.fastq_gz_list_current) == 2:
                raise RuntimeError('Odd number of fastq files')
            # Here we have only two files for a given barcode id and we can process accorrdingly
            self._handle_one_pair_fastq_files()

    def _handle_one_pair_fastq_files(self, read_tup=None, use=True, is_rep=False, cat=np.nan, col=np.nan, size_dict_passed=None):
        """ Columns for the information dict
        Either the fwd and rev can be supplied. Or if not supplied they will be found from the fastq_gz list
        if use is supplied as false, then we will not submit to the symportal datasheet
        """
        # barcode_id, read_set, fastq_fwd_file_name, fastq_rev_file_name,
        # will_use, genoscope_directory, is_replicate, access_date_time, replicate_class, replicate_color_code,
        # fwd_file_size_compressed, rev_file_size_compressed
        if read_tup:
            (fwd_read, rev_read) = read_tup
        else:
            fwd_read = [_ for _ in self.fastq_gz_list_current if 'R1' in _][0]
            rev_read = [_ for _ in self.fastq_gz_list_current if 'R1' in _][0]
        barcode_id = fwd_read.split('_')[1]
        # There may be more than one set of reads in the readset df for a given barcode id even
        # if there is only one set of fastq in the genoscope site
        if 'BID' in fwd_read:
            element_one = fwd_read.split('-')[-3]
            element_two = fwd_read.split('-')[-4].split('_')[-1]
            read_int = int(fwd_read.split('_')[-2][-1])
            bid_element = fwd_read.split('-')[-2]
            readset_str = f'{read_int}_{element_two}.{element_one}-{bid_element}'
        else:
            element_one = fwd_read.split('-')[-2]
            element_two = fwd_read.split('-')[-3].split('_')[-1]
            read_int = int(fwd_read.split('_')[-2][-1])
            readset_str = f'{read_int}_{element_two}.{element_one}'
        # now we need to look for this string in the indices
        # if we find more than one then raise error
        index_list = [_ for _ in self.readset_df[self.readset_df['barcode_id'] == barcode_id].index if readset_str in _]
        if len(index_list) != 1:
            raise RuntimeError
        else:
            readset = index_list[0]

        if use:
            size_dict = self._download_file_if_necessary(fwd_read, rev_read)
        else:
            if size_dict_passed:
                size_dict = size_dict_passed
            else:
                size_dict = self._get_sizes_trough_head_request(fwd_read, rev_read)

        # Now we can populate the output information dict and the sp_datasheet
        self._populate_output_information_list(barcode_id, fwd_read, readset, rev_read, size_dict, use, is_rep, cat, col)

        if use:
            self._populate_sp_datasheet_dict_item(barcode_id, fwd_read, rev_read)

    def _get_sizes_trough_head_request(self, fwd_read, rev_read):
        # then we need to get the file size without downloading the file
        # we can do this by sending a head request
        # https://stackoverflow.com/questions/14270698/get-file-size-using-python-requests-while-only-getting-the-header
        size_dict = {}
        for read in [fwd_read, rev_read]:
            response = self.s.head(os.path.join(self.current_remote_dir, read))
            size_dict[read] = int(response.headers['content-length'])
        return size_dict

    def _download_file_if_necessary(self, fwd_read, rev_read):
        # Check to see if the files already exist in the download directory
        # If they don't, then download them
        size_dict = {}
        for read in [fwd_read, rev_read]:
            local_file_path = os.path.join(self.seq_file_download_directory, read)
            # Download the file if it is not already downloaded
            # https://stackoverflow.com/questions/14270698/get-file-size-using-python-requests-while-only-getting-the-header
            if not os.path.isfile(local_file_path):
                self._download_fastq_gz(local_file_path, read)
                size_local = Path(local_file_path).stat().st_size
            else:
                #Check to see that the file on disk is as it should be
                response = self.s.head(os.path.join(self.current_remote_dir, read))
                size_remote = int(response.headers['content-length'])
                size_local = Path(local_file_path).stat().st_size
                if size_remote != size_local:
                    self._download_fastq_gz(local_file_path, read)
                    size_local = Path(local_file_path).stat().st_size
                    if size_remote != size_local:
                        raise RuntimeError('there seems to be size compatability issue')
            size_dict[read] = size_local

        return size_dict

    def _download_fastq_gz(self, local_file_path, read):
        print(f'downloading {read}')
        r = self.s.get(os.path.join(self.current_remote_dir, read), stream=True)
        with open(local_file_path, 'wb') as f:
            f.write(r.raw.data)

    def _populate_output_information_list(self, barcode_id, fwd_read, readset, rev_read, size_dict, use, is_rep, cat, col):
        # First the output information dict
        dat_string = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f%Z")
        temp_list = [
            barcode_id, readset, fwd_read, rev_read, use, self.current_remote_dir,
            is_rep, dat_string, cat, col, size_dict[fwd_read], size_dict[rev_read]
        ]
        self.output_information_list.append(temp_list)

    def _populate_sp_datasheet_dict_item(self, barcode_id, fwd_read, rev_read):
        """Now populate the sp datasheet
        Check to see if sample is coral or not
        headers are:
        'sample_name','fastq_fwd_file_name','fastq_rev_file_name','sample_type','host_phylum',
        'host_class','host_order','host_family','host_genus','host_species','collection_latitude','collection_longitude',
        'collection_date','collection_depth'"""
        barcode_prov_series = self.sample_provenance_df.loc[barcode_id]
        collection_date = barcode_prov_series['EVENT DATETIME START (UTC)']
        collection_depth = np.nan
        latitude = barcode_prov_series['EVENT LATITUDE START']
        longitude = barcode_prov_series['EVENT LONGITUDE START']
        if barcode_prov_series['SAMPLE ENVIRONMENT, short'] == 'C-CORAL':
            nominal_tax_name = barcode_prov_series["Sample Material label, organismal system level, taxonomic, nominal"]
            host_phylum = 'cnidaria'
            if nominal_tax_name == 'Porites':
                host_class = 'anthozoa'
                host_order = 'scleractinia'
                host_family = 'poritidae'
                host_genus = 'porites'
                host_species = 'lobata'
            elif nominal_tax_name == 'Porites panamensis':
                host_class = 'anthozoa'
                host_order = 'scleractinia'
                host_family = 'poritidae'
                host_genus = 'porites'
                host_species = 'panamensis'
            elif nominal_tax_name == 'Pocillopora':
                host_class = 'anthozoa'
                host_order = 'scleractinia'
                host_family = 'pocilloporidae'
                host_genus = 'pocillopora'
                host_species = 'meandrina'
            elif nominal_tax_name == 'Millepora':
                host_class = 'hydrozoa'
                host_order = 'anthoathecata'
                host_family = 'milleporidae'
                host_genus = 'millepora'
                host_species = 'dichotoma'
            elif nominal_tax_name == 'Pocillopora eydouxi':
                host_class = 'anthozoa'
                host_order = 'scleractinia'
                host_family = 'pocilloporidae'
                host_genus = 'pocillopora'
                host_species = 'eydouxi'
            elif nominal_tax_name == 'Millepora platyphylla':
                host_class = 'hydrozoa'
                host_order = 'anthoathecata'
                host_family = 'milleporidae'
                host_genus = 'millepora'
                host_species = 'platyphylla'
            elif nominal_tax_name == 'Heliopora':
                host_class = 'anthozoa'
                host_order = 'heliopracea'
                host_family = 'helioporidae'
                host_genus = 'heliopora'
                host_species = 'unknown'
            else:
                raise NotImplementedError

            self.coral_sp_datasheet_df_dict[barcode_id] = [fwd_read, rev_read, 'coral', host_phylum, host_class,
                                                           host_order, host_family, host_genus, host_species, latitude,
                                                           longitude, collection_date, collection_depth]
        else:
            # Non coral sample
            self.non_coral_sp_datasheet_df_dict[barcode_id] = [fwd_read, rev_read, 'seawater', np.nan, np.nan, np.nan,
                                                               np.nan, np.nan, np.nan, latitude, longitude,
                                                               collection_date, collection_depth]

    def _handle_multiple_fastq_files(self):
        # For the sediments there are multiple sample IDs worth of fastq
        # files in a single directory so we need to check for the ratio
        # of fastq files to sample IDs.
        barcode_id_set = set([_.split('_')[1] for _ in self.fastq_gz_list_current])
        if len(self.fastq_gz_list_current) > len(barcode_id_set) * 2:
            # This is where we need to categorise the replication as
            # one of three different categories,
            # When all seq files have the same base_name, then we will call
            # this 'sequencing_replication'.
            # Where the base_name is different, but the pcr and dna names are the
            # same, we will call this 'unknown_replication'
            # Finally where there are any differences in pcr or dna names
            # I will call these 'method_replication'
            base_names = {'-'.join(_.split('-')[:-1]) for _ in self.fastq_gz_list_current}
            if '' in base_names:
                print(self.fastq_gz_list_current)
                foo = 'bar'
            if len(base_names) > 1:
                self._process_unkn_method_replication(base_names, barcode_id_set)
            else:
                self._process_seq_replication()
        else:
            # here we have a case of multiple barcode ids worth of fastq pairs
            # The easiest way to treat this situation is to is to split it up into pairs of fastq
            # files and send these into the
            for barcode_id in barcode_id_set:
                barcode_reads = [_ for _ in self.fastq_gz_list_current if barcode_id in _]
                fwd_read = [_ for _ in barcode_reads if 'R1' in _][0]
                rev_read = [_ for _ in barcode_reads if 'R2' in _][0]
                self._handle_one_pair_fastq_files(read_tup=(fwd_read, rev_read))

    def _process_unkn_method_replication(self, base_names, barcode_id_set):
        # Then we need to check their PCR_sample_name and DNA_sample_name
        # To see if they are different.
        # This will produce two classes of seq difference
        # One where ReadSet names are different, but PCR and DNA are the same
        # Have to ask why?
        # The other set will be where ReadSet names are different, and one
        # of the PCR or DNA names are different.
        pcr_names_set = set()
        dna_names_set = set()
        for fastq_fwd in [_ for _ in self.fastq_gz_list_current if 'R1' in _]:
            barcode_id = fastq_fwd.split('_')[1]
            if 'BID' in fastq_fwd:
                element_one = fastq_fwd.split('-')[-3]
                element_two = fastq_fwd.split('-')[-4].split('_')[-1]
                read_int = int(fastq_fwd.split('_')[-2][-1])
                bid_element = fastq_fwd.split('-')[-2]
                readset_str = f'{read_int}_{element_two}.{element_one}-{bid_element}'
            else:
                element_one = fastq_fwd.split('-')[-2]
                element_two = fastq_fwd.split('-')[-3].split('_')[-1]
                read_int = int(fastq_fwd.split('_')[-2][-1])
                readset_str = f'{read_int}_{element_two}.{element_one}'
            readset_list = [_ for _ in self.readset_df[self.readset_df['barcode_id'] == barcode_id].index if
                          readset_str in _]
            if len(readset_list) == 1:
                readset = readset_list[0]
            else:
                raise RuntimeError
            pcr_sample_name = self.readset_df.at[readset, 'pcr_sample_name']
            dna_sample_name = self.readset_df.at[readset, 'dna_sample_name']


            pcr_names_set.add(pcr_sample_name)
            dna_names_set.add(dna_sample_name)

        # print(fastq_gz_list)
        if (len(pcr_names_set) != len(dna_names_set)) or (len(pcr_names_set) > len(barcode_id_set)):
            self._log_method_replication()
        else:
            self._log_unknown_replication()

    def _log_unknown_replication(self):
        # Then this a unknown_replication
        # Here we want to do the same as when we were doing _process_seq_replication
        # however, there may be more than two sets of sequences. We're looking to keep the biggest pair

        size_dict = {}
        for fastq_fwd in [_ for _ in self.fastq_gz_list_current if 'R1' in _]:
            size_dict[fastq_fwd] = sum(self._get_sizes_trough_head_request(
                fwd_read=fastq_fwd, rev_read=fastq_fwd.replace('R1', 'R2')
            ).values())
        # now simply sort to get the largest fwd read and submit that as keep, submit all others as no keep
        fwd_fastq_to_keep = sorted(size_dict, key=size_dict.get, reverse=True)[0]
        for fastq_fwd in [_ for _ in self.fastq_gz_list_current if 'R1' in _]:
            if fastq_fwd == fwd_fastq_to_keep:
                # This is the keep
                self._handle_one_pair_fastq_files(
                    read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=True, is_rep=True,
                    cat='unknown_replication', col='yellow'
                )
            else:
                # this is a no keep
                self._handle_one_pair_fastq_files(
                    read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=False, is_rep=True,
                    cat='unknown_replication', col='yellow'
                )

    def _log_method_replication(self):
        # then this is a 'method_replication'
        # We can use the elements in the temp_list to submit each of these
        # to the single fastq pair handler. But before we do that we need to identify readset (as this is the
        # UID essentially) that corresponds to the largest set of the files.
        # Actually TODO we need to check the do not use lists and check to see if we have identified the correct
        # TWO samples that Julie said we can check for.s

        # First check whether we are dealing with one of the barcodes that we don't have keep information for
        set_of_barcodes = set([_.split('_')[1] for _ in self.fastq_gz_list_current])
        if not len(set_of_barcodes) == 1:
            raise RuntimeError
        if list(set_of_barcodes)[0] in self.no_keep_info_red_barcodes_list:
            # Then we are working with one of the barcodes that we don't have keep info for.
            # We need to caculate the size and keep the biggest. Same as the sediment
            # We need to take into account that there may also be seq reps.
            size_dict = {}
            for fastq_fwd in [_ for _ in self.fastq_gz_list_current if 'R1' in _]:
                size_dict[fastq_fwd] = sum(self._get_sizes_trough_head_request(fastq_fwd, fastq_fwd.replace('R1', 'R2')).values())
            # Now we just want to keep the largest fastq_fwd in the size dict.
            fwd_fastq_to_keep = sorted(size_dict, key=size_dict.get, reverse=True)[0]
            for fastq_fwd in [_ for _ in self.fastq_gz_list_current if 'R1' in _]:
                if fastq_fwd == fwd_fastq_to_keep:
                    # This is the keep
                    self._handle_one_pair_fastq_files(
                        read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=True, is_rep=True,
                        cat='method_replication', col='red'
                    )
                else:
                    # this is a no keep
                    self._handle_one_pair_fastq_files(
                        read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=False, is_rep=True,
                        cat='method_replication', col='red'
                    )
        else:
            # Then one of the fwd fastq files should be in the self.fwd_readset_strings_to_keep list
            # For each of the read pairs we need to grab the readset strings and see if one of the keeps
            # is in that list
            readset_string_list = []
            for fastq_fwd in [_ for _ in self.fastq_gz_list_current if 'R1' in _]:
                if 'BID' in fastq_fwd:
                    element_one = fastq_fwd.split('-')[-3]
                    element_two = fastq_fwd.split('-')[-4].split('_')[-1]

                else:
                    element_one = fastq_fwd.split('-')[-2]
                    element_two = fastq_fwd.split('-')[-3].split('_')[-1]

                readset_str = f'{element_two}.{element_one}'
                readset_string_list.append(readset_str)
            # There should be exactly one match
            if not len(set(readset_string_list).intersection(set(self.fwd_readset_strings_to_keep))) == 1:
                # Then we are not finding a read that matches one of the keepers given by Juli
                raise RuntimeError
            # If we get here, then go back through and get the sizes of those seq pairs that have a readset that
            # matches the one given by Julies
            size_dict = {}
            for fastq_fwd in [_ for _ in self.fastq_gz_list_current if 'R1' in _]:
                if 'BID' in fastq_fwd:
                    element_one = fastq_fwd.split('-')[-3]
                    element_two = fastq_fwd.split('-')[-4].split('_')[-1]
                else:
                    element_one = fastq_fwd.split('-')[-2]
                    element_two = fastq_fwd.split('-')[-3].split('_')[-1]
                readset_str = f'{element_two}.{element_one}'
                if readset_str in self.fwd_readset_strings_to_keep:
                    size_dict[fastq_fwd] = sum(self._get_sizes_trough_head_request(fastq_fwd, fastq_fwd.replace('R1', 'R2')).values())
            #TODO we are here.
            # Now we just want to keep the largest fastq_fwd in the size dict.
            fwd_fastq_to_keep = sorted(size_dict, key=size_dict.get, reverse=True)[0]
            for fastq_fwd in [_ for _ in self.fastq_gz_list_current if 'R1' in _]:
                if fastq_fwd == fwd_fastq_to_keep:
                    # This is the keep
                    self._handle_one_pair_fastq_files(
                        read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=True, is_rep=True,
                        cat='method_replication', col='red'
                    )
                else:
                    # this is a no keep
                    self._handle_one_pair_fastq_files(
                        read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=False, is_rep=True,
                        cat='method_replication', col='red'
                    )

    def _process_seq_replication(self):
        # Then this is a case of sequence_replication
        # We should be able to find a readset per fastq.
        # The readset should contain two bits of information
        # in the fastq and it should containing the -1 or -2
        # This is a pain in the arse!
        size_dict = {}
        for fastq_fwd in [_ for _ in self.fastq_gz_list_current if 'R1' in _]:
            size_dict[fastq_fwd] = sum(self._get_sizes_trough_head_request(fastq_fwd, fastq_fwd.replace('R1', 'R2')).values())
        # Now we just want to keep the largest fastq_fwd in the size dict.
        fwd_fastq_to_keep = sorted(size_dict, key=size_dict.get, reverse=True)[0]
        for fastq_fwd in [_ for _ in self.fastq_gz_list_current if 'R1' in _]:
            if fastq_fwd == fwd_fastq_to_keep:
                # This is the keep
                self._handle_one_pair_fastq_files(
                    read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=True, is_rep=True,
                    cat='sequencing_replicate', col='green'
                )
            else:
                # this is a no keep
                self._handle_one_pair_fastq_files(
                    read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=False, is_rep=True,
                    cat='sequencing_replicate', col='green'
                )

    def _make_auth_tup(self):
        auth_path = os.path.join(self.exe_path, 'auth.txt')
        with open(auth_path, 'r') as f:
            auth_lines = [line.rstrip() for line in f]
        return (auth_lines[0], auth_lines[1])

    def _make_readset_info_dir(self):
        # read in the three sepearate csv files
        coral_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "coral_readset_info.csv"), skiprows=[0],
                                       names=['readset', 'primers', 'barcode_id', 'pcr_sample_name',
                                              'dna_sample_name'])
        sed_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "ssed_readset_info.csv"), skiprows=[0],
                                     names=['readset', 'primers', 'barcode_id', 'pcr_sample_name',
                                            'dna_sample_name'])
        fish_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "fish_readset_info.csv"))
        # fish_readset_df.drop(columns='PCR FL sample name', inplace=True)
        fish_readset_df.columns = ['readset', 'primers', 'barcode_id', 'pcr_sample_name', 'pcr_fl_sample_name',
                                   'dna_sample_name']
        plankton_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "plankton_readset_info.csv"))
        # plankton_readset_df.drop(columns='PCR FL sample name', inplace=True)
        plankton_readset_df.columns = ['readset', 'primers', 'barcode_id', 'pcr_sample_name', 'pcr_fl_sample_name',
                                       'dna_sample_name']
        df = pd.concat([coral_readset_df, sed_readset_df, fish_readset_df, plankton_readset_df])
        return df.set_index('readset', drop=True)



    def _make_sample_provenance_df(self):
        # The BARCODE ID is in the sequencing file name so we will be able to use this to get the latitute
        df = pd.read_csv(self.sample_provenance_path, skiprows=[1])
        df.set_index(keys='BARCODE ID', drop=True, inplace=True)
        return df


    def _make_sp_data_sheet(self):
        """Run this method to create a .csv file that can then be copy and pasted into a SymPortal
        datasheet template. We will then use this to start the symportal data. We are going to implement
        a new functionality to SymPortal that will allow for full paths to be provided in the datasheet in which
        case the the sequencing file pair will be looked for at that directory else the file name will be looked
        for in the directory that was provided as the argument to --load.
        As such, for the datasheet that we will want to create here we will want to have the full paths to the 
        seq files so that we don't need to create a single directory in which the files.
        
        THe columns for the SP datasheet are:
        sample_name	fastq_fwd_file_name	fastq_rev_file_name	sample_type	host_phylum	host_class	
        host_order	host_family	host_genus	host_species	collection_latitude	collection_longitude	collection_date	collection_depth
        """

        spsh = SPDataSheet(info_df=self.info_df, out_dir=self.output_dir)
        return spsh.create_sp_df()

class SequenceAndProfilePlotterCoralOnly:
    def __init__(self, info_df, sp_datasheet, seq_meta_info_df, seq_color_dict, prof_color_dict,
            sample_uid_to_name_dict,
            seq_absolute_abundance_df,
            seq_relative_abundance_df,
            prof_meta_info_df,
            prof_uid_to_name_dict,
            prof_absolute_abundance_df,
            prof_relative_abundance_df, output_dir, fig_type):
        self.info_df = info_df
        self.sp_datasheet = sp_datasheet
        self.seq_meta_info_df = seq_meta_info_df
        self.sample_uid_to_name_dict = sample_uid_to_name_dict
        self.seq_absolute_abundance_df = seq_absolute_abundance_df
        self.seq_relative_abundance_df = seq_relative_abundance_df
        self.prof_meta_info_df = prof_meta_info_df
        self.prof_uid_to_name_dict = prof_uid_to_name_dict
        self.prof_absolute_abundance_df = prof_absolute_abundance_df
        self.prof_relative_abundance_df = prof_relative_abundance_df
        self.output_dir = output_dir
        self.fig_type = fig_type
        self.species = ['PORITES', 'POCILLOPORA', 'MILLEPORA']
        # Get list of islands
        # Dict that will have island names as key and a list of sites as value
        self.island_site_dict = self._determine_sites_and_island()
        self.islands = sorted(list(self.island_site_dict.keys()))
        # Get the colour dicts that we will use for plotting
        self.seq_color_dict = seq_color_dict
        self.prof_color_dict = prof_color_dict
        if self.fig_type == 'genera':
            # Convert the seq_color_dict to genera colours
            temp_dict = {}
            for k, v in self.seq_color_dict.items():
                if k.startswith('A') or k.endswith('A'):
                    temp_dict[k] = "#DCDCDC"
                elif k.startswith('C') or k.endswith('C'):
                    temp_dict[k] = "#A9A9A9"
                elif k.startswith('D') or k.endswith('D'):
                    temp_dict[k] = "#696969"
                else:
                    temp_dict[k] = "#FF0000"
            self.seq_color_dict = temp_dict
        # Setup the plot
        self.fig = plt.figure(figsize=(14, 10))
        self.gs = self.fig.add_gridspec(18, 18, figure=self.fig, height_ratios=[1 for _ in range(18)],
                                        width_ratios=[1 for _ in range(18)])
        # The axis list that has been linearised
        # we will go in order of: for island, for site, for species
        # we will linearize the axes and go in order of for island, for site, for species
        for island in self.islands:
            site_list = sorted(list(self.island_site_dict[island]))[:3]
            for site in site_list:
                for species in self.species:
                    ax_row_index = int((int(self.islands.index(island)/6)*3) + site_list.index(site))
                    ax_col_index = int(((self.islands.index(island)%6)*3) + self.species.index(species))
                    # In here we can do the actual plotting
                    single_sp_time_start = time.time()
                    ax = self.fig.add_subplot(self.gs[ax_row_index, ax_col_index])
                    spip = SeqProfIndiPlot(parent=self, ax=ax, island=island, site=site, species=species, fig_type=self.fig_type)
                    spip.do_plotting()
                    single_sp_time_tot = time.time() - single_sp_time_start
                    print(f'The whole single took {single_sp_time_tot}')
        print('Writing .svg')
        plt.savefig(os.path.join(self.output_dir, f'seq_and_profile_fig_{self.fig_type}.svg'))
        print('Writing .png')
        plt.savefig(os.path.join(self.output_dir, f'seq_and_profile_fig_{self.fig_type}.png'), dpi=1200)
        foo = 'bar'


    def _determine_sites_and_island(self):
        # For each sample pull out the island and site from the info_df
        # Dict that will have island names as key and a list of sites as value
        island_site_dict = defaultdict(set)
        for sample_index in self.sp_datasheet.index:
            island = self.info_df.at[sample_index, 'location']
            site = self.info_df.at[sample_index, 'site']
            island_site_dict[island].add(site)
        return island_site_dict
        
class SeqProfIndiPlot:
    def __init__(self, parent, ax, island, site, species, fig_type):
        self.parent = parent
        self.ax = ax
        self.island = island
        self.site = site
        self.species = species
        self.fig_type = fig_type
        self.samples = self.parent.info_df[
            (self.parent.info_df['location'] == self.island) &
            (self.parent.info_df['site'] == self.site) &
            (self.parent.info_df['spp_water'] == self.species) &
            (self.parent.info_df['coral_plankton'] == 'CORAL') &
            (self.parent.info_df['spp_water'] != 'HELIOPORA') &
            (self.parent.info_df['spp_water'] != 'PORITES_PANAMENSIS')
        ].index.values.tolist()
        foo = 'bar'
        self.patches_list = []
        self.ind = 0
        self.color_list = []
        self.num_smp_in_this_subplot = len(self.samples)

    def do_plotting(self):
        div_over_total = 0
        type_under_total = 0
        for sample_to_plot in self.samples:
            sys.stdout.write(f'\rPlotting sample: {self.island} {self.site} {self.species} {sample_to_plot}')

            # PLOT DIVs
            div_over_start = time.time()
            self._plot_div_over_type(sample_to_plot)
            div_over_stop = time.time()
            div_over_total += (div_over_stop - div_over_start)

            # PLOT type
            if self.fig_type == 'all':
                type_under_start = time.time()
                self._plot_type_under_div(sample_to_plot)
                type_under_stop = time.time()
                type_under_total += (type_under_stop - type_under_start)
            self.ind += 1
        paint_start = time.time()
        self._paint_rect_to_axes_div_and_type()
        paint_stop = time.time()
        paint_total = paint_stop - paint_start

        print(f'\n\ndiv_over took {div_over_total}')
        print(f'type_under took {type_under_total}')
        print(f'paint took {paint_total}')



    def _plot_div_over_type(self, sample_to_plot):
        bottom_div = 0
        # In order that the sequences are listed in the seq_relative_abundance_df for those that are
        # present in the sample, plot a rectangle.
        smp_series = self.parent.seq_relative_abundance_df.loc[sample_to_plot]
        non_zero_series = smp_series.iloc[smp_series.to_numpy().nonzero()[0]]
        for seq_name, seq_rel_abund in non_zero_series.iteritems():
            # class matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
            self.patches_list.append(Rectangle((self.ind - 0.5, bottom_div), 1, seq_rel_abund, color=self.parent.seq_color_dict[seq_name]))
            self.color_list.append(self.parent.seq_color_dict[seq_name])
            bottom_div += seq_rel_abund

    def _plot_type_under_div(self, sample_to_plot):
        # the idea of the type is to put it as a reflection below the y=0 line
        # as such we should just want to make everything negative
        bottom_prof = 0
        # for each sequence, create a rect patch
        # the rect will be 1 in width and centered about the ind value.
        # we want to plot the rects so that they add to 1. As such we want to divide
        # each value by the total for that sample.
        non_zero_indices = self.parent.prof_absolute_abundance_df.loc[sample_to_plot].to_numpy().nonzero()[0]
        non_zero_series = self.parent.prof_absolute_abundance_df.loc[sample_to_plot].iloc[non_zero_indices]
        total = sum(non_zero_series.values)
        non_zero_series_relative = non_zero_series / total
        for profile_uid, profile_rel_abund in non_zero_series_relative.iteritems():
            # We will scale the profile so that it is 0.2 of the length of the seq info
            depth = -0.2 * profile_rel_abund
            self.patches_list.append(
                Rectangle((self.ind - 0.5, bottom_prof), 1, depth,
                            color=self.parent.prof_color_dict[profile_uid]))
            # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
            self.color_list.append(self.parent.prof_color_dict[profile_uid])
            bottom_prof += depth

    def _paint_rect_to_axes_div_and_type(self, max_num_smpls_in_subplot=10):
        # We can try making a custom colour map
        # https://matplotlib.org/api/_as_gen/matplotlib.colors.ListedColormap.html
        this_cmap = ListedColormap(self.color_list)
        
        # here we should have a list of Rectangle patches
        # now create the PatchCollection object from the patches_list
        patches_collection = PatchCollection(self.patches_list, cmap=this_cmap)
        patches_collection.set_array(np.arange(len(self.patches_list)))
        # if n_subplots is only 1 then we can refer directly to the axarr object
        # else we will need ot reference the correct set of axes with i
        # Add the pathces to the axes
        self.ax.add_collection(patches_collection)
        self.ax.autoscale_view()
        # self.ax.figure.canvas.draw()
        # also format the axes.
        # make it so that the x axes is constant length
        self.ax.set_xlim(0 - 0.5, max_num_smpls_in_subplot - 0.5)
        if self.fig_type == 'all':
            self.ax.set_ylim(-0.2, 1)
        else:
            self.ax.set_ylim(0,1)
        self._remove_axes_but_allow_labels()

        # as well as getting rid of the top and right axis splines
        # I'd also like to restrict the bottom spine to where there are samples plotted but also
        # maintain the width of the samples
        # I think the easiest way to do this is to hack a bit by setting the x axis spines to invisible
        # and then drawing on a line at y = 0 between the smallest and largest ind (+- 0.5)
        # ax.spines['bottom'].set_visible(False)
        if self.fig_type == 'all':
            self.ax.add_line(Line2D((0 - 0.5, self.num_smp_in_this_subplot - 0.5), (0, 0), linewidth=0.5, color='black'))

    def _remove_axes_but_allow_labels(self):
        self.ax.set_frame_on(False)
        self.ax.set_yticks([])
        self.ax.set_xticks([])

class QCPlotterCoralOnly:
    def __init__(self, info_df, sp_datasheet, seq_meta_info_df, output_dir):
        self.info_df = info_df
        self.sp_datasheet = sp_datasheet
        self.seq_meta_info_df = seq_meta_info_df
        self.output_dir = output_dir
        # Dict that will have island names as key and a list of sites as value
        self.island_site_dict = self._determine_sites_and_island()
        self.islands = sorted(list(self.island_site_dict.keys()))
        self.num_coral_islands = len(self.island_site_dict.keys())
        self.num_coral_sites = 0
        for k, v in self.island_site_dict.items():
            self.num_coral_sites += len(v)
        self.fig = plt.figure(figsize=(14, 10))
        self.gs = self.fig.add_gridspec(3, 32, figure=self.fig, height_ratios=[1, 1, 1], width_ratios=[1 for _ in range(32)])

        # setup each of the axes in separate lists for raw_contgs, non_symbiodiniaceae and symbiodiniaceae
        self.raw_contig_ax_list = []
        for i in range(32):
            self.raw_contig_ax_list.append(self.fig.add_subplot(self.gs[0,i]))
        self.non_sym_ax_list = []
        for i in range(32):
            self.non_sym_ax_list.append(self.fig.add_subplot(self.gs[1, i]))
        self.sym_ax_list = []
        for i in range(32):
            self.sym_ax_list.append(self.fig.add_subplot(self.gs[2, i]))

        # Here we have each of the axes setup. Now we just need to plot in them
        
        foo = 'bar'


    def _determine_sites_and_island(self):
        # For each sample pull out the island and site from the info_df
        # Dict that will have island names as key and a list of sites as value
        island_site_dict = defaultdict(set)
        for sample_index in self.sp_datasheet.index:
            island = self.info_df.at[sample_index, 'location']
            site = self.info_df.at[sample_index, 'site']
            island_site_dict[island].add(site)
        return island_site_dict

    def plot_qc_data(self):
        # We can go island by island.
        # For each island get the sites and sort
        # For each site of island get the samples
        # Then for each sample plot

        # We will need to keep track of some grand maximums so that we can scale
        # each of the subplot axes to the same values
        raw_contig_max = int(self.seq_meta_info_df['raw_contigs'].max())
        raw_contig_min = int(self.seq_meta_info_df['raw_contigs'].min())
        non_sym_total_max = int(self.seq_meta_info_df['post_taxa_id_absolute_non_symbiodinium_seqs'].max())
        non_sym_total_min = int(self.seq_meta_info_df['post_taxa_id_absolute_non_symbiodinium_seqs'].min())
        non_sym_distinct_max = int(self.seq_meta_info_df['post_taxa_id_unique_non_symbiodinium_seqs'].max())
        non_sym_distinct_min = int(self.seq_meta_info_df['post_taxa_id_unique_non_symbiodinium_seqs'].min())
        sym_total_max = int(self.seq_meta_info_df['post_taxa_id_absolute_symbiodinium_seqs'].max())
        sym_total_min = int(self.seq_meta_info_df['post_taxa_id_absolute_symbiodinium_seqs'].min())
        sym_distinct_max = int(self.seq_meta_info_df['post_taxa_id_unique_symbiodinium_seqs'].max())
        sym_distinct_min = int(self.seq_meta_info_df['post_taxa_id_unique_symbiodinium_seqs'].min())

        # We want to output a few stats to go with this figure. It would also be good to annotate these stats on
        # the figure. I think it would be good to do this using a horizonatl line
        print(f'Number of samples = {len(self.seq_meta_info_df.index)}')
        print(f'\tPorites lobata: {len(self.sp_datasheet[self.sp_datasheet["host_species"] == "lobata"].index)}')
        print(f'\tPocillopora meandrina: {len(self.sp_datasheet[self.sp_datasheet["host_species"] == "meandrina"].index)}')
        print(f'\tMillepora dichotoma: {len(self.sp_datasheet[self.sp_datasheet["host_species"] == "dichotoma"].index)}')
        average_raw_contigs_per_sample = int(self.seq_meta_info_df['raw_contigs'].mean())
        print(f'average_raw_contigs_per_sample: {average_raw_contigs_per_sample}')
        average_non_sym_total = int(self.seq_meta_info_df['post_taxa_id_absolute_non_symbiodinium_seqs'].mean())
        print(f'average_non_sym_total: {average_non_sym_total}')
        average_non_sym_distinct = int(self.seq_meta_info_df['post_taxa_id_unique_non_symbiodinium_seqs'].mean())
        print(f'average_non_sym_distinct: {average_non_sym_distinct}')
        average_sym_total = int(self.seq_meta_info_df['post_taxa_id_absolute_symbiodinium_seqs'].mean())
        print(f'average_sym_total: {average_sym_total}')
        average_sym_distinct = int(self.seq_meta_info_df['post_taxa_id_unique_symbiodinium_seqs'].mean())
        print(f'average_sym_distinct: {average_sym_distinct}')
        print(f'\n\ntotal_raw_contigs: {self.seq_meta_info_df["raw_contigs"].sum()}')

        for island_to_plot in self.islands:
            sites = sorted(list(self.island_site_dict[island_to_plot]))
            # For each of the subplots we will treat them as a single scatter plot
            # We will work with the x axis being limited to 0-->1
            # x axis coordinates will depend on the number of sites
            # we will want four plotting points per site,
            # one for the individual points and one for the avergae for total seqs
            # and then we'll want one for distnct sequences
            num_sites = len(sites)
            x_space = 1 / ((num_sites * 4) + 1)
            x_coord_vals = [i * x_space for i in range(1, num_sites*4 + 1, 1)]
            # Pair up the x_coord_vals so that they are easier to index
            # one pair per site. 0 will be for individual datapoints, 1 for the mean
            x_coord_vals = [(x_coord_vals[i], x_coord_vals[i + 1], x_coord_vals[i + 2], x_coord_vals[i + 3]) for i in range(0, num_sites*4, 4)]

            raw_contigs_ax = self.raw_contig_ax_list[self.islands.index(island_to_plot)]
            raw_contigs_ax.set_xlim(0, 1)
            raw_contigs_ax.set_ylim(10000, raw_contig_max)
            # set all x axes labels off
            raw_contigs_ax.set_xticks([])

            non_sym_ax = self.non_sym_ax_list[self.islands.index(island_to_plot)]
            non_sym_ax.set_xlim(0, 1)
            non_sym_ax.set_ylim(non_sym_total_min, non_sym_total_max)
            non_sym_ax.set_xticks([])
            non_sym_ax_distinct = non_sym_ax.twinx()
            non_sym_ax_distinct.set_ylim(non_sym_distinct_min, non_sym_distinct_max)

            sym_ax = self.sym_ax_list[self.islands.index(island_to_plot)]
            sym_ax.set_xlim(0, 1)
            sym_ax.set_ylim(sym_total_min, sym_total_max)
            sym_ax.set_xticks([])
            sym_ax_distinct = sym_ax.twinx()
            sym_ax_distinct.set_ylim(100, sym_distinct_max)
            # place the island label on the x axis rotated
            sym_ax.set_xlabel(island_to_plot, rotation='vertical')

            if self.islands.index(island_to_plot) == 0:
                # If first plot then only remove top and right
                # and we only want to remove the right y axis
                raw_contigs_ax.spines['right'].set_visible(False)
                raw_contigs_ax.spines['top'].set_visible(False)
                raw_contigs_ax.set_yscale('symlog')
                raw_contigs_ax.spines['left'].set_color(c='blue')
                raw_contigs_ax.tick_params('y', colors='blue')

                # plot average line
                raw_contigs_ax.axhline(y=average_raw_contigs_per_sample, xmin=0, xmax=1, c='blue')

                non_sym_ax.spines['right'].set_visible(False)
                non_sym_ax.spines['top'].set_visible(False)
                non_sym_ax.set_yscale('symlog')
                non_sym_ax.spines['left'].set_color(c='blue')
                non_sym_ax.tick_params('y', colors='blue')

                non_sym_ax_distinct.spines['right'].set_visible(False)
                non_sym_ax_distinct.spines['top'].set_visible(False)
                non_sym_ax_distinct.set_yscale('symlog')
                non_sym_ax_distinct.set_yticks([])
                non_sym_ax_distinct.minorticks_off()
                non_sym_ax_distinct.spines['left'].set_color(c='blue')

                sym_ax.spines['right'].set_visible(False)
                sym_ax.spines['top'].set_visible(False)
                sym_ax.set_yscale('symlog')
                sym_ax.spines['left'].set_color(c='blue')
                sym_ax.tick_params('y', colors='blue')

                sym_ax_distinct.spines['right'].set_visible(False)
                sym_ax_distinct.spines['top'].set_visible(False)
                sym_ax_distinct.set_yscale('symlog')
                sym_ax_distinct.set_yticks([])
                sym_ax_distinct.minorticks_off()
                sym_ax_distinct.spines['left'].set_color(c='blue')

            elif self.islands.index(island_to_plot) == 31:
                # If the last plot then we want to annotate the right y axis
                raw_contigs_ax.spines['left'].set_visible(False)
                raw_contigs_ax.spines['top'].set_visible(False)
                raw_contigs_ax.spines['right'].set_visible(False)
                raw_contigs_ax.set_yscale('symlog')
                raw_contigs_ax.set_yticks([])
                raw_contigs_ax.minorticks_off()


                non_sym_ax.spines['left'].set_visible(False)
                non_sym_ax.spines['top'].set_visible(False)
                non_sym_ax.set_yscale('symlog')
                non_sym_ax.set_yticks([])
                non_sym_ax.minorticks_off()
                non_sym_ax.spines['right'].set_color(c='red')


                non_sym_ax_distinct.spines['left'].set_visible(False)
                non_sym_ax_distinct.spines['top'].set_visible(False)
                non_sym_ax_distinct.set_yscale('symlog')
                non_sym_ax_distinct.spines['right'].set_color(c='red')
                non_sym_ax_distinct.tick_params('y', colors='red')


                sym_ax.spines['left'].set_visible(False)
                sym_ax.spines['top'].set_visible(False)
                sym_ax.set_yscale('symlog')
                sym_ax.set_yticks([])
                sym_ax.minorticks_off()
                sym_ax.spines['right'].set_color(c='red')

                sym_ax_distinct.spines['left'].set_visible(False)
                sym_ax_distinct.spines['top'].set_visible(False)
                sym_ax_distinct.set_yscale('symlog')
                sym_ax_distinct.spines['right'].set_color(c='red')
                sym_ax_distinct.tick_params('y', colors='red')
            else:
                raw_contigs_ax.spines['left'].set_visible(False)
                raw_contigs_ax.spines['top'].set_visible(False)
                raw_contigs_ax.spines['right'].set_visible(False)
                raw_contigs_ax.set_yscale('symlog')
                raw_contigs_ax.set_yticks([])
                raw_contigs_ax.minorticks_off()

                non_sym_ax.spines['left'].set_visible(False)
                non_sym_ax.spines['top'].set_visible(False)
                non_sym_ax.spines['right'].set_visible(False)
                non_sym_ax.set_yscale('symlog')
                non_sym_ax.set_yticks([])
                non_sym_ax.minorticks_off()

                non_sym_ax_distinct.spines['left'].set_visible(False)
                non_sym_ax_distinct.spines['top'].set_visible(False)
                non_sym_ax_distinct.spines['right'].set_visible(False)
                non_sym_ax_distinct.set_yscale('symlog')
                non_sym_ax_distinct.set_yticks([])
                non_sym_ax_distinct.minorticks_off()

                sym_ax.spines['left'].set_visible(False)
                sym_ax.spines['top'].set_visible(False)
                sym_ax.spines['right'].set_visible(False)
                sym_ax.set_yscale('symlog')
                sym_ax.set_yticks([])
                sym_ax.minorticks_off()

                sym_ax_distinct.spines['left'].set_visible(False)
                sym_ax_distinct.spines['top'].set_visible(False)
                sym_ax_distinct.spines['right'].set_visible(False)
                sym_ax_distinct.set_yscale('symlog')
                sym_ax_distinct.set_yticks([])
                sym_ax_distinct.minorticks_off()


            for site in sites:
                samples = self.info_df[(self.info_df['location'] == island_to_plot) & (self.info_df['site'] == site) & (self.info_df['coral_plankton'] == 'CORAL') & (self.info_df['spp_water'] != 'HELIOPORA') & (self.info_df['spp_water'] != 'PORITES_PANAMENSIS')]
                # Now we can populate the sample information for each of the plots
                # lists that we will create the avergae and stdev points from
                raw_contigs_vals = []
                non_sym_vals_total = []
                non_sym_vals_distinct = []
                sym_vals_total = []
                sym_vals_distinct = []



                # the x_coordinates
                indi_data_point_x_total = x_coord_vals[sites.index(site)][0]
                mean_data_point_x_total = x_coord_vals[sites.index(site)][1]
                indi_data_point_x_distinct = x_coord_vals[sites.index(site)][2]
                mean_data_point_x_distinct = x_coord_vals[sites.index(site)][3]

                # Now plot up the individual points and collect for calculating the mean
                for sample_name, sample_ser in samples.iterrows():
                    # raw_contigs
                    y_raw_contig = self.seq_meta_info_df.at[sample_name, 'raw_contigs']
                    raw_contigs_vals.append(y_raw_contig)
                    raw_contigs_ax.scatter(x=indi_data_point_x_total, y=y_raw_contig, marker='.', s=1, c='b')

                    # non_sym_vals total and distinct
                    y_non_sym_total = self.seq_meta_info_df.at[sample_name, 'post_taxa_id_absolute_non_symbiodinium_seqs']
                    non_sym_vals_total.append(y_non_sym_total)
                    non_sym_ax.scatter(x=indi_data_point_x_total, y=y_non_sym_total, marker='.', s=1, c='b')
                    y_non_sym_distinct = self.seq_meta_info_df.at[sample_name, 'post_taxa_id_unique_non_symbiodinium_seqs']
                    non_sym_vals_distinct.append(y_non_sym_distinct)
                    non_sym_ax_distinct.scatter(x=indi_data_point_x_distinct, y=y_non_sym_distinct, marker='.', s=1, c='r')

                    # sym counts total and distinct
                    y_sym_total = self.seq_meta_info_df.at[
                        sample_name, 'post_taxa_id_absolute_symbiodinium_seqs']
                    sym_vals_total.append(y_sym_total)
                    sym_ax.scatter(x=indi_data_point_x_total, y=y_sym_total, marker='.', s=1, c='b')
                    y_sym_distinct = self.seq_meta_info_df.at[
                        sample_name, 'post_taxa_id_unique_symbiodinium_seqs']
                    sym_vals_distinct.append(y_sym_distinct)
                    sym_ax_distinct.scatter(x=indi_data_point_x_distinct, y=y_sym_distinct, marker='.', s=1, c='r')

                # Here we should have all of the individual data points plotted up
                # We should also have collected all of the opints in a list so that we can now
                # calculate the mean and the standard deviations
                # For the time being just plot the mean and worry about stdev bars later
                mean_point_size = 40
                raw_contigs_ax.scatter(x=mean_data_point_x_total, y=sum(raw_contigs_vals)/len(raw_contigs_vals), marker='.', s=mean_point_size, c='b')
                non_sym_ax.scatter(x=mean_data_point_x_total, y=sum(non_sym_vals_total)/len(non_sym_vals_total), marker='.', s=mean_point_size, c='b')
                non_sym_ax_distinct.scatter(x=mean_data_point_x_distinct, y=sum(non_sym_vals_distinct)/len(non_sym_vals_distinct), marker='.', s=mean_point_size, c='r')
                sym_ax.scatter(x=mean_data_point_x_total, y=sum(sym_vals_total) / len(sym_vals_total), marker='.', s=mean_point_size, c='b')
                sym_ax_distinct.scatter(x=mean_data_point_x_distinct,
                                   y=sum(sym_vals_distinct) / len(sym_vals_distinct), marker='.', s=mean_point_size, c='r')

            # plot up the total average lines that will be extended by hand
            # plot on the 0 and 1 so that we can still see both lines despite overlap.
            if self.islands.index(island_to_plot) == 0:
                raw_contigs_ax.axhline(y=average_raw_contigs_per_sample, xmin=0, xmax=1, c='blue')
                non_sym_ax.axhline(y=average_non_sym_total, xmin=0, xmax=1, c='blue')
                sym_ax.axhline(y=average_sym_total, xmin=0, xmax=1, c='blue')
            if self.islands.index(island_to_plot) == 1:
                non_sym_ax_distinct.axhline(y=average_non_sym_distinct, xmin=0, xmax=1, c='red')
                sym_ax_distinct.axhline(y=average_sym_distinct, xmin=0, xmax=1, c='red')

        plt.savefig(os.path.join(self.output_dir, 'coral_only_qc_fig.svg'))
        plt.savefig(os.path.join(self.output_dir, 'coral_only_qc_fig.png'), dpi=1200)

class SPDataSheet:
    """Class to make the SPDataSheet from the sample info dataframe."""
    def __init__(self, info_df, out_dir):
        self.datasheet_cols = ['sample_name','fastq_fwd_file_name','fastq_rev_file_name','sample_type','host_phylum',
        'host_class','host_order','host_family','host_genus','host_species','collection_latitude','collection_longitude',
        'collection_date','collection_depth']
        self.rows = []
        self.info_df = info_df
        self.out_dir = out_dir
    
    def create_sp_df(self):
        # For each sample that is a coral, populate a row that will be entered
        # into the new datasheet dataframe.
        # NB I think we should be able to use None values where we want blank values,
        # i.e. for the collection dates
        non_three_list = defaultdict(int)
        for sample_name, ser in self.info_df[self.info_df['coral_plankton'] == 'CORAL'].iterrows():
            fastq_fwd_file_path = ser['fastq_fwd_file_path']
            fastq_rev_file_path = ser['fastq_rev_file_path']
            sample_type = 'coral'
            species = ser['spp_water']
            host_phylum = 'cnidaria'
            if species == 'PORITES':
                host_class = 'anthozoa'
                host_order = 'scleractinia'
                host_family = 'poritidae'
                host_genus = 'porites'
                host_species = 'lobata'
            elif species == 'POCILLOPORA':
                host_class = 'anthozoa'
                host_order = 'scleractinia'
                host_family = 'pocilloporidae'
                host_genus = 'pocillopora'
                host_species = 'meandrina'
            elif species == 'MILLEPORA':
                host_class = 'hydrozoa'
                host_order = 'anthoathecata'
                host_family = 'milleporidae'
                host_genus = 'millepora'
                host_species = 'dichotoma'
            else:
                non_three_list[species] += 1
                continue
            latitude = ser['lat']
            longitude = ser['lon']
            collection_date = None
            collection_depth = None

            self.rows.append([sample_name, fastq_fwd_file_path, fastq_rev_file_path, sample_type, host_phylum, host_class, host_order, host_family, host_genus, host_species, latitude, longitude, collection_date, collection_depth])

        print(f'There were {sum(non_three_list.values())} coral samples collected that were not of the normal three species.')
        print(non_three_list)
        spds_df = pd.DataFrame(self.rows, columns=self.datasheet_cols)
        df_output_path = os.path.join(self.out_dir, 'coral_only_sp_data_sheet.csv')
        spds_df.to_csv(df_output_path, index=False)
        
        # Read in the csv file and add the top two rows that contain the master headers and headers
        with open(df_output_path, 'r') as f:
            lines = [line.rstrip() for line in f]
        # insert in reverse order
        header_row_one = ',,,,host_info if applicable,,,,,,sampling info if applicable,,,'
        lines.insert(0, header_row_one)

        # Write out the csv with the new headers added
        with open(df_output_path, 'w') as f:
            for line in lines:
                f.write(f'{line}\n')

        lines_for_df = [line.split(',') for line in lines]
        df = pd.DataFrame(lines_for_df[2:], columns=lines_for_df[1])
        df.set_index(keys='sample_name', drop=True, inplace=True)
        return df

class GenerateInfoCollectionDF:
    """Class concerned with creating the information dataframe"""
    def __init__(self, base_dir, provenance_df):
        self.base_dir = base_dir
        self.provenance_df = provenance_df
        self.rows = []

    def generate_df(self):
        self._create_data_rows()
        return pd.DataFrame(data=self.rows, columns=['sample_name', 'fastq_fwd_file_path', 'fastq_rev_file_path', 'coral_plankton', 'spp_water', 'location', 'site', 'lat', 'lon']).set_index(keys='sample_name')
    
    def _create_data_rows(self):
        """ Parse through the seq data file structure gathering 
        inferring the info for each sequecning file pair as we descend the structure/
        Collect a single row of data for each sample"""
        for location in os.listdir(self.base_dir):
            if 'ISLAND' in location:
                parsing_dir_loc = os.path.join(self.base_dir, location)
                for site in os.listdir(parsing_dir_loc):
                    parsing_dir_site = os.path.join(parsing_dir_loc, site)
                    for sample_type in os.listdir(parsing_dir_site):
                        parsing_dir_sample_type = os.path.join(parsing_dir_site, sample_type)
                        if sample_type == 'CORAL':
                            for species in os.listdir(parsing_dir_sample_type):
                                parsing_dir_species = os.path.join(parsing_dir_sample_type, species)
                                for individual in os.listdir(parsing_dir_species):
                                    parsing_dir_indi = os.path.join(parsing_dir_species, individual, 'CS4L')
                                    # now we are in the directory that contains the actual paired fastq.gz files for a
                                    # given coral individual
                                    # collect the information we need
                                    srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_indi, sample_type=sample_type, site=site, water_species=species, provenance_df=self.provenance_df)
                                    self.rows.append(srg.create_sample_row())
                                    print(f'Processed: {parsing_dir_indi}')

                        elif sample_type == 'PLANKTON':
                            for water_type in os.listdir(parsing_dir_sample_type):
                                if water_type == 'CSW':
                                    parsing_dir_water_type = os.path.join(parsing_dir_sample_type, water_type)
                                    for individual in os.listdir(parsing_dir_water_type):
                                        parsing_dir_indi = os.path.join(parsing_dir_water_type, individual, 'S320')
                                        # now we are in the directory that contains the actual paired fastq.gz files for a
                                        # given water sample
                                        # collect the information we need
                                        srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_indi, sample_type=sample_type, site=site, water_species=water_type, provenance_df=self.provenance_df)
                                        self.rows.append(srg.create_sample_row())
                                        print(f'Processed: {parsing_dir_indi}')

                                elif water_type == 'SURFACE':
                                    # then this is a SURFACE sample and there are no individuals
                                    parsing_dir_water_type = os.path.join(parsing_dir_sample_type, water_type, 'S320')
                                    # collect the information we need
                                    srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_water_type, sample_type=sample_type, site=site, water_species=water_type, provenance_df=self.provenance_df)
                                    self.rows.append(srg.create_sample_row())
                                    print(f'Processed: {parsing_dir_water_type}')

            elif 'OA' in location:
                parsing_dir_loc = os.path.join(self.base_dir, location, 'PLANKTON', 'SURFACE', 'S320')
                # NB the arguments are a little messed up here on purpose as we don't have seperate site and location
                srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_loc, sample_type='OA', site=location, water_species='PLANKTON', provenance_df=self.provenance_df)
                self.rows.append(srg.create_sample_row())
                print(f'Processed: {parsing_dir_loc}')


class SampleRowGenerator:
    def __init__(self, location, parsing_dir, sample_type, site, water_species, provenance_df):
        self.parsing_dir = parsing_dir
        self.location = location
        self.sample_type = sample_type
        self.site = site
        self.water_type = water_species
        self.provenance_df=provenance_df
        self.sample_name = '_'.join(os.listdir(self.parsing_dir)[0].split('_')[:2])
        self.lat = self.provenance_df.at[self.sample_name, 'lat']
        self.lon = self.provenance_df.at[self.sample_name, 'lon']
        self.files = os.listdir(self.parsing_dir)
        if len(self.files) != 2:
            self._concatenate_seq_files()
            self.files = os.listdir(self.parsing_dir)
            if len(self.files) != 2:
                raise RuntimeError(f'Error in concatenation of seq files in {self.parsing_dir}')
        self.fwd_found = False
        self.rev_found = False
        self.fwd_path = None
        self.rev_path = None
        self._get_seq_paths()

    def _concatenate_seq_files(self):
        """ Sometime there is more than 1 set of sequencing files. In this case
        we want to concatenate all of the R1 reads together and all of the R2 reads
        together. Importantly we want to merge the respective memebers of the sequencing
        pairs in the same order.
        We will get a list of R1 files in the director (arbitrary order)
        We will then get a corresponding list of R2 files in the same realtive order
        by doing a replace function on the R1 list.
        We will then concatenate each list of files. Importantly, we will do this
        directly on the .gz files as this is apparently OK.
        https://www.biostars.org/p/81924/
        """
        # get files
        r1_files_list = [file_name for file_name in self.files if 'R1' in file_name]
        r2_files_list = [file_name.replace("R1", "R2") for file_name in r1_files_list]
        
        # Create the list that holds the cat command and the arguments
        exe_1 = [os.path.join(self.parsing_dir, _) for _ in r1_files_list]
        exe_1.insert(0, 'cat')
        exe_2 = [os.path.join(self.parsing_dir, _) for _ in r2_files_list]
        exe_2.insert(0, 'cat')
        
        # outpaths that will be used to create a stream that will be passed as stdout
        out_path_1 = os.path.join(self.parsing_dir, r1_files_list[0].replace("R1.fastq.gz", "merged.R1.fastq.gz"))
        out_path_2 = os.path.join(self.parsing_dir, r2_files_list[0].replace("R2.fastq.gz", "merged.R2.fastq.gz"))

        # do cat
        with open(out_path_1, 'wb') as f:
            subprocess.run(exe_1, stdout=f)
        with open(out_path_2, 'wb') as f:
            subprocess.run(exe_2, stdout=f)
        
        # rm the old files that have now been concatenated
        for file_path in exe_1[1:]:
            os.remove(file_path)
        for file_path in exe_2[1:]:
            os.remove(file_path)

    def _get_seq_paths(self):
        for file_name in self.files:
                if 'R1' in file_name:
                    self.fwd_path = os.path.join(self.parsing_dir, file_name)
                    self.fwd_found = True
                elif 'R2' in file_name:
                    self.rev_path = os.path.join(self.parsing_dir, file_name)
                    self.rev_found = True

    def create_sample_row(self):
        if not self.fwd_found or not self.rev_found:
            print('fwd or rev read not found')
            sys.exit(1)
        else:
            return [self.sample_name, self.fwd_path, self.rev_path, self.sample_type, self.water_type, self.location, self.site, self.lat, self.lon]

def human_readable_size(size, decimal_places=3):
    for unit in ['B','KiB','MiB','GiB','TiB']:
        if size < 1024.0:
            break
        size /= 1024.0
    return f"{size:.{decimal_places}f}{unit}"

its2processing = ITS2Processing().start_walking()

