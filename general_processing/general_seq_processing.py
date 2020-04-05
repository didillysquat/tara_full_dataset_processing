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
import compress_pickle

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

    # We need to update so that we are using the new tara_sample_provenance table.
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


    def __init__(self, marker, seq_file_download_directory=None, date_string=None, download=False):
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
        if date_string:
            dat_string = date_string
        else:
            dat_string = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f%Z").replace(':', '_')
        self.df_dict = {}
        self.sp_coral_datasheet_path = os.path.join(self.output_dir, f'sp_coral_datasheet_{dat_string}.csv')
        self.non_coral_sp_datasheet_df_dict = {}
        self.sp_non_coral_datasheet_path = os.path.join(self.output_dir, f'sp_non_coral_datasheet_{dat_string}.csv')
        # The dict that will hold the info for the information df to document which files we used and which we did not
        self.output_information_df = {}
        self.marker = marker
        self.output_information_df_path = os.path.join(
            self.output_dir, f'output_information_df_all_fastqs_{self.marker}_{dat_string}.csv')
        self. output_information_df_cols = [
            'barcode_id', 'readset', 'fwd_read_name', 'rev_read_name', 'use', 'URL', 'is_replicate',
            'access_time', 'replication_category', 'replication_color', 'fwd_read_size_compressed_bytes',
            'rev_read_size_compressed_bytes', 'pcr_code', 'dna_extraction_code', 'genescope_comment_1', 'genescope_comment_2']
        # We can download the files that we are going to keep while we're at it
        # We should save them to a single directory

        if self.marker == 'its2':
            self.remote_base_dir = "https://www.genoscope.cns.fr/sadc/tarapacific/METABARCODING/ITS2/ITS2_SYM_VAR_5.8S2_SYM_VAR_REV/"
            self.seq_file_download_directory = seq_file_download_directory
            if not download:
                self.download_data = False
            else:
                self.download_data = True
            self.no_keep_info_red_barcodes_list = ['CO-0002385', 'CO-0004425']
        elif self.marker == '18s':
            self.remote_base_dir = "https://www.genoscope.cns.fr/sadc/tarapacific/METABARCODING/18S_V9/18S_V9_1389F_1510R/"
            if download:
                self.download_data = True
                self.seq_file_download_directory = seq_file_download_directory
            else:
                self.download_data = False
            self.no_keep_info_red_barcodes_list = [
                'AW-0000001', 'AW-0000015', 'AW-0000029', 'AW-0000032', 'AW-0000042', 'AW-0000048', 'AW-0000136',
                'CO-0001018', 'CO-0001954', 'CO-0001958', 'CO-0003055', 'CO-0003361', 'CO-0005681', 'CO-0005914',
                'CO-1014469', 'CO-1014900', 'CO-1018002', 'CO-1018836', 'CO-1018837',
                'IW-0000009', 'IW-0000034', 'IW-0000039', 'IW-0000045', 'IW-0000047', 'IW-0000072', 'IW-0000077',
                'OA-0000011', 'OA-0000027', 'OA-0002159', 'OA-0002546'
            ]
        elif self.marker == '16s_45':
            self.remote_base_dir = "https://www.genoscope.cns.fr/sadc/tarapacific/METABARCODING/16S_V4V5/Fuhrman_primers/"
            self.download_data = False
            self.no_keep_info_red_barcodes_list = []
        elif self.marker == '16s_full_45':
            self.remote_base_dir = "https://www.genoscope.cns.fr/sadc/tarapacific/METABARCODING/16S_Full_Length_plus_16S_V4V5/16S_FL_27F_1492R_plus_Fuhrman_primers/"
            self.download_data = False
            self.no_keep_info_red_barcodes_list = [
                'CO-0002037', 'CO-0002038', 'CO-0002039', 'CO-0002040', 'CO-0002043', 'CO-0002046',
                'CO-0002048', 'CO-0002049', 'CO-0002050', 'CO-0002051', 'CO-0002052', 'CO-0002053',
                'CO-0002059', 'CO-0002060', 'CO-0002061', 'CO-0002062', 'CO-0002063', 'CO-0002064', 'CO-0002065',
                'FH-0000297', 'FH-0000303', 'FH-0000308', 'FH-0000313', 'FH-0000318', 'FH-0000323',
                'FH-0000328', 'FH-0000333', 'FH-0000338', 'FH-0000343', 'FH-0000348', 'FH-0000353',
            ]
        self.exe_path = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.authorisation_tup = self._make_auth_tup()
        self.headers = {'User-Agent': 'Benjamin Hume', 'From': 'benjamin.hume@kaust.edu.sa'}
        # Var to collect the output of the MP processing
        self.mp_output_list_of_tups = None
        self.sp_datasheet_cols = ['fastq_fwd_file_name', 'fastq_rev_file_name', 'sample_type', 'host_phylum',
         'host_class', 'host_order', 'host_family', 'host_genus', 'host_species', 'collection_latitude',
         'collection_longitude',
         'collection_date', 'collection_depth']

    def _make_readset_info_dir(self):
        # read in the three sepearate csv files
        coral_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "coral_readset_info.csv"), skiprows=[0],
                                       names=['readset', 'primers', 'barcode_id', 'pcr_code',
                                              'dna_extraction_code'], dtype={'readset': str, 'primers': str, 'barcode_id': str, 'pcr_code': str, 'dna_extraction_code': str})
        sed_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "ssed_readset_info.csv"), skiprows=[0],
                                     names=['readset', 'primers', 'barcode_id', 'pcr_code',
                                            'dna_extraction_code'], dtype={'readset': str, 'primers': str, 'barcode_id': str, 'pcr_code': str, 'dna_extraction_code': str})
        fish_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "fish_readset_info.csv"), skiprows=[0], names=['readset', 'primers', 'barcode_id', 'pcr_code', 'pcr_fl_sample_name',
                                   'dna_extraction_code'], dtype={'readset': str, 'primers': str, 'barcode_id': str, 'pcr_code': str, 'dna_extraction_code': str, 'pcr_fl_sample_name': str})
        # fish_readset_df.drop(columns='PCR FL sample name', inplace=True)
        # fish_readset_df.columns = ['readset', 'primers', 'barcode_id', 'pcr_code', 'pcr_fl_sample_name',
        #                            'dna_extraction_code']
        plankton_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "plankton_readset_info.csv"), skiprows=[0], names=['readset', 'primers', 'barcode_id', 'pcr_code', 'pcr_fl_sample_name',
                                   'dna_extraction_code'], dtype={'readset': str, 'primers': str, 'barcode_id': str, 'pcr_code': str, 'dna_extraction_code': str, 'pcr_fl_sample_name': str})
        # # plankton_readset_df.drop(columns='PCR FL sample name', inplace=True)
        # plankton_readset_df.columns = ['readset', 'primers', 'barcode_id', 'pcr_code', 'pcr_fl_sample_name',
        #                                'dna_extraction_code']
        df = pd.concat([coral_readset_df, sed_readset_df, fish_readset_df, plankton_readset_df])
        df = df.set_index('readset', drop=True)
        # Here add in the comments for Julie.
        juli_comment_df = pd.read_csv('comments_from_julie.csv', index_col=0)
        new_df = pd.concat([df, juli_comment_df], axis=1)
        return new_df

    def _make_auth_tup(self):
        auth_path = os.path.join(self.exe_path, 'auth.txt')
        with open(auth_path, 'r') as f:
            auth_lines = [line.rstrip() for line in f]
        return (auth_lines[0], auth_lines[1])

    def start_walking(self):
        if os.path.isfile(os.path.join(self.cache_dir, f'mp_output_list_of_tups_{self.marker}X.p.bz')):
            self.mp_output_list_of_tups = compress_pickle.load(os.path.join(self.cache_dir, f'mp_output_list_of_tups_{self.marker}.p.bz'))
        else:
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
            # worker_base_dirs = [_ for _ in worker_base_dirs if 'ISLAND06' in _]

            # NB We were originally mapping the rep_walker_list directly to the ReplicationWalkerWork class and running
            # its _walk function from within the __init__. However this was causing problems when running
            # for some of the markers and giving us an issue about errors returning the error and recursion (Pool error).
            # Instead we now pass in instances of the ReplicationWalkerWorker class and then run its _walk method, returning
            # only the list.
            # Create a ReplicationWalker for every worker_base_dir
            rep_walker_list = []
            for w_dir in worker_base_dirs:
                rep_walker_list.append(ReplicationWalkerWorker(w_dir, marker=self.marker, prov_df=self.sample_provenance_df, readset_df=self.readset_df))
            with Pool(20) as p:
                self.mp_output_list_of_tups = p.map(self._run_walk_on_rep_walker_item, rep_walker_list)
            compress_pickle.dump(self.mp_output_list_of_tups, os.path.join(self.cache_dir, f'mp_output_list_of_tups_{self.marker}.p.bz'))

        # At this point we will have the info required to make the sp_data sheets and output the
        # information dataframe
        # The information comes out as a single list, each element in the list is a tuple of three elements
        # The three elements are in the order of
        # self.coral_sp_datasheet_df_dict, self.non_coral_sp_datasheet_df_dict, self.output_information_list

        # Firstly, do sanity check to  make sure that all of the barcodes that we don't have keep info for
        # have been visited
        visited = [list_item for tup in self.mp_output_list_of_tups for list_item in tup[3]]
        if set(visited) != set(self.no_keep_info_red_barcodes_list):
            raise RuntimeError('no_keep lists do not match')
        if self.marker == 'its2':
            # SP coral datasheet
            self._sp_df_to_csv(tup_index=0, csv_path=self.sp_coral_datasheet_path)

            # SP non-coral datasheet
            self._sp_df_to_csv(tup_index=1, csv_path=self.sp_non_coral_datasheet_path)

        # output information dataframe
        self._output_output_information_df()

        # Finally output a list of the barcodes that are on the server but not in Julies tables
        problem_b_codes = [list_item for tup in self.mp_output_list_of_tups for list_item in tup[4]]
        if problem_b_codes:
            print('barcode_id s that were on server but not in Julie tables:')
            for b_code in problem_b_codes:
                print(b_code)

    def _output_output_information_df(self):
        self.output_information_df = pd.DataFrame(
            [_ for t in self.mp_output_list_of_tups for _ in t[2]], columns=self.output_information_df_cols
        )
        self.output_information_df.to_csv(self.output_information_df_path, index=False)

    def _sp_df_to_csv(self, tup_index, csv_path):
        self.df_dict = {k: v for t in self.mp_output_list_of_tups for k, v in t[tup_index].items()}
        sp_datasheet_df = pd.DataFrame.from_dict(
            self.df_dict, columns=self.sp_datasheet_cols, orient='index'
        )
        sp_datasheet_df.index.name = 'sample_name'
        sp_datasheet_df.to_csv(csv_path, index=True)
        self._add_master_headers_to_sp_datasheet(csv_path)

    def _add_master_headers_to_sp_datasheet(self, output_path):
        # Read in the csv file and add the top two rows that contain the master headers and headers
        with open(output_path, 'r') as f:
            lines = [line.rstrip() for line in f]
        # insert in reverse order
        header_row_one = ',,,,host_info if applicable,,,,,,sampling info if applicable,,,'
        lines.insert(0, header_row_one)
        # Write out the csv with the new headers added
        with open(output_path, 'w') as f:
            for line in lines:
                f.write(f'{line}\n')

    @staticmethod
    def _run_walk_on_rep_walker_item(rep_walker_class_instance):
        return rep_walker_class_instance._walk()

    def _make_sample_provenance_df(self):
        # The BARCODE ID is in the sequencing file name so we will be able to use this to get the latitute
        df = pd.read_csv(self.sample_provenance_path, skiprows=[1])
        df.set_index(keys='BARCODE ID', drop=True, inplace=True)
        return df

class ReplicationWalkerWorker:
    def __init__(self, remote_base_dir, marker, prov_df=None, readset_df=None):
        self.marker = marker
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
        # Two versions of the fastq_gz_list_current
        # This one refers to the list in the directory
        self.fastq_gz_list_current_dir = []
        # This one refers to the list that belong to a specific barcode within the directory (updated dynamically)
        self.fastq_gz_list_current_barcode = []
        self.links_to_visit_current = []
        self.last_fork_location = None
        self.home_dir_reached = None
        self.s = requests.Session()
        self.s.auth = self.authorisation_tup
        self.s.headers = self.headers
        self.no_keep_info_red_barcodes_list_visited = []
        self.coral_sp_datasheet_df_dict = {}
        self.non_coral_sp_datasheet_df_dict = {}
        if self.marker == 'its2':
            # Two dictionaries that will hold the information for creating the dataframes form that will become the
            # symportal datasheets for doing the loading

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
            self.download_data = True
        elif self.marker == '18s':
            self.download_data = False
            # For the 18S things are a little tricker than the ITS2. For some barcodes,
            # there are multiple not to use and multiple that we can use.
            # As such we'll have to adapt the code to look for multiple readsets
            self.fwd_readset_strings_to_keep = [
                'HMNC2BBXX.12BA193',
                'BWFDM.12BA159-BID05',
                'BWFDM.12BA159-BID09',
                'HMLTWBBXX.12BA258',
                'HCYFJBBXX.12BA111',
                'HCYFJBBXX.12BA135',
                'BWFDM.12BA171-BID08',
                'BWFDM.12BA171-BID06',
                'BWFDM.12BA171-BID09',
                'HGWFYBCX2.12BA001',
                'H2NNHBCX3.12BA289-BID01',
                'HGWFYBCX2.12BA002',
                'H2TTMBCX3.12BA001-BID01',
                'H2NNHBCX3.12BA295-BID01',
                'HGWFYBCX2.12BA003',
                'BWFDM.12BA159-BID08'
            ]
            # She gave us the reads to keep for all but two of the barcodes. The two barcodes that she didn't give us
            # keep reads for are:
            self.no_keep_info_red_barcodes_list = [
                'AW-0000001', 'AW-0000015', 'AW-0000029', 'AW-0000032', 'AW-0000042', 'AW-0000048', 'AW-0000136',
                'CO-0001018', 'CO-0001954', 'CO-0001958', 'CO-0003055', 'CO-0003361', 'CO-0005681', 'CO-0005914',
                'CO-1014469', 'CO-1014900', 'CO-1018002', 'CO-1018836', 'CO-1018837',
                'IW-0000009', 'IW-0000034', 'IW-0000039', 'IW-0000045', 'IW-0000047', 'IW-0000072', 'IW-0000077',
                'OA-0000011', 'OA-0000027', 'OA-0002159', 'OA-0002546'
            ]
        elif self.marker == '16s_45':
            self.download_data = False
            self.fwd_readset_strings_to_keep = [
                'HHHYMDRXX.12BA104-BID11',
                'HHKFGDRXX.12BA313-BID08',
                'H3J7FBCX2.12BA194',
                'H3J7FBCX2.12BA193',
                'H3J7FBCX2.12BA253',
                'H3J7FBCX2.12BA277',
                'H3J7FBCX2.12BA265',
                'H3J7FBCX2.12BA217',
                'H3J7FBCX2.12BA229'
            ]
            self.no_keep_info_red_barcodes_list = []
        elif self.marker == '16s_full_45':
            self.download_data = False
            self.fwd_readset_strings_to_keep = [
                'HFK2KBCX2.12BA172',
                'H7H7CBCX2.12BA106',
                'H7H7CBCX2.12BA118',
                'HFK2KBCX2.12BA101',
                'H7H7CBCX2.12BA130',
                'H7H7CBCX2.12BA142',
                'C8W52.12BA124',
                'C8W52.12BA123'
            ]
            self.no_keep_info_red_barcodes_list = [
                'CO-0002037', 'CO-0002038', 'CO-0002039', 'CO-0002040', 'CO-0002043', 'CO-0002046',
                'CO-0002048', 'CO-0002049', 'CO-0002050', 'CO-0002051', 'CO-0002052', 'CO-0002053',
                'CO-0002059', 'CO-0002060', 'CO-0002061', 'CO-0002062', 'CO-0002063', 'CO-0002064', 'CO-0002065',
                'FH-0000297', 'FH-0000303', 'FH-0000308', 'FH-0000313', 'FH-0000318', 'FH-0000323',
                'FH-0000328', 'FH-0000333', 'FH-0000338', 'FH-0000343', 'FH-0000348', 'FH-0000353',
            ]
        # TODO write code to sanity check to make sure that all of the no_keep info barcodes have been visited
        # The list that will hold the info for the information df to document which files we used and which we did not
        self.output_information_list = []
        # Keep track of the barcode ids that are on the genoscope webstie but are not in the Julie tables
        # and will not be considered in the outputs of this script.
        self.not_in_julie_tables_barcode_list = []


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
            self.fastq_gz_list_current_dir = [link.string for link in self.links_to_visit_current if
                                          'fastq.gz' in link.string]
            self.links_to_visit_current = [link for link in self.links_to_visit_current if
                                           link not in self.fastq_gz_list_current_dir]
            if len(self.links_to_visit_current) > 1:
                self.last_fork_location = self.current_remote_dir
            else:
                self.done_and_empty_list.add(self.current_remote_dir)
            if self.current_remote_dir == self.remote_base_dir and not self.links_to_visit_current:
                break

            if self.fastq_gz_list_current_dir:
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
        return (self.coral_sp_datasheet_df_dict, self.non_coral_sp_datasheet_df_dict, self.output_information_list, self.no_keep_info_red_barcodes_list_visited, self.not_in_julie_tables_barcode_list)

    def _document_fastq_files(self):
        # Then we can count how many there are, add current dir to done
        # and continue walking the directories
        if len(self.fastq_gz_list_current_dir) > 2:
            try:
                self._handle_multiple_fastq_files()
            except IndexError as e:
                foo = 'bar'
        else:
            if not len(self.fastq_gz_list_current_dir) == 2:
                raise RuntimeError('Odd number of fastq files')
            # Here we have only two files for a given barcode id and we can process accorrdingly
            try:
                self._handle_one_pair_fastq_files()
            except IndexError as e:
                foo = 'bar'

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
            fwd_read = [_ for _ in self.fastq_gz_list_current_dir if 'R1' in _][0]
            rev_read = [_ for _ in self.fastq_gz_list_current_dir if 'R2' in _][0]
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
            # For the 18S ther are some atmosphere samples eg. barcode G-0000185, that are not listed in Julie's
            # tables. I will write into the code to ignore these. If we fail to find an good index_list I will
            # check to see if the barcode_id is found in the Julies tables. If it is found, then we have a logic
            # problem. But if its not find, then we know that this is simply a sample that isn't in Julies
            # tables and so we should simply return from this method
            if len(self.readset_df[self.readset_df['barcode_id'] == barcode_id].index) == 0:
                print(f'barcode_id {barcode_id} does not appear to be in the Julie tables but is on the genoscope server')
                self.not_in_julie_tables_barcode_list.append(barcode_id)
                return
            else:
                raise RuntimeError
        else:
            readset = index_list[0]

        if use and self.download_data:
            size_dict = self._download_file_if_necessary(fwd_read, rev_read)
        else:
            if size_dict_passed:
                size_dict = size_dict_passed
            else:
                size_dict = self._get_sizes_trough_head_request(fwd_read, rev_read)
        #TODO add the pcr name and extraction nme
        pcr_code = self.readset_df.at[readset, 'pcr_code']
        dna_extraction_code = self.readset_df.at[readset, 'dna_extraction_code']
        # Now we can populate the output information dict and the sp_datasheet
        self._populate_output_information_list(barcode_id, fwd_read, readset, rev_read, size_dict, use, is_rep, cat, col, pcr_code, dna_extraction_code)

        if use and self.marker == 'its2':
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

    def _populate_output_information_list(self, barcode_id, fwd_read, readset, rev_read, size_dict, use, is_rep, cat, col, pcr_code, dna_extraction_code):
        # First the output information dict
        dat_string = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f%Z")
        genescope_comment_1 = self.readset_df.at[readset, 'genescope_comment_1']
        genescope_comment_2 = self.readset_df.at[readset, 'genescope_comment_2']
        temp_list = [
            barcode_id, readset, fwd_read, rev_read, use, self.current_remote_dir,
            is_rep, dat_string, cat, col, size_dict[fwd_read], size_dict[rev_read], pcr_code, dna_extraction_code,
            genescope_comment_1, genescope_comment_2
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
        barcode_id_set = set([_.split('_')[1] for _ in self.fastq_gz_list_current_dir])
        if 'CO-0004552' in barcode_id_set:
            foo = 'bar'
        if len(self.fastq_gz_list_current_dir) > len(barcode_id_set) * 2:
            # We should be careful to run this on a barcode_id basis
            # As there may be more than one barcode in a directory
            # And there may be different replication statuses per barcode
            for barcode_id in barcode_id_set:
                if len(barcode_id_set) > 1:
                    foo = 'bar'
                # Hopefully it is enough to simply change the fastq_gz list as this is what we
                # infre all other information from
                self.fastq_gz_list_current_barcode = [_ for _ in self.fastq_gz_list_current_dir if barcode_id in _]
                # This is where we need to categorise the replication as
                # one of three different categories,
                # When all seq files have the same base_name, then we will call
                # this 'sequencing_replicate_same_run'.
                # Where the base_name is different, but the pcr and dna names are the
                # same, we will call this 'sequencing_replicate_different_run'
                # Finally where there are any differences in pcr or dna names
                # I will call these 'method_replicate'
                base_names = {'-'.join(_.split('-')[:-1]) for _ in self.fastq_gz_list_current_barcode}
                if '' in base_names:
                    print(self.fastq_gz_list_current_barcode)
                    foo = 'bar'
                if len(base_names) > 1:
                    try:
                        self._process_unkn_method_replicate(barcode_id_set)
                    except IndexError as e:
                        foo = 'bar'
                else:
                    try:
                        self._process_sequencing_replication_same_run()
                    except IndexError as e:
                        foo = 'bar'
        else:
            # here we have a case of multiple barcode ids worth of fastq pairs
            # The easiest way to treat this situation is to is to split it up into pairs of fastq
            # files and send these into the
            for barcode_id in barcode_id_set:
                barcode_reads = [_ for _ in self.fastq_gz_list_current_dir if barcode_id in _]
                fwd_read = [_ for _ in barcode_reads if 'R1' in _][0]
                rev_read = [_ for _ in barcode_reads if 'R2' in _][0]
                self._handle_one_pair_fastq_files(read_tup=(fwd_read, rev_read))

    def _process_unkn_method_replicate(self, barcode_id_set):
        # Then we need to check their PCR_sample_name and DNA_sample_name
        # To see if they are different.
        # This will produce two classes of seq difference
        # One where ReadSet names are different, but PCR and DNA are the same
        # Have to ask why?
        # The other set will be where ReadSet names are different, and one
        # of the PCR or DNA names are different.
        pcr_names_set = set()
        dna_names_set = set()
        for fastq_fwd in [_ for _ in self.fastq_gz_list_current_barcode if 'R1' in _]:
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
            pcr_code = self.readset_df.at[readset, 'pcr_code']
            dna_extraction_code = self.readset_df.at[readset, 'dna_extraction_code']


            pcr_names_set.add(pcr_code)
            dna_names_set.add(dna_extraction_code)

        # print(fastq_gz_list)
        if (len(pcr_names_set) != len(dna_names_set)) or (len(pcr_names_set) > len(barcode_id_set)):
            self._log_method_replicate()
        else:
            self._log_sequencing_replication_different_run()

    def _log_sequencing_replication_different_run(self):
        # Then this a sequencing_replicate_different_run
        # Here we want to do the same as when we were doing _process_seq_replication
        # however, there may be more than two sets of sequences. We're looking to keep the biggest pair
        # TODO, there can only be a case of yellow/green if there are more than 3 sets of fastqfiles
        # so we will check for this

        size_dict = {}
        for fastq_fwd in [_ for _ in self.fastq_gz_list_current_barcode if 'R1' in _]:
            size_dict[fastq_fwd] = sum(self._get_sizes_trough_head_request(
                fwd_read=fastq_fwd, rev_read=fastq_fwd.replace('R1', 'R2')
            ).values())
        # now simply sort to get the largest fwd read and submit that as keep, submit all others as no keep
        fwd_fastq_to_keep = sorted(size_dict, key=size_dict.get, reverse=True)[0]
        list_of_fwd_reads = [_ for _ in self.fastq_gz_list_current_barcode if 'R1' in _]
        list_of_fwd_base_names = ['-'.join(_.split('-')[:-1]) for _ in list_of_fwd_reads]
        for fastq_fwd in list_of_fwd_reads:
            individual_fwd_base_name = '-'.join(fastq_fwd.split('-')[:-1])
            if list_of_fwd_base_names.count(individual_fwd_base_name) > 1:
                # Then this is a sequencing replicate same run
                if fastq_fwd == fwd_fastq_to_keep:
                    # This is the keep
                    self._handle_one_pair_fastq_files(
                        read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=True, is_rep=True,
                        cat='sequencing_replicate_different_run/sequencing_replicate_same_run', col='yellow/green'
                    )
                else:
                    # this is a no keep
                    self._handle_one_pair_fastq_files(
                        read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=False, is_rep=True,
                        cat='sequencing_replicate_different_run/sequencing_replicate_same_run', col='yellow/green'
                    )
            elif list_of_fwd_base_names.count(individual_fwd_base_name) == 1:
                # Then this is not a sequencing replicate same run
                if fastq_fwd == fwd_fastq_to_keep:
                    # This is the keep
                    self._handle_one_pair_fastq_files(
                        read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=True, is_rep=True,
                        cat='sequencing_replicate_different_run', col='yellow'
                    )
                else:
                    # this is a no keep
                    self._handle_one_pair_fastq_files(
                        read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=False, is_rep=True,
                        cat='sequencing_replicate_different_run', col='yellow'
                    )

    def _log_method_replicate(self):
        # then this is a 'method_replicate'
        # We can use the elements in the temp_list to submit each of these
        # to the single fastq pair handler. But before we do that we need to identify readset (as this is the
        # UID essentially) that corresponds to the largest set of the files.
        # Actually we need to check the do not use lists and check to see if we have identified the correct
        # TWO samples that Julie said we can check for.s
        # TODO in the case where there is sequence replication as well as method replication
        # we will mark as red/green and method_replicate/sequencing_replicate
        # TODO check to see if there are anycases of red/yellow. and yellow/green

        # First check whether we are dealing with one of the barcodes that we don't have keep information for
        set_of_barcodes = set([_.split('_')[1] for _ in self.fastq_gz_list_current_barcode])
        if not len(set_of_barcodes) == 1:
            raise RuntimeError
        barcode_id = list(set_of_barcodes)[0]
        if barcode_id in self.no_keep_info_red_barcodes_list:
            # Then we are working with one of the barcodes that we don't have keep info for.
            self.no_keep_info_red_barcodes_list_visited.append(barcode_id)
            # We need to caculate the size and keep the biggest. Same as the sediment
            # We need to take into account that there may also be seq reps.
            size_dict = {}
            for fastq_fwd in [_ for _ in self.fastq_gz_list_current_barcode if 'R1' in _]:
                size_dict[fastq_fwd] = sum(self._get_sizes_trough_head_request(fastq_fwd, fastq_fwd.replace('R1', 'R2')).values())

            # Now we just want to keep the largest fastq_fwd in the size dict.
            # We also want to infer whether this is a sequencing replicate as well
            # This is not so easy as we may have several method sets that Julie has said we can keep.
            # Perhaps the best way to check whether there is a seq replicate
            # is to get the base name for the fastq_fwd in quetsion and then see how many of these are
            # in a list of basenames generated from the fastq_gz_list_current.
            # We will likely need to take into account BID. Let me check.
            # No we will not need to take into account BID.
            self._process_red_keep_no_keep(size_dict)

        else:
            # Then one of the fwd fastq files should be in the self.fwd_readset_strings_to_keep list
            # For each of the read pairs we need to grab the readset strings and see if one of the keeps
            # is in that list.
            # For the 18S things are a little more complicated. It is possible that there are a choice of serveral
            # readsets that can be kept. As such, we need to look for AT LEAST one matching readset.
            # Otherwise, the code actually doesn't need to change.
            readset_string_list = []
            for fastq_fwd in [_ for _ in self.fastq_gz_list_current_barcode if 'R1' in _]:
                if 'BID' in fastq_fwd:
                    element_one = fastq_fwd.split('-')[-3]
                    element_two = fastq_fwd.split('-')[-4].split('_')[-1]
                    bid_element = fastq_fwd.split('-')[-2]
                    readset_str = f'{element_two}.{element_one}-{bid_element}'
                else:
                    element_one = fastq_fwd.split('-')[-2]
                    element_two = fastq_fwd.split('-')[-3].split('_')[-1]
                    readset_str = f'{element_two}.{element_one}'
                readset_string_list.append(readset_str)
            # There should be at least one match
            if not len(set(readset_string_list).intersection(set(self.fwd_readset_strings_to_keep))) >= 1:
                # Then we are not finding a read that matches one of the keepers given by Juli
                raise RuntimeError
            # If we get here, then go back through and get the sizes of those seq pairs that have a readset that
            # matches the one given by Julie's
            size_dict = {}
            for fastq_fwd in [_ for _ in self.fastq_gz_list_current_barcode if 'R1' in _]:
                if 'BID' in fastq_fwd:
                    element_one = fastq_fwd.split('-')[-3]
                    element_two = fastq_fwd.split('-')[-4].split('_')[-1]
                    bid_element = fastq_fwd.split('-')[-2]
                    readset_str = f'{element_two}.{element_one}-{bid_element}'
                else:
                    element_one = fastq_fwd.split('-')[-2]
                    element_two = fastq_fwd.split('-')[-3].split('_')[-1]
                    readset_str = f'{element_two}.{element_one}'
                if readset_str in self.fwd_readset_strings_to_keep:
                    size_dict[fastq_fwd] = sum(self._get_sizes_trough_head_request(fastq_fwd, fastq_fwd.replace('R1', 'R2')).values())
            # Now we just want to keep the largest fastq_fwd in the size dict.
            self._process_red_keep_no_keep(size_dict)

    def _process_red_keep_no_keep(self, size_dict):
        fwd_fastq_to_keep = sorted(size_dict, key=size_dict.get, reverse=True)[0]
        list_of_fwd_reads = [_ for _ in self.fastq_gz_list_current_barcode if 'R1' in _]
        list_of_fwd_base_names = ['-'.join(_.split('-')[:-1]) for _ in list_of_fwd_reads]
        for fastq_fwd in list_of_fwd_reads:
            individual_fwd_base_name = '-'.join(fastq_fwd.split('-')[:-1])
            if list_of_fwd_base_names.count(individual_fwd_base_name) > 1:
                # Then this is a sequencing replicate same run
                if fastq_fwd == fwd_fastq_to_keep:
                    # This is the keep
                    self._handle_one_pair_fastq_files(
                        read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=True, is_rep=True,
                        cat='method_replicate/sequencing_replicate_same_run', col='red/green'
                    )
                else:
                    # this is a no keep
                    self._handle_one_pair_fastq_files(
                        read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=False, is_rep=True,
                        cat='method_replicate/sequencing_replicate_same_run', col='red/green'
                    )
            elif list_of_fwd_base_names.count(individual_fwd_base_name) == 1:
                # Then this is not a sequencing replicate same run
                if fastq_fwd == fwd_fastq_to_keep:
                    # This is the keep
                    self._handle_one_pair_fastq_files(
                        read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=True, is_rep=True,
                        cat='method_replicate', col='red'
                    )
                else:
                    # this is a no keep
                    self._handle_one_pair_fastq_files(
                        read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=False, is_rep=True,
                        cat='method_replicate', col='red'
                    )
            else:
                # Then something has gone wrong as we are not finding the base name in the list
                # of base names
                raise RuntimeError('inidividual base name not found in list of base names')

    def _process_sequencing_replication_same_run(self):
        # Then this is a case of sequencing_replicate_same_run
        # We should be able to find a readset per fastq.
        # The readset should contain two bits of information
        # in the fastq and it should containing the -1 or -2
        # This is a pain in the arse!
        size_dict = {}
        for fastq_fwd in [_ for _ in self.fastq_gz_list_current_barcode if 'R1' in _]:
            size_dict[fastq_fwd] = sum(self._get_sizes_trough_head_request(fastq_fwd, fastq_fwd.replace('R1', 'R2')).values())
        # Now we just want to keep the largest fastq_fwd in the size dict.
        fwd_fastq_to_keep = sorted(size_dict, key=size_dict.get, reverse=True)[0]
        for fastq_fwd in [_ for _ in self.fastq_gz_list_current_barcode if 'R1' in _]:
            if fastq_fwd == fwd_fastq_to_keep:
                # This is the keep
                self._handle_one_pair_fastq_files(
                    read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=True, is_rep=True,
                    cat='sequencing_replicate_same_run', col='green'
                )
            else:
                # this is a no keep
                self._handle_one_pair_fastq_files(
                    read_tup=(fastq_fwd, fastq_fwd.replace('R1', 'R2')), use=False, is_rep=True,
                    cat='sequencing_replicate_same_run', col='green'
                )

    def _make_auth_tup(self):
        auth_path = os.path.join(self.exe_path, 'auth.txt')
        with open(auth_path, 'r') as f:
            auth_lines = [line.rstrip() for line in f]
        return (auth_lines[0], auth_lines[1])

    def _make_readset_info_dir(self):
        # read in the three sepearate csv files
        coral_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "coral_readset_info.csv"), skiprows=[0],
                                       names=['readset', 'primers', 'barcode_id', 'pcr_code',
                                              'dna_extraction_code'])
        sed_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "ssed_readset_info.csv"), skiprows=[0],
                                     names=['readset', 'primers', 'barcode_id', 'pcr_code',
                                            'dna_extraction_code'])
        fish_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "fish_readset_info.csv"))
        # fish_readset_df.drop(columns='PCR FL sample name', inplace=True)
        fish_readset_df.columns = ['readset', 'primers', 'barcode_id', 'pcr_code', 'pcr_fl_sample_name',
                                   'dna_extraction_code']
        plankton_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "plankton_readset_info.csv"))
        # plankton_readset_df.drop(columns='PCR FL sample name', inplace=True)
        plankton_readset_df.columns = ['readset', 'primers', 'barcode_id', 'pcr_code', 'pcr_fl_sample_name',
                                       'dna_extraction_code']
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


def human_readable_size(size, decimal_places=3):
    for unit in ['B','KiB','MiB','GiB','TiB']:
        if size < 1024.0:
            break
        size /= 1024.0
    return f"{size:.{decimal_places}f}{unit}"

dat_string = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f%Z").replace(':', '_')
ITS2Processing(marker='its2', seq_file_download_directory="/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20200326_tara_its2_data", date_string=dat_string, download=True).start_walking()
ITS2Processing(marker='18s', seq_file_download_directory="/home/humebc/projects/tara/18s_data", date_string=dat_string, download=False).start_walking()
ITS2Processing(marker='16s_45', date_string=dat_string).start_walking()
ITS2Processing(marker='16s_full_45', date_string=dat_string).start_walking()


