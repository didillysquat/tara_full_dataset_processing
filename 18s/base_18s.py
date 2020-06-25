"""This is the base class for the other major 18S classes that are involved with 'doing'
the tara 18S analysis. This takes an output_information table that is produced by the
general_seq_processing.py script.
The general_seq_processing.py needs to be run before this can be used."""
import os
import sys
import compress_pickle
import pandas as pd
from collections import defaultdict
import subprocess
import re

class EighteenSBase:
    def __init__(self):
        self.root_dir = '/home/humebc/projects/tara/tara_full_dataset_processing'
        self.eighteens_dir = os.path.join(self.root_dir, '18s')
        self.output_dir_18s = os.path.join(self.eighteens_dir, 'output')
        self.input_dir = os.path.join(self.root_dir, 'input')
        # Directory where the fastq.gz 18s files are
        self.seq_dir = os.path.join(self.eighteens_dir, '18s_data')
        self.output_dir = os.path.join(self.root_dir, 'output')
        
        self.input_dir_18s = os.path.join(self.eighteens_dir, 'input')
        os.makedirs(self.output_dir, exist_ok=True)
        self.fig_output_dir = os.path.join(self.root_dir, 'figures')
        self.fig_output_dir_18s = os.path.join(self.eighteens_dir, 'figures')
        os.makedirs(self.fig_output_dir, exist_ok=True)
        # The directory where the finalised post qc and post taxa screening files will be written
        self.qc_dir = os.path.join(self.eighteens_dir, 'seq_qc')
        os.makedirs(self.qc_dir, exist_ok=True)
        self.cache_dir = os.path.join(self.root_dir, 'cache')
        self.cache_dir_18s = os.path.join(self.eighteens_dir, 'cache')
        os.makedirs(self.cache_dir, exist_ok=True)
        self.temp_dir_18s = os.path.join(self.eighteens_dir, 'temp')
        os.makedirs(self.temp_dir_18s, exist_ok=True)
        # This is the file that should be used if possible to get meta
        # data for samples
        self.sample_provenance_path = os.path.join(self.input_dir, "sample_provenance_20200201.csv")
        self.sample_provenance_df = self._make_sample_provenance_df()
        # This is the dataframe that was produced using the general_seq_processing.py
        # file (for the 18s). We will work through this to make sure that
        # we process the relevant directories
        self.fastq_info_df_path = os.path.join(self.input_dir, 'output_information_df_all_fastqs_18s_2020-04-09T07_11_24.231392UTC.csv')
        self.fastq_info_df = self._make_fastq_info_df()
        # List of the readsets that are from coral samples only
        self.coral_readsets = [readset for readset, ser in self.fastq_info_df.iterrows() if self.sample_provenance_df.at[ser['sample-id'] ,'SAMPLE ENVIRONMENT, short'] == 'C-CORAL']
        self.genera = ['Pocillopora', 'Millepora', 'Porites']

        # read in the classifications for the SNP
        self.poc_snp_classifications = self._get_snp_classifications(genus='Pocillopora')
        self.por_snp_classifications = self._get_snp_classifications(genus='Porites')
        
    def _get_snp_classifications(self, genus):
        if genus == 'Pocillopora':
            # First check to see if the cached version exists
            snp_cache_dir = os.path.join(self.input_dir_18s, 'snp_classifications', f'poc_snp_class_df.p.bz')
        elif genus == 'Porites':
            snp_cache_dir = os.path.join(self.input_dir_18s, 'snp_classifications', f'por_snp_class_df.p.bz')
        
        if os.path.exists(snp_cache_dir):
            return compress_pickle.load(snp_cache_dir)
        else:
            # Need to create it from scratch
            if genus == 'Pocillopora':
                raw_snp_class_path = os.path.join(self.input_dir_18s, 'snp_classifications', f'POC_SNP_classifications.csv')
            elif genus == 'Porites':
                raw_snp_class_path = os.path.join(self.input_dir_18s, 'snp_classifications', f'POR_SNP_classifications.csv')
        
            snp_class_df = pd.read_csv(raw_snp_class_path, index_col=0)
            snp_class_df.index = self._convert_index_to_sample_ids(snp_class_df.index)
            snp_class_df.dropna(inplace=True)
            snp_class_df.columns = ['label']
            compress_pickle.dump(snp_class_df, snp_cache_dir)
            return snp_class_df
            
    def _convert_index_to_sample_ids(self, index):
        # We want to convert these indices to samplie-id
        island_re = re.compile('I\d+')
        site_re = re.compile('S\d+')
        co_re = re.compile('C\d+')
        # A dict that converts from the current sample name (e.g. I01S01C011POR) to the proper sample-id
        # (e.g. TARA_CO-1016606)
        sample_name_dict = {}
        for ind in index:
            island = island_re.search(ind).group()
            site = site_re.search(ind).group()
            colony = co_re.search(ind).group()
            sample_id = self.sample_provenance_df[
                (self.sample_provenance_df['ISLAND#'] == island) & 
                (self.sample_provenance_df['SITE#'] == site) & 
                (self.sample_provenance_df['COLONY# (C000) FISH# (F000) MACROALGAE# (MA00)'] == colony) & 
                (self.sample_provenance_df['SAMPLE PROTOCOL LABEL, level 2'] == 'CS4L')].index.values.tolist()
            if len(sample_id) != 1:
                raise RunTimeError('More than one matching sample-id')
            else:
                sample_id = sample_id[0]
            sample_name_dict[ind] = sample_id
        
        # Convert the index to sample-id
        return [sample_name_dict[ind] for ind in index]

    @staticmethod
    def decompress_read_compress(path_to_read_in):
        """Because we now have all of the seq_qc files compressed we will make this utility
        class to decompress a given file, read it in to a list and then compress it again.
        Then return the list"""
        if '.gz' in path_to_read_in:
            # Then we have been handed a compressed file
            if os.path.isfile(path_to_read_in):
                subprocess.run(['gzip', '-d', path_to_read_in])
                path_to_compress = path_to_read_in.replace('.gz', '')
            else:
                raise RuntimeError(f'{path_to_read_in} does not exist')
        else:
            # Then we have been handed a non_compressed
            if os.path.isfile(path_to_read_in):
                # great it already exists, nothing to do
                path_to_compress = path_to_read_in
                pass
            elif os.path.isfile(path_to_read_in + '.gz'):
                # Then the target path is currently compressed
                subprocess.run(['gzip', '-d', path_to_read_in + '.gz'])
                path_to_compress = path_to_read_in
            else:
                # The file of interest doesn't seem to exist
                raise RuntimeError(f'could not find a path related to {path_to_read_in}')
        # now read it into a list
        with open(path_to_compress, 'r') as f:
            list_to_return = [line.rstrip() for line in f]
        
        # now compress the file
        subprocess.run(['gzip', path_to_compress])

        return list_to_return

    def _make_fastq_info_df(self):
        df = pd.read_csv(self.fastq_info_df_path)
        return df.set_index(keys='readset', drop=False)

    def _make_sample_provenance_df(self):
        # The SAMPLE ID is in the sequencing file name so we will be able to use this to get the latitute
        df = pd.read_csv(self.sample_provenance_path, skiprows=[1])
        df.set_index(keys='SAMPLE ID', drop=True, inplace=True)
        # df.rename(columns={'EVENT LATITUDE START': 'lat', 'EVENT LONGITUDE START': 'lon'}, inplace=True)
        return df

if __name__ == "__main__":
    EighteenSBase()