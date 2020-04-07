import os
import sys
import compress_pickle
import pandas as pd
from collections import defaultdict
class EighteenSBase:
    def __init__(self):
        self.root_dir = '/home/humebc/projects/tara/tara_full_dataset_processing'
        self.eighteens_dir = os.path.join(self.root_dir, '18s')
        self.input_dir = os.path.join(self.root_dir, 'input')
        # Directory where the fastq.gz 18s files are
        self.seq_dir = os.path.join(self.eighteens_dir, '18s_data')
        self.output_dir = os.path.join(self.root_dir, 'output')
        os.makedirs(self.output_dir, exist_ok=True)
        self.fig_output_dir = os.path.join(self.root_dir, 'figures')
        os.makedirs(self.fig_output_dir, exist_ok=True)
        # The directory where the finalised post qc and post taxa screening files will be written
        self.qc_dir = os.path.join(self.eighteens_dir, 'seq_qc')
        os.makedirs(self.qc_dir, exist_ok=True)
        self.cache_dir = os.path.join(self.root_dir, 'cache')
        os.makedirs(self.cache_dir, exist_ok=True)
        # This is the file that should be used if possible to get meta
        # data for samples
        self.sample_provenance_path = os.path.join(self.input_dir, "sample_provenance_20200201.csv")
        self.sample_provenance_df = self._make_sample_provenance_df()
        # This is the dataframe that was produced using the general_seq_processing.py
        # file (for the 18s). We will work through this to make sure that
        # we process the relevant directories
        self.fastq_info_df_path = os.path.join(self.input_dir, 'output_information_df_all_fastqs_18s_2020-04-06T08_43_59.156006UTC.csv')
        self.fastq_info_df = self._make_fastq_info_df()
        foo = 'bar'
        # The main info df that we will use
        # Sample name as key, fwd and rev path to seq files, coral species
        # TODO, we will see if we can walk away from relying on this and rather use the above
        # 
        # self.info_df = self._make_info_df()

    def _make_fastq_info_df(self):
        df = pd.read_csv(self.fastq_info_df_path)
        return df.set_index(keys='readset', drop=False)

    # def _make_info_df(self):
    #     try:
    #         return compress_pickle.load(os.path.join(self.cache_dir, 'info_df.p.bz'))
    #     except FileNotFoundError:
    #         info_df = MakeInfoDF(seq_dir=self.seq_dir, sample_provenance_df=self.sample_provenance_df).make_info_df()
    #         compress_pickle.dump(info_df, os.path.join(self.cache_dir, 'info_df.p.bz'))
    #         return info_df

    def _make_sample_provenance_df(self):
        # The SAMPLE ID is in the sequencing file name so we will be able to use this to get the latitute
        df = pd.read_csv(self.sample_provenance_path, skiprows=[1])
        df.set_index(keys='BARCODE ID', drop=True, inplace=True)
        # df.rename(columns={'EVENT LATITUDE START': 'lat', 'EVENT LONGITUDE START': 'lon'}, inplace=True)
        return df

# class MakeInfoDF:
#     def __init__(self, seq_dir, sample_provenance_df):
#         self.paired_files_dict = defaultdict(list)
#         self.seq_dir = seq_dir
#         self.seq_files = [file for file in os.listdir(self.seq_dir) if 'fastq' in file]
#         self.sample_provenance_df = sample_provenance_df
#         # Get a list of the samples that we have seq files for
#         self.sample_names = self._get_sample_names_list()

#     def make_info_df(self):
#         # for each of the sample names, find the respective seq files, pair them up and populate the paired_files dict
#         self._find_sample_respective_seq_files()

#         # here we have the dictionary that tells us which sample names relate to which seqfiles populated
#         # Now populate the info dict splitting up those samples that have multiple pairs of seq files
#         # Dict that we will create info df from eventually
#         return self._return_info_df()

#     def _return_info_df(self):
#         info_df_dict = {}
#         for sample_id in self.paired_files_dict.keys():
#             species = self.sample_provenance_df.at[sample_id, 'Sample Material label, organismal system level, taxonomic, nominal']
#             path_lists = self.paired_files_dict[sample_id]
#             island = self.sample_provenance_df.at[sample_id, 'ISLAND#']
#             site = self.sample_provenance_df.at[sample_id, 'SITE#']
#             if len(path_lists) > 1:
#                 for i in range(len(path_lists)):
#                     fwd_path = os.path.join(self.seq_dir, path_lists[i][0])
#                     rev_path = os.path.join(self.seq_dir, path_lists[i][1])
#                     info_df_dict[f'{sample_id}_{i}'] = [fwd_path, rev_path, species, island, site]
#             else:
#                 info_df_dict[sample_id] = [os.path.join(self.seq_dir, path_lists[0][0]), os.path.join(self.seq_dir, path_lists[0][1]), species, island, site]
#         return pd.DataFrame.from_dict(info_df_dict, orient='index', columns=['fwd_path', 'rev_path', 'species', 'island', 'site'])

#     def _find_sample_respective_seq_files(self):
#         for sample_name in self.sample_names:
#             matching_seq_files = []
#             for seq_file in self.seq_files:
#                 if sample_name in seq_file:
#                     matching_seq_files.append(seq_file)
#             # here we have a list of all of the seq files
#             # if there are only two, then this is simple
#             if not matching_seq_files:
#                 raise RuntimeError(f'No sequencing files found for sample {sample_name}')
#             if len(matching_seq_files) == 2:
#                 if '_R1.' in matching_seq_files[0]:
#                     self.paired_files_dict[sample_name] = [[matching_seq_files[0], matching_seq_files[1]]]
#                 else:
#                     self.paired_files_dict[sample_name] = [[matching_seq_files[1], matching_seq_files[0]]]
#             else:
#                 # if there are more than two then we need to pair them up
#                 for seq_file in [seq_file for seq_file in matching_seq_files if '_R1.' in seq_file]:
#                     self.paired_files_dict[sample_name].append([seq_file, seq_file.replace('_R1.', '_R2.')])

#     def _get_sample_names_list(self):
#         sample_names = set()
#         for seq_file in self.seq_files:
#             sample_names.add('_'.join(seq_file.split('_')[:2]))
#         sample_names = list(sample_names)
#         return sample_names

if __name__ == "__main__":
    EighteenSBase()