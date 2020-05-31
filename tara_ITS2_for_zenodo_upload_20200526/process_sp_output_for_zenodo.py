#!/usr/bin/env python3

"""
This script will take in the paths to the output directories related to the
coral and non-coral SymPortal outputs. It will process them in a way to end up with
a set of files ready for upload to zenodo.
"""
import os
import pandas as pd
import shutil
import subprocess

class ProcessSPOutputForZenodo:
    def __init__(self):
        self.base_dir = os.path.dirname(os.path.realpath(__file__))
        self.coral_output_dir_path = os.path.join(self.base_dir, "20200528_tara_corals")
        self.non_coral_output_dir_path = os.path.join(self.base_dir, "20200528_tara_non_corals")
        self.negative_output_dir_path = os.path.join(self.base_dir, "20200528_tara_negatives")
        # A dynamic directory path
        self.sub_dir = None

    def process(self):
        for p_dir in [self.coral_output_dir_path, self.non_coral_output_dir_path, self.negative_output_dir_path]:
            # Process the betweeen_sample_distances
            # This directory will exist in both the coral and non_coral
            self.sub_dir = os.path.join(p_dir, 'between_sample_distances')
            for clade_dir in os.listdir(self.sub_dir): # For each cladal directory
                self._p_between_dir(os.path.join(self.sub_dir, clade_dir))
            
            # Process the between_profile_distances
            # Only for coral samples
            if os.path.exists(os.path.join(p_dir, 'between_profile_distances')):
                self.sub_dir = os.path.join(p_dir, 'between_profile_distances')
                for clade_dir in os.listdir(self.sub_dir): # For each cladal directory
                    self._p_between_dir(os.path.join(self.sub_dir, clade_dir))

            # Process the post_med directory
            self.sub_dir = os.path.join(p_dir, 'post_med_seqs')
            self._p_post_med_dir()
    
            # Process the pre_med directory
            self.sub_dir = os.path.join(p_dir, 'pre_med_seqs')
            self._p_pre_med_dir()

            # Process the profile directory
            if os.path.exists(os.path.join(p_dir, 'its2_type_profiles')):
                self.sub_dir = os.path.join(p_dir, 'its2_type_profiles')
                self._p_profile_dir()

            # Delete the html dir
            if os.path.exists(os.path.join(p_dir, 'html')):
                print('Deleting html dir')
                shutil.rmtree(os.path.join(p_dir, 'html'))

            # Delete the non_sym_and_size_violation_sequences dir
            if os.path.exists(os.path.join(p_dir, 'non_sym_and_size_violation_sequences')):
                print('Deleting non_sym_and_size_violation_sequences dir')
                shutil.rmtree(os.path.join(p_dir, 'non_sym_and_size_violation_sequences'))

            # Delete the study_data.js file
            if os.path.exists(os.path.join(p_dir, 'study_data.js')):
                print('Deleting study_data.js')
                os.remove(os.path.join(p_dir, 'study_data.js'))
    
    def _p_profile_dir(self):
        """
        The count table for the profiles is already in a suitable format.
        As before we will keep only the absolute count table
        We will also keep the additional info .txt file
        """
        # Read in the absolute abund and meta file
        # write back out as csv
        df_path = False
        for file_ in os.listdir(self.sub_dir):
            if 'profiles.absolute.abund_and_meta.txt' in file_:
                df_path = os.path.join(self.sub_dir, file_)
        
        if df_path:
            df = pd.read_table(df_path)
            df.to_csv(df_path.replace('.txt', '.csv'), sep=',', index=False)

        # Now delete all other files except for the addtional info file
        for file_ in [_ for _ in os.listdir(self.sub_dir) if not any(substring in _ for substring in ['additional_info.txt', 'profiles.absolute.abund_and_meta.csv'])]:
            print(f'Deleting {file_}')
            os.remove(os.path.join(self.sub_dir, file_))


    def _p_pre_med_dir(self):
        # rename the headers
        # sample_uid --> symportal_datasetsample_uid
        # sample_name --> sample-id
        # Reading int the df takes too long so do using sed.
        new_out_path = os.path.join(self.sub_dir, 'new_out.csv')
        new_out = open(new_out_path, 'w')
        pre_med_abund_path = os.path.join(self.sub_dir, 'pre_med_absolute_abundance_df.csv')
        subprocess.run(['sed', '-e', 's/sample_uid/symportal_datasetsample_uid/g', '-e', 's/sample_name/sample-id/g', pre_med_abund_path], stdout=new_out)
        new_out.close()
        os.rename(new_out_path, pre_med_abund_path)

        # Here just keep the absolute count table and the pre_med_fasta
        for file_ in os.listdir(self.sub_dir):
            if file_ not in ['pre_med_absolute_abundance_df.csv', 'pre_med_master_seqs.fasta']:
                print(f'Deleting {file_}')
                if os.path.isdir(os.path.join(self.sub_dir, file_)):
                    shutil.rmtree(os.path.join(self.sub_dir, file_))
                else:
                    os.remove(os.path.join(self.sub_dir, file_))
        
    def _p_post_med_dir(self):
        # We will aim for only a single count table
        # This count table will be based on the ....seqs.absolute.abund_and_meta.txt
        # table, but it will have some of the columns removed.
        # We will remove all of the unnamed sequence total columns, we will remove all of the host
        # depth and location columns too as this should come from the TARA data
        # We will delete the meta only and abund only tables and all of the
        # relative abund tables as these can easily be calculated from the one count table
        # We will also keep the addtional info file
        # as this tells us about the dataset (i.e. which ID etc.)
        # And we will also keep the ...seqs.fasta so that we know the exact sequences
        # referred to in the count table
        df_path = False
        for file_ in os.listdir(self.sub_dir):
            if '.seqs.absolute.abund_and_meta.txt' in file_:
                df_path = os.path.join(self.sub_dir, file_)
        
        if df_path:
            df = pd.read_table(df_path)
            cols_to_drop = [_ for _ in list(df) if any(substring in _ for substring in ['host', 'collection', 'noName', 'sample_type'])]
            df.drop(columns=cols_to_drop, inplace=True)

            # rename the sample_uid col to symportal_datasetsample_uid
            # rename the sample_name to sample-id
            temp_dict = {'sample_uid':'symportal_datasetsample_uid', 'sample_name':'sample-id'}
            df.rename(columns=temp_dict, inplace=True)

            # now write out the df to its original path only with a csv extension instead of .txt
            df.to_csv(df_path.replace('.txt', '.csv'), sep=',', index=False)

        # then delete all unwanted files in the directory
        for file_ in [_ for _ in os.listdir(self.sub_dir) if not _.endswith(('.csv', '.fasta', '.additional_info.txt'))]:
            print(f'Deleting {file_}')
            os.remove(os.path.join(self.sub_dir, file_))

    def _p_between_dir(self, dir_):
            # There are 11 files we want to keep
            ## DISTANCE FILES ##
            # The 4 .dist files that are the distance matrices that have been computed
            # with Bray-Curtis and Unifrac. Each is provided with source abundance
            # sqrt transformed and not.

            ## PCOA ##
            # The pcoa coordinates that have been computed respective to each of the distance files
            
            ## FASTA FILES ##
            # Aligned and unaligned .fasta files of the sequences used in the computation
            # of the trees 

            ## TREE AND MODEL OUTPUT ##
            # The tree and the model used to compute the tree in iqtree

            # Delete all other files
            for file_ in [_ for _ in os.listdir(dir_) if not _.endswith(('.dist', '.csv', '.fasta', '.iqtree'))]:
                print(f"deleting {file_}")
                os.remove(os.path.join(dir_, file_))

            # TODO the .dist files do not need renaming but the .csv headers do
            if 'profile' in dir_:
                # For the between profile distances
                # sample --> profile
                for file_ in [_ for _ in os.listdir(dir_) if _.endswith('.csv')]:
                    df = pd.read_csv(os.path.join(dir_, file_))
                    # enforce int for uid
                    df['analysis_type_uid'] = df['analysis_type_uid'].astype(int)
                    temp_dict = {'sample':'profile'}
                    df.rename(columns=temp_dict, inplace=True)
                    
                    df.to_csv(os.path.join(dir_, file_), index=False)
            else:
                # For the between sample distances
                # sample --> sample-id
                # sample_uid --> symportal_datasetsample_uid
                for file_ in [_ for _ in os.listdir(dir_) if _.endswith('.csv')]:
                    df = pd.read_csv(os.path.join(dir_, file_))
                    if 'sample' in list(df) or 'sample_uid' in list(df):
                        temp_dict = {'sample':'sample-id', 'sample_uid':'symportal_datasetsample_uid'}
                        df['sample_uid'] = df['sample_uid'].astype(int)
                        df.rename(columns=temp_dict, inplace=True)
                    df.to_csv(os.path.join(dir_, file_), index=False)

ProcessSPOutputForZenodo().process()