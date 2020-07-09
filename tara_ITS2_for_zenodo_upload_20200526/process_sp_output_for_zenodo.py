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
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import re

class ProcessSPOutputForZenodo:
    def __init__(self):
        self.base_dir = os.path.dirname(os.path.realpath(__file__))
        self.coral_output_dir_path = os.path.join(self.base_dir, "coral")
        self.non_coral_output_dir_path = os.path.join(self.base_dir, "non_coral")
        self.negative_output_dir_path = os.path.join(self.base_dir, "negative")
        # The file that Stephane produced that maps negative readset names
        # to TARA sample-id names
        self.negative_mapping_file_path = os.path.join(os.path.dirname(self.base_dir), 'TARA-PACIFIC_list-of-negatives-ITS2_for-Caroline.csv')
        self.neg_mapping_dict = self._remove_tech_rep_negatives()
        # A dynamic directory path
        self.sub_dir = None
        # This will be need to be updated if we create new versions
        self.version_string = 'v1'

    def _remove_tech_rep_negatives(self):
        """
        The samples that we have processed for the negative controls contains technical
        replicates. We want to get rid of these for the zenodo uploads.
        To do this we will use Stephanes table to identify those readsets that share
        a TARA id (these will be tech reps) and we will only keep the largest sample
        (i.e highest number of contigs). We can get the contig information from the
        abund_meta file in the post_med directory.
        We can then process the post_med and pre_med count files here.
        """

        # For each row in the post_med_df, find the mapping key that is a substring
        # Should be only one, check this.
        # Then once you have found the one, check all samples in the post_med df to see if it matches any other
        # if you return multiple matches, then keep only the one with the biggest number of contigs,
        # and all others to a drop list. Keep a checked list so that we don't have to check readsets twice.
        # Also up date a dictionary as you go that is the full readset to the sample-id that it needs to become.
        # Once this has been done for the post-med do it for the pre-med.
        # For the pre-med, use the dictionary we created while doing the post-med

        # Get the post med df. Read it in with index as false and set index manually without dropping
        # this way we can work with the index, but then we can not write it out later so as not
        # to disturb the column orders.
        post_med_count_path = os.path.join(self.negative_output_dir_path, 'post_med_seqs', [_ for _ in os.listdir(
            os.path.join(self.negative_output_dir_path, 'post_med_seqs')) if 'abund' in _][0])
        post_med_df = pd.read_csv(post_med_count_path, index_col=False)
        post_med_df = post_med_df.set_index('sample-id', drop=False)

        # Same for the pre_med
        pre_med_count_path = os.path.join(self.negative_output_dir_path, 'pre_med_seqs', [_ for _ in os.listdir(
            os.path.join(self.negative_output_dir_path, 'pre_med_seqs')) if 'abund' in _][0])
        pre_med_df = pd.read_csv(pre_med_count_path, index_col=False)
        pre_med_df = pre_med_df.set_index('sample-id', drop=False)

        # First check to see if the sample-ids have already been fixed
        if 'TARA' in pre_med_df.index[0] and 'TARA' in post_med_df.index[0]:
            return
        if 'TARA' in pre_med_df.index[0] and 'TARA' not in post_med_df.index[0]:
            raise RuntimeError
        if 'TARA' not in pre_med_df.index[0] and 'TARA' in post_med_df.index[0]:
            raise RuntimeError

        # The dictionary df that Setphane produced
        mapping_df = pd.read_csv(self.negative_mapping_file_path, index_col=0)
        # Make the mapping dictionary from the Stephane df
        raw_mapping_dict = {}
        for df_ind in mapping_df.index:
            raw_mapping_dict[df_ind] = mapping_df.at[df_ind, 'sample-id_source']

        # This is the dictionary we are going to populate that had the full genoscope readset
        # as the key and the equivalent TARA sample-id as the value
        curated_mapping_dict = {}

        # Check that the assumption holds that both of the indeces are identifcal except for order.
        # NB the post med df has an annoying row at the end.
        assert(set(post_med_df.index[:-1]) == set(pre_med_df.index))
        contig_dict = {readset: contig for readset, contig in zip(post_med_df['sample-id'][:-1], post_med_df['raw_contigs'][:-1])}

        to_drop_list = []
        checked_list = []
        for pm_ind in post_med_df.index[:-1]:
            if pm_ind in checked_list:
                continue
            match = []
            for map_ind in mapping_df.index:
                if map_ind in pm_ind:
                    match.append(map_ind)
            if len(match) == 0:
                print(f'pm_ind: {pm_ind} found 0 matches. This sample will be dropped.')
                to_drop_list.append(pm_ind)
                continue
            elif len(match) > 1:
                raise RuntimeError

            # Now we have the mapping indice that matches
            match = match[0]
            pm_matches = []
            for pm_ind_again in post_med_df.index[:-1]:
                if match in pm_ind_again:
                    pm_matches.append(pm_ind_again)
            assert(len(pm_matches) > 0)
            if len(pm_matches) > 1:
                # Then we have technical replicates and we only want to keep the largest
                contig_match_dict = {pm_match: contig_dict[pm_match] for pm_match in pm_matches}
                sorted_keys = sorted(contig_match_dict, key=contig_match_dict.get, reverse=True)
                # Add all of the matches to the check_list
                checked_list.extend(sorted_keys)
                curated_mapping_dict[sorted_keys[0]] = raw_mapping_dict[match]
                to_drop_list.extend(sorted_keys[1:])
            else:
                checked_list.append(pm_matches[0])
                curated_mapping_dict[pm_matches[0]] = raw_mapping_dict[match]

        # drop the rows
        post_med_df.drop(index=to_drop_list, inplace=True)
        # We now need to get rid of any sequence count columns that only have 0s after dropping the samples
        # The last meta column is post_med_unique
        cols = list(post_med_df)
        c_ind = cols.index('post_med_unique') + 1
        cols_to_check = cols[c_ind:]
        cols_to_drop = []
        for col in cols_to_check:
            if (post_med_df[col][:-1] == 0).all():
                cols_to_drop.append(col)

        # drop the cols
        post_med_df.drop(columns=cols_to_drop, inplace=True)

        # rename
        for ind in post_med_df.index[:-1]:
            current = post_med_df.at[ind, 'sample-id']
            post_med_df.at[ind, 'sample-id'] = curated_mapping_dict[current]

        # Here we have the curated mapping dict popualted and we can now use this to
        # process the pre_med df
        pre_med_df.drop(index=to_drop_list, inplace=True)
        # We now need to get rid of any sequence count columns that only have 0s after dropping the samples
        # The last meta column is post_med_unique
        cols = list(pre_med_df)
        c_ind = cols.index('sample-id') + 1
        cols_to_check = cols[c_ind:]
        cols_to_drop = []
        for col in cols_to_check:
            if (pre_med_df[col][:-1] == 0).all():
                cols_to_drop.append(col)

        # drop the cols
        pre_med_df.drop(columns=cols_to_drop, inplace=True)

        # rename
        for ind in pre_med_df.index:
            current = pre_med_df.at[ind, 'sample-id']
            pre_med_df.at[ind, 'sample-id'] = curated_mapping_dict[current]

        # Now convert the columns to int32
        d_type_dict = {col_name : pd.Int32Dtype() for col_name in list(post_med_df)[2:]}
        post_med_df = post_med_df.astype(d_type_dict)
        d_type_dict = {col_name : pd.Int32Dtype() for col_name in list(pre_med_df)[2:]}
        pre_med_df = pre_med_df.astype(d_type_dict)

        # Important to write out with index as false
        post_med_df.to_csv(post_med_count_path, index=False, header=True)
        pre_med_df.to_csv(pre_med_count_path, index=False, header=True)

    def process(self):
        for p_dir in [self.coral_output_dir_path, self.non_coral_output_dir_path, self.negative_output_dir_path]:
            # Process the betweeen_sample_distances
            # If we are working with the negatives, delete the between sample directory
            self.sub_dir = os.path.join(p_dir, 'between_sample_distances')
            if p_dir == self.negative_output_dir_path:
                if os.path.exists(self.sub_dir):
                    shutil.rmtree(self.sub_dir)
            else:
                # This directory will exist in both the coral and non_coral
                for clade_dir in [_ for _ in os.listdir(self.sub_dir) if _ in list('ABCDEFGHI')]: # For each cladal directory
                    self._p_between_dir(os.path.join(self.sub_dir, clade_dir))
            
            # Process the between_profile_distances
            # Only for coral samples
            self.sub_dir = os.path.join(p_dir, 'between_profile_distances')
            if os.path.exists(self.sub_dir):
                for clade_dir in [_ for _ in os.listdir(self.sub_dir) if _ in list('ABCDEFGHI')]: # For each cladal directory
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
            if not os.path.basename(self.base_dir) in file_:
                print(f'Deleting {file_}')
                os.remove(os.path.join(self.sub_dir, file_))

        # Now rename the files
        # Create a rename dict
        # Key is old full path, value is new full path
        rename_dict = {}
        parent_directories = self._get_parent_dir_list(self.sub_dir)

        # Rename the count table
        count_file = [_ for _ in os.listdir(self.sub_dir) if 'abund_and_meta' in _]
        assert (len(count_file) == 1)
        count_file_old = count_file[0]
        if os.path.basename(self.base_dir) not in count_file_old:
            match = re.search("profiles", count_file_old)
            if match is None:
                raise RuntimeError
            # we need to discard from the match onwards
            new_count_file = count_file_old[match.span()[0]:]
            ext = new_count_file.split('.')[-1]
            new_count_file = new_count_file.replace('.', '_')
            new_count_file = new_count_file.replace('_profiles', '')
            # insert the version number
            new_count_file = new_count_file.replace(f'_{ext}', f'_{self.version_string}.{ext}')
            # add the directory structure info
            new_count_file = '_'.join(parent_directories) + f'_{new_count_file}'
            # add to rename dict
            rename_dict[os.path.join(self.sub_dir, count_file_old)] = os.path.join(self.sub_dir, new_count_file)

        # Rename the additional output file
        info_file = [_ for _ in os.listdir(self.sub_dir) if 'additional' in _]
        assert (len(info_file) == 1)
        info_file_old = info_file[0]
        if os.path.basename(self.base_dir) not in info_file_old:
            match = re.search("additional", info_file_old)
            if match is None:
                raise RuntimeError
            # we need to discard from the match onwards
            new_info_file = info_file_old[match.span()[0]:]
            # add the directory structure info
            new_info_file = '_'.join(parent_directories) + '_symportal_output_' + new_info_file
            # insert the version number
            new_info_file = new_info_file.replace('.txt', f'_{self.version_string}.txt')
            # add to rename dict
            rename_dict[os.path.join(self.sub_dir, info_file_old)] = os.path.join(self.sub_dir, new_info_file)

        # At this point we can do the renaming
        self._rename_from_dict(rename_dict)

    def _rename_from_dict(self, r_dict):
        for k, v in r_dict.items():
            os.rename(k, v)

    def _p_pre_med_dir(self):
        # rename the headers
        # sample_uid --> symportal_datasetsample_uid
        # sample_name --> sample-id
        # Reading int the df takes too long so do using sed.
        pre_med_abund_path = os.path.join(self.sub_dir, 'pre_med_absolute_abundance_df.csv')
        if os.path.exists(pre_med_abund_path):
            new_out_path = os.path.join(self.sub_dir, 'new_out.csv')
            new_out = open(new_out_path, 'w')
            pre_med_abund_path = os.path.join(self.sub_dir, 'pre_med_absolute_abundance_df.csv')
            subprocess.run(['sed', '-e', 's/sample_uid/symportal_datasetsample_uid/g', '-e', 's/sample_name/sample-id/g', pre_med_abund_path], stdout=new_out)
            new_out.close()
            os.rename(new_out_path, pre_med_abund_path)

        # Here just keep the absolute count table and the pre_med_fasta
        for file_ in [_ for _ in os.listdir(self.sub_dir) if os.path.basename(self.base_dir) not in _]:
            if file_ not in ['pre_med_absolute_abundance_df.csv', 'pre_med_master_seqs.fasta']:
                print(f'Deleting {file_}')
                if os.path.isdir(os.path.join(self.sub_dir, file_)):
                    shutil.rmtree(os.path.join(self.sub_dir, file_))
                else:
                    os.remove(os.path.join(self.sub_dir, file_))

        # Now rename the files
        # Create a rename dict
        # Key is old full path, value is new full path
        rename_dict = {}
        parent_directories = self._get_parent_dir_list(self.sub_dir)

        # Rename the pre_med count
        count_file = [_ for _ in os.listdir(self.sub_dir) if 'abund' in _]
        assert (len(count_file) == 1)
        count_file_old = count_file[0]
        if os.path.basename(self.base_dir) not in count_file_old:
            new_count_file = f'sequences_absolute_abund_{self.version_string}.csv'
            # add the directory structure info
            new_count_file = '_'.join(parent_directories) + f'_{new_count_file}'
            # add to rename dict
            rename_dict[os.path.join(self.sub_dir, count_file_old)] = os.path.join(self.sub_dir, new_count_file)

        # Rename the fasta
        fasta_file = [_ for _ in os.listdir(self.sub_dir) if _.endswith('.fasta')]
        assert (len(fasta_file) == 1)
        fasta_file_old = fasta_file[0]
        if os.path.basename(self.base_dir) not in fasta_file_old:
            new_fasta_file = f'sequences_{self.version_string}.fasta'
            # add the directory structure info
            new_fasta_file = '_'.join(parent_directories) + f'_{new_fasta_file}'
            # add to rename dict
            rename_dict[os.path.join(self.sub_dir, fasta_file_old)] = os.path.join(self.sub_dir, new_fasta_file)

        # At this point we can do the renaming
        self._rename_from_dict(rename_dict)

    def _p_post_med_dir(self, rename_negatives=False):
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
        # If rename_negatives is true, then we are working with the negatives
        # and we need to check that the technical replicates have been removed
        # and the that the genoscope readset IDs have been replaced with sample-ids

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
        for file_ in [_ for _ in os.listdir(self.sub_dir) if not _.endswith(('.csv', '.fasta')) and 'additional_info' not in _]:
            print(f'Deleting {file_}')
            os.remove(os.path.join(self.sub_dir, file_))

        # Now rename the files
        # Create a rename dict
        # Key is old full path, value is new full path
        rename_dict = {}
        parent_directories = self._get_parent_dir_list(self.sub_dir)

        # Rename the additional output file
        info_file = [_ for _ in os.listdir(self.sub_dir) if 'additional' in _]
        assert (len(info_file) == 1)
        info_file_old = info_file[0]
        if os.path.basename(self.base_dir) not in info_file_old:
            match = re.search("additional", info_file_old)
            if match is None:
                raise RuntimeError
            # we need to discard from the match onwards
            new_info_file = info_file_old[match.span()[0]:]
            # add the directory structure info
            new_info_file = '_'.join(parent_directories) + '_symportal_output_' + new_info_file
            # insert the version number
            new_info_file = new_info_file.replace('.txt', f'_{self.version_string}.txt')
            # add to rename dict
            rename_dict[os.path.join(self.sub_dir, info_file_old)] = os.path.join(self.sub_dir, new_info_file)

        # Rename the count table
        count_file = [_ for _ in os.listdir(self.sub_dir) if 'abund_and_meta' in _]
        assert (len(count_file) == 1)
        count_file_old = count_file[0]
        if os.path.basename(self.base_dir) not in count_file_old:
            match = re.search("seqs", count_file_old)
            if match is None:
                raise RuntimeError
            # we need to discard from the match onwards
            new_count_file = count_file_old[match.span()[0]:]
            new_count_file = new_count_file.replace('seqs', 'sequences')
            ext = new_count_file.split('.')[-1]
            new_count_file = new_count_file.replace('.', '_')
            # insert the version number
            new_count_file = new_count_file.replace(f'_{ext}', f'_{self.version_string}.{ext}')
            # add the directory structure info
            new_count_file = '_'.join(parent_directories) + f'_{new_count_file}'
            # add to rename dict
            rename_dict[os.path.join(self.sub_dir, count_file_old)] = os.path.join(self.sub_dir, new_count_file)

        # Rename the fasta
        fasta_file = [_ for _ in os.listdir(self.sub_dir) if _.endswith('.fasta')]
        assert (len(fasta_file) == 1)
        fasta_file_old = fasta_file[0]
        if os.path.basename(self.base_dir) not in fasta_file_old:
            match = re.search("seqs", fasta_file_old)
            if match is None:
                raise RuntimeError
            # we need to discard from the match onwards
            new_fasta_file = fasta_file_old[match.span()[0]:]
            new_fasta_file = new_fasta_file.replace('seqs', 'sequences')
            ext = new_fasta_file.split('.')[-1]
            new_fasta_file = new_fasta_file.replace('.', '_')
            # insert the version number
            new_fasta_file = new_fasta_file.replace(f'_{ext}', f'_{self.version_string}.{ext}')
            # add the directory structure info
            new_fasta_file = '_'.join(parent_directories) + f'_{new_fasta_file}'
            # add to rename dict
            rename_dict[os.path.join(self.sub_dir, fasta_file_old)] = os.path.join(self.sub_dir, new_fasta_file)

        # At this point we can do the renaming
        self._rename_from_dict(rename_dict)

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
            # The tree and the model used to compute the tree in iqtree format

            ## RENAMING ##
            # We will also do some renaming.
            # The name will incorportate the directory structure in order
            # It will incorporate the version identifier.

            # Delete all other files
            for file_ in [_ for _ in os.listdir(dir_) if not _.endswith(('.dist', '.csv', '.fasta', '.iqtree'))]:
                print(f"deleting {file_}")
                os.remove(os.path.join(dir_, file_))

            # The braycurtis files have an additional line at the top that should be got rid of
            # This line is a single int that is the number of samples in the distance file
            for file_ in [_ for _ in os.listdir(dir_) if _.endswith('.dist') and 'braycurtis' in _]:
                # Remove the top line
                with open(os.path.join(dir_, file_), 'r') as f:
                    line_list = [_.rstrip() for _ in f]
                try:
                    # This should raise a value error if the line in question doesn't exist
                    foo = int(line_list[0])
                    # If no error raised, then remove the first line and write out
                    new_list = line_list[1:]
                    with open(os.path.join(dir_, file_), 'w') as f:
                        for _ in new_list:
                            f.write(f'{_}\n')
                except ValueError:
                    # Then the line has already been deleted or no long exists
                    pass

            # The .dist files do not need renaming but the .csv headers do
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

            # Now we should rename the files
            # First check to see if the basename dir is incorporated into filename
            # as a check to see whether the renaming has already been completed
            # We want to remove the time date component of the file name
            # We want to remove the clade component of the file name.
            # We want to replace the period (if it exists before the 'unifrac' and 'braycurtis')
            # We want to replace the periods with underscores

            # Create a rename dict
            # Key is old full path, value is new full path
            rename_dict = {}
            parent_directories = self._get_parent_dir_list(dir_)

            # Rename the .dist and .csv files
            for file_ in [_ for _ in os.listdir(dir_) if _.endswith('.dist') or _.endswith('.csv')]:
                if os.path.basename(self.base_dir) in file_:
                    # File has already been renamed
                    continue
                # First process the file_ name
                match = re.search("braycurtis|unifrac", file_)
                if match is None:
                    raise RuntimeError
                # we need to discard from the match onwards
                new_file = file_[match.span()[0]:]

                # For the pcoa we can replace the 'sample_PCoA_' with just 'PCoA_'
                new_file = new_file.replace('sample_PCoA_', 'PCoA_')
                new_file = new_file.replace('samples_PCoA_', 'PCoA_')
                new_file = new_file.replace('_btwn', '')
                new_file = new_file.replace('sample_distances_', 'distances_')
                new_file = new_file.replace('_within_clade_profile', '')
                new_file = new_file.replace('_profiles', '')
                new_file = new_file.replace('_profile', '')

                # now remove the clade designation
                # the clade designation is the last parent directory
                new_file = new_file.replace(f'_{parent_directories[-1]}', '')

                # Insert the version number
                ext = new_file.split('.')[-1]
                new_file = new_file.replace(f'.{ext}', f'_{self.version_string}.{ext}')

                # finally, concatenate the parent directories to generate the name
                # Then add to rename_dict
                new_file = '_'.join(parent_directories) + f'_{new_file}'
                rename_dict[os.path.join(dir_, file_)] = os.path.join(dir_, new_file)

            # Rename the .fasta and .iqtree files
            for file_ in [_ for _ in os.listdir(dir_) if _.endswith('.fasta') or _.endswith('.iqtree')]:
                if os.path.basename(self.base_dir) in file_:
                    # File has already been renamed
                    continue
                match = re.search("seqs", file_)
                if match is None:
                    raise RuntimeError
                # we need to discard from the match onwards
                new_file = file_[match.span()[0]:]
                new_file = new_file.replace('seqs', 'sequences')
                new_file = new_file.replace('.fasta.iqtree', '.iqtree')

                # Remove periods and insert the version number
                ext = new_file.split('.')[-1]
                new_file = new_file.replace('.', '_')
                new_file = new_file.replace(f'_{ext}', f'_{self.version_string}.{ext}')

                # finally, concatenate the parent directories to generate the name
                # Then add to rename_dict
                new_file = '_'.join(parent_directories) + f'_{new_file}'
                rename_dict[os.path.join(dir_, file_)] = os.path.join(dir_, new_file)

            # At this point we can do the renaming
            # TODO we are here.
            self._rename_from_dict(rename_dict)

    def _get_parent_dir_list(self, dir_):
        # Get a list of the parent directories to make the new file name from
        parent_directories = []
        parent_path = os.path.dirname(os.path.join(dir_, os.listdir(dir_)[0]))
        parent_dir = os.path.basename(parent_path)
        parent_directories.append(parent_dir)
        while parent_path != self.base_dir:
            parent_path = os.path.dirname(parent_path)
            parent_dir = os.path.basename(parent_path)
            parent_directories.append(parent_dir)
        parent_directories.reverse()
        return parent_directories


ProcessSPOutputForZenodo().process()

class AssessSeqDepth:
    """
    Quick utility class to investigate whether we need to be using any of the tech replicates to increase
    sequencing depth.
    Load in the post_med dfs for coral and non coral and look to see that all samples are above 10000 contigs.
    Turns out that the average is >60000 for both. And the smallest raw_contigs are about 15 000. So no need
    for tech reps.
    """
    def __init__(self):
        self.path_to_post_med_seq_tab_coral = "/home/humebc/projects/tara/tara_full_dataset_processing/tara_ITS2_for_zenodo_upload_20200526/TARA_PACIFIC_METAB_ITS2_CORAL/post_med_seqs/107_20200527_2020-05-27_22-53-40.050470.seqs.absolute.abund_and_meta.csv"
        self.path_to_post_med_seq_tab_non_coral = "/home/humebc/projects/tara/tara_full_dataset_processing/tara_ITS2_for_zenodo_upload_20200526/TARA_PACIFIC_METAB_ITS2_NON_CORAL/post_med_seqs/2020-05-27_22-56-23.172444.seqs.absolute.abund_and_meta.csv"
        self.path_to_its2_replication_tab = "/home/humebc/projects/tara/tara_full_dataset_processing/output/output_information_df_all_fastqs_its2_2020-05-31T11_06_03.996934UTC.csv"
        # We will use the two above tables to append a column to the replciation
        # table that is the num of contigs. Then we will work with this one table.
        self.readset_df = self._make_readset_df()

    def _make_readset_df(self):
        post_med_coral_df = pd.read_csv(self.path_to_post_med_seq_tab_coral, index_col=1)
        post_med_non_coral_df = pd.read_csv(self.path_to_post_med_seq_tab_non_coral, index_col=1)
        
        # look for coral
        for ind, val in post_med_coral_df['raw_contigs'].items():
            if str(ind) != 'nan':
                if int(val) < 10000:
                    print(ind, val)
        print(f"The average for coral contigs was: {np.average(post_med_coral_df['raw_contigs'][:-1])}")

        # look for non-coral
        for ind, val in post_med_non_coral_df['raw_contigs'].items():
            if str(ind) != 'nan':
                if int(val) < 10000:
                    print(ind, val)
        print(f"The average for coral contigs was: {np.average(post_med_non_coral_df['raw_contigs'][:-1])}")


        
# AssessSeqDepth()._make_readset_df()



