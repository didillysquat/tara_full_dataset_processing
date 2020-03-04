from base_18s import EighteenSBase
import os
import compress_pickle
import pandas as pd
import hashlib
import sys
from plumbum import local
import subprocess
from scipy.spatial.distance import braycurtis
from skbio.stats.ordination import pcoa
from skbio.tree import TreeNode
from skbio.diversity import beta_diversity
import numpy as np
import itertools
from functools import partial

class EighteenSDistance(EighteenSBase):
    """
    Class concerned with the calculation and presentation of the between sample
    distances calculated based on the low level sequence diversity found within each of the samples.
    We will work with UniFrac and BrayCurtis versions of all of the distances we calculate.
    Eventually we will want to be comparing these distances to the clusters that Didier is going
    to produce. However, for the time being we will just look to see if we can find any natural clustering patterns.
    I think we will likely have to work in a somewhat nested manner of fashion.
    Fistly we should split by host species. From this we should verify that the samples that have
    the secondary sequences (i.e. different most abundant sequences to the 'normal' most abundant sequence for the
    genus) can be differentiated in the PCoA. It should be quite obvious.
    If this is the case then the second set of distances that we should compute should be seperated according to 
    the most abundant sequences (i.e. seperate out those samples that have the secondary sequences as their most abunant sequence).
    Once we have done this we should look at the plots and see how the structuring is looking.
    """
    def __init__(self):
        super().__init__()
        # We will need to which coral species the consolidated coral dict belongs
        # to for each of the samples. This differs from the TARA logsheet annotation
        # sometimes. We will add this information to the info_df and then pickle out
        # a new info_df so that we only have to do this once. We will also need to know
        # what the most abundant sequence is in the consolidated host dictionary for 
        # each sample so that we can later differentiate according to this and also
        # check for separation accorrding to this on the first PCoA that we will produce
        self._add_additional_info_to_info_df()
        
    def _add_additional_info_to_info_df(self):
        if os.path.isfile(os.path.join(self.cache_dir, 'info_df_with_additional_info.p.bz')):
            self.info_df = compress_pickle.load(os.path.join(self.cache_dir, 'info_df_with_additional_info.p.bz'))
        else:
            print('updating info_df')
            most_abund_coral_genus_df_list = []
            most_abund_seq_of_coral_genus_df_list = []
            for sample_name in self.info_df.index:
                sys.stdout.write(f'\r{sample_name}')
                sample_qc_dir = os.path.join(self.qc_dir, sample_name)
                rel_all_seq_abundance_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'rel_all_seq_abundance_dict.p.bz'))
                coral_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'coral_annotation_dict.p.bz'))
                most_abund_coral_genus_df_list.append(self._identify_most_abund_coral_genus(rel_all_seq_abundance_dict, coral_annotation_dict))
                consolidated_host_seqs_abund_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
                most_abund_seq_of_coral_genus_df_list.append(sorted([_ for _ in consolidated_host_seqs_abund_dict.items()], key=lambda x: x[0], reverse=True)[0][0])
            self.info_df['most_abund_coral_genus'] = most_abund_coral_genus_df_list
            self.info_df['most_abund_seq_of_coral_genus'] = most_abund_seq_of_coral_genus_df_list
            compress_pickle.dump(self.info_df, os.path.join(self.cache_dir, 'info_df_with_additional_info.p.bz'))
            print()

    def _identify_most_abund_coral_genus(self, rel_all_seq_abundance_dict, coral_annotation_dict):
        for sorted_tup in sorted(
            [(seq_name, rel_abund) for seq_name, rel_abund in rel_all_seq_abundance_dict.items()], 
            key=lambda x: x[1], 
            reverse=True
            ):
                try:
                    genus = coral_annotation_dict[sorted_tup[0]]
                    if genus == 'Porites':
                        return 'Porites'
                    elif genus == 'Pocillopora':
                        return 'Pocillopora'
                    elif genus == 'Millepora':
                        return 'Millepora'
                except KeyError:
                    continue
    
    def make_and_plot_dist_and_pcoa(self, resolution_type):
        # Call the class that will be responsible for calculating the distances
        # Calculating the PCoA and then plotting
        da = DistanceAnlyses(resolution_type=resolution_type, info_df=self.info_df, qc_dir=self.qc_dir, output_dir=self.output_dir)
        da.compute_distances()
        # TODO we can do the plotting here

class DistanceAnlyses:
    def __init__(self, resolution_type, info_df, qc_dir, output_dir):
        self.resolution_type = resolution_type
        self.qc_dir = qc_dir
        self.info_df = info_df
        self.output_dir = output_dir
        # get the category names and list of samples that will be used in calculating
        # the distances for each of those categories. For example when
        # resolution_type == 'host_only', categories will be Porites, Millepora, Pocillopora
        # When host_secondary_seq == e.g. Porites_first, Porites_second etc.
        self.categories, self.sample_lists = self._generate_distance_categories()
        self.dist_methods = ['unifrac', 'braycurtis']
    
    def compute_distances(self):
        for dist_method in self.dist_methods:
            for distance_cat, sample_list in zip(self.categories, self.sample_lists):
                indi_dist = IndiDistanceAnalysis(
                    distance_cat=distance_cat, sample_list=sample_list, dist_method=dist_method, qc_dir=self.qc_dir, output_dir=self.output_dir)
                indi_dist.do_analysis()

    def _generate_distance_categories(self):
        if self.resolution_type == 'host_only':
            sample_lists = []
            categories = ['Pocillopora', 'Millepora', 'Porites']
            for category in categories:
                sample_lists.append(list(self.info_df[self.info_df['most_abund_coral_genus'] == category].index))
            return categories, sample_lists

class IndiDistanceAnalysis:
    def __init__(self, distance_cat, sample_list, dist_method, qc_dir, output_dir):
        """This is the class that will take care of creating a single between sample
        distance method
        
        TODO move the tree files to the cache file rather than the temp file so that we can delete the temp file
        Also use compression when writing out and then reading in the pcoa output files and the distance files
        """
        self.qc_dir = qc_dir
        self.category = distance_cat
        self.samples = sample_list
        self.dist_method = dist_method
        self.output_dir = output_dir
        
        # A temp directory where we can write out the unaligned fastas
        # and aligned fasta and any other intermediate files
        self.temp_dir = os.path.join(os.path.dirname(self.qc_dir), 'temp')
        os.makedirs(self.temp_dir, exist_ok=True)
        # In the initial 18s we were using a normalisation of 10000 seqs for the braycurtis
        # but 1000 for the unifrac due to the amount of time it was taking to create
        # the trees. We will start here with 1000 seqs for the unifrac and see
        # how long it takes to make the tree. If its not too bad then we can up the number
        if self.dist_method == 'unifrac':
            self.num_seqs_to_normalise_to = 1000
        else:
            self.num_seqs_to_normalise_to = 10000
        self.abundance_df = self._create_abundance_df()

        # Generic variables shared between braycurtis and unifrac
        self.dist_out_path = os.path.join(self.output_dir, f'{self.category}_{self.dist_method}.dist')
        self.pcoa_out_path = os.path.join(self.output_dir, f'{self.category}_{self.dist_method}.csv')
        self.pcoa_df = None

        # Variables concerned with unifrac
        self.unaligned_fasta_path = os.path.join(self.temp_dir, 'unaligned_fasta.fasta')
        self.aligned_fasta_path = os.path.join(self.temp_dir, 'aligned_fasta.fasta')
        self.tree_path = self.aligned_fasta_path + '.treefile'
        self.wu = None

        # Variables concerned with braycurtis
        self.braycurtis_btwn_sample_distance_dictionary = None
        
        # The pseudo code for this is
        # Create an abundance dataframe (we need to choose whether to normalise this or not)
        # We probably should normalise.
        # Do this by creating a dict of dicts and using the dicts methods that exist
        # in the initial 18S code
        # master fasta file of the sequences from the samples,
        # If we are doing unifrac then we will need to compute a tree
        # then we can do either the braycurtis calculations or
        # the unifrac calculations.
        # then from this produce pcoa
        # then plot this
        # The tree creation will likely be one of the most expensive parts
        # so we should aim to pickle out the trees if possible, same with alignments.
        # We can possibly pickle these items out as their hashes.

    def _create_abundance_df(self):
        """For each sample in self.samples, load up the consolidated_host_seqs_abund_dict
        and create a dict that is sequence to normalised abundance. Then add this dictionary
        to a master dictionary where sample_name is key. Then finally create a df from this
        dict of dicts"""
        dict_to_create_df_from = {}
        print('Creating abundance df')
        for sample_name in self.samples:
            sys.stdout.write(f'\r{sample_name}')
            temp_sample_dict = {}
            sample_qc_dir = os.path.join(self.qc_dir, sample_name)
            consolidated_host_seqs_abund_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
            
            # We need to remove the most abundant sequence from the equation
            most_abund_sequence = max(consolidated_host_seqs_abund_dict.keys(), key=(lambda key: consolidated_host_seqs_abund_dict[key]))
            # remove the most abund seq
            del consolidated_host_seqs_abund_dict[most_abund_sequence]
            # renormalise
            tot = sum(consolidated_host_seqs_abund_dict.values())
            consolidated_host_seqs_abund_dict = {k: v/tot for k, v in consolidated_host_seqs_abund_dict.items()}

            for sequence, rel_abund in consolidated_host_seqs_abund_dict.items():
                normalised_abund = int(rel_abund*self.num_seqs_to_normalise_to)
                if normalised_abund:
                    temp_sample_dict[sequence] = normalised_abund
            dict_to_create_df_from[sample_name] = temp_sample_dict
        
        df = pd.DataFrame.from_dict(dict_to_create_df_from, orient='index')
        df[pd.isna(df)] = 0
        print('\ndf creation complete\n')
        return df

    def do_analysis(self):
        if self.dist_method == 'unifrac':
            self._do_unifrac_analysis()
        else:
            self._do_braycurtis_analysis()

    # BRAYCURTIS METHODS
    def _do_braycurtis_analysis(self):
        self.braycurtis_btwn_sample_distance_dictionary = self._compute_braycurtis_distance_dict()
        self.braycurtis_btwn_sample_distance_file = self._make_and_write_braycurtis_distance_file()
        self.pcoa_df = self._make_pcoa_df_braycurtis()

    def _make_pcoa_df_braycurtis(self):
        # simultaneously grab the sample names in the order of the distance matrix and put the matrix into
        # a twoD list and then convert to a numpy array
        temp_two_D_list = []
        sample_names_from_dist_matrix = []
        for line in self.braycurtis_btwn_sample_distance_file[1:]:
            temp_elements = line.split('\t')
            sample_names_from_dist_matrix.append(temp_elements[0])
            temp_two_D_list.append([float(a) for a in temp_elements[1:]])
        uni_frac_dist_array = np.array(temp_two_D_list)
        sys.stdout.write('\rcalculating PCoA coordinates')

        pcoa_df = pcoa(uni_frac_dist_array)

        # rename the dataframe index as the sample names
        pcoa_df.samples['sample'] = sample_names_from_dist_matrix
        renamed_dataframe = pcoa_df.samples.set_index('sample')

        # now add the variance explained as a final row to the renamed_dataframe
        renamed_dataframe = renamed_dataframe.append(pcoa_df.proportion_explained.rename('proportion_explained'))

        return renamed_dataframe

    def _make_and_write_braycurtis_distance_file(self):
        # Generate the distance out file from this dictionary
        # from this dict we can produce the distance file that can be passed into the generate_PCoA_coords method
        distance_out_file = [len(self.samples)]
        for sample_outer in self.samples:
            # The list that will hold the line of distance. This line starts with the name of the sample
            temp_sample_dist_string = [sample_outer]

            for sample_inner in self.samples:
                if sample_outer == sample_inner:
                    temp_sample_dist_string.append(0)
                else:
                    temp_sample_dist_string.append(self.braycurtis_btwn_sample_distance_dictionary[
                                                        f'{sample_outer}_{sample_inner}'])
            distance_out_file.append(
                '\t'.join([str(distance_item) for distance_item in temp_sample_dist_string]))
        
        with open(self.dist_out_path, 'w') as f:
            for line in distance_out_file:
                f.write(f'{line}\n')
        return distance_out_file

    def _compute_braycurtis_distance_dict(self):
        # Create a dictionary that will hold the distance between the two samples
        spp_distance_dict = {}
        # For pairwise comparison of each of these sequences
        for smp_one, smp_two in itertools.combinations(self.samples, 2):
            print('Calculating distance for {}_{}'.format(smp_one, smp_two))
            # Get a set of the sequences found in either one of the samples
            smp_one_abund_dict = self.abundance_df[smp_one]
            smp_two_abund_dict = self.abundance_df[smp_two]
            list_of_seqs_of_pair = []
            list_of_seqs_of_pair.extend(list(smp_one_abund_dict.keys()))
            list_of_seqs_of_pair.extend(list(smp_two_abund_dict.keys()))
            list_of_seqs_of_pair = list(set(list_of_seqs_of_pair))

            # then create a list of abundances for sample one by going through the above list and checking
            sample_one_abundance_list = []
            for seq_name in list_of_seqs_of_pair:
                if seq_name in smp_one_abund_dict.keys():
                    sample_one_abundance_list.append(smp_one_abund_dict[seq_name])
                else:
                    sample_one_abundance_list.append(0)

            # then create a list of abundances for sample two by going through the above list and checking
            sample_two_abundance_list = []
            for seq_name in list_of_seqs_of_pair:
                if seq_name in smp_two_abund_dict.keys():
                    sample_two_abundance_list.append(smp_two_abund_dict[seq_name])
                else:
                    sample_two_abundance_list.append(0)

            # Do the Bray Curtis.
            distance = braycurtis(sample_one_abundance_list, sample_two_abundance_list)

            # Add the distance to the dictionary using both combinations of the sample names
            spp_distance_dict[f'{smp_one}_{smp_two}'] = distance
            spp_distance_dict[f'{smp_two}_{smp_one}'] = distance
        
        return spp_distance_dict

    # UNIFRAC METHODS
    def _do_unifrac_analysis(self):
        if os.path.isfile(self.pcoa_out_path):
            print(f'pcoa output file {self.pcoa_out_path} already exists. Skipping calculation.')
            return
        self._create_tree()
        self._compute_weighted_unifrac()
        self._write_out_unifrac_dist_file()
        self._make_pcoa_df_unifrac()

    def _create_tree(self):
        """Unifrac requires a tree."""
        # Create a fasta file that has the hashed sequence as the sequence names
        hash_cols, fasta_as_list = self._make_fasta_and_hash_names()
        # Set the df column names to the hashes
        self._set_df_cols_to_hashes(hash_cols)
        # Align the fasta
        self._align_seqs(fasta_as_list)
        # Make and root the tree
        self._make_and_root_tree()

    def _make_fasta_and_hash_names(self):
            # here we have a set of all of the sequences
            seq_fasta_list = []
            # we will change the df columns so that they match the seq names in the tree
            columns = []
            for seq in list(self.abundance_df):
                hash_of_seq = hashlib.md5(seq.encode()).hexdigest()
                if hash_of_seq in columns:
                    sys.exit('non-unique hash')
                seq_fasta_list.extend([f'>{hash_of_seq}', seq])
                columns.append(hash_of_seq)
            return columns, seq_fasta_list
    
    def _set_df_cols_to_hashes(self, columns):
            # change the columns of the df
            self.abundance_df.columns = columns

    def _align_seqs(self, seq_fasta_list):
            # Write out the unaligned fasta
            with open(self.unaligned_fasta_path, 'w') as f:
                for line in seq_fasta_list:
                    f.write(f'{line}\n')
            
            self._mafft_align_fasta(input_path=self.unaligned_fasta_path, output_path=self.aligned_fasta_path, method='auto', num_proc=6)

            # read in the aligned fasta
            with open(self.aligned_fasta_path, 'r') as f:
                aligned_fasta = [line.rstrip() for line in f]


            sequential_fasta = self._convert_interleaved_to_sequencial_fasta(aligned_fasta)
            
            # Write out the sequential aligned fasta to the same aligned fasta path
            with open(self.aligned_fasta_path, 'w') as f:
                for line in sequential_fasta:
                    f.write(f'{line}\n')

    def _make_and_root_tree(self):
        # make the tree
        print('Testing models and making phylogenetic tree')
        print('This could take some time...')
        # making the tree is very computationally expensive.
        # To see if we have a computed a tree of the aligned fasta in question,
        # we will name the output tree the md5sum of the aligned fasta file
        # This way we can check to see if there is already a tree that we can use
        # by looking for a file called <md5sum_of_fasta>.treefile
        # Fist get the md5sum of the aligned fasta
        # https://stackoverflow.com/questions/7829499/using-hashlib-to-compute-md5-digest-of-a-file-in-python-3
        hash_of_aligned_fasta = self._md5sum(self.aligned_fasta_path)
        if os.path.isfile(os.path.join(self.temp_dir, f'{hash_of_aligned_fasta}.treefile')):
            # Then we have already computed the tree and we can use this tree
            self.tree_path = os.path.join(self.temp_dir, f'{hash_of_aligned_fasta}.treefile')
            self.rooted_tree = TreeNode.read(self.tree_path)
        else:
            # Then we need to do the tree from scratch
            subprocess.run(
                ['iqtree', '-nt', 'AUTO', '-s', f'{self.aligned_fasta_path}', '-redo', '-mredo'])
            print('Tree creation complete')
            print('Rooting the tree at midpoint')
            tree = TreeNode.read(self.tree_path)
            self.rooted_tree = tree.root_at_midpoint()
            self.rooted_tree.write(self.tree_path)
            # And then rename the tree so that it is the md5sum of the aligned fasta
            os.rename(self.tree_path, os.path.join(self.temp_dir, f'{hash_of_aligned_fasta}.treefile'))
            self.tree_path = os.path.join(self.temp_dir, f'{hash_of_aligned_fasta}.treefile')

    @staticmethod
    def _md5sum(filename):
        with open(filename, mode='rb') as f:
            d = hashlib.md5()
            for buf in iter(partial(f.read, 128), b''):
                d.update(buf)
        return d.hexdigest()

    def _compute_weighted_unifrac(self):
        print('Performing unifrac calculations')
        self.wu = beta_diversity(
            metric='weighted_unifrac', counts=self.abundance_df.to_numpy(),
            ids=[str(_) for _ in list(self.abundance_df.index)],
            tree=self.rooted_tree, otu_ids=[str(_) for _ in list(self.abundance_df.columns)])

    def _write_out_unifrac_dist_file(self):
        wu_df = self.wu.to_data_frame()
        wu_df.to_csv(path_or_buf=self.dist_out_path, index=True, header=False)
        
    def _make_pcoa_df_unifrac(self):
        self.pcoa_df = self._do_pcoa_unifrac(self.wu)
        self.pcoa_df.to_csv(self.pcoa_out_path, index=True, header=True)

    def _do_pcoa_unifrac(self, wu):
        # compute the pcoa
        pcoa_output = pcoa(wu.data)
        self._rescale_pcoa(pcoa_output)
        pcoa_output.samples['sample'] = self.samples
        spp_unifrac_pcoa_df = pcoa_output.samples.set_index('sample')
        # now add the variance explained as a final row to the renamed_dataframe
        spp_unifrac_pcoa_df = spp_unifrac_pcoa_df.append(
            pcoa_output.proportion_explained.rename('proportion_explained'))
        return spp_unifrac_pcoa_df

    def _rescale_pcoa(self, pcoa_output):
        # work through the magnitudes of order and see what the bigest scaler we can work with is
        # whilst still remaining below 1
        query = 0.1
        scaler = 10
        while 1:
            if pcoa_output.samples.max().max() > query:
                # then we cannot multiply by the scaler
                # revert back and break
                scaler /= 10
                break
            else:
                # then we can safely multiply by the scaler
                # increase by order of magnitude and test again
                # we also need to test the negative if it is negative
                min_val = pcoa_output.samples.min().min()
                if min_val < 0:
                    if min_val > (-1 * query):
                        scaler *= 10
                        query /= 10
                    else:
                        scaler /= 10
                        break
                else:
                    scaler *= 10
                    query /= 10
        # now scale the df by the scaler unless it is 1
        if scaler != 1:
            pcoa_output.samples = pcoa_output.samples * scaler

    @staticmethod
    def _mafft_align_fasta(input_path, output_path, method='auto', mafft_exec_string='mafft', num_proc=1, iterations=1000):
        # http://plumbum.readthedocs.io/en/latest/local_commands.html#pipelining
        print(f'Aligning {input_path}')
        if method == 'auto':
            mafft = local[f'{mafft_exec_string}']
            (mafft['--auto', '--thread', f'{num_proc}', input_path] > output_path)()
        elif method == 'linsi':
            mafft = local[f'{mafft_exec_string}']
            (mafft['--localpair', '--maxiterate', f'{iterations}', '--thread', f'{num_proc}', input_path] > output_path)()
        elif method == 'unifrac':  # These are the alignment settings specifically for doing the unifrac alignments
            mafft = local[f'{mafft_exec_string}']
            (mafft[
                '--thread', f'{num_proc}', '--maxiterate', f'{iterations}',
                '--ep', '0', '--genafpair', input_path] > output_path)()
        print(f'Writing to {output_path}')

    @staticmethod
    def _convert_interleaved_to_sequencial_fasta(fasta_as_list):
        new_fasta = []
        temp_seq_string_list = []
        for fasta_line in fasta_as_list:
            if fasta_line.startswith('>'):
                if temp_seq_string_list:
                    new_fasta.append(''.join(temp_seq_string_list))
                    temp_seq_string_list = []
                    new_fasta.append(fasta_line)
                else:
                    new_fasta.append(fasta_line)
            else:
                temp_seq_string_list.append(fasta_line)
        new_fasta.append(''.join(temp_seq_string_list))
        return new_fasta

if __name__ == "__main__":
    dist = EighteenSDistance()
    dist.make_and_plot_dist_and_pcoa(resolution_type='host_only')