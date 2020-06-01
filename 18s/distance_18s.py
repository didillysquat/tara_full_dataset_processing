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
import shutil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import defaultdict, Counter
import re
import json

class EighteenSDistance(EighteenSBase):
    """
    This script can only be run once we have made the meta_info_table from output_tables.py
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

    #TODO we have looked at the distance plots containing the samples with different secondary sequences and we know that they are
    considerably different and cluster together. We should output a dist plot where these samples are included
    for supp materials.
    #TODO we now have Didier's distance matrices so we can read these in a process and cache them out. These will be our comparisons.
    We should also visualise these, and maybe create some kmers groupings from them so that we can try some categoristaion approaches
    although, Didier did say that he is also going to send through some categorisations.
    """
    def __init__(
        self, exclude_secondary_seq_samples=True, samples_at_least_threshold=0.0, 
        remove_majority_sequence=True, mafft_num_proc=6, braycurtis_normalisation_abundance=10000, 
        unifrac_normalisation_abundance=1000, normalisation_method='pwr', approach='dist', 
        only_snp_samples=False, use_replicates=False, snp_distance_type='biallelic', min_num_distinct_seqs_per_sample=3, most_abund_seq_cutoff=0):
        super().__init__()
        # Overwrite self.genera to only include Pocillopora and Porites as we are not currently working with Millepora
        self.genera = ['Pocillopora', 'Porites']
        # Use the meta_info_table produced by output_tables.py
        self.meta_info_table_path = os.path.join(self.output_dir, 'coral_18S_meta_info_table_2020-04-15T12_53_51.189071UTC.csv')
        self.meta_info_df = self._get_meta_info_df()
        # Dict where genus if key and primary seq (nucleotide seq) is value.
        self.primary_seq_dict = self._make_primary_seq_dict()
        # Here we can make use of the meta_info_table that we created in output_tables.py
        self.exclude_secondary_seq_samples = exclude_secondary_seq_samples
        self.only_snp_samples = only_snp_samples
        self.use_replicates = use_replicates
        self.samples_at_least_threshold = samples_at_least_threshold
        self.most_abund_seq_cutoff = most_abund_seq_cutoff
        if self.most_abund_seq_cutoff > 0 and self.samples_at_least_threshold != 0:
            print('Both samples_at_least_threshold and most_abund_seq_cutoff have been passed.')
            print('most_abund_seq_cutoff will be used and samples_at_least_threshold will be ignored.')
            print(f'uing most_abund_seq_cutoff = {self.most_abund_seq_cutoff}')
            self.samples_at_least_threshold = 0
        self.min_num_distinct_seqs_per_sample = min_num_distinct_seqs_per_sample
        if self.most_abund_seq_cutoff !=0 and self.most_abund_seq_cutoff < self.min_num_distinct_seqs_per_sample:
            raise RuntimeError(f'most_abund_seq_cutoff ({self.most_abund_seq_cutoff}) <' 
                                f' min_num_distinct_seqs_per_sample ({self.min_num_distinct_seqs_per_sample}).'
                                '\nPlease adjust one of the prarmeters accordingly.')
        
        try:
            assert(0 <= self.samples_at_least_threshold <= 1)
        except AssertionError:
            raise AssertionError('samples_at_least_threshold must be between 0 and 1')
        if snp_distance_type == 'biallelic':
            self.snp_dist_df_dict, self.snp_sample_list_dict = self._generate_biallelic_snp_dist_dfs()
        # A dict of genus to the list of representative readsets to include
        self.genus_to_representative_readset_dict = self._generate_genera_readset_dict()
        self.dist_methods = ['unifrac', 'braycurtis']
        self.remove_maj_seq = remove_majority_sequence
        self.mafft_num_proc = mafft_num_proc
        self.unifrac_norm_abund = unifrac_normalisation_abundance
        self.braycurtis_norm_abund = braycurtis_normalisation_abundance
        self.normalisation_method = normalisation_method
        self.approach = approach
        
        foo = 'bar'
    
    def _generate_biallelic_snp_dist_dfs(self):
        poc_snp_dist_path = os.path.join(self.input_dir_18s, 'snp_dist_matrices', 'MP_PocG111_biallelic_gaps.tree.distances.txt')
        por_snp_dist_path = os.path.join(self.input_dir_18s, 'snp_dist_matrices', 'MP_PorG111_biallelic_gaps.tree.distances.txt')
        snp_dist_df_dict = {}
        snp_sample_list_dict = {}
        for dist_path, genus in zip([poc_snp_dist_path, por_snp_dist_path], self.genera):
            # Read in the distance file. There is no header so we will read in as list and then process
            with open(dist_path, 'r') as f:
                dat = [line.rstrip().split('\t') for line in f]
            # Indices are currently in format (e.g. I01S01C011POR). Need to convert to sample-id.
            index = [_[0] for _ in dat]
            index = self._convert_index_to_sample_ids(index)
            dat = np.array([_[1:] for _ in dat], dtype=float).astype(int)
            # make the df
            snp_dist_df_dict[genus] = pd.DataFrame(data=dat, columns=index, index=index).astype(int)
            snp_sample_list_dict[genus] = index
        return snp_dist_df_dict, snp_sample_list_dict

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

    def _get_meta_info_df(self):
        return pd.read_csv(self.meta_info_table_path).set_index('readset')

    def _make_primary_seq_dict(self):
        # Get the primary sequences for each of the genera
        # This is the sequence that is found as the most abundant for the largest number of readsets
        # for a given genus.
        primary_seq_dict = {}  # genus key to primary seq value
        for genus in self.genera:
            primary_seq_dict[genus] = self.meta_info_df[
                self.meta_info_df['genetic_18S_genus_taxonomic_annotation'] == genus
                ]['primary_sequence'].mode().values.tolist()[0]
        return primary_seq_dict

    def _get_abundance_info_df(self):
        if os.path.isfile(os.path.join(self.cache_dir, 'abundance_info_df.p.bz')):
            return compress_pickle.load(os.path.join(self.cache_dir, 'abundance_info_df.p.bz'))
        else:
            raise RuntimeError('abundance_info_df.p.bz should have been created as part of processing_18s.py')

    def compute_distances(self):
        for dist_method in self.dist_methods:
            for genus, readset_list in self.genus_to_representative_readset_dict.items():
                indi_dist = IndiDistanceAnalysis(parent=self,
                                                 genus=genus, sample_list=readset_list,
                                                 dist_method=dist_method
                                                 )
                indi_dist.do_analysis()

    def plot_computed_distances(self):
        """Plot up the PCoA outputs from the computed distances.
        To start with we will work with a simple 2x3 plot
        Options for color by are 'none', 'island', 'maj_seq', 'intra_diversity'"""
        for color_by in ['island', 'tech_reps', None,  'maj_seq', 'intra_diversity']:
            dist_plotter = DistancePlotter(
                parent=self, color_by=color_by, plot_unifrac=True
            )
            dist_plotter.plot()

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
    
    def make_and_plot_dist_and_pcoa(self):
        # Call the class that will be responsible for calculating the distances
        # Calculating the PCoA and then plotting
        self.compute_distances()
        self.plot_computed_distances()

    def _generate_genera_readset_dict(self):
        """
        For each target coral (i.e. Pocillopora and Porites) get the list of readsets that we should be working with
        We will work with readsets as there are sometimes > 1 readset per sample-id.
        Start with all sequences and exclude according to:
        Target coral
        exclude_secondary_seq_samples
        use_replicates
        only_snp_samples
        """
        temp_dict = {}
        # First screen by coral host
        for genus in self.genera:
            temp_dict[genus] = self.meta_info_df[self.meta_info_df['genetic_18S_genus_taxonomic_annotation'] == genus]
        # Then by replicates
        if not self.use_replicates:
            for genus in self.genera:
                temp_dict[genus] = temp_dict[genus][temp_dict[genus]['is_representative_for_sample'] == True]
        # Then by exclude_secondary_seq_samples
        if self.exclude_secondary_seq_samples:
            for genus in self.genera:
                temp_dict[genus] = temp_dict[genus][temp_dict[genus]['is_different_primary_sequence'] == False]
        # Then by only_snp_samples
        if self.only_snp_samples:
            for genus in self.genera:
                temp_dict[genus] = temp_dict[genus][temp_dict[genus]['sample-id'].isin(self.snp_sample_list_dict[genus])]
        # Finally, we only want the list of the readsets, so get the indices from the df
        for genus in self.genera:
            temp_dict[genus] = list(temp_dict[genus].index)
        return temp_dict

class DistancePlotter:
    def __init__(self, parent, color_by=None, plot_unifrac=True):
        self.parent = parent
        if plot_unifrac:
            self.fig = plt.figure(figsize=(10, 7))
        else:
            self.fig = plt.figure(figsize=(12, 4))
        self.color_by = color_by
        self.color_dict = self._init_color_dict()
        self.plot_unifrac = plot_unifrac

    def _init_color_dict(self):
        color_list = get_colour_list()

        if self.color_by == 'island':
            island_c_dict = {island: color for island, color in zip(self.parent.info_df['island'].unique(), color_list[:len(self.parent.info_df['island'].unique())])}
            return {sample_name: island_c_dict[self.parent.info_df.at[sample_name, 'island']] for sample_name in self.parent.info_df.index}
        elif self.color_by is None:
            return {sample_name: 'black' for sample_name in self.parent.info_df.index}
        elif self.color_by == 'maj_seq':
            maj_seq_c_dict = {maj_seq: color for maj_seq, color in
                             zip(self.parent.info_df['most_abund_seq_of_coral_genus'].unique(), color_list[:len(self.parent.info_df['most_abund_seq_of_coral_genus'].unique())])}
            return {sample_name: maj_seq_c_dict[self.parent.info_df.at[sample_name, 'most_abund_seq_of_coral_genus']] for sample_name in
                    self.parent.info_df.index}
        elif self.color_by == 'intra_diversity':
            # If we are colouring by intra div then we will use the return a dict of
            # sample_name to number of unique sequences
            return self._create_intra_diversity_dict()
        elif self.color_by == 'tech_reps':
            # We are implementing this coloring as a sanity check for the unifrac
            # find all the tech reps and give them the same color
            # Keep the same color dict
            if os.path.isfile(os.path.join(self.parent.cache_dir, 'tech_reps_color_dict.p.bz')):
                return compress_pickle.load(os.path.join(self.parent.cache_dir, 'tech_reps_color_dict.p.bz'))
            tech_rep_base = set()
            color_dict = {}
            for sample_name in self.parent.info_df.index:
                if sample_name[-2] == '_':
                    # Then this is a tech rep sample name
                    tech_rep_base.add('_'.join(sample_name.split('_')[:-1]))
            # Just do the first 10 tech reps
            tech_rep_color_dict = {tech_rep: color for tech_rep, color in zip(list(tech_rep_base)[:25], color_list[:25])}
            for sample_name in self.parent.info_df.index:
                if sample_name[-2] == '_':
                    try:
                        # Then this is tech rep
                        color_dict[sample_name] = tech_rep_color_dict['_'.join(sample_name.split('_')[:-1])]
                        foo = 'asdf'
                    except KeyError:
                        color_dict[sample_name] = (1,1,1,0)
                else:
                    color_dict[sample_name] = (1,1,1,0)
            compress_pickle.dump(color_dict, os.path.join(self.parent.cache_dir, 'tech_reps_color_dict.p.bz'))
            return color_dict
        else:
            raise NotImplementedError


    def _create_intra_diversity_dict(self):
        diversity_dict = {}
        print('Creating intra diversity dict')
        for sample_name in self.parent.info_df.index:
            sys.stdout.write(f'\r{sample_name}')
            sample_qc_dir = os.path.join(self.parent.qc_dir, sample_name)
            consolidated_host_seqs_abund_dict = compress_pickle.load(
                os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))

            # We need to remove the most abundant sequence from the equation
            most_abund_sequence = max(consolidated_host_seqs_abund_dict.keys(),
                                      key=(lambda key: consolidated_host_seqs_abund_dict[key]))
            # remove the most abund seq
            del consolidated_host_seqs_abund_dict[most_abund_sequence]
            diversity_dict[sample_name] = len(consolidated_host_seqs_abund_dict.keys())
        return diversity_dict

    def plot(self):
        if self.plot_unifrac:
            self.gs = gridspec.GridSpec(3, 3, height_ratios=[0.2, 1, 1])
        else:
            self.gs = gridspec.GridSpec(2, 3, height_ratios=[0.2,1])
        if self.plot_unifrac:
            dist_methods = ['braycurtis', 'unifrac']
        else:
            dist_methods = ['braycurtis']

        if self.parent.exclude_secondary_seq_samples:
            secondary_seq_string = 'exluded'
        else:
            secondary_seq_string = 'included'


        for i, dist_method in enumerate(dist_methods):
            for j, genus in enumerate(['Pocillopora', 'Millepora', 'Porites']):
                ax = self.fig.add_subplot(self.gs[i + 1, j])
                try:
                    indi_dist_plotter = IndiDistancePlotter(
                        eightsdistparent=self.parent, fig=self.fig, ax=ax, genus=genus, dist_method=dist_method,
                        path_to_csv=os.path.join(
                            self.parent.output_dir_18s,
                            f'secondary_seq_{secondary_seq_string}_'
                            f'{self.parent.samples_at_least_threshold}_{genus}_{dist_method}.csv'),
                        color_dict=self.color_dict)
                    indi_dist_plotter.plot()
                except FileNotFoundError as e:
                    # Then the .csv hasn't been output yet
                    print(e)
                    continue

        if self.plot_unifrac:
            with_unifrac_string = '_w_unifrac'
        else:
            with_unifrac_string = ''
        title_ax = self.fig.add_subplot(self.gs[0,:])
        if self.parent.exclude_secondary_seq_samples:
            secondary_seq_string = 'exluded'
        else:
            secondary_seq_string = 'included'

        title_ax.text(x=0.5, y=0.5, s=f'secondary_seq_{secondary_seq_string}_{self.parent.samples_at_least_threshold}_{self.color_by}{with_unifrac_string}', ha='center', va='center')
        remove_axes_but_allow_labels(title_ax)
        plt.tight_layout()
        print(f'saving secondary_seq_{secondary_seq_string}_{self.parent.samples_at_least_threshold}_{self.color_by}_pcoa{with_unifrac_string}.png')
        plt.savefig(os.path.join(self.parent.figure_dir, f'secondary_seq_{secondary_seq_string}_{self.parent.samples_at_least_threshold}_{self.color_by}_pcoa{with_unifrac_string}.png'), dpi=1200)
        print(f'saving secondary_seq_{secondary_seq_string}_{self.parent.samples_at_least_threshold}_{self.color_by}_pcoa{with_unifrac_string}.svg')
        plt.savefig(os.path.join(self.parent.figure_dir, f'secondary_seq_{secondary_seq_string}_{self.parent.samples_at_least_threshold}_{self.color_by}_pcoa{with_unifrac_string}.svg'))


class IndiDistancePlotter:
    """Class for plotting an individual axis"""
    def __init__(self, eightsdistparent, fig, ax, genus, dist_method, path_to_csv, color_dict):
        self.eightsdistparent = eightsdistparent
        self.fig = fig
        self.ax = ax
        self.category = genus
        self.dist_method = dist_method
        self.pcoa_df = pd.read_csv(path_to_csv)
        self.color_dict = color_dict


        # if self.category == 'Porites':
        #     if self.dist_method == 'unifrac':
        #         # First get the samples that are in the bottom right and get the seqs that are found in common
        #         master_set = None
        #         for i, sample_name in enumerate(self.pcoa_df[self.pcoa_df['PC1'] > 0.18]['sample'].values.tolist()[:-1]):
        #             master_set = self.get_seqs_in_common(i, master_set, sample_name)
        #         self.print_master_set(master_set=master_set, prefix='br')
        #         # second get the samples that are at the top center and get the seqs that are found in common
        #         master_set = None
        #         for i, sample_name in enumerate(
        #                 self.pcoa_df[self.pcoa_df['PC2'] > 0.05]['sample'].values.tolist()[:-1]):
        #             master_set = self.get_seqs_in_common(i, master_set, sample_name)
        #         self.print_master_set(master_set=master_set, prefix='top')
        #
        #         # Go through all samples and get an abundance dict of the sequences
        #         # We could even plot this is a hist
        #         master_count = defaultdict(int)
        #         for sample_name in self.info_df.index:
        #             sample_qc_dir = os.path.join(self.qc_dir, sample_name)
        #             consolidated_host_seqs_abund_dict = compress_pickle.load(
        #                 os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
        #             for seq in consolidated_host_seqs_abund_dict.keys():
        #                 master_count[seq] += 1
        #         # now sort the list
        #         hist_list = [val for val in master_count.values()]
        #         # plt.hist(hist_list, bins=100)
        #         sorted_rev = sorted([_ for _ in master_count.items()], key=lambda x: x[1], reverse=False)
        #         # Want to do a scatter plot of abundance cutoff effect on the total number of seqs we are working with
        #         cutoff = 0
        #         cutoff_list = [_ for _ in master_count.items() if _[1] > cutoff]
        #         x = []
        #         y = []
        #         while cutoff_list:
        #             print(cutoff)
        #             x.append(cutoff)
        #             y.append(len(cutoff_list))
        #             cutoff += 1
        #             cutoff_list = [_ for _ in master_count.items() if _[1] > cutoff]
        #         plt.scatter(x=x, y=y, c='black', s=2)
        # self.opt = 'foo'

    def print_master_set(self, master_set, prefix):
        for i, seq in enumerate(list(master_set)):
            print(f'>{prefix}_{i}')
            print(seq)

    def get_seqs_in_common(self, i, master_set, sample_name):
        sample_qc_dir = os.path.join(self.qc_dir, sample_name)
        consolidated_host_seqs_abund_dict = compress_pickle.load(
            os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
        # We need to remove the most abundant sequence from the equation
        most_abund_sequence = max(consolidated_host_seqs_abund_dict.keys(),
                                  key=(lambda key: consolidated_host_seqs_abund_dict[key]))
        # remove the most abund seq
        del consolidated_host_seqs_abund_dict[most_abund_sequence]
        if i == 0:
            master_set = set(consolidated_host_seqs_abund_dict.keys())
        else:
            master_set = master_set.intersection(set(consolidated_host_seqs_abund_dict.keys()))
            print(f'now {len(master_set)} elements')
        island = self.info_df.at[sample_name, "island"]
        site = self.info_df.at[sample_name, "site"]
        print(f'{sample_name}: {island} {site}')
        return master_set

    def plot(self):
        if self.color_by == 'intra_diversity':
        # Then we need to use a color map
        # https://medium.com/better-programming/how-to-use-colormaps-with-matplotlib-to-create-colorful-plots-in-python-969b5a892f0c
            im = self.ax.scatter(x=self.pcoa_df['PC1'].values[:-1], y=self.pcoa_df['PC2'].values[:-1],
                            c=[self.color_dict[sample_name] for sample_name in self.pcoa_df['sample'][:-1].tolist()], cmap='cool', s=2)
            cbar = self.ax.figure.colorbar(im, ax=self.ax)
            cbar.ax.set_ylabel('unique seqs', rotation=-90, va='bottom')
        elif self.color_by=='tech_reps':
            self.ax.scatter(x=self.pcoa_df['PC1'].values[:-1], y=self.pcoa_df['PC2'].values[:-1],
                            c=[self.color_dict[sample_name] for sample_name in self.pcoa_df['sample'][:-1].tolist()],
                            s=8)
        else:
            self.ax.scatter(x=self.pcoa_df['PC1'].values[:-1], y=self.pcoa_df['PC2'].values[:-1], c=[self.color_dict[sample_name] for sample_name in self.pcoa_df['sample'][:-1].tolist()], s=2)
        self.ax.set_title(f'{self.category}_{self.dist_method}')
        pc1_var = format(self.pcoa_df['PC1'].values[-1], '.3f')
        pc2_var = format(self.pcoa_df['PC2'].values[-1], '.3f')
        self.ax.set_xlabel(f'PC1 var explained: {pc1_var}')
        self.ax.set_ylabel(f'PC2 var explained: {pc2_var}')
        foo = 'bar'


class IndiDistanceAnalysis():
    def __init__(self, parent, genus, sample_list, dist_method):
        """This is the class that will take care of creating a single between sample
        distance method

        Also use compression when writing out and then reading in the pcoa output files and the distance files
        """
        self.parent = parent
        self.genus = genus
        self.samples = sample_list
        self.dist_method = dist_method
        # A temp directory where we can write out the unaligned fastas
        # and aligned fasta and any other intermediate files
        self.temp_dir = os.path.join(os.path.dirname(self.parent.qc_dir), 'temp')
        os.makedirs(self.temp_dir, exist_ok=True)
        # In the initial 18s we were using a normalisation of 10000 seqs for the braycurtis
        # but 1000 for the unifrac due to the amount of time it was taking to create
        # the trees. We will start here with 1000 seqs for the unifrac and see
        # how long it takes to make the tree. If its not too bad then we can up the number
        if self.dist_method == 'unifrac':
            self.num_seqs_to_normalise_to = self.parent.unifrac_norm_abund
            self.abundance_df = None
        else:
            self.num_seqs_to_normalise_to = self.parent.braycurtis_norm_abund
            self.abundance_dict = None


        # Generic variables shared between braycurtis and unifrac
        if self.parent.exclude_secondary_seq_samples:
            secondary_seq_string = 'exluded'
        else:
            secondary_seq_string = 'included'
        self.dist_out_path = os.path.join(
            self.parent.output_dir_18s,
            f'secondary_seq_{secondary_seq_string}_maj_seq_remove_{self.parent.remove_maj_seq}_min_samp_{self.parent.samples_at_least_threshold}_{self.genus}_{self.parent.normalisation_method}_{self.dist_method}.dist')
        self.pcoa_out_path = os.path.join(
            self.parent.output_dir_18s,
            f'secondary_seq_{secondary_seq_string}_maj_seq_remove_{self.parent.remove_maj_seq}_min_smp_{self.parent.samples_at_least_threshold}_{self.genus}_{self.parent.normalisation_method}_{self.dist_method}.csv')

        self.pcoa_df = None

        # Variables concerned with unifrac
        self.unaligned_fasta_path = os.path.join(self.temp_dir, 'unaligned_fasta.fasta')
        self.aligned_fasta_path = os.path.join(self.temp_dir, 'aligned_fasta.fasta')
        self.tree_path = self.aligned_fasta_path + '.treefile'
        self.wu = None

        # Variables concerned with braycurtis
        self.braycurtis_btwn_sample_distance_dictionary = None

    def _create_abundance_df(self, obj_to_return):
        """For each sample in self.samples, load up the consolidated_host_seqs_abund_dict
        and create a dict that is sequence to normalised abundance. Then add this dictionary
        to a master dictionary where sample_name is key. Then finally create a df from this
        dict of dicts.
        
        If obj_to_return=='abund_dict', we will return the master dictionary
        else obj_to_return=='abund_df', we will return a df made from the dict
        """
        if self.parent.samples_at_least_threshold > 0: # TODO plot cutoff against mean seqs remaining and std.
            seq_to_sample_occurence_dict = self._get_seq_to_sample_occurence_dict()
            threshold_set = {k for k, v in seq_to_sample_occurence_dict.items() if v > self.parent.samples_at_least_threshold}
        dict_to_create_df_from = {}
        print(f'Creating {obj_to_return}')
        # NB there is one sample that only had one sequence in the consolidated_host_seqs_abund_dict
        # and this sample is causing errors. We will check for empty conolidated dict,
        # add the sample to a list and then remove it from the self.samples list
        samples_to_remove_list = []
        for sample_name in self.samples:
            sys.stdout.write(f'\r{sample_name}')
            
            sample_qc_dir = os.path.join(self.parent.qc_dir, sample_name)
            consolidated_host_seqs_abund_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
            
            # remove most abundant sequences if this option is set
            if self.parent.remove_maj_seq:
                most_abund_sequence = max(consolidated_host_seqs_abund_dict.keys(), 
                                        key=(lambda key: consolidated_host_seqs_abund_dict[key]))
                # remove the most abund seq
                del consolidated_host_seqs_abund_dict[most_abund_sequence]
            
            # if working with samples_at_least_threshold, screen out rare seqs here
            if self.parent.samples_at_least_threshold > 0:
                consolidated_host_seqs_abund_dict = {k: v for k, v in consolidated_host_seqs_abund_dict.items() if k in threshold_set}
            elif self.parent.most_abund_seq_cutoff > 0:
                sorted_dict_keys = sorted(consolidated_host_seqs_abund_dict, key=consolidated_host_seqs_abund_dict.get, reverse=True)[:self.parent.most_abund_seq_cutoff]
                consolidated_host_seqs_abund_dict =  {k: v for k, v in consolidated_host_seqs_abund_dict.items() if k in sorted_dict_keys}
            
            # normalise the consolidated_host_seqs_abund_dict back to 1
            tot = sum(consolidated_host_seqs_abund_dict.values())
            consolidated_host_seqs_abund_dict = {k: v/tot for k, v in consolidated_host_seqs_abund_dict.items()}

            # To prevent downstream problems we will insist on there always being at least
            # three sequences. Later we can put this to much higher to look at the effect
            if len(consolidated_host_seqs_abund_dict.keys()) < self.parent.min_num_distinct_seqs_per_sample:
                samples_to_remove_list.append(sample_name)
                continue

            # Here normalise acorrding to the given normalisation method
            if self.parent.normalisation_method == 'rai':
                # Relative abundance integer conversion
                temp_sample_dict = {}
                for sequence, rel_abund in consolidated_host_seqs_abund_dict.items():
                    normalised_abund = int(rel_abund*self.num_seqs_to_normalise_to)
                    if normalised_abund:
                        temp_sample_dict[sequence] = normalised_abund
            else:
                # pwr
                # pick without replacement.
                # If subsample, we will produce a list using np.random.choice
                # https://docs.scipy.org/doc/numpy-1.16.0/reference/generated/numpy.random.choice.html
                # We will then convert this to a dict using counter
                # This dict will have absolute abundances. These will be converted to relaive abundances
                # outside of this script.
                seqs, probs = zip(*consolidated_host_seqs_abund_dict.items())
                seqs_list = np.random.choice(seqs, self.num_seqs_to_normalise_to, p=probs)
                temp_sample_dict = dict(Counter(seqs_list))

            dict_to_create_df_from[sample_name] = temp_sample_dict
        
        for sample_name in samples_to_remove_list:
            self.samples.remove(sample_name)

        # # It may also be very helpful to look at the distribution of the number of minor sequences
        # # a given sample has
        # hist_list = [len(sub_dict.keys()) for sub_dict in dict_to_create_df_from.values()]
        # print('making and writing histogram of sequence diversity')
        # plt.hist(hist_list, bins=30)
        # plt.savefig(os.path.join(self.parent.output_dir_18s, f'seq_diversity_hist_3_cutoff_{self.category}_{self.dist_method}.png'))
        # plt.close()

        if obj_to_return == 'abund_dict':
            return dict_to_create_df_from
        else:
            df = pd.DataFrame.from_dict(dict_to_create_df_from, orient='index')
            df[pd.isna(df)] = 0
            print('\ndf creation complete\n')
            return df

    def _get_seq_to_sample_occurence_dict(self):
        """
        Return a dictionary with key as a given sequence, and the number of readsets that sequence was found
        in as the the value. We only want to consider the readsets that are to be used in this given ordination.
        To save recompution time, we will sort the name of readsets and get a md5sum of the resulting jsoned
        string. We will then save then pickle out the resultant seq_to_sample_occurence_dict using
        this hash.
        """
        print('Creating sample occurence dictionary')
        print('Checking to see if cache available')
        sample_list_hash = self._md5sum_from_python_object(self.samples)
        pickle_path_to_find = os.path.join(self.parent.cache_dir, f'{sample_list_hash}_seq_to_sample_occurence_dict.p.bz')
        if os.path.isfile(pickle_path_to_find):
            print('Cache is available. Loading from cache.')
            seq_to_sample_occurence_dict = compress_pickle.load(pickle_path_to_find)
        else:
            print('Cache is not available. Generating from scratch.')
            seq_to_sample_occurence_dict = defaultdict(int)
            for read_set in self.samples:
                sys.stdout.write('f\r{read_set}')
                sample_qc_dir = os.path.join(self.parent.qc_dir, read_set)
                consolidated_host_seqs_abund_dict = compress_pickle.load(
                    os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
                for seq_name in consolidated_host_seqs_abund_dict.keys():
                    seq_to_sample_occurence_dict[seq_name] += 1
            print('Complete')
            # Now convert the absolute abundances to relative abudances for the given set of readsets
            tot = len(self.samples)
            seq_to_sample_occurence_dict = {k: v/tot for k, v in seq_to_sample_occurence_dict.items()}
            compress_pickle.dump(seq_to_sample_occurence_dict, pickle_path_to_find)
        return seq_to_sample_occurence_dict

    def _create_abundance_dicts(self):
        """
        Method that produces the abundance dictionary that will hold a set of sub dictionaries
        one for each readset we are concerned with.
        """
        if self.parent.samples_at_least_threshold:
            seq_to_sample_occurence_dict = self._get_seq_to_sample_occurence_dict()
            threshold_set = {k for k, v in seq_to_sample_occurence_dict.items() if v > self.parent.samples_at_least_threshold}

        abundance_dict = {}
        print('Creating abundance dict')
        
        samples_to_remove_list = []
        for sample_name in self.samples:
            sys.stdout.write(f'\r{sample_name}')
            
            
            sample_qc_dir = os.path.join(self.parent.qc_dir, sample_name)
            consolidated_host_seqs_abund_dict = compress_pickle.load(
                os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))

            if self.parent.remove_maj_seq:
                # We need to remove the most abundant sequence from the equation
                most_abund_sequence = max(consolidated_host_seqs_abund_dict.keys(),
                                        key=(lambda key: consolidated_host_seqs_abund_dict[key]))
                # remove the most abund seq
                del consolidated_host_seqs_abund_dict[most_abund_sequence]
            
            # renormalise
            if self.parent.samples_at_least_threshold > 0:
                # screen out those sequences that are not found in X or more samples
                consolidated_host_seqs_abund_dict = {
                    k: v for k, v in consolidated_host_seqs_abund_dict.items() if k in threshold_set}
            
            # normalise the consolidated_host_seqs_abund_dict back to 1
            tot = sum(consolidated_host_seqs_abund_dict.values())
            consolidated_host_seqs_abund_dict = {k: v / tot for k, v in consolidated_host_seqs_abund_dict.items()}

            # To prevent downstream problems we will insist on there always being at least
            # three sequences. Later we can put this to much higher to look at the effect
            if len(consolidated_host_seqs_abund_dict.keys()) < self.parent.min_num_distinct_seqs_per_sample:
                samples_to_remove_list.append(sample_name)
                continue

            if self.parent.normalisation_method == 'rai':
                # Relative abundance integer conversion
                temp_sample_dict = {}
                for sequence, rel_abund in consolidated_host_seqs_abund_dict.items():
                    normalised_abund = int(rel_abund*self.num_seqs_to_normalise_to)
                    if normalised_abund:
                        temp_sample_dict[sequence] = normalised_abund
            else:
                # pwr
                # pick without replacement.
                # If subsample, we will produce a list using np.random.choice
                # https://docs.scipy.org/doc/numpy-1.16.0/reference/generated/numpy.random.choice.html
                # We will then convert this to a dict using counter
                # This dict will have absolute abundances. These will be converted to relaive abundances
                # outside of this script.
                seqs, probs = zip(*consolidated_host_seqs_abund_dict.items())
                seqs_list = np.random.choice(seqs, self.num_seqs_to_normalise_to, p=probs)
                temp_sample_dict = dict(Counter(seqs_list))
            
            abundance_dict[sample_name] = temp_sample_dict
        
        for sample_name in samples_to_remove_list:
            self.samples.remove(sample_name)

        # It may also be very helpful to look at the distribution of the number of minor sequences
        # a given sample has
        # hist_list = [len(sub_dict.keys()) for sub_dict in abundance_dict.values()]
        # print('making and writing histogram of sequence diversity')
        # plt.hist(hist_list, bins=30)
        # plt.savefig(
        #     os.path.join(self.parent.output_dir_18s, f'seq_diversity_hist_3_cutoff_{self.genus}_{self.dist_method}.png'))
        # plt.close()

        return abundance_dict

    def do_analysis(self):
        if self.dist_method == 'unifrac':
            self._do_unifrac_analysis()
        else:
            self._do_braycurtis_analysis()

    # BRAYCURTIS METHODS
    def _do_braycurtis_analysis(self):
        if self._check_if_pcoa_already_computed():
            return
        # Key is readset, value is dict
        # For the sub dicts, the key is nucleotide sequences and value is 
        # abundance of that sequence normalised to a given number of seuences that depend
        # on the distance calculation method. For braycurtis, this is currently set to 10000
        # or UniFrac this is currently set to 1000.
        self.abundance_dict = self._create_abundance_df(obj_to_return='abund_dict')
        self.braycurtis_btwn_sample_distance_dictionary = self._compute_braycurtis_distance_dict()
        self.braycurtis_btwn_sample_distance_file = self._make_and_write_braycurtis_distance_file()
        self.pcoa_df = self._make_pcoa_df_braycurtis()
        self.pcoa_df.to_csv(self.pcoa_out_path, index=True, header=True)

    def _make_pcoa_df_braycurtis(self):
        # simultaneously grab the sample names in the order of the distance matrix and put the matrix into
        # a twoD list and then convert to a numpy array
        temp_two_D_list = []
        sample_names_from_dist_matrix = []
        for line in self.braycurtis_btwn_sample_distance_file[1:]:
            temp_elements = line.split('\t')
            sample_names_from_dist_matrix.append(temp_elements[0])
            temp_two_D_list.append([float(a) for a in temp_elements[1:]])
        bc_dist_array = np.array(temp_two_D_list)
        sys.stdout.write('\rcalculating PCoA coordinates')

        pcoa_df = pcoa(bc_dist_array)

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
        print('Making braycurtis distance file')
        for sample_outer in self.samples:
            sys.stdout.write(f'\r{sample_outer}')
            # The list that will hold the line of distance. This line starts with the name of the sample
            temp_sample_dist_string = [sample_outer]
            temp_dist_list = [
                self.braycurtis_btwn_sample_distance_dictionary[frozenset([sample_outer, sample_inner])] if
                sample_outer != sample_inner else
                0 for sample_inner in self.samples
            ]
            temp_sample_dist_string.extend(temp_dist_list)

            distance_out_file.append(
                '\t'.join([str(distance_item) for distance_item in temp_sample_dist_string]))
        print('\nwriting dist file')
        with open(self.dist_out_path, 'w') as f:
            for line in distance_out_file:
                f.write(f'{line}\n')
        print('complete')
        return distance_out_file

    def _compute_braycurtis_distance_dict(self):
        # NB I did re-write this function using the abundance df that I made for the unifrac
        # clculations. I have left this method in below this one. It was much slower though so we
        # will stick with this. I have re-written this a bit to speed it up.
        # Create a dictionary that will hold the distance between the two samples
        distance_dict = {}
        # For pairwise comparison of each of these sequences
        print('\nComputing braycurtis distance dictionary\n')
        tot = len([_ for _ in itertools.combinations(self.samples, 2)])
        count = 0
        for smp_one, smp_two in itertools.combinations(self.samples, 2):
            count += 1
            sys.stdout.write(f'\r{count} out of {tot}: {smp_one}_{smp_two}')

            # Get a set of the sequences found in either one of the samples
            smp_one_abund_dict = self.abundance_dict[smp_one]
            smp_two_abund_dict = self.abundance_dict[smp_two]
            list_of_seqs_of_pair = set(smp_one_abund_dict.keys())
            list_of_seqs_of_pair.update(smp_two_abund_dict.keys())
            list_of_seqs_of_pair =  list(list_of_seqs_of_pair)

            # Create lists of abundances in order of the seq in list_of_seqs_of_pairs
            sample_one_abundance_list = [
                smp_one_abund_dict[seq_name] if seq_name in smp_one_abund_dict else
                0 for seq_name in list_of_seqs_of_pair]
            sample_two_abundance_list = \
                [smp_two_abund_dict[seq_name] if seq_name in smp_two_abund_dict else
                 0 for seq_name in list_of_seqs_of_pair]

            # Do the Bray Curtis.
            distance = braycurtis(sample_one_abundance_list, sample_two_abundance_list)

            # Add the distance to the dictionary
            distance_dict[frozenset([smp_one, smp_two])] = distance

        print('\nComplete')
        return distance_dict

    def _compute_braycurtis_distance_dict_using_df(self):
        # NB Although this code looks great. It is very slow. The second sub_df assignment
        # takes a second or two to do. The above messy code using dicts is much faster.
        # As such I will stick with that.
        # Create a dictionary that will hold the distance between the two samples
        spp_distance_dict = {}
        # For pairwise comparison of each of these sequences
        for smp_one, smp_two in itertools.combinations(self.samples, 2):
            print('Calculating distance for {}_{}'.format(smp_one, smp_two))
            # df containing only the sequences found in either smp_one or smp_two
            sub_df = self.abundance_df.loc[[smp_one, smp_two]]
            sub_df = sub_df.loc[:, (sub_df == 0).all(axis=0)]

            # Do the Bray Curtis.
            distance = braycurtis(sub_df.loc[smp_one].values.tolist(), sub_df.loc[smp_two].values.tolist())

            # Add the distance to the dictionary
            spp_distance_dict[frozenset([smp_one, smp_two])] = distance

        return spp_distance_dict

    # UNIFRAC METHODS
    def _do_unifrac_analysis(self):
        if self._check_if_pcoa_already_computed():
            return
        self.abundance_df = self._create_abundance_df(obj_to_return='abund_df')
        self._create_tree()
        self._compute_weighted_unifrac()
        self._write_out_unifrac_dist_file()
        self._make_pcoa_df_unifrac()
        self._clean_temp_dir()

    def _clean_temp_dir(self):
        shutil.rmtree(self.temp_dir)
        os.makedirs(self.temp_dir, exist_ok=True)

    def _check_if_pcoa_already_computed(self):
        if os.path.isfile(self.pcoa_out_path):
            print(f'pcoa output file {self.pcoa_out_path} already exists. Skipping calculation.')
            return True

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
            
            self._mafft_align_fasta(input_path=self.unaligned_fasta_path, output_path=self.aligned_fasta_path, method='auto', num_proc=self.parent.mafft_num_proc)

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
        if os.path.isfile(os.path.join(self.parent.cache_dir, f'{hash_of_aligned_fasta}.treefile')):
            # Then we have already computed the tree and we can use this tree
            self.tree_path = os.path.join(self.parent.cache_dir, f'{hash_of_aligned_fasta}.treefile')
            self.rooted_tree = TreeNode.read(self.tree_path)
        else:
            # Then we need to do the tree from scratch
            subprocess.run(
                ['iqtree', '-T', 'AUTO', '--threads-max', '2', '-s', f'{self.aligned_fasta_path}', '-redo', '-mredo'])
            print('Tree creation complete')
            print('Rooting the tree at midpoint')
            tree = TreeNode.read(self.tree_path)
            self.rooted_tree = tree.root_at_midpoint()
            self.rooted_tree.write(self.tree_path)
            # And then rename the tree so that it is the md5sum of the aligned fasta
            os.rename(self.tree_path, os.path.join(self.parent.cache_dir, f'{hash_of_aligned_fasta}.treefile'))
            self.tree_path = os.path.join(self.parent.cache_dir, f'{hash_of_aligned_fasta}.treefile')

    def _md5sum_from_python_object(self, p_obj):
        """ 
        A wrapper function around the self._md5sum function, that will take in a python object,
        sort it, write it out to temp, md5sum the written out file, delete the file,
        and return the md5sum hash.
        """
        sorted_p_obj = sorted(p_obj)
        temp_out_path = os.path.join(self.temp_dir, 'p_obj.out')
        with open(temp_out_path, 'w') as f:
            json.dump(sorted_p_obj, f)
        md5sum_hash = self._md5sum(temp_out_path)
        os.remove(temp_out_path)
        return md5sum_hash
    
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
        self.wu_df = self.wu.to_data_frame()
        self.wu_df.to_csv(path_or_buf=self.dist_out_path, index=True, header=False)
        
    def _make_pcoa_df_unifrac(self):
        self.pcoa_df = self._do_pcoa_unifrac()
        self.pcoa_df.to_csv(self.pcoa_out_path, index=True, header=True)

    def _do_pcoa_unifrac(self):
        # compute the pcoa
        pcoa_output = pcoa(self.wu.data)
        self._rescale_pcoa(pcoa_output)

        # NB We were originally assigning samples as self.samples
        # This was wrong and was causing a bad bug. The order of self.samples
        # is not the same order as the order used in the PCoA index.
        # Rather we need to use the index used in any one of the dataframes we created
        pcoa_output.samples['sample'] = self.wu_df.index.values.tolist()

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

def get_colour_list():
    colour_list = ["#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", "#FFDBE5",
                "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693",
                "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900",
                "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
                "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744",
                "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68",
                "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED",
                "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", "#7900D7",
                "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F",
                "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55",
                "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C", "#83AB58", "#001C1E", "#D1F7CE", "#004B28",
                "#C8D0F6", "#A3A489", "#806C66", "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59",
                "#8ADBB4", "#1E0200", "#5B4E51", "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94",
                "#7ED379", "#012C58", "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393",
                "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
                "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5", "#E773CE",
                "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4", "#00005F", "#A97399",
                "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01", "#6B94AA", "#51A058", "#A45B02",
                "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966", "#64547B", "#97979E", "#006A66", "#391406",
                "#F4D749", "#0045D2", "#006C31", "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9",
                "#FFFFFE", "#C6DC99", "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527",
                "#8BB400", "#797868", "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C",
                "#B88183", "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
                "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F", "#003109",
                "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E", "#1A3A2A", "#494B5A",
                "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F", "#BDC9D2", "#9FA064", "#BE4700",
                "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00", "#061203", "#DFFB71", "#868E7E", "#98D058",
                "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66", "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F",
                "#545C46", "#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"]
    return colour_list

def remove_axes_but_allow_labels(ax):
    ax.set_frame_on(False)
    ax.set_yticks([])
    ax.set_xticks([])

if __name__ == "__main__":
    """
    samples_at_least_threshold = The minimum number of samples (as a realtive abundance of the total
    number of readsets that we will be working with for the given configuratin of samples) that a given sequence
    must be found in for it to be considered part of the calculation of the distance matrices. E.g. at 0.5, a
    sequence must be found in half of the readsets we are working with for the ordination. [0.0]

    remove_majority_sequence = Whether to remove the most abundant sequence from each of the readsets before
    calculating the distance ordinations.

    exclude_secondary_seq_samples = Whether to exclude samples that returned a different most abundant sequence
    compared to the most commonly abundant sequence for the coral taxa in question.

    braycurtis_normalisation_abundance = The absolute abundance that normalisation will be conducted to for braycurtis [10000].

    unifrac_normalisation_abundance = The absolute abundance that the normalistaion will be conducted to for unifrac [1000].

    normalisation method = either 'pwr' (pick without replacement), 'rai' (realtive abundance integer convertion) [pwr]

    mafft_num_proc = the numper of processes to pass to MAFFT for doing the alignment required for tree creation for
    unifrac-based distance calculations [6].

    approach = 'cat' (category), 'dist' (distance) [dist]. Whether we are doing categorical or distance comparisons.

    only_snp_samples = whether we are only working with the samples that have snp data available for them [False]

    use_replicates = whether to use all available technical replicates for samples (i.e. multiple readsets per
    sample) [False]. If False we will use the representative readset noted in the fastq_info_df

    TODO
    exclude_no_use_samples = whether to exlude samples (readsets) that are listed as use==False in the fastq_info_df

    snp_distance_type = The type of distance matrix we will use for the snp data [biallelic].

    min_num_distinct_seqs_per_sample = The minimum number of distinct sequences a given readset must have in its
    abundance dictionary after normalisation. Else the readset will be discarded [3].

    TODO
    most_abund_seq_cutoff = If this is 0. It will effectively be turned off. If this  is any other value, then
    a cutoff will be applied so that only the i most abundant sequences in a given sample will be considered. This
    value must be greaer than the min_num_distinct_seqs_per_sample else an error will be thrown. If a non-zero
    value is provided to this argument. It will override the samples_at_lesast_threshold. [0]
    """
    dist = EighteenSDistance(
        exclude_secondary_seq_samples=True, samples_at_least_threshold=0.5, remove_majority_sequence=True, most_abund_seq_cutoff=4)
    # Options for resolution_type are:
    # 'host_only' = 0 filtering of the seqs in the samples and they are only separated by
    # the gentically identified majority genus.

    dist.make_and_plot_dist_and_pcoa()