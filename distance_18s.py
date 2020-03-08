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
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import defaultdict

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
                most_abund_seq_of_coral_genus_df_list.append(sorted([_ for _ in consolidated_host_seqs_abund_dict.items()], key=lambda x: x[1], reverse=True)[0][0])
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
        da = DistanceAnlyses(
            resolution_type=resolution_type, info_df=self.info_df,
            qc_dir=self.qc_dir, output_dir=self.output_dir, cache_dir=self.cache_dir,
            figure_dir=self.fig_output_dir
        )
        da.compute_distances()
        # TODO we can do the plotting here
        da.plot_computed_distances()
        foo = 'bar'

class DistanceAnlyses:
    def __init__(self, resolution_type, info_df, qc_dir, output_dir, cache_dir, figure_dir):
        self.resolution_type = resolution_type
        self.qc_dir = qc_dir
        self.info_df = info_df
        self.output_dir = output_dir
        self.cache_dir = cache_dir
        self.figure_dir = figure_dir
        # get the category names and list of samples that will be used in calculating
        # the distances for each of those categories. For example when
        # resolution_type == 'host_only', categories will be Porites, Millepora, Pocillopora
        # When host_secondary_seq == e.g. Porites_first, Porites_second etc.
        self.categories, self.sample_lists = self._generate_distance_categories()
        self.dist_methods = ['braycurtis', 'unifrac']
    
    def compute_distances(self):
        for dist_method in self.dist_methods:
            for distance_cat, sample_list in zip(self.categories, self.sample_lists):
                indi_dist = IndiDistanceAnalysis(
                    distance_cat=distance_cat, sample_list=sample_list,
                    dist_method=dist_method, qc_dir=self.qc_dir,
                    output_dir=self.output_dir, cache_dir=self.cache_dir, resolution_type=self.resolution_type,
                    info_df=self.info_df
                )
                indi_dist.do_analysis()

    def plot_computed_distances(self):
        """Plot up the PCoA outputs from the computed distances.
        To start with we will work with a simple 2x3 plot
        Options for color by are 'none', 'island', 'maj_seq', 'intra_diversity'"""
        for color_by in ['island', 'tech_reps', None,  'maj_seq', 'intra_diversity']:
            dist_plotter = DistancePlotter(
                resolution_type=self.resolution_type, output_dir=self.output_dir,
                info_df=self.info_df, qc_dir=self.qc_dir, cache_dir=self.cache_dir, color_by=color_by,
                figure_dir=self.figure_dir, plot_unifrac=True
            )
            dist_plotter.plot()

    def _generate_distance_categories(self):
        if self.resolution_type == 'host_only':
            sample_lists = []
            categories = ['Pocillopora', 'Millepora', 'Porites']
            for category in categories:
                sample_lists.append(list(self.info_df[self.info_df['most_abund_coral_genus'] == category].index))
            return categories, sample_lists
        elif self.resolution_type in ['host_primary_seq_only', 'host_primary_seq_only_5_samples_at_least', 'host_primary_seq_only_50_samples_at_least']:
            sample_lists = []
            categories = ['Pocillopora', 'Millepora', 'Porites']
            for category in categories:
                # First for each species we need to find what the primary sequence is
                cat_df = self.info_df[self.info_df['most_abund_coral_genus'] == category]
                primary_seq = sorted([_ for _ in cat_df["most_abund_seq_of_coral_genus"].value_counts().items()], key=lambda x: x[1], reverse=True)[0][0]
                sample_lists.append(list(cat_df[cat_df["most_abund_seq_of_coral_genus"] == primary_seq].index))
            return categories, sample_lists
        else:
            raise NotImplementedError

class DistancePlotter:
    def __init__(self, resolution_type, output_dir, info_df, figure_dir, qc_dir, cache_dir, color_by=None, plot_unifrac=True):
        if plot_unifrac:
            self.fig = plt.figure(figsize=(10, 7))
        else:
            self.fig = plt.figure(figsize=(12, 4))
        self.resolution_type = resolution_type
        self.output_dir = output_dir
        self.color_by = color_by
        self.cache_dir = cache_dir
        self.info_df = info_df
        self.figure_dir = figure_dir
        self.qc_dir = qc_dir
        self.color_dict = self._init_color_dict()
        self.plot_unifrac = plot_unifrac

    def _init_color_dict(self):
        color_list = get_colour_list()

        if self.color_by == 'island':
            island_c_dict = {island: color for island, color in zip(self.info_df['island'].unique(), color_list[:len(self.info_df['island'].unique())])}
            return {sample_name: island_c_dict[self.info_df.at[sample_name, 'island']] for sample_name in self.info_df.index}
        elif self.color_by is None:
            return {sample_name: 'black' for sample_name in self.info_df.index}
        elif self.color_by == 'maj_seq':
            maj_seq_c_dict = {maj_seq: color for maj_seq, color in
                             zip(self.info_df['most_abund_seq_of_coral_genus'].unique(), color_list[:len(self.info_df['most_abund_seq_of_coral_genus'].unique())])}
            return {sample_name: maj_seq_c_dict[self.info_df.at[sample_name, 'most_abund_seq_of_coral_genus']] for sample_name in
                    self.info_df.index}
        elif self.color_by == 'intra_diversity':
            # If we are colouring by intra div then we will use the return a dict of
            # sample_name to number of unique sequences
            return self._create_intra_diversity_dict()
        elif self.color_by == 'tech_reps':
            # We are implementing this coloring as a sanity check for the unifrac
            # find all the tech reps and give them the same color
            # Keep the same color dict
            if os.path.isfile(os.path.join(self.cache_dir, 'tech_reps_color_dict.p.bz')):
                return compress_pickle.load(os.path.join(self.cache_dir, 'tech_reps_color_dict.p.bz'))
            tech_rep_base = set()
            color_dict = {}
            for sample_name in self.info_df.index:
                if sample_name[-2] == '_':
                    # Then this is a tech rep sample name
                    tech_rep_base.add('_'.join(sample_name.split('_')[:-1]))
            # Just do the first 10 tech reps
            tech_rep_color_dict = {tech_rep: color for tech_rep, color in zip(list(tech_rep_base)[:25], color_list[:25])}
            for sample_name in self.info_df.index:
                if sample_name[-2] == '_':
                    try:
                        # Then this is tech rep
                        color_dict[sample_name] = tech_rep_color_dict['_'.join(sample_name.split('_')[:-1])]
                        foo = 'asdf'
                    except KeyError:
                        color_dict[sample_name] = (1,1,1,0)
                else:
                    color_dict[sample_name] = (1,1,1,0)
            compress_pickle.dump(color_dict, os.path.join(self.cache_dir, 'tech_reps_color_dict.p.bz'))
            return color_dict
        else:
            raise NotImplementedError


    def _create_intra_diversity_dict(self):
        diversity_dict = {}
        print('Creating intra diversity dict')
        for sample_name in self.info_df.index:
            sys.stdout.write(f'\r{sample_name}')
            sample_qc_dir = os.path.join(self.qc_dir, sample_name)
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
        if self.resolution_type in ['host_only', 'host_primary_seq_only', 'host_primary_seq_only_5_samples_at_least', 'host_primary_seq_only_50_samples_at_least']:
            if self.plot_unifrac:
                self.gs = gridspec.GridSpec(3, 3, height_ratios=[0.2, 1, 1])
            else:
                self.gs = gridspec.GridSpec(2, 3, height_ratios=[0.2,1])
            if self.plot_unifrac:
                dist_methods = ['braycurtis', 'unifrac']
            else:
                dist_methods = ['braycurtis']
            for i, dist_method in enumerate(dist_methods):
                for j, cat in enumerate(['Pocillopora', 'Millepora', 'Porites']):
                    ax = self.fig.add_subplot(self.gs[i + 1, j])
                    try:
                        indi_dist_plotter = IndiDistancePlotter(
                            fig=self.fig, ax=ax, category=cat, dist_method=dist_method,
                            path_to_csv=os.path.join(self.output_dir, f'{self.resolution_type}_{cat}_{dist_method}.csv'),
                            color_dict=self.color_dict, color_by=self.color_by, info_df=self.info_df, qc_dir=self.qc_dir)
                        indi_dist_plotter.plot()
                    except FileNotFoundError as e:
                        # Then the .csv hasn't been output yet
                        print(e)
                        continue
        else:
            raise NotImplementedError
        if self.plot_unifrac:
            with_unifrac_string = '_w_unifrac'
        else:
            with_unifrac_string = ''
        title_ax = self.fig.add_subplot(self.gs[0,:])
        title_ax.text(x=0.5, y=0.5, s=f'{self.resolution_type}_{self.color_by}{with_unifrac_string}', ha='center', va='center')
        remove_axes_but_allow_labels(title_ax)
        plt.tight_layout()
        print(f'saving {self.resolution_type}_{self.color_by}_pcoa{with_unifrac_string}.png')
        plt.savefig(os.path.join(self.figure_dir, f'{self.resolution_type}_{self.color_by}_pcoa{with_unifrac_string}.png'), dpi=1200)
        print(f'saving {self.resolution_type}_{self.color_by}_pcoa{with_unifrac_string}.svg')
        plt.savefig(os.path.join(self.figure_dir, f'{self.resolution_type}_{self.color_by}_pcoa{with_unifrac_string}.svg'))

class IndiDistancePlotter:
    """Class for plotting an individual axis"""
    def __init__(self, fig, ax, category, dist_method, path_to_csv, color_dict, color_by, info_df, qc_dir):
        self.fig = fig
        self.ax = ax
        self.category = category
        self.dist_method = dist_method
        self.pcoa_df = pd.read_csv(path_to_csv)
        self.color_dict = color_dict
        self.color_by = color_by
        self.info_df = info_df
        self.qc_dir = qc_dir

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


class IndiDistanceAnalysis:
    def __init__(self, distance_cat, sample_list, dist_method, qc_dir, output_dir, cache_dir, resolution_type, info_df):
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
        self.resolution_type = resolution_type
        self.info_df = info_df
        # A temp directory where we can write out the unaligned fastas
        # and aligned fasta and any other intermediate files
        self.temp_dir = os.path.join(os.path.dirname(self.qc_dir), 'temp')
        self.cache_dir = cache_dir
        os.makedirs(self.temp_dir, exist_ok=True)
        # In the initial 18s we were using a normalisation of 10000 seqs for the braycurtis
        # but 1000 for the unifrac due to the amount of time it was taking to create
        # the trees. We will start here with 1000 seqs for the unifrac and see
        # how long it takes to make the tree. If its not too bad then we can up the number
        if self.dist_method == 'unifrac':
            self.num_seqs_to_normalise_to = 1000
            self.abundance_df = None
        else:
            self.num_seqs_to_normalise_to = 10000
            self.abundance_dict = None


        # Generic variables shared between braycurtis and unifrac
        self.dist_out_path = os.path.join(self.output_dir, f'{self.resolution_type}_{self.category}_{self.dist_method}.dist')
        self.pcoa_out_path = os.path.join(self.output_dir, f'{self.resolution_type}_{self.category}_{self.dist_method}.csv')
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

        # Then we need to screen any sequences we use to make sure that they are found in at least
        # 5 samples. We need to build a dict to allow us to do this.
        if self.resolution_type in ['host_primary_seq_only_5_samples_at_least', 'host_primary_seq_only_50_samples_at_least']:
            seq_to_total_abund_dict = self._get_seq_to_total_abund_dict()
            if self.resolution_type == 'host_primary_seq_only_5_samples_at_least':
                threshold_set = {k for k, v in seq_to_total_abund_dict.items() if v > 4}
            elif self.resolution_type == 'host_primary_seq_only_50_samples_at_least':
                threshold_set = {k for k, v in seq_to_total_abund_dict.items() if v > 49}
            else:
                raise NotImplementedError
        dict_to_create_df_from = {}
        print('Creating abundance df')

        # NB there is one sample that only had one sequence in the consolidated_host_seqs_abund_dict
        # and this sample is causing errors. We will check for empty conolidated dict,
        # add the sample to a list and then remove it from the self.samples list
        samples_to_remove_list = []
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
            if self.resolution_type in ['host_primary_seq_only_5_samples_at_least',
                                        'host_primary_seq_only_50_samples_at_least']:
                # screen out those sequences that are not found in 5 or more samples
                consolidated_host_seqs_abund_dict = {k: v for k, v in consolidated_host_seqs_abund_dict.items() if k in threshold_set}
            tot = sum(consolidated_host_seqs_abund_dict.values())
            consolidated_host_seqs_abund_dict = {k: v/tot for k, v in consolidated_host_seqs_abund_dict.items()}

            # To prevent downstream problems we will insist on there always being at least
            # three sequences. Later we can put this to much higher to look at the effect
            if len(consolidated_host_seqs_abund_dict.keys()) < 3:
                samples_to_remove_list.append(sample_name)
                continue
            for sequence, rel_abund in consolidated_host_seqs_abund_dict.items():
                normalised_abund = int(rel_abund*self.num_seqs_to_normalise_to)
                if normalised_abund:
                    temp_sample_dict[sequence] = normalised_abund
            dict_to_create_df_from[sample_name] = temp_sample_dict
        for sample_name in samples_to_remove_list:
            self.samples.remove(sample_name)

        # It may also be very helpful to look at the distribution of the number of minor sequences
        # a given sample has
        hist_list = [len(sub_dict.keys()) for sub_dict in dict_to_create_df_from.values()]
        print('making and writing histogram of sequence diversity')
        plt.hist(hist_list, bins=30)
        plt.savefig(os.path.join(self.output_dir, f'seq_diversity_hist_3_cutoff_{self.category}_{self.dist_method}.png'))
        plt.close()

        df = pd.DataFrame.from_dict(dict_to_create_df_from, orient='index')
        df[pd.isna(df)] = 0
        print('\ndf creation complete\n')
        return df

    def _get_seq_to_total_abund_dict(self):
        if os.path.isfile(os.path.join(self.cache_dir, 'seq_to_total_abund_dict.p.bz')):
            seq_to_total_abund_dict = compress_pickle.load(os.path.join(self.cache_dir, 'seq_to_total_abund_dict.p.bz'))
        else:
            seq_to_total_abund_dict = defaultdict(int)
            for sample_name in self.info_df.index:
                sample_qc_dir = os.path.join(self.qc_dir, sample_name)
                consolidated_host_seqs_abund_dict = compress_pickle.load(
                    os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
                for seq_name in consolidated_host_seqs_abund_dict.keys():
                    seq_to_total_abund_dict[seq_name] += 1
            compress_pickle.dump(seq_to_total_abund_dict, os.path.join(self.cache_dir, 'seq_to_total_abund_dict.p.bz'))
        return seq_to_total_abund_dict

    def _create_abundance_dicts(self):
        # If resolution_type == host_primary_seq_only_5_samples_at_least
        # Then we need to screen any sequences we use to make sure that they are found in at least
        # 5 samples. We need to build a dict to allow us to do this.
        if self.resolution_type in ['host_primary_seq_only_5_samples_at_least', 'host_primary_seq_only_50_samples_at_least']:
            seq_to_total_abund_dict = self._get_seq_to_total_abund_dict()
            if self.resolution_type == 'host_primary_seq_only_5_samples_at_least':
                threshold_set = {k for k, v in seq_to_total_abund_dict.items() if v > 4}
            elif self.resolution_type == 'host_primary_seq_only_50_samples_at_least':
                threshold_set = {k for k, v in seq_to_total_abund_dict.items() if v > 49}
            else:
                raise NotImplementedError

        abundance_dict = {}
        print('Creating abundance df')
        samples_to_remove_list = []
        for sample_name in self.samples:
            sys.stdout.write(f'\r{sample_name}')
            temp_sample_dict = {}
            sample_qc_dir = os.path.join(self.qc_dir, sample_name)
            consolidated_host_seqs_abund_dict = compress_pickle.load(
                os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))

            # We need to remove the most abundant sequence from the equation
            most_abund_sequence = max(consolidated_host_seqs_abund_dict.keys(),
                                      key=(lambda key: consolidated_host_seqs_abund_dict[key]))
            # remove the most abund seq
            del consolidated_host_seqs_abund_dict[most_abund_sequence]
            # renormalise
            if self.resolution_type in ['host_primary_seq_only_5_samples_at_least',
                                        'host_primary_seq_only_50_samples_at_least']:
                # screen out those sequences that are not found in 5 or more samples
                consolidated_host_seqs_abund_dict = {
                    k: v for k, v in consolidated_host_seqs_abund_dict.items() if k in threshold_set}
            tot = sum(consolidated_host_seqs_abund_dict.values())
            consolidated_host_seqs_abund_dict = {k: v / tot for k, v in consolidated_host_seqs_abund_dict.items()}

            # To prevent downstream problems we will insist on there always being at least
            # three sequences. Later we can put this to much higher to look at the effect
            if len(consolidated_host_seqs_abund_dict.keys()) < 3:
                samples_to_remove_list.append(sample_name)
                continue

            for sequence, rel_abund in consolidated_host_seqs_abund_dict.items():
                normalised_abund = int(rel_abund * self.num_seqs_to_normalise_to)
                if normalised_abund:
                    temp_sample_dict[sequence] = normalised_abund
            abundance_dict[sample_name] = temp_sample_dict
        for sample_name in samples_to_remove_list:
            self.samples.remove(sample_name)

        # It may also be very helpful to look at the distribution of the number of minor sequences
        # a given sample has
        hist_list = [len(sub_dict.keys()) for sub_dict in abundance_dict.values()]
        print('making and writing histogram of sequence diversity')
        plt.hist(hist_list, bins=30)
        plt.savefig(
            os.path.join(self.output_dir, f'seq_diversity_hist_3_cutoff_{self.category}_{self.dist_method}.png'))
        plt.close()

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
        self.abundance_dict = self._create_abundance_dicts()
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
        compress_pickle.dump(distance_dict, os.path.join(self.cache_dir, f'{self.resolution_type}_{self.category}_distance_dict.p.bz'))
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
        self.abundance_df = self._create_abundance_df()
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
        if os.path.isfile(os.path.join(self.cache_dir, f'{hash_of_aligned_fasta}.treefile')):
            # Then we have already computed the tree and we can use this tree
            self.tree_path = os.path.join(self.cache_dir, f'{hash_of_aligned_fasta}.treefile')
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
            os.rename(self.tree_path, os.path.join(self.cache_dir, f'{hash_of_aligned_fasta}.treefile'))
            self.tree_path = os.path.join(self.cache_dir, f'{hash_of_aligned_fasta}.treefile')

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
    dist = EighteenSDistance()
    # Options for resolution_type are:
    # 'host_only' = 0 filtering of the seqs in the samples and they are only separated by
    # the gentically identified majority genus.

    # 'host_primary_seq_only' = As above, but only the samples hosting the primary sequence are used

    # host_primary_seq_only_5_samples_at_least = As above but the sequences are only used if they
    # have been found in 5 or more samples.

    # host_primary_seq_only_50_samples_at_least = As above but the sequences are only used if they
    # have been found in 50 or more samples.

    # TODO we should seperate out the samples_at_least_threshold
    # and the fact that we are exluding secondary seqs

    # TODO, There are a couple of extra things that I think we should explore.
    # Firstly, I want to do screening for only using sequences found in at least 5 samples or above.
    # Firstly, I want to screening for the number of sequences in a sample

    dist.make_and_plot_dist_and_pcoa(resolution_type='host_primary_seq_only_100_samples_at_least')