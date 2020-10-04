"""
18s_processing must have been run before this script.

This script is concerned with plotting stacked bar plots of 18S data.
It is used to produce the plots for the 18S release 1 that correspond and support
the output release tables.

The script 18s_processing.py takes care of all of the processing of the samples.
From doing that processing we end up with a directory called seq_qc that has a directory for
each readset in it. In each of these directories we have three dictionaries pickled out as well as
a fasta and names file. The fasta gives us all of the sequences in the readset after mothur processing
(i.e. no taxonomic exclusion) and the names file gives us the abundances of those samples. Have a look
at the 18s_processing.py script to get exactaly what the three dicts are but essentially, one is 
all sequences taxonomically annotated, one is just Symbiodiniaceae sequences and one is just the 
coral sequnces.

I was worried about misannoation of sequences but I don't think we have to worry about this so much
becauase the only contaminating sequnces should be from the other corals that were sampled
i.e. the Millepora, Pocillopora or porites and these should be really quite distinct from each
other that the differences should be obvious. 
"""
import os
import sys
import pandas as pd
from collections import defaultdict, Counter
from multiprocessing import Pool
import subprocess
import compress_pickle
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import time
import numpy as np
import operator
import matplotlib.gridspec as gridspec
from base_18s import EighteenSBase

class EighteenSAnalysis(EighteenSBase):
    def __init__(self):
        super().__init__()
        # This will hold the additional variables used only in this analysis class
        self.island_site_dict = self._determine_sites_and_island()
        self.islands = sorted(list(self.island_site_dict.keys()))
        self.host_species = ["Pocillopora", "Porites", "Millepora"]

    def _determine_sites_and_island(self):
        # For each readset pull out the island and site from the fastq_info_df
        # Dict that will have island names as key and a list of sites as value
        island_site_dict = defaultdict(set)
        for readset in self.coral_readsets:
            sample_id = self.fastq_info_df.at[readset, 'sample-id']
            island = self.sample_provenance_df.at[sample_id, 'ISLAND#']
            site = self.sample_provenance_df.at[sample_id, 'SITE#']
            island_site_dict[island].add(site)
        return island_site_dict

    def do_stacked_bar_plots(self, plot_type, in_sample_cutoff=None, norm_abund=None, norm_method=None, labels='colony'):
        sbp = StackedBarPlotter(
            plot_type=plot_type, islands=self.islands,
            island_site_dict=self.island_site_dict, host_species=self.host_species, 
            fig_output_dir=self.fig_output_dir, qc_dir=self.qc_dir, fastq_info_df=self.fastq_info_df, sample_provenance_df=self.sample_provenance_df, cache_dir=self.cache_dir,
            in_sample_cutoff=in_sample_cutoff, norm_abund=norm_abund, norm_method=norm_method, label_type=labels)
        sbp.plot()


class StackedBarPlotter:
    """
    This class will be responsible for plotting the stacked bar plots
    that will be arranged in a large 18 x 18 matrice of the islands sites and coral species.
    We will use this plot to get an overview of the sequencing results.
    """
    def __init__(
            self, plot_type, islands, island_site_dict, host_species,
            fig_output_dir, qc_dir, fastq_info_df, sample_provenance_df, 
            cache_dir, in_sample_cutoff, norm_abund, norm_method, label_type):
        # We will use this plot type variable to change between the different types of plots being produced
        # We will start with 'all_taxa' that can be all of the taxa in a single plot
        # We will now do 'all_coral_genus'. This will be a plot of only the coral sequences and to a genus
        # level resolution. I.e. only the three colour and then a light grey.
        # We will also 'all_coral_sequence'. This will colour the sequences according to overall abundance
        # of the sequences. This will show us whether it is exactly the same sequence within the genera
        # that have been amplified. This will require a bit more work as we will need to do sequence consolidation
        # (reuse the code from SymPortal DataSetSampleSequencePM ref seq consolidation. We will also need to generate
        # a suite of colors to use.
        # After that we should do 'minor_coral'. This will be only the genus of coral that the most abundant
        # sequence belongs to, but the most abundant sequnce wil not have been included.
        self.plot_type = plot_type
        self.fastq_info_df = fastq_info_df
        self.sample_provenance_df = sample_provenance_df
        self.qc_dir = qc_dir
        self.cache_dir = cache_dir
        self.islands = islands
        # This will be updated every island, and passed in so that we know
        # how many sites we are working with for a given indi plot as this effects
        # how big the label needs to be.
        self.site_list = None
        self.island_site_dict = island_site_dict
        self.host_species = host_species
        self.fig_output_dir = fig_output_dir
        self.in_sample_cutoff = in_sample_cutoff
        self.norm_abund = norm_abund
        self.norm_method = norm_method
        self.label_type = label_type
        self.plotting_categories, self.color_dict = self._init_color_dict()
        if self.plot_type in ['all_coral_sequence','minor_coral_sequence']:
            self.ordered_seq_name_list = self._get_ordered_seq_name_list()
        # Setup the plot
        self.fig = plt.figure(figsize=(14, 10))
        
        self.gs = gridspec.GridSpec(7, 6, figure=self.fig, height_ratios=([0.2, 1, 1, 1, 1, 1, 1]))
        self._plot_species_headers()
        self._do_legend()

    def _plot_species_headers(self):
        for i in range(6):
            sub_gs = self.gs[0, i].subgridspec(1, 3)
            for j, lab in enumerate(['POC', 'POR', 'MIL']):
                ax = self.fig.add_subplot(sub_gs[j])
                ax.text(x=0.5, y=0.5, s=lab, ha='center', va='center')
                self._remove_axes_but_allow_labels(ax)

    def _get_ordered_seq_name_list(self):
        coral_seq_abund_dict = compress_pickle.load(os.path.join(self.cache_dir, 'final_consolidated_host_seqs_rel_abundance_dict.p.bz'))
        return [tup[0] for tup in sorted([(k, v) for k, v in coral_seq_abund_dict.items()], key=lambda x:x[1], reverse=True)]
    
    def _init_color_dict(self):
        if self.plot_type == 'all_taxa':
            col_dict = {'Porites': '#FFFF00', 'Pocillopora': '#87CEFA', 'Millepora': '#FF6347',
                            'other_coral': '#C0C0C0', 'Symbiodiniaceae': '#00FF00', 'other_taxa': '#696969',
                        'not_annotated': '#282828'}
            return ['Pocillopora', 'Porites', 'Millepora', 'other_coral',
                    'Symbiodiniaceae', 'other_taxa', 'not_annotated'], col_dict
        elif self.plot_type == 'all_coral_genus':
            col_dict = {'Porites': '#FFFF00', 'Pocillopora': '#87CEFA', 'Millepora': '#FF6347',
                            'other_coral': '#C0C0C0'}
            return ['Pocillopora', 'Porites', 'Millepora', 'other_coral'], col_dict
        elif self.plot_type in ['all_coral_sequence', 'minor_coral_sequence']:
            col_dict = compress_pickle.load(os.path.join(self.cache_dir, 'all_coral_sequence_color_dict.p.bz'))
            return None, col_dict
        else:
            raise NotImplementedError()

    def plot(self):
        # we will go in order of: for island, for site, for species
        for island in self.islands:
            self.site_list = sorted(list(self.island_site_dict[island]))
            # Make the subgridspec here
            ax_row = int(self.islands.index(island)/6) + 1
            ax_col = self.islands.index(island)%6
            # Number of rows will be the number of sites + 1 for the title
            # Number of columns will be constant and the number of hosts
            sub_gs = self.gs[ax_row, ax_col].subgridspec(len(self.site_list) + 1, 3)
            # Put the island name in the title plot.
            title_ax = self.fig.add_subplot(sub_gs[0, :])
            self._do_island_title_plot(island=island, ax=title_ax)
            for site in self.site_list:
                for species in self.host_species:
                    # The last part to the row index is to incorporate the plot that we will do the naming in.
                    # The row will be the site
                    ax_row_index = self.site_list.index(site) + 1
                    # The col will always be the species
                    ax_col_index = self.host_species.index(species)
                    # In here we can do the actual plotting
                    ax = self.fig.add_subplot(sub_gs[ax_row_index, ax_col_index])
                    # Do the plotting for a given island, site, species set of samples
                    sbip = StackedBarIndiPlot(parent=self, ax=ax, island=island, 
                    site=site, species=species)
                    if sbip.readsets:
                        sbip.do_plotting()
                    else:
                        self._remove_axes_but_allow_labels(ax)

        self.fig.suptitle(f'18S {self.plot_type}', fontsize=16)
        if self.plot_type == 'minor_coral_sequence':
            if self.in_sample_cutoff and self.norm_abund:
                svg_path = os.path.join(self.fig_output_dir,
                                        f'stacked_bar_18S_{self.plot_type}_in_sample_cutoff_{self.in_sample_cutoff}_norm_{self.norm_method}_{self.norm_abund}_{self.label_type}.svg')
                png_path = os.path.join(self.fig_output_dir,
                                        f'stacked_bar_18S_{self.plot_type}_in_sample_cutoff_{self.in_sample_cutoff}_norm_{self.norm_method}_{self.norm_abund}_{self.label_type}.png')
            elif self.in_sample_cutoff:
                svg_path = os.path.join(self.fig_output_dir,
                                        f'stacked_bar_18S_{self.plot_type}_in_sample_cutoff_{self.in_sample_cutoff}_{self.label_type}.svg')
                png_path = os.path.join(self.fig_output_dir,
                                        f'stacked_bar_18S_{self.plot_type}_in_sample_cutoff_{self.in_sample_cutoff}_{self.label_type}.png')
            elif self.norm_abund:
                svg_path = os.path.join(self.fig_output_dir,
                                        f'stacked_bar_18S_{self.plot_type}_norm_{self.norm_method}_{self.norm_abund}_{self.label_type}.svg')
                png_path = os.path.join(self.fig_output_dir,
                                        f'stacked_bar_18S_{self.plot_type}_norm_{self.norm_method}_{self.norm_abund}_{self.label_type}.png')
            else:
                # if no normalisation and no sample cutoff used
                svg_path = os.path.join(self.fig_output_dir, f'stacked_bar_18S_{self.plot_type}_{self.label_type}.svg')
                png_path = os.path.join(self.fig_output_dir, f'stacked_bar_18S_{self.plot_type}_{self.label_type}.png')
        else:
            svg_path = os.path.join(self.fig_output_dir, f'stacked_bar_18S_{self.plot_type}_{self.label_type}.svg')
            png_path = os.path.join(self.fig_output_dir, f'stacked_bar_18S_{self.plot_type}_{self.label_type}.png')

        print(f'\nWriting .svg to {svg_path}')
        plt.savefig(svg_path)
        print(f'Writing .png to {png_path}')
        plt.savefig(png_path, dpi=1200)

    def _do_legend(self):
        """Plot a legend at the bottom of the figure in the remaining
        space. For the all_coral_sequence and minor_coral_sequence plots
        This should be a lgend of the most abundant sequences. For the
        all_taxa plot and the all_coral_genus plot this should be the
        categories of the color dictionaries"""
        # We will use the remaing space of the figure to make the axis
        ax = self.fig.add_subplot(self.gs[6, 2:])

        if self.plot_type in ['all_coral_sequence', 'minor_coral_sequence']:
            # To start with let's attempt to have the top 9 sequence
            # We need to pass in artists, i.e. pathches, i.e. rectangles and labels
            # We will also want to output the top 24 seqs as a fasta
            fasta = []
            import math
            ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(math.floor(n/10)%10!=1)*(n%10<4)*n%10::4])
            legend_tups = []
            num_top_seqs = 24
            for i, top_seq in enumerate(self.ordered_seq_name_list[:num_top_seqs]):
                legend_tups.append((
                    Rectangle((0 - 0.5, 0), 1, 0.5, color=self.color_dict[top_seq]),
                    f'{ordinal(i+1)} most abund'
                    ))
                fasta.extend([f'>{i+1}', top_seq])
            print(f'Writing out top_{num_top_seqs}_seqs fasta')
            with open(os.path.join(self.fig_output_dir, f'top_{num_top_seqs}_seqs.fasta'), 'w') as f:
                for line in fasta:
                    f.write(f'{line}\n')
            ax.legend([tup[0] for tup in legend_tups], [tup[1] for tup in legend_tups], loc=10, ncol=5,
                      fontsize='x-small')
        elif self.plot_type in ['all_taxa', 'all_coral_genus']:
            legend_tups = []
            for plotting_cat in self.plotting_categories:
                legend_tups.append((
                    Rectangle((0 - 0.5, 0), 1, 0.5, color=self.color_dict[plotting_cat]),
                    plotting_cat
                ))
            if self.plot_type == 'all_taxa':
                ax.legend([tup[0] for tup in legend_tups], [tup[1] for tup in legend_tups], loc=10, ncol=4,
                          fontsize='x-small')
            else:
                ax.legend([tup[0] for tup in legend_tups], [tup[1] for tup in legend_tups], loc=10, ncol=6,
                          fontsize='x-small')

        self._remove_axes_but_allow_labels(ax)
        foo = 'bar'

    def _do_island_title_plot(self, island, ax):
        # The coordinates of the title plot will always be the same
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.text(x=0.5, y=0.5, s=f'Island {island}', ha='center', va='center')
        # ax.set_frame_on(False)
        ax.set_yticks([])
        ax.set_xticks([])

    @staticmethod
    def _remove_axes_but_allow_labels(ax):
        ax.set_frame_on(False)
        ax.set_yticks([])
        ax.set_xticks([])


class StackedBarIndiPlot():
    def __init__(self, parent, ax, island, site, species):
        #TODO see if we can get rid of this base reference. Seems messy
        
        self.parent = parent
        self.ax = ax
        self.island = island
        self.site = site
        self.species = species
        # As well as the sample names, get the sample ids (i.e. C001).
        # We will plot these very small underneath the plots so that
        # we can see exactly which individual we are working with
        self.sample_ids, self.indi_annotations, self.readsets = self._get_readset_list()
        # self.readsets, self.individuals = self._get_readset_list()
        self.patches_list = []
        self.ind = 0
        self.color_list = []
        self.num_smp_in_this_subplot = len(self.readsets)
        # We need to create an abundance dictionary for each of the samples
        # We will use the pickled out dictionaries to do this
        # TODO it might be a good idea to do this for all samples at once so that it can be pickled out
        # rather than for a variable collection of samples at one time
        if self.parent.plot_type in ['all_taxa', 'all_coral_genus']:
            # If doing these plots then we are working with set categories
            # and we can work with a DataFrame
            self.abundance_df = self._make_abundance_df()
        elif self.parent.plot_type in ['all_coral_sequence', 'minor_coral_sequence']:
            # Then we are working with unknown sequences and we need to work with a dictionary
            # That we will then plot in the order of self.parent.ordered_seq_name_list
            # This is a dict of dicts where first key is sample name,
            # second key is consolidated seq_name, and value is rel abund in sample
            self.abundance_dicts = self._make_abundance_dicts()
        else:
            raise NotImplementedError
        
    def _get_readset_list(self):
        # We want to return a list of the readsets, sample_ids and individual numbers
        # that match the island, site and species specs.
        # We want these to be ordered the same and be ordered according to the individual number i.e. C001, C002 etc.
        sample_id_matches = self.parent.sample_provenance_df[
            (self.parent.sample_provenance_df['ISLAND#'] == self.island) &
            (self.parent.sample_provenance_df['SITE#'] == self.site) &
            (self.parent.sample_provenance_df['Sample Material label, organismal system level, taxonomic, nominal'].str.contains(self.species, na=False)) & 
            (self.parent.sample_provenance_df['SAMPLE ENVIRONMENT, short'] == 'C-CORAL')
        ].index.values.tolist()
        # The above contains samples for things like photos etc. We need to filter it down further to only the sample-ids
        # that are in the fastq_info_df. We have already filtered for only coral samples above so we shuoldn't pick up
        # the OA sample-ids in the fastq_info_df.
        sample_id_matches = [_ for _ in sample_id_matches if _ in self.parent.fastq_info_df['sample-id'].values.tolist()]
        # At this point we have the sample_ids in a random order
        # Now get the sample individual numbers that match the sample_ids in the same order
        indi_list = []
        for sample_id in sample_id_matches:
            indi_list.append(self.parent.sample_provenance_df.at[sample_id, 'COLONY# (C000) FISH# (F000) MACROALGAE# (MA00)'])
        # Now make indi to sample_id dict
        indi_to_sample_dict = {indi:sample_id for indi, sample_id in zip(indi_list, sample_id_matches)}
        # Now order the indi_list and re_order the sample_idslist
        indi_list = sorted(indi_list)
        sample_id_matches = [indi_to_sample_dict[_] for _ in indi_list]
        
        # Now for each indi sample_id, and therefore indi, get the list of readsets
        # Then extend the various lists with the number of readsets or the actual readset values
        readset_list = []
        sample_ids = []
        individual_ids = []
        for sample_id, indi_val in zip(sample_id_matches, indi_list):
            readsets = self.parent.fastq_info_df[
                    self.parent.fastq_info_df['sample-id'] == sample_id
                    ].index.values.tolist()
            readset_list.extend(readsets)
            sample_ids.extend([sample_id for _ in readsets])
            individual_ids.extend([indi_val for _ in readsets])
            
        
        # Here we now have the sample_id, indi and readsets in order of the indi values.
        return sample_ids, individual_ids, readset_list


    def _make_numeric(self, id):
        # remove the CO_0 part from the name
        return id.replace('C0', '')

    def _annotation_dicts_present(self, sample_name):
        sample_qc_dir = os.path.join(self.parent.qc_dir, sample_name)
        if os.path.isfile(os.path.join(sample_qc_dir, 'sample_annotation_dict.p.bz')):
            if os.path.isfile(os.path.join(sample_qc_dir, 'coral_annotation_dict.p.bz')):
                return True
        return False

    def _get_seq_to_total_abund_dict(self):
        return compress_pickle.load(os.path.join(self.cache_dir, 'seq_to_total_abund_dict.p.bz'))
        
    def _make_abundance_dicts(self):
        # If we are using an in_sample_cutoff (minimum number of samples a sequence must be found in)
        # load up the master abundance dictionary and use this to screen the sample_consolidated_abund_dict
        if self.parent.in_sample_cutoff:
            seq_to_total_abund_dict = self._get_seq_to_total_abund_dict()
            threshold_set = {k for k, v in seq_to_total_abund_dict.items() if v > self.parent.in_sample_cutoff}

        df_dict = {}
        for readset in self.readsets:
            sample_qc_dir = os.path.join(self.parent.qc_dir, readset)
            sample_consolidated_abund_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
            # When we do the minor seqs only then we will want to get rid of the most abundant seq and work with this
            if self.parent.plot_type == 'all_coral_sequence':
                df_dict[readset] = sample_consolidated_abund_dict
            elif self.parent.plot_type == 'minor_coral_sequence':
                # Remove the most abundant sequence from the dict
                del sample_consolidated_abund_dict[max(sample_consolidated_abund_dict, key=sample_consolidated_abund_dict.get)]
                # Now do normalisation if requested
                if self.parent.norm_method:
                    sample_consolidated_abund_dict = self._normalise_dict_abund(sample_consolidated_abund_dict, self.parent.norm_method)
                if self.parent.in_sample_cutoff:
                    sample_consolidated_abund_dict = {
                        k: v for k, v in sample_consolidated_abund_dict.items() if k in threshold_set}
                tot = sum(sample_consolidated_abund_dict.values())
                sample_consolidated_abund_dict = {k: v/tot for k, v in sample_consolidated_abund_dict.items()}
                df_dict[readset] = sample_consolidated_abund_dict
            else:
                raise NotImplementedError
        return df_dict

    def _normalise_dict_abund(self, dict_to_norm, method):
        # First normalise the dict so that it adds up to 1
        tot = sum(dict_to_norm.values())
        dict_to_norm = {k: v/tot for k, v in dict_to_norm.items()}

        if method == 'relative':
            # If relative then we simply want to multiply the relative abunds by the norm_abund and take the int
            return {k: int(v*self.parent.norm_abund) for k, v in dict_to_norm.items() if int(v*self.parent.norm_abund) > 0}

        elif method == 'subsample':
            # If subsample, we will produce a list using np.random.choice
            # https://docs.scipy.org/doc/numpy-1.16.0/reference/generated/numpy.random.choice.html
            # We will then convert this to a dict using counter
            # This dict will have absolute abundances. These will be converted to relaive abundances
            # outside of this script.
            seqs, probs = zip(*dict_to_norm.items())
            seqs_list = np.random.choice(seqs, self.parent.norm_abund, p=probs)
            return dict(Counter(seqs_list))

        elif method == 'num_seqs':
            # If num seqs then we want to only keep the nth most abundant sequence
            sorted_dict_keys = sorted(dict_to_norm, key=dict_to_norm.get, reverse=True)[:self.parent.norm_abund]
            return {k: dict_to_norm[k] for k in sorted_dict_keys}

    def _make_abundance_df(self):
        # Dict that we will populate and then use to make the abundance_df
        df_dict = {}
        print('making abundance df')
        for readset in self.readsets:
            sys.stdout.write(f'\r{readset}')
            sample_qc_dir = os.path.join(self.parent.qc_dir, readset)
            # make a seq_name to abundance dict from the fasta and .names pair
            seq_abund_dict = self._make_abund_dict_from_names_path(readset=readset)
            # For the all_taxa, we will go sequence by sequence through the fasta file
            fasta_path = os.path.join(sample_qc_dir, 'stability.trim.contigs.good.unique.abund.pcr.unique.fasta')
            fasta_file_as_list = EighteenSBase.decompress_read_compress(fasta_path)
            fasta_names = [line.split('\t')[0][1:] for line in fasta_file_as_list if line[0] == '>']
            
            # then load the three dictionaries
            sample_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'sample_annotation_dict.p.bz'))
            coral_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'coral_annotation_dict.p.bz'))
            sample_count_dict = {cat: 0 for cat in self.parent.plotting_categories}

            if self.parent.plot_type == 'all_taxa':
                self._log_abundances_all_taxa(sample_annotation_dict, sample_count_dict, seq_abund_dict, coral_annotation_dict, fasta_names)
            elif self.parent.plot_type == 'all_coral_genus':
                self._log_abundances_all_coral_genus(sample_annotation_dict, sample_count_dict, seq_abund_dict, coral_annotation_dict)
            else:
                raise NotImplementedError

            # Now add the collected abundances to the sample df_dict
            # Making them relative by dividing by the total of the sample_count_dd
            df_dict[readset] = [
                sample_count_dict[cat_key]/sum(sample_count_dict.values()) for
                cat_key in self.parent.plotting_categories
            ]

        # Now create the df from the df_dict
        return pd.DataFrame.from_dict(data=df_dict, orient='index', columns=self.parent.plotting_categories)

    def _log_abundances_all_taxa(
            self, sample_annotation_dict, sample_count_dict, seq_abund_dict, coral_annotation_dict, fasta_names):
        for fasta_name in fasta_names:
            try:
                annotation = sample_annotation_dict[fasta_name]
                if annotation[2] in ['Scleractinia', 'Anthoathecata']:
                    # Then this is a coral seq and we should add the count to either one of the target genera
                    # or to an other coral count
                    coral_genus = coral_annotation_dict[fasta_name]
                    if coral_genus == 'Porites':
                        key = 'Porites'
                    elif coral_genus == 'Pocillopora':
                        key = 'Pocillopora'
                    elif coral_genus == 'Millepora':
                        key = 'Millepora'
                    else:
                        key = 'other_coral'
                elif annotation[1] == 'Symbiodiniaceae':
                    key = 'Symbiodiniaceae'
                else:
                    key = 'other_taxa'
            except KeyError:
                key = 'not_annotated'
            # now log the abundance
            sample_count_dict[key] += seq_abund_dict[fasta_name]

    def _log_abundances_all_coral_genus(
            self, sample_annotation_dict, sample_count_dict, seq_abund_dict, coral_annotation_dict):
        for blasted_seq, annotation in sample_annotation_dict.items():
            if annotation[2] in ['Scleractinia', 'Anthoathecata']:
                # Then this is a coral seq and we should add the count to either one of the target genera
                # or to an other coral count
                coral_genus = coral_annotation_dict[blasted_seq]
                if coral_genus == 'Porites':
                    key = 'Porites'
                elif coral_genus == 'Pocillopora':
                    key = 'Pocillopora'
                elif coral_genus == 'Millepora':
                    key = 'Millepora'
                else:
                    key = 'other_coral'

                # now log the abundance
                sample_count_dict[key] += seq_abund_dict[blasted_seq]

    def _make_abund_dict_from_names_path(self, readset):
        name_path = os.path.join(
                self.parent.qc_dir, readset, 'stability.trim.contigs.good.unique.abund.pcr.names')
        name_file = EighteenSBase.decompress_read_compress(name_path)
        return {line.split('\t')[0]: len(line.split('\t')[1].split(',')) for line in name_file}

    def do_plotting(self):
        for readset in self.readsets:
            sys.stdout.write(f'\rPlotting sample: {self.island} {self.site} {self.species} {readset}')
            if self.parent.plot_type in ['all_taxa', 'all_coral_genus']:
                self._plot_bars_from_df(readset)
            else:
                self._plot_bars_from_dicts(readset)
            self.ind += 1
        # If <= 10 readsets, we want to set the x lim to 10 so that bars are constant width
        # Unless there are more than 10 readsets then we want to use this number
        if len(self.readsets) > 10:
            self._paint_rect_to_axes(max_num_smpls_in_subplot=len(self.readsets))
        else:
            self._paint_rect_to_axes()

    def _plot_bars_from_dicts(self, readset):
        bottom_div = 0
        readset_abund_dict = self.abundance_dicts[readset]
        order_to_plot = [seq_name for seq_name in self.parent.ordered_seq_name_list if seq_name in readset_abund_dict]
        # In order of the master consolidated seqs
        for seq_name in order_to_plot:
            seq_rel_abund = readset_abund_dict[seq_name] 
            self.patches_list.append(Rectangle((self.ind - 0.5, bottom_div), 1, seq_rel_abund, color=self.parent.color_dict[seq_name]))
            self.color_list.append(self.parent.color_dict[seq_name])
            bottom_div += seq_rel_abund

    def _plot_bars_from_df(self, readset):
        bottom_div = 0
        # In the order of the plotting categories
        for plot_cat in self.parent.plotting_categories:
            cat_rel_abund = self.abundance_df.at[readset, plot_cat]
            self.patches_list.append(Rectangle((self.ind - 0.5, bottom_div), 1, cat_rel_abund, color=self.parent.color_dict[plot_cat]))
            self.color_list.append(self.parent.color_dict[plot_cat])
            bottom_div += cat_rel_abund

    def _paint_rect_to_axes(self, max_num_smpls_in_subplot=10):
        # Makie a custom colour map
        # https://matplotlib.org/api/_as_gen/matplotlib.colors.ListedColormap.html
        this_cmap = ListedColormap(self.color_list)
        
        # Here we have a list of Rectangle patches
        # Create the PatchCollection object from the patches_list
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
        self.ax.set_ylim(0,1)
        # This is where we chose what to annotate the x axis
        if self.parent.label_type == 'colony':
            x_labels = self.indi_annotations
            fontsize=2
        elif self.parent.label_type == 'sample-id':
            if len(self.parent.site_list) > 3:
                fontsize=1
            else:
                fontsize=2
            x_labels = self.sample_ids
        elif self.parent.label_type == 'readset':
            x_labels = []
            for readset in self.readsets:
                one_el = readset.split('_')[1].replace('OSTA', '')
                two_el = readset.split('_')[3].split('.')[0]
                x_labels.append(f'{one_el}_{two_el}')
            if len(self.parent.site_list) > 4:
                fontsize=1
            else:
                fontsize=1.5
        self._remove_axes_but_allow_labels(x_labels=x_labels, fontsize=fontsize)
        self.ax.set_ylabel(self.site, fontsize=5, labelpad=0)

    def _remove_axes_but_allow_labels(self, ax=None, x_labels=None, fontsize=2):
        if ax is None:
            self.ax.set_frame_on(False)
            self.ax.set_yticks([])
            if x_labels:
                #TODO Here we want to write the labels directly onto the bars if possible
                self.ax.set_xticks([])
                for i, lab in enumerate(x_labels):
                    self.ax.text(x=i-0.5, y=0, s=lab, va='bottom', ha='left', rotation=90, fontsize=fontsize)
                    # plt.savefig('/home/humebc/draft.png', dpi=1200)
                    # self.ax.set_xticks([_-0.5 for _ in range(len(self.indi_annotations))])
                    # self.ax.set_xticklabels(self.indi_annotations, fontsize=2, rotation=-90)
                    # self.ax.xaxis.set_tick_params(which='major',  length=0, pad=0)
            else:
                self.ax.set_xticks([])
        else:
            ax.set_frame_on(False)
            ax.set_yticks([])
            if x_labels:
                ax.set_xticks([x_labels])
            else:
                ax.set_xticks([])


if __name__ == "__main__":
    # Plot types that can be provided to do_stacked_bar_plots are:
    # all_taxa
    # all_coral_genus
    # all_coral_sequences
    # minor_coral_sequence
    # You can also use the in_sample_cutoff argument to apply a cutoff that will be used when plotting
    # It is the minimum number of samples a given sample must be found in, else it will not be used.
    # This cutoff appears to have a big effect on the unifrac method in particular.
    # norm_abund is the value that sequences should be normalised to
    # norm_method is the method for normilisation.
    # norm_method can be 'hard', 'subsample' or 'num_seqs'.
    # Normalisation will only be applied when plotting minor_coral_sequences
    # The Normalistaion will be applied AFTER the most abundant sequence has been removed.
    # A hard normalisation will use the relative abundances of the sequences and do a pick without replacement
    # A rel normalisation will multiply the normalisation abund by the relative abundance of each sequences and take
    # an int() function of this result.
    # The num_seqs normalisation will only use the ith most abundant sequences where the ith degree is defined by
    # normalisation_abund
    # A 'labels' option may also be passed to the do_stacked_bar_plots function. This will denote
    # What labels are plotted on the each of the columns. By default 'colony' is used. This uses the individual
    # annotations for example C001. 'sample-id' - will use the sample-id. 'readset' will use readset.

    # Do all the plots
    # For the release 1, I think it makes sense to just provide the three plots
    # 'all_taxa', 'all_coral_genus', 'all_coral_sequence'
    for plot_type in ['all_taxa', 'all_coral_genus', 'all_coral_sequence']:
        for label in ['sample-id', 'readset']:
            EighteenSAnalysis().do_stacked_bar_plots(
                    plot_type=plot_type, in_sample_cutoff=None,
                    norm_abund=None, norm_method=None, labels=label)
    
    # for plot_type in ['all_taxa', 'all_coral_genus', 'all_coral_sequence', 'minor_coral_sequence']:
    #     if plot_type == 'minor_coral_sequence':
    #         for cutoff in [300, 500, 700, 1000]:
    #             EighteenSAnalysis().do_stacked_bar_plots(
    #                 plot_type=plot_type, in_sample_cutoff=cutoff,
    #                 norm_abund=None, norm_method=None)
    #         for norm_abund in [1, 5, 10, 20, 25, 50, 75, 100]:
    #             EighteenSAnalysis().do_stacked_bar_plots(
    #                 plot_type=plot_type, in_sample_cutoff=None,
    #                 norm_abund=norm_abund, norm_method='num_seqs')
    #     else:
    #         EighteenSAnalysis().do_stacked_bar_plots(
    #             plot_type=plot_type, in_sample_cutoff=None,
    #             norm_abund=None, norm_method=None)
