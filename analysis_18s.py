"""The script 18s_processing.py takes care of all of the processing of the samples.
From doing that processing we end up with a directory called seq_qc that has a directory for
each sample in it. In each of these directories we have three dictionaries pickled out as well as
a fasta and names file. The fasta gives us all of the sequences in the sample after mothur processing
(i.e. no taxonomic dislusion) and the names file gives us the abundances of those samples. Have a look
at the 18s_processing.py script to get exactaly what the three dicts are but essentially, one is 
all sequences taxonomically annotated, one is just Symbiodiniaceae sequences and one is just the 
coral sequnces.

In this script we will hold all methods concerned with the further processing of these samples.
We have cached out some of the utility objects like information dataframes from the 18s_processing.py
script and we will make use of these here in this script. If the caches don't exist then we will call
the 18s_processing.py script to make them. The 18s_processing.py script is fully cahce enabled and 
so will only redo parts of the processing that are required.

I think a useful thing to do will be the bar plots similar to what we did with the ITS2. This will give us
an overview of what we are working with. As with before we can use this format to plot, all coral seqs,
just the minor coral seqs and all seqs, regardless of taxa.

I was worried about misannoation of sequences but I don't think we have to worry about this so much
becauase the only contaminating sequnces should be from the other corals that were sampled
i.e. the Millepora, Pocillopora or porites and these should be really quite distinct from each
other that the differences should be obvious. 
"""
from processing_18s import MakeInfoDF
import os
import sys
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool
import subprocess
import compress_pickle
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import time
import numpy as np

# TODO later we will be able to write this as a subclass of EighteenSProcessing
# But for the time being we don't want to interfere with any of that code because 
# it is currently running to create the 18S taxa annotations.
class EighteenSAnalysis:
    def __init__(self):
        self.root_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.seq_dir = os.path.join(self.root_dir, 'seq_files')
        # The directory where the finalised post qc and post taxa screening files will be written
        self.qc_dir = os.path.join(self.root_dir, 'seq_qc')
        self.cache_dir = os.path.join(self.root_dir, 'cache')
        os.makedirs(self.cache_dir, exist_ok=True)
        self.sample_provenance_path = os.path.join(self.root_dir, "tara_samples_provenance.csv")
        self.sample_provenance_df = self._make_sample_provenance_df()
        # The main info df that we will use
        # Sample name as key, fwd and rev path to seq files, coral species
        self.info_df = self._make_info_df()

        # This will hold the additional variables used only in this analysis class
        self.island_site_dict = self._determine_sites_and_island()
        self.islands = sorted(list(self.island_site_dict.keys()))
        self.fig_output_dir = os.path.join(self.root_dir, 'figures')
        self.host_species = ["Porites", "Millepora", "Pocillopora"]

    def _determine_sites_and_island(self):
        # For each sample pull out the island and site from the info_df
        # Dict that will have island names as key and a list of sites as value
        island_site_dict = defaultdict(set)
        for sample_index in self.info_df.index:
            island = self.info_df.at[sample_index, 'island']
            site = self.info_df.at[sample_index, 'site']
            island_site_dict[island].add(site)
        return island_site_dict

    def do_stacked_bar_plots(self):
        sbp = StackedBarPlotter(
            plot_type='all_taxa', islands=self.islands, 
            island_site_dict=self.island_site_dict, host_species=self.host_species, 
            fig_output_dir=self.fig_output_dir, qc_dir=self.qc_dir)
        sbp.plot()

    def _make_info_df(self):
        try:
            return compress_pickle.load(os.path.join(self.cache_dir, 'info_df.p.bz'))
        except FileNotFoundError:
            info_df = MakeInfoDF(seq_dir=self.seq_dir, sample_provenance_df=self.sample_provenance_df).make_info_df()
            compress_pickle.dump(info_df, os.path.join(self.cache_dir, 'info_df.p.bz'))
            return info_df

    def _make_sample_provenance_df(self):
        # The SAMPLE ID is in the sequencing file name so we will be able to use this to get the latitute
        df = pd.read_csv(self.sample_provenance_path)
        df.set_index(keys='SAMPLE ID', drop=True, inplace=True)
        df.rename(columns={'EVENT latitude start (North)': 'lat', 'EVENT longitude start (East)': 'lon'}, inplace=True)
        return df

class StackedBarPlotter:
    """
    This class will be responsible for plotting the stacked bar plots
    that will be arranged in a large 16 x 16 matrice of the islands sites and coral species.
    We will use this plot to get an overview of the sequencing results.
    """
    def __init__(self, plot_type, islands, island_site_dict, host_species, fig_output_dir, qc_dir):
        # We will use this plot type variable to change between the different types of plots being produced
        # We will start with 'all_taxa' that can be all of the taxa in a single plot
        self.plot_type = plot_type
        self.qc_dir = qc_dir
        self.islands = islands
        self.island_site_dict = island_site_dict
        self.host_species = host_species
        self.fig_output_dir = fig_output_dir
        self.plotting_categories, self.color_dict = self._init_color_dict()
        # Setup the plot
        self.fig = plt.figure(figsize=(14, 10))
        self.gs = self.fig.add_gridspec(18, 18, figure=self.fig, height_ratios=[1 for _ in range(18)],
                                        width_ratios=[1 for _ in range(18)])
    
    def _init_color_dict(self):
        if self.plot_type == 'all_taxa':
            col_dict = {'Porites': '#FFFF00', 'Pocillopora': '#87CEFA', 'Millepora': '#FF6347',
                            'other_coral': '#C0C0C0', 'Symbiodiniaceae': '#00FF00', 'other_taxa': '#696969'}
            return col_dict, ['Porites', 'Pocillopora', 'Millepora', 'other_coral', 'Symbiodiniaceae', 'other_taxa'] 
        else:
            raise NotImplementedError()

    def plot(self):
        # we will go in order of: for island, for site, for species
        for island in self.islands:
            site_list = sorted(list(self.island_site_dict[island]))[:3]
            for site in site_list:
                for species in self.host_species:
                    ax_row_index = int((int(self.islands.index(island)/6)*3) + site_list.index(site))
                    ax_col_index = int(((self.islands.index(island)%6)*3) + self.host_species.index(species))
                    # In here we can do the actual plotting
                    ax = self.fig.add_subplot(self.gs[ax_row_index, ax_col_index])
                    
                    # Do the plotting for a given island, site, species set of samples
                    sbip = StackedBarIndiPlot(parent=self, ax=ax, island=island, 
                    site=site, species=species)
                    if sbip.samples:
                        sbip.do_plotting()
        
        print('Writing .svg')
        plt.savefig(os.path.join(self.fig_output_dir, f'stacked_bar_18s_{self.plot_type}.svg'))
        print('Writing .png')
        plt.savefig(os.path.join(self.fig_output_dir, f'stacked_bar_18s_{self.plot_type}.png'), dpi=1200)

class StackedBarIndiPlot:
    def __init__(self, parent, ax, island, site, species):
        self.parent = parent
        self.ax = ax
        self.island = island
        self.site = site
        self.species = species
        self.samples = self._get_sample_name_list()
        self.patches_list = []
        self.ind = 0
        self.color_list = []
        self.num_smp_in_this_subplot = len(self.samples)
        # We need to create an abundance dictionary for each of the samples
        # We will use the pickled out dictionaries to do this
        # TODO it might be a good idea to do this for all samples at once so that it can be pickled out
        # rather than for a variable collection of samples at one time
        self.abundance_df = self._make_abundance_df()

    def _get_sample_name_list(self):
        """TODO because we still have the taxa annotation running
        here, we will screen the samples to be plotted to only plot those that
        have already had the annotation analysis completed."""
        init_sample_list = self.parent.info_df[
            (self.parent.info_df['island'] == self.island) &
            (self.parent.info_df['site'] == self.site) &
            (self.parent.info_df['species'] == self.species)
        ].index.values.tolist()

        # Filter out those samples that do not have the annotation dicts already created
        return [sample_name for sample_name in init_sample_list if self._annotation_dicts_present(sample_name)]

    def _annotation_dicts_present(self, sample_name):
        sample_qc_dir = os.path.join(self.parent.qc_dir, sample_name)
        if os.path.isfile(os.path.join(sample_qc_dir, 'sample_annotation_dict.p.bz')):
            if os.path.isfile(os.path.join(sample_qc_dir, 'coral_annotation_dict.p.bz')):
                return True
        return False

    def _make_abundance_df(self):
        # Dict that we will populate and then use to make the abundance_df
        df_dict = {}
        for sample_name in self.samples:
            sample_qc_dir = os.path.join(self.parent.qc_dir, sample_name)
            # make a seq_name to abundance dict from the fasta and .names pair
            sample_abund_dict = self._make_abund_dict_from_names_path(sample_name=sample_name)

            # then load the three dictionaries
            sample_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'sample_annotation_dict.p.bz'))
            coral_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'coral_annotation_dict.p.bz'))
            # symbiodiniaceae_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'symbiodiniaceae_annotation_dict.p.bz'))

            # the order of the categories is ['Porites', 'Pocillopora', 'Millepora', 'other_coral', 'Symbiodiniaceae', 'other_taxa']
            # We will use this list order for the abundances too
            sample_count_dd = defaultdict(int)
            if self.parent.plot_type == 'all_taxa':
                self._log_abundances_all_taxa(sample_annotation_dict, sample_count_dd, sample_abund_dict, coral_annotation_dict)
            else:
                raise NotImplementedError

            # Now add the collected abundances to the sample df_dict
            # Making them relative by dividing by the total of the sample_count_dd
            df_dict[sample_name] = [sample_count_dd[cat_key]/sum(sample_count_dd.values()) for cat_key in self.parent.plotting_categories]

        # Now create the df from the df_dict
        return pd.DataFrame.from_dict(data=df_dict, orient='index', columns=self.parent.plotting_categories)

    def _log_abundances_all_taxa(self, sample_annotation_dict, sample_count_dd, sample_abund_dict, coral_annotation_dict):
        for blasted_seq, annotation in sample_annotation_dict.items():
                if annotation == 'Scleractinia_Anthoathecata':
                    # Then this is a coral seq and we should add the count to either one of the target genera
                    # or to an other coral count
                    # TODO this is where we can change our logic to according to what type of plot we are doing
                    coral_genus = coral_annotation_dict[blasted_seq]
                    key = None
                    if coral_genus == 'Porites':
                        key = 'Porites'
                    elif coral_genus == 'Pocillopora':
                        key = 'Pocillopora'
                    elif coral_genus == 'Millepora':
                        key = 'Millepora'
                    else:
                        key = 'other_coral'
                elif annotation == 'Symbiodiniaceae':
                    key = 'Symbiodiniaceae'
                else:
                    key = 'other_taxa'
                
                # now log the abundance
                sample_count_dd[key] += sample_abund_dict[blasted_seq]

    def _make_abund_dict_from_names_path(self, sample_name):
        with open(os.path.join(self.parent.qc_dir, sample_name, 'stability.trim.contigs.good.unique.abund.pcr.names'), 'r') as f:
            name_file = [line.rstrip() for line in f]
        return {line.split('\t')[0]: len(line.split('\t')[1].split(',')) for line in name_file}

    def do_plotting(self):
        for sample_to_plot in self.samples:
            sys.stdout.write(f'\rPlotting sample: {self.island} {self.site} {self.species} {sample_to_plot}')
            
            self._plot_bars(sample_to_plot)
            self.ind += 1
        self._paint_rect_to_axes()

    def _plot_bars(self, sample_to_plot):
        bottom_div = 0
        # In order that the sequences are listed in the seq_relative_abundance_df for those that are
        # present in the sample, plot a rectangle.
        for plot_cat in self.parent.plotting_categories:
            cat_rel_abund = self.abundance_df.at[sample_to_plot, plot_cat]
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
        self._remove_axes_but_allow_labels()

    def _remove_axes_but_allow_labels(self):
        self.ax.set_frame_on(False)
        self.ax.set_yticks([])
        self.ax.set_xticks([])