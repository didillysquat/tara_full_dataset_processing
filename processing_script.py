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

class ITS2Processing:
    def __init__(self):
        self.base_directory_of_sequencing_data = "/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20200116_tara_pacific_its2/"
        self.input_dir = os.path.abspath(os.path.join('.', 'input'))
        self.output_dir = os.path.abspath(os.path.join('.', 'output'))
        os.makedirs(self.input_dir, exist_ok=True)
        os.makedirs(self.output_dir, exist_ok=True)
        self.sample_provenance_path = os.path.join(self.input_dir, 'tara_samples_provenance.csv')
        self.sample_provenance_df = self._make_sample_provenance_df()
        self.cache_dir = os.path.abspath(os.path.join('.', 'cache'))
        self.info_df = self._make_info_df()
        self.sp_datasheet = self._make_sp_data_sheet()
        foo = 'bar'
        # At this point we have an information dataframe that can be used to create datasheets.
        # This information dataframe contains the information for all of the samples
        # It does not contain information for the negative controls
        # Ideally we want to run the negative samples that correspond to the coral samples
        # in the same dataloading and analysis
        # However, it appears that a given set of samples (i.e. multiple samples) have 2 negative controls.
        # As such I'm not sure how much use it will be to run the negatives with the samples.
        # The negative control sequencing files are all located in a single directory
        # and there is a mapping file that is space delimited. In this file the sequencing file name
        # i.e. TARA_AW-0000064_METAB-ITS2_CP8RR-12BA001-1
        # and then in the second and third column there is an abbreviated version
        # of the two sets of negative control sequencing files that map to the samples
        # i.e. CEB_AWN CEB_AWO. It seems that a very large number of samples
        # map to two negative controls. As such I'm guessing that the negative controls are
        # for the master mixes or something.
        # I think a good approach from here will be to first run the negatives and see what they contain.
        # If they all fail (i.e. good news) then we don't need to consider them futher.
        # If some of them have a considerable amount of ITS2 Symbiodiniaceae in them 
        # Then we can run all of the negative controls together and see where they fall out in the analysis.

        # At this point we have used the generated datasheet to load the coral samples and run them in an
        # analysis.


        self.sp_output_base_dir = os.path.join(self.input_dir, 'tara_coral_output_20200126')
        # Now create count and meta info df for the sequencing data
        self.seq_meta_info_df, self.sample_uid_to_name_dict = self._make_seq_meta_info_df()
        self.seq_absolute_abundance_df, self.seq_relative_abundance_df = self._make_seq_abundance_info_dfs()


        # Now create count and meta info df for the profile data
        self.prof_meta_info_df, self.prof_uid_to_name_dict = self._make_profile_meta_info_df()
        self.prof_absolute_abundance_df, self.prof_relative_abundance_df = self._make_prof_abundance_info_dfs()

        # Color dicts for plotting
        self.seq_color_dict, self.prof_color_dict = self._get_json_color_dicts()

    def _get_json_color_dicts(self):
        with open(os.path.join(self.sp_output_base_dir, 'html', 'color_dict_post_med.json'), 'r') as f:
            text = f.read()
            seq_color_dict = json.loads(text)
        with open(os.path.join(self.sp_output_base_dir, 'html', 'prof_color_dict.json'), 'r') as f:
            text = f.read()
            prof_color_dict = json.loads(text)
        return seq_color_dict, prof_color_dict

    def _make_profile_meta_info_df(self):
        prof_meta_abs_info_path = os.path.join(
            self.sp_output_base_dir,
            'its2_type_profiles',
            '90_20200125_DBV_2020-01-27_03-21-24.288417.profiles.absolute.meta_only.txt'
        )
        prof_meta_info_df = pd.read_csv(prof_meta_abs_info_path, sep='\t')
        # Make our uid to sample name dict
        profile_uid_to_sample_name_dict = {uid: name for uid, name in zip(prof_meta_info_df['ITS2 type profile UID'].values,
                                                                         prof_meta_info_df['ITS2 type profile'].values)}
        prof_meta_info_df.set_index('ITS2 type profile UID', inplace=True, drop=True)
        return prof_meta_info_df, profile_uid_to_sample_name_dict

    def _make_prof_abundance_info_dfs(self):
        abs_path = os.path.join(
            self.sp_output_base_dir,
            'its2_type_profiles',
            '90_20200125_DBV_2020-01-27_03-21-24.288417.profiles.absolute.abund_only.txt')

        rel_path = os.path.join(
            self.sp_output_base_dir,
            'its2_type_profiles',
            '90_20200125_DBV_2020-01-27_03-21-24.288417.profiles.relative.abund_only.txt')

        abs_df = pd.read_csv(abs_path, sep='\t')
        rel_df = pd.read_csv(rel_path, sep='\t')

        abs_df['sample_name'] = [self.sample_uid_to_name_dict[sample_uid] for sample_uid in abs_df['sample_uid'].values]
        rel_df['sample_name'] = [self.sample_uid_to_name_dict[sample_uid] for sample_uid in rel_df['sample_uid'].values]

        abs_df.set_index('sample_name', drop=True, inplace=True)
        rel_df.set_index('sample_name', drop=True, inplace=True)

        abs_df.drop(columns='sample_uid', inplace=True)
        rel_df.drop(columns='sample_uid', inplace=True)

        return abs_df, rel_df

    def _make_seq_abundance_info_dfs(self):
        abs_path = os.path.join(
            self.sp_output_base_dir,
            'post_med_seqs',
            '90_20200125_DBV_2020-01-27_03-21-24.288417.seqs.absolute.abund_only.txt')

        rel_path = os.path.join(
            self.sp_output_base_dir,
            'post_med_seqs',
            '90_20200125_DBV_2020-01-27_03-21-24.288417.seqs.relative.abund_only.txt')

        abs_df = pd.read_csv(abs_path, sep='\t')
        rel_df = pd.read_csv(rel_path, sep='\t')

        abs_df['sample_name'] = [self.sample_uid_to_name_dict[sample_uid] for sample_uid in abs_df['sample_uid'].values]
        rel_df['sample_name'] = [self.sample_uid_to_name_dict[sample_uid] for sample_uid in rel_df['sample_uid'].values]

        abs_df.set_index('sample_name', drop=True, inplace=True)
        rel_df.set_index('sample_name', drop=True, inplace=True)

        abs_df.drop(columns='sample_uid', inplace=True)
        rel_df.drop(columns='sample_uid', inplace=True)

        return abs_df, rel_df

    def _make_seq_meta_info_df(self):
        seq_meta_abs_info_path = os.path.join(
            self.sp_output_base_dir,
            'post_med_seqs',
            '90_20200125_DBV_2020-01-27_03-21-24.288417.seqs.absolute.meta_only.txt'
        )
        seq_meta_info_df = pd.read_csv(seq_meta_abs_info_path, sep='\t')
        # Make our uid to sample name dict
        sample_uid_to_sample_name_dict = {uid:name for uid, name in zip(seq_meta_info_df['sample_uid'].values, seq_meta_info_df['sample_name'].values)}
        seq_meta_info_df.drop(columns='sample_uid', inplace=True)
        seq_meta_info_df.set_index('sample_name', inplace=True, drop=True)
        return seq_meta_info_df, sample_uid_to_sample_name_dict

    def qc_overview_figure(self):
        qcof = QCPlotterCoralOnly(info_df=self.info_df, sp_datasheet=self.sp_datasheet,
                                  seq_meta_info_df=self.seq_meta_info_df, output_dir=self.output_dir)
        qcof.plot_qc_data()

    def seq_and_profile_results_figure(self, fig_type='all'):
        """I envisage that we are going to have a quite a few variations of this figure. I think we can
        work with a fig_type variable that can toggle between them. To start with lets start with 'all'.
        This will be a 3x3 plot for each island (there will be some islands with more than one site but
        we will deal with that hurdle when we get to it). It is going to have information for all of the species.
        It is going to be pretty massive so it will likely make a lot of sense to have the different plot types
        that can have only one species or some other variation, like clade colouring."""
        sap = SequenceAndProfilePlotterCoralOnly(
            info_df=self.info_df, sp_datasheet=self.sp_datasheet,
            seq_meta_info_df=self.seq_meta_info_df,
            seq_color_dict=self.seq_color_dict,
            prof_color_dict=self.prof_color_dict,
            sample_uid_to_name_dict=self.sample_uid_to_name_dict,
            seq_absolute_abundance_df=self.seq_absolute_abundance_df,
            seq_relative_abundance_df=self.seq_relative_abundance_df,
            prof_meta_info_df=self.prof_meta_info_df,
            prof_uid_to_name_dict=self.prof_uid_to_name_dict,
            prof_absolute_abundance_df=self.prof_absolute_abundance_df,
            prof_relative_abundance_df=self.prof_relative_abundance_df, output_dir=self.output_dir)


    def _make_sample_provenance_df(self):
        # The SAMPLE ID is in the sequencing file name so we will be able to use this to get the latitute
        df = pd.read_csv(self.sample_provenance_path)
        df.set_index(keys='SAMPLE ID', drop=True, inplace=True)
        df.rename(columns={'EVENT latitude start (North)':'lat', 'EVENT longitude start (East)': 'lon'}, inplace=True)
        return df

    def _make_info_df(self):
        if os.path.isfile(os.path.join(self.cache_dir, 'info_df.pickle')):
            info_df = pickle.load(open(os.path.join(self.cache_dir, 'info_df.pickle'), 'rb'))
        else:
            gicdf = GenerateInfoCollectionDF(base_dir=self.base_directory_of_sequencing_data, provenance_df=self.sample_provenance_df)
            info_df = gicdf.generate_df()

            pickle.dump(info_df, open(os.path.join(self.cache_dir, 'info_df.pickle'), 'wb'))

        return info_df

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

class SequenceAndProfilePlotterCoralOnly:
    def __init__(self, info_df, sp_datasheet, seq_meta_info_df, seq_color_dict, prof_color_dict,
            sample_uid_to_name_dict,
            seq_absolute_abundance_df,
            seq_relative_abundance_df,
            prof_meta_info_df,
            prof_uid_to_name_dict,
            prof_absolute_abundance_df,
            prof_relative_abundance_df, output_dir):
        self.info_df = info_df
        self.sp_datasheet = sp_datasheet
        self.seq_meta_info_df = seq_meta_info_df
        self.sample_uid_to_name_dict = sample_uid_to_name_dict
        self.seq_absolute_abundance_df = seq_absolute_abundance_df
        self.seq_relative_abundance_df = seq_relative_abundance_df
        self.prof_meta_info_df = prof_meta_info_df
        self.prof_uid_to_name_dict = prof_uid_to_name_dict
        self.prof_absolute_abundance_df = prof_absolute_abundance_df
        self.prof_relative_abundance_df = prof_relative_abundance_df
        self.output_dir = output_dir
        self.species = ['PORITES', 'POCILLOPORA', 'MILLEPORA']
        # Get list of islands
        # Dict that will have island names as key and a list of sites as value
        self.island_site_dict = self._determine_sites_and_island()
        self.islands = sorted(list(self.island_site_dict.keys()))
        # Get the colour dicts that we will use for plotting
        self.seq_color_dict = seq_color_dict
        self.prof_color_dict = prof_color_dict
        # Setup the plot
        self.fig = plt.figure(figsize=(14, 10))
        self.gs = self.fig.add_gridspec(18, 18, figure=self.fig, height_ratios=[1 for _ in range(18)],
                                        width_ratios=[1 for _ in range(18)])
        # The axis list that has been linearised
        # we will go in order of: for island, for site, for species
        # we will linearize the axes and go in order of for island, for site, for species
        for island in self.islands:
            site_list = sorted(list(self.island_site_dict[island]))[:3]
            for site in site_list:
                for species in self.species:
                    ax_row_index = int((int(self.islands.index(island)/6)*3) + site_list.index(site))
                    ax_col_index = int(((self.islands.index(island)%6)*3) + self.species.index(species))
                    # In here we can do the actual plotting
                    single_sp_time_start = time.time()
                    ax = self.fig.add_subplot(self.gs[ax_row_index, ax_col_index])
                    spip = SeqProfIndiPlot(parent=self, ax=ax, island=island, site=site, species=species)
                    spip.do_plotting()
                    single_sp_time_tot = time.time() - single_sp_time_start
                    print(f'The whole single took {single_sp_time_tot}')
        plt.savefig(os.path.join(self.output_dir, 'seq_and_profile_fig_all.svg'))
        plt.savefig(os.path.join(self.output_dir, 'seq_and_profile_fig_all.png'), dpi=1200)
        foo = 'bar'


    def _determine_sites_and_island(self):
        # For each sample pull out the island and site from the info_df
        # Dict that will have island names as key and a list of sites as value
        island_site_dict = defaultdict(set)
        for sample_index in self.sp_datasheet.index:
            island = self.info_df.at[sample_index, 'location']
            site = self.info_df.at[sample_index, 'site']
            island_site_dict[island].add(site)
        return island_site_dict
        
class SeqProfIndiPlot:
    def __init__(self, parent, ax, island, site, species):
        self.parent = parent
        self.ax = ax
        self.island = island
        self.site = site
        self.species = species
        self.samples = self.parent.info_df[
            (self.parent.info_df['location'] == self.island) &
            (self.parent.info_df['site'] == self.site) &
            (self.parent.info_df['spp_water'] == self.species) &
            (self.parent.info_df['coral_plankton'] == 'CORAL') &
            (self.parent.info_df['spp_water'] != 'HELIOPORA') &
            (self.parent.info_df['spp_water'] != 'PORITES_PANAMENSIS')
        ].index.values.tolist()
        foo = 'bar'
        self.patches_list = []
        self.ind = 0
        self.color_list = []
        self.num_smp_in_this_subplot = len(self.samples)

    def do_plotting(self):
        div_over_total = 0
        type_under_total = 0
        for sample_to_plot in self.samples:
            sys.stdout.write(f'\rPlotting sample: {self.island} {self.site} {self.species} {sample_to_plot}')

            # PLOT DIVs
            div_over_start = time.time()
            self._plot_div_over_type(sample_to_plot)
            div_over_stop = time.time()
            div_over_total += (div_over_stop - div_over_start)

            # PLOT type
            type_under_start = time.time()
            self._plot_type_under_div(sample_to_plot)
            type_under_stop = time.time()
            type_under_total += (type_under_stop - type_under_start)
            self.ind += 1
        paint_start = time.time()
        self._paint_rect_to_axes_div_and_type()
        paint_stop = time.time()
        paint_total = paint_stop - paint_start

        print(f'\n\ndiv_over took {div_over_total}')
        print(f'type_under took {type_under_total}')
        print(f'paint took {paint_total}')



    def _plot_div_over_type(self, sample_to_plot):
        bottom_div = 0
        # In order that the sequences are listed in the seq_relative_abundance_df for those that are
        # present in the sample, plot a rectangle.
        smp_series = self.parent.seq_relative_abundance_df.loc[sample_to_plot]
        non_zero_series = smp_series.iloc[smp_series.to_numpy().nonzero()[0]]
        for seq_name, seq_rel_abund in non_zero_series.iteritems():
            # class matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
            self.patches_list.append(Rectangle((self.ind - 0.5, bottom_div), 1, seq_rel_abund, color=self.parent.seq_color_dict[seq_name]))
            self.color_list.append(self.parent.seq_color_dict[seq_name])
            bottom_div += seq_rel_abund

    def _plot_type_under_div(self, sample_to_plot):
        # the idea of the type is to put it as a reflection below the y=0 line
        # as such we should just want to make everything negative
        bottom_prof = 0
        # for each sequence, create a rect patch
        # the rect will be 1 in width and centered about the ind value.
        # we want to plot the rects so that they add to 1. As such we want to divide
        # each value by the total for that sample.
        non_zero_indices = self.parent.prof_absolute_abundance_df.loc[sample_to_plot].to_numpy().nonzero()[0]
        non_zero_series = self.parent.prof_absolute_abundance_df.loc[sample_to_plot].iloc[non_zero_indices]
        total = sum(non_zero_series.values)
        non_zero_series_relative = non_zero_series / total
        for profile_uid, profile_rel_abund in non_zero_series_relative.iteritems():
            # We will scale the profile so that it is 0.2 of the length of the seq info
            depth = -0.2 * profile_rel_abund
            self.patches_list.append(
                Rectangle((self.ind - 0.5, bottom_prof), 1, depth,
                            color=self.parent.prof_color_dict[profile_uid]))
            # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
            self.color_list.append(self.parent.prof_color_dict[profile_uid])
            bottom_prof += depth

    def _paint_rect_to_axes_div_and_type(self, max_num_smpls_in_subplot=10):
        # We can try making a custom colour map
        # https://matplotlib.org/api/_as_gen/matplotlib.colors.ListedColormap.html
        this_cmap = ListedColormap(self.color_list)
        
        # here we should have a list of Rectangle patches
        # now create the PatchCollection object from the patches_list
        patches_collection = PatchCollection(self.patches_list, cmap=this_cmap)
        patches_collection.set_array(np.arange(len(self.patches_list)))
        # if n_subplots is only 1 then we can refer directly to the axarr object
        # else we will need ot reference the correct set of axes with i
        # Add the pathces to the axes
        self.ax.add_collection(patches_collection)
        self.ax.autoscale_view()
        self.ax.figure.canvas.draw()
        # also format the axes.
        # make it so that the x axes is constant length
        self.ax.set_xlim(0 - 0.5, max_num_smpls_in_subplot - 0.5)
        self.ax.set_ylim(-0.2, 1)
        self._remove_axes_but_allow_labels()

        # as well as getting rid of the top and right axis splines
        # I'd also like to restrict the bottom spine to where there are samples plotted but also
        # maintain the width of the samples
        # I think the easiest way to do this is to hack a bit by setting the x axis spines to invisible
        # and then drawing on a line at y = 0 between the smallest and largest ind (+- 0.5)
        # ax.spines['bottom'].set_visible(False)
        self.ax.add_line(Line2D((0 - 0.5, self.num_smp_in_this_subplot - 0.5), (0, 0), linewidth=0.5, color='black'))

    def _remove_axes_but_allow_labels(self):
        self.ax.set_frame_on(False)
        self.ax.set_yticks([])
        self.ax.set_xticks([])

class QCPlotterCoralOnly:
    def __init__(self, info_df, sp_datasheet, seq_meta_info_df, output_dir):
        self.info_df = info_df
        self.sp_datasheet = sp_datasheet
        self.seq_meta_info_df = seq_meta_info_df
        self.output_dir = output_dir
        # Dict that will have island names as key and a list of sites as value
        self.island_site_dict = self._determine_sites_and_island()
        self.islands = sorted(list(self.island_site_dict.keys()))
        self.num_coral_islands = len(self.island_site_dict.keys())
        self.num_coral_sites = 0
        for k, v in self.island_site_dict.items():
            self.num_coral_sites += len(v)
        self.fig = plt.figure(figsize=(14, 10))
        self.gs = self.fig.add_gridspec(3, 32, figure=self.fig, height_ratios=[1, 1, 1], width_ratios=[1 for _ in range(32)])

        # setup each of the axes in separate lists for raw_contgs, non_symbiodiniaceae and symbiodiniaceae
        self.raw_contig_ax_list = []
        for i in range(32):
            self.raw_contig_ax_list.append(self.fig.add_subplot(self.gs[0,i]))
        self.non_sym_ax_list = []
        for i in range(32):
            self.non_sym_ax_list.append(self.fig.add_subplot(self.gs[1, i]))
        self.sym_ax_list = []
        for i in range(32):
            self.sym_ax_list.append(self.fig.add_subplot(self.gs[2, i]))

        # Here we have each of the axes setup. Now we just need to plot in them
        
        foo = 'bar'


    def _determine_sites_and_island(self):
        # For each sample pull out the island and site from the info_df
        # Dict that will have island names as key and a list of sites as value
        island_site_dict = defaultdict(set)
        for sample_index in self.sp_datasheet.index:
            island = self.info_df.at[sample_index, 'location']
            site = self.info_df.at[sample_index, 'site']
            island_site_dict[island].add(site)
        return island_site_dict

    def plot_qc_data(self):
        # We can go island by island.
        # For each island get the sites and sort
        # For each site of island get the samples
        # Then for each sample plot

        # We will need to keep track of some grand maximums so that we can scale
        # each of the subplot axes to the same values
        raw_contig_max = int(self.seq_meta_info_df['raw_contigs'].max())
        raw_contig_min = int(self.seq_meta_info_df['raw_contigs'].min())
        non_sym_total_max = int(self.seq_meta_info_df['post_taxa_id_absolute_non_symbiodinium_seqs'].max())
        non_sym_total_min = int(self.seq_meta_info_df['post_taxa_id_absolute_non_symbiodinium_seqs'].min())
        non_sym_distinct_max = int(self.seq_meta_info_df['post_taxa_id_unique_non_symbiodinium_seqs'].max())
        non_sym_distinct_min = int(self.seq_meta_info_df['post_taxa_id_unique_non_symbiodinium_seqs'].min())
        sym_total_max = int(self.seq_meta_info_df['post_taxa_id_absolute_symbiodinium_seqs'].max())
        sym_total_min = int(self.seq_meta_info_df['post_taxa_id_absolute_symbiodinium_seqs'].min())
        sym_distinct_max = int(self.seq_meta_info_df['post_taxa_id_unique_symbiodinium_seqs'].max())
        sym_distinct_min = int(self.seq_meta_info_df['post_taxa_id_unique_symbiodinium_seqs'].min())

        # We want to output a few stats to go with this figure. It would also be good to annotate these stats on
        # the figure. I think it would be good to do this using a horizonatl line
        print(f'Number of samples = {len(self.seq_meta_info_df.index)}')
        print(f'\tPorites lobata: {len(self.sp_datasheet[self.sp_datasheet["host_species"] == "lobata"].index)}')
        print(f'\tPocillopora meandrina: {len(self.sp_datasheet[self.sp_datasheet["host_species"] == "meandrina"].index)}')
        print(f'\tMillepora dichotoma: {len(self.sp_datasheet[self.sp_datasheet["host_species"] == "dichotoma"].index)}')
        average_raw_contigs_per_sample = int(self.seq_meta_info_df['raw_contigs'].mean())
        print(f'average_raw_contigs_per_sample: {average_raw_contigs_per_sample}')
        average_non_sym_total = int(self.seq_meta_info_df['post_taxa_id_absolute_non_symbiodinium_seqs'].mean())
        print(f'average_non_sym_total: {average_non_sym_total}')
        average_non_sym_distinct = int(self.seq_meta_info_df['post_taxa_id_unique_non_symbiodinium_seqs'].mean())
        print(f'average_non_sym_distinct: {average_non_sym_distinct}')
        average_sym_total = int(self.seq_meta_info_df['post_taxa_id_absolute_symbiodinium_seqs'].mean())
        print(f'average_sym_total: {average_sym_total}')
        average_sym_distinct = int(self.seq_meta_info_df['post_taxa_id_unique_symbiodinium_seqs'].mean())
        print(f'average_sym_distinct: {average_sym_distinct}')
        print(f'\n\ntotal_raw_contigs: {self.seq_meta_info_df["raw_contigs"].sum()}')

        for island_to_plot in self.islands:
            sites = sorted(list(self.island_site_dict[island_to_plot]))
            # For each of the subplots we will treat them as a single scatter plot
            # We will work with the x axis being limited to 0-->1
            # x axis coordinates will depend on the number of sites
            # we will want four plotting points per site,
            # one for the individual points and one for the avergae for total seqs
            # and then we'll want one for distnct sequences
            num_sites = len(sites)
            x_space = 1 / ((num_sites * 4) + 1)
            x_coord_vals = [i * x_space for i in range(1, num_sites*4 + 1, 1)]
            # Pair up the x_coord_vals so that they are easier to index
            # one pair per site. 0 will be for individual datapoints, 1 for the mean
            x_coord_vals = [(x_coord_vals[i], x_coord_vals[i + 1], x_coord_vals[i + 2], x_coord_vals[i + 3]) for i in range(0, num_sites*4, 4)]

            raw_contigs_ax = self.raw_contig_ax_list[self.islands.index(island_to_plot)]
            raw_contigs_ax.set_xlim(0, 1)
            raw_contigs_ax.set_ylim(10000, raw_contig_max)
            # set all x axes labels off
            raw_contigs_ax.set_xticks([])

            non_sym_ax = self.non_sym_ax_list[self.islands.index(island_to_plot)]
            non_sym_ax.set_xlim(0, 1)
            non_sym_ax.set_ylim(non_sym_total_min, non_sym_total_max)
            non_sym_ax.set_xticks([])
            non_sym_ax_distinct = non_sym_ax.twinx()
            non_sym_ax_distinct.set_ylim(non_sym_distinct_min, non_sym_distinct_max)

            sym_ax = self.sym_ax_list[self.islands.index(island_to_plot)]
            sym_ax.set_xlim(0, 1)
            sym_ax.set_ylim(sym_total_min, sym_total_max)
            sym_ax.set_xticks([])
            sym_ax_distinct = sym_ax.twinx()
            sym_ax_distinct.set_ylim(100, sym_distinct_max)
            # place the island label on the x axis rotated
            sym_ax.set_xlabel(island_to_plot, rotation='vertical')

            if self.islands.index(island_to_plot) == 0:
                # If first plot then only remove top and right
                # and we only want to remove the right y axis
                raw_contigs_ax.spines['right'].set_visible(False)
                raw_contigs_ax.spines['top'].set_visible(False)
                raw_contigs_ax.set_yscale('symlog')
                raw_contigs_ax.spines['left'].set_color(c='blue')
                raw_contigs_ax.tick_params('y', colors='blue')

                # plot average line
                raw_contigs_ax.axhline(y=average_raw_contigs_per_sample, xmin=0, xmax=1, c='blue')

                non_sym_ax.spines['right'].set_visible(False)
                non_sym_ax.spines['top'].set_visible(False)
                non_sym_ax.set_yscale('symlog')
                non_sym_ax.spines['left'].set_color(c='blue')
                non_sym_ax.tick_params('y', colors='blue')

                non_sym_ax_distinct.spines['right'].set_visible(False)
                non_sym_ax_distinct.spines['top'].set_visible(False)
                non_sym_ax_distinct.set_yscale('symlog')
                non_sym_ax_distinct.set_yticks([])
                non_sym_ax_distinct.minorticks_off()
                non_sym_ax_distinct.spines['left'].set_color(c='blue')

                sym_ax.spines['right'].set_visible(False)
                sym_ax.spines['top'].set_visible(False)
                sym_ax.set_yscale('symlog')
                sym_ax.spines['left'].set_color(c='blue')
                sym_ax.tick_params('y', colors='blue')

                sym_ax_distinct.spines['right'].set_visible(False)
                sym_ax_distinct.spines['top'].set_visible(False)
                sym_ax_distinct.set_yscale('symlog')
                sym_ax_distinct.set_yticks([])
                sym_ax_distinct.minorticks_off()
                sym_ax_distinct.spines['left'].set_color(c='blue')

            elif self.islands.index(island_to_plot) == 31:
                # If the last plot then we want to annotate the right y axis
                raw_contigs_ax.spines['left'].set_visible(False)
                raw_contigs_ax.spines['top'].set_visible(False)
                raw_contigs_ax.spines['right'].set_visible(False)
                raw_contigs_ax.set_yscale('symlog')
                raw_contigs_ax.set_yticks([])
                raw_contigs_ax.minorticks_off()


                non_sym_ax.spines['left'].set_visible(False)
                non_sym_ax.spines['top'].set_visible(False)
                non_sym_ax.set_yscale('symlog')
                non_sym_ax.set_yticks([])
                non_sym_ax.minorticks_off()
                non_sym_ax.spines['right'].set_color(c='red')


                non_sym_ax_distinct.spines['left'].set_visible(False)
                non_sym_ax_distinct.spines['top'].set_visible(False)
                non_sym_ax_distinct.set_yscale('symlog')
                non_sym_ax_distinct.spines['right'].set_color(c='red')
                non_sym_ax_distinct.tick_params('y', colors='red')


                sym_ax.spines['left'].set_visible(False)
                sym_ax.spines['top'].set_visible(False)
                sym_ax.set_yscale('symlog')
                sym_ax.set_yticks([])
                sym_ax.minorticks_off()
                sym_ax.spines['right'].set_color(c='red')

                sym_ax_distinct.spines['left'].set_visible(False)
                sym_ax_distinct.spines['top'].set_visible(False)
                sym_ax_distinct.set_yscale('symlog')
                sym_ax_distinct.spines['right'].set_color(c='red')
                sym_ax_distinct.tick_params('y', colors='red')
            else:
                raw_contigs_ax.spines['left'].set_visible(False)
                raw_contigs_ax.spines['top'].set_visible(False)
                raw_contigs_ax.spines['right'].set_visible(False)
                raw_contigs_ax.set_yscale('symlog')
                raw_contigs_ax.set_yticks([])
                raw_contigs_ax.minorticks_off()

                non_sym_ax.spines['left'].set_visible(False)
                non_sym_ax.spines['top'].set_visible(False)
                non_sym_ax.spines['right'].set_visible(False)
                non_sym_ax.set_yscale('symlog')
                non_sym_ax.set_yticks([])
                non_sym_ax.minorticks_off()

                non_sym_ax_distinct.spines['left'].set_visible(False)
                non_sym_ax_distinct.spines['top'].set_visible(False)
                non_sym_ax_distinct.spines['right'].set_visible(False)
                non_sym_ax_distinct.set_yscale('symlog')
                non_sym_ax_distinct.set_yticks([])
                non_sym_ax_distinct.minorticks_off()

                sym_ax.spines['left'].set_visible(False)
                sym_ax.spines['top'].set_visible(False)
                sym_ax.spines['right'].set_visible(False)
                sym_ax.set_yscale('symlog')
                sym_ax.set_yticks([])
                sym_ax.minorticks_off()

                sym_ax_distinct.spines['left'].set_visible(False)
                sym_ax_distinct.spines['top'].set_visible(False)
                sym_ax_distinct.spines['right'].set_visible(False)
                sym_ax_distinct.set_yscale('symlog')
                sym_ax_distinct.set_yticks([])
                sym_ax_distinct.minorticks_off()


            for site in sites:
                samples = self.info_df[(self.info_df['location'] == island_to_plot) & (self.info_df['site'] == site) & (self.info_df['coral_plankton'] == 'CORAL') & (self.info_df['spp_water'] != 'HELIOPORA') & (self.info_df['spp_water'] != 'PORITES_PANAMENSIS')]
                # Now we can populate the sample information for each of the plots
                # lists that we will create the avergae and stdev points from
                raw_contigs_vals = []
                non_sym_vals_total = []
                non_sym_vals_distinct = []
                sym_vals_total = []
                sym_vals_distinct = []



                # the x_coordinates
                indi_data_point_x_total = x_coord_vals[sites.index(site)][0]
                mean_data_point_x_total = x_coord_vals[sites.index(site)][1]
                indi_data_point_x_distinct = x_coord_vals[sites.index(site)][2]
                mean_data_point_x_distinct = x_coord_vals[sites.index(site)][3]

                # Now plot up the individual points and collect for calculating the mean
                for sample_name, sample_ser in samples.iterrows():
                    # raw_contigs
                    y_raw_contig = self.seq_meta_info_df.at[sample_name, 'raw_contigs']
                    raw_contigs_vals.append(y_raw_contig)
                    raw_contigs_ax.scatter(x=indi_data_point_x_total, y=y_raw_contig, marker='.', s=1, c='b')

                    # non_sym_vals total and distinct
                    y_non_sym_total = self.seq_meta_info_df.at[sample_name, 'post_taxa_id_absolute_non_symbiodinium_seqs']
                    non_sym_vals_total.append(y_non_sym_total)
                    non_sym_ax.scatter(x=indi_data_point_x_total, y=y_non_sym_total, marker='.', s=1, c='b')
                    y_non_sym_distinct = self.seq_meta_info_df.at[sample_name, 'post_taxa_id_unique_non_symbiodinium_seqs']
                    non_sym_vals_distinct.append(y_non_sym_distinct)
                    non_sym_ax_distinct.scatter(x=indi_data_point_x_distinct, y=y_non_sym_distinct, marker='.', s=1, c='r')

                    # sym counts total and distinct
                    y_sym_total = self.seq_meta_info_df.at[
                        sample_name, 'post_taxa_id_absolute_symbiodinium_seqs']
                    sym_vals_total.append(y_sym_total)
                    sym_ax.scatter(x=indi_data_point_x_total, y=y_sym_total, marker='.', s=1, c='b')
                    y_sym_distinct = self.seq_meta_info_df.at[
                        sample_name, 'post_taxa_id_unique_symbiodinium_seqs']
                    sym_vals_distinct.append(y_sym_distinct)
                    sym_ax_distinct.scatter(x=indi_data_point_x_distinct, y=y_sym_distinct, marker='.', s=1, c='r')

                # Here we should have all of the individual data points plotted up
                # We should also have collected all of the opints in a list so that we can now
                # calculate the mean and the standard deviations
                # For the time being just plot the mean and worry about stdev bars later
                mean_point_size = 40
                raw_contigs_ax.scatter(x=mean_data_point_x_total, y=sum(raw_contigs_vals)/len(raw_contigs_vals), marker='.', s=mean_point_size, c='b')
                non_sym_ax.scatter(x=mean_data_point_x_total, y=sum(non_sym_vals_total)/len(non_sym_vals_total), marker='.', s=mean_point_size, c='b')
                non_sym_ax_distinct.scatter(x=mean_data_point_x_distinct, y=sum(non_sym_vals_distinct)/len(non_sym_vals_distinct), marker='.', s=mean_point_size, c='r')
                sym_ax.scatter(x=mean_data_point_x_total, y=sum(sym_vals_total) / len(sym_vals_total), marker='.', s=mean_point_size, c='b')
                sym_ax_distinct.scatter(x=mean_data_point_x_distinct,
                                   y=sum(sym_vals_distinct) / len(sym_vals_distinct), marker='.', s=mean_point_size, c='r')

            # plot up the total average lines that will be extended by hand
            # plot on the 0 and 1 so that we can still see both lines despite overlap.
            if self.islands.index(island_to_plot) == 0:
                raw_contigs_ax.axhline(y=average_raw_contigs_per_sample, xmin=0, xmax=1, c='blue')
                non_sym_ax.axhline(y=average_non_sym_total, xmin=0, xmax=1, c='blue')
                sym_ax.axhline(y=average_sym_total, xmin=0, xmax=1, c='blue')
            if self.islands.index(island_to_plot) == 1:
                non_sym_ax_distinct.axhline(y=average_non_sym_distinct, xmin=0, xmax=1, c='red')
                sym_ax_distinct.axhline(y=average_sym_distinct, xmin=0, xmax=1, c='red')

        plt.savefig(os.path.join(self.output_dir, 'coral_only_qc_fig.svg'))
        plt.savefig(os.path.join(self.output_dir, 'coral_only_qc_fig.png'), dpi=1200)

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

class GenerateInfoCollectionDF:
    """Class concerned with creating the information dataframe"""
    def __init__(self, base_dir, provenance_df):
        self.base_dir = base_dir
        self.provenance_df = provenance_df
        self.rows = []

    def generate_df(self):
        self._create_data_rows()
        return pd.DataFrame(data=self.rows, columns=['sample_name', 'fastq_fwd_file_path', 'fastq_rev_file_path', 'coral_plankton', 'spp_water', 'location', 'site', 'lat', 'lon']).set_index(keys='sample_name')
    
    def _create_data_rows(self):
        """ Parse through the seq data file structure gathering 
        inferring the info for each sequecning file pair as we descend the structure/
        Collect a single row of data for each sample"""
        for location in os.listdir(self.base_dir):
            if 'ISLAND' in location:
                parsing_dir_loc = os.path.join(self.base_dir, location)
                for site in os.listdir(parsing_dir_loc):
                    parsing_dir_site = os.path.join(parsing_dir_loc, site)
                    for sample_type in os.listdir(parsing_dir_site):
                        parsing_dir_sample_type = os.path.join(parsing_dir_site, sample_type)
                        if sample_type == 'CORAL':
                            for species in os.listdir(parsing_dir_sample_type):
                                parsing_dir_species = os.path.join(parsing_dir_sample_type, species)
                                for individual in os.listdir(parsing_dir_species):
                                    parsing_dir_indi = os.path.join(parsing_dir_species, individual, 'CS4L')
                                    # now we are in the directory that contains the actual paired fastq.gz files for a
                                    # given coral individual
                                    # collect the information we need
                                    srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_indi, sample_type=sample_type, site=site, water_species=species, provenance_df=self.provenance_df)
                                    self.rows.append(srg.create_sample_row())
                                    print(f'Processed: {parsing_dir_indi}')

                        elif sample_type == 'PLANKTON':
                            for water_type in os.listdir(parsing_dir_sample_type):
                                if water_type == 'CSW':
                                    parsing_dir_water_type = os.path.join(parsing_dir_sample_type, water_type)
                                    for individual in os.listdir(parsing_dir_water_type):
                                        parsing_dir_indi = os.path.join(parsing_dir_water_type, individual, 'S320')
                                        # now we are in the directory that contains the actual paired fastq.gz files for a
                                        # given water sample
                                        # collect the information we need
                                        srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_indi, sample_type=sample_type, site=site, water_species=water_type, provenance_df=self.provenance_df)
                                        self.rows.append(srg.create_sample_row())
                                        print(f'Processed: {parsing_dir_indi}')

                                elif water_type == 'SURFACE':
                                    # then this is a SURFACE sample and there are no individuals
                                    parsing_dir_water_type = os.path.join(parsing_dir_sample_type, water_type, 'S320')
                                    # collect the information we need
                                    srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_water_type, sample_type=sample_type, site=site, water_species=water_type, provenance_df=self.provenance_df)
                                    self.rows.append(srg.create_sample_row())
                                    print(f'Processed: {parsing_dir_water_type}')

            elif 'OA' in location:
                parsing_dir_loc = os.path.join(self.base_dir, location, 'PLANKTON', 'SURFACE', 'S320')
                # NB the arguments are a little messed up here on purpose as we don't have seperate site and location
                srg = SampleRowGenerator(location=location, parsing_dir=parsing_dir_loc, sample_type='OA', site=location, water_species='PLANKTON', provenance_df=self.provenance_df)
                self.rows.append(srg.create_sample_row())
                print(f'Processed: {parsing_dir_loc}')


class SampleRowGenerator:
    def __init__(self, location, parsing_dir, sample_type, site, water_species, provenance_df):
        self.parsing_dir = parsing_dir
        self.location = location
        self.sample_type = sample_type
        self.site = site
        self.water_type = water_species
        self.provenance_df=provenance_df
        self.sample_name = '_'.join(os.listdir(self.parsing_dir)[0].split('_')[:2])
        self.lat = self.provenance_df.at[self.sample_name, 'lat']
        self.lon = self.provenance_df.at[self.sample_name, 'lon']
        self.files = os.listdir(self.parsing_dir)
        if len(self.files) != 2:
            self._concatenate_seq_files()
            self.files = os.listdir(self.parsing_dir)
            if len(self.files) != 2:
                raise RuntimeError(f'Error in concatenation of seq files in {self.parsing_dir}')
        self.fwd_found = False
        self.rev_found = False
        self.fwd_path = None
        self.rev_path = None
        self._get_seq_paths()

    def _concatenate_seq_files(self):
        """ Sometime there is more than 1 set of sequencing files. In this case
        we want to concatenate all of the R1 reads together and all of the R2 reads
        together. Importantly we want to merge the respective memebers of the sequencing
        pairs in the same order.
        We will get a list of R1 files in the director (arbitrary order)
        We will then get a corresponding list of R2 files in the same realtive order
        by doing a replace function on the R1 list.
        We will then concatenate each list of files. Importantly, we will do this
        directly on the .gz files as this is apparently OK.
        https://www.biostars.org/p/81924/
        """
        # get files
        r1_files_list = [file_name for file_name in self.files if 'R1' in file_name]
        r2_files_list = [file_name.replace("R1", "R2") for file_name in r1_files_list]
        
        # Create the list that holds the cat command and the arguments
        exe_1 = [os.path.join(self.parsing_dir, _) for _ in r1_files_list]
        exe_1.insert(0, 'cat')
        exe_2 = [os.path.join(self.parsing_dir, _) for _ in r2_files_list]
        exe_2.insert(0, 'cat')
        
        # outpaths that will be used to create a stream that will be passed as stdout
        out_path_1 = os.path.join(self.parsing_dir, r1_files_list[0].replace("R1.fastq.gz", "merged.R1.fastq.gz"))
        out_path_2 = os.path.join(self.parsing_dir, r2_files_list[0].replace("R2.fastq.gz", "merged.R2.fastq.gz"))

        # do cat
        with open(out_path_1, 'wb') as f:
            subprocess.run(exe_1, stdout=f)
        with open(out_path_2, 'wb') as f:
            subprocess.run(exe_2, stdout=f)
        
        # rm the old files that have now been concatenated
        for file_path in exe_1[1:]:
            os.remove(file_path)
        for file_path in exe_2[1:]:
            os.remove(file_path)

    def _get_seq_paths(self):
        for file_name in self.files:
                if 'R1' in file_name:
                    self.fwd_path = os.path.join(self.parsing_dir, file_name)
                    self.fwd_found = True
                elif 'R2' in file_name:
                    self.rev_path = os.path.join(self.parsing_dir, file_name)
                    self.rev_found = True

    def create_sample_row(self):
        if not self.fwd_found or not self.rev_found:
            print('fwd or rev read not found')
            sys.exit(1)
        else:
            return [self.sample_name, self.fwd_path, self.rev_path, self.sample_type, self.water_type, self.location, self.site, self.lat, self.lon]

its2processing = ITS2Processing()
its2processing.seq_and_profile_results_figure()
