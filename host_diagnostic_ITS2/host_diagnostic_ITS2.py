"""
Script for investigating whether we can find diagnostic ITS2 sequences that can 'predict' host genotypes as determined by Didier
First approach is simply to look for single diagnostic sequences. This means that for each ITS2 sequences we look to
see if it is present in all samples of a given host grouping but not present in all other host groupings.
Perhaps easiest to produce a multiway venn for this.
"""

import pandas as pd
import re
from collections import defaultdict
# NB this had to be installed with pip install venn. Couldn't find a conda install
from venn import venn
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from sputils.sphierarchical import SPHierarchical
from sputils.spbars import SPBars
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import pathlib
import os

class HostDiagnosticZooxs:
    def __init__(self, pre_post='post', genus='Cladocopium', coral='Pocillopora'):
        self.base_dir = pathlib.Path(__file__).parent.absolute()
        self.input_dir = os.path.join(self.base_dir, 'inputs')
        self.figure_dir = os.path.join(self.base_dir, "figures")
        self.pre_post = pre_post
        self.genus = genus
        self.coral = coral
        # Read in the host group data
        # 'sampling-design_label' 'sample-id_source'
        if coral == 'Pocillopora':
            self.host_groupings = pd.read_table(os.path.join(self.input_dir, "host_groupings_POC.txt"), names=["INDV", "clus", "svd_cluster"])
            self.host_groupings.index = self.host_groupings["INDV"]
            self.host_groupings = self.host_groupings["svd_cluster"]
        elif coral == 'Porites':
            self.host_groupings = pd.read_table(os.path.join(self.input_dir, "host_groupings_POR.txt"),
                                                names=["INDV", "cluster_snmf", "sub_cluster_snmf", "SVDQuartet"])
            self.host_groupings.index = self.host_groupings["INDV"]
            self.host_groupings = self.host_groupings["sub_cluster_snmf"]
            self.host_groupings = self.host_groupings[~pd.isnull(self.host_groupings)]
        
        # Read in Guillems golden table
        if coral == 'Pocillopora':
            self.golden_table = pd.read_table(os.path.join(self.input_dir, "data_available_gold_dataset_11_islands_POC.tsv"))
            self.golden_table.index = self.golden_table["INDV"]
            # We want to work with the subset of samples for which there is host grouping info and ITS2 sequences
            self.golden_table = self.golden_table[~pd.isnull(self.golden_table["dataset_ITS2_sample_id"])]
            # Also remove the one sample that there isn't an SVDQ value for
            self.golden_table = self.golden_table[~pd.isnull(self.golden_table["SVDQ"])]
        elif coral == 'Porites':
            self.golden_table = pd.read_table(
                os.path.join(self.input_dir, "data_available_gold_dataset_11_islands_POR.tsv"))
            self.golden_table.index = self.golden_table["Ind"]
            # We want to work with the subset of samples for which there is host grouping info and ITS2 sequences
            self.golden_table = self.golden_table[~pd.isnull(self.golden_table["dataset_ITS2_sample_id"])]
            # Also remove any samples that there isn't an sub_cluster_snmf value for
            self.golden_table = self.golden_table[~pd.isnull(self.golden_table["sub_cluster_snmf"])]

        if self.pre_post == 'post':
            self.path_to_zooxs_counts = os.path.join(self.input_dir, "TARA_PACIFIC_METAB_ITS2_coral_post_med_seqs_sequences_absolute_abund_and_meta_v1.csv")
            # Get zooxs count data
            self.counts_df = pd.read_csv(self.path_to_zooxs_counts)
            self.counts_df = self.counts_df.iloc[:-1,]
            self.counts_df.index = self.counts_df["sample-id"]
            self.counts_df_with_host = self.counts_df.loc[[_ for _ in self.counts_df.index if _ in self.golden_table["dataset_ITS2_sample_id"].values],]
            # For the count table this will likely leave us with a lot of columns where it is only 0s. we will want to get rid of these
            self.counts_df_with_host = self.counts_df_with_host.loc[:, (self.counts_df_with_host != 0).any(axis=0)]
        
        else:
            self.path_to_zooxs_counts = os.path.join(self.input_dir, "TARA_PACIFIC_METAB_ITS2_coral_pre_med_seqs_sequences_absolute_abund_v1.csv.zip")
            with open(self.path_to_zooxs_counts, 'r') as f:
                header = None
                rows = []
                for i, line in enumerate(f):
                    if i == 0:
                        header = line.rstrip().split(',')
                    else:
                        split_line_list = line.rstrip().split(',')
                        if split_line_list[1] in self.golden_table["dataset_ITS2_sample_id"].values:
                            rows.append(split_line_list)
            self.counts_df_with_host = pd.DataFrame(rows, columns=header)
            self.counts_df_with_host.index = self.counts_df_with_host["sample-id"]
            self.counts_df_with_host = self.counts_df_with_host.loc[:, (self.counts_df_with_host != 0).any(axis=0)]
            foo = 'bar'

        
        # Now strip this down to only the abundances
        # Get index of first sequence
        drop_index = 0
        for i, col in enumerate(list(self.counts_df_with_host)):
            if col[0] in list('ABCDEFGHI') or col[-1] in list('ABCDEFGHI'):
                drop_index = i
                break
        self.counts_df_with_host = self.counts_df_with_host.iloc[:,drop_index:]

        # Need to convert between didiers naming system and the coral sample-id format
        didier_to_sp_name_dict = {k:v for k, v in zip(self.golden_table.index.values, self.golden_table["dataset_ITS2_sample_id"]) if not pd.isnull(v)}
        self.host_groupings = self.host_groupings.loc[[_ for _ in self.host_groupings.index.values if _ in didier_to_sp_name_dict],]
        self.sample_to_island_dict = {}
        for did_name in self.host_groupings.index:
            sp_name = didier_to_sp_name_dict[did_name]
            island = f"Island {did_name[1:3]}"
            self.sample_to_island_dict[sp_name] = island
        self.host_groupings.index = [didier_to_sp_name_dict[_] for _ in self.host_groupings.index.values]
        self.group_to_sample_list_dict = {k: list(self.host_groupings[self.host_groupings==k].index.values) for k in self.host_groupings.unique() if not pd.isnull(k)}
        self.sample_to_host_group_dict = dict(self.host_groupings)
        
    def look_for_diagnostic_sqeuences(self):
        """ 
        For every sequence in the ITS2 count table, look to see if it is found in every individual
            of a given group if it is, then count this sequence. Keep track of which groups a given
            sequence is found in and then we will output this info as a Venn or table. We will obviously
            be particularly interested in those sequences that are only found in one of the groups.

        Result:
        In the end there are very few sequences that are found in common with all samples of a host group.
        This approach will not work. Single sequences cannot be diagnostic.
        """
        host_group_to_seq_dd = defaultdict(set)
        for svd_group, sample_list in self.group_to_sample_list_dict.items():
            tot = len(list(self.counts_df_with_host))
            for i, seq in enumerate(list(self.counts_df_with_host)):
                print(f"{svd_group}:{i}/{tot}")
                ser = self.counts_df_with_host[seq]
                ser = ser[ser != 0]
                if set(set(sample_list)).issubset(ser.index.values):
                    # Then this sequence is found in all samples of the given host group
                    host_group_to_seq_dd[svd_group].add(seq)
        # At this point we know which sequences are found in all samples of a given group
        # now we can plot this up as a venn
        venn_obj = venn(host_group_to_seq_dd)
        venn_obj.set_title("Venn of sequnces found in all\nsamples of a given host group")
        plt.savefig('/home/humebc/projects/tara/tara_full_dataset_processing/host_diagnostic_ITS2/venn_plot.png' )
        
        # host_group_to_seq_dd = defaultdict(set)
        # for seq in list(self.counts_df_with_host):
        #     ser = self.counts_df_with_host[seq]
        #     ser = ser[ser != 0]
        #     # Check to see if, of the samples this seq is found in, if at least
        #     # one of those samples if from one of the host groups
        #     for svd_group, sample_list in self.group_to_sample_list_dict.items():
        #         if len(set(ser.index.values).intersection(set(sample_list))) > 1:
        #             # Then at least one of the samples that this seq is found in is of the host group
        #             host_group_to_seq_dd[svd_group].add(seq)
        # venn_obj = venn(host_group_to_seq_dd)
        # venn_obj.set_title("Venn of sequnces found in all\nsamples of a given host group")
        # plt.savefig('/home/humebc/projects/tara/tara_full_dataset_processing/host_diagnostic_ITS2/venn_plot.png' )
        # foo = 'this'

        # Venn is not really right for what we want to show here.
        # let's just straight up search for what we're after
        host_group_to_seq_dd = defaultdict(set)
        for svd_group, sample_list in self.group_to_sample_list_dict.items():
            for seq in list(self.counts_df_with_host):
                ser = self.counts_df_with_host[seq]
                ser = ser[ser != 0]
                # Check to see if this sequences is found in all samples of this group
                # and also none of the samples of the other groups
                if set(set(sample_list)).issubset(ser.index.values):
                    # Seq is found in all samples of this host group
                    found_in_other = False
                    for svd_group_other, sample_list_other in [(k, v) for k, v in self.group_to_sample_list_dict.items() if k != svd_group]:
                        # For all of the other svd_groups
                        if len(set(ser.index.values).intersection(set(sample_list_other))) > 0:
                            found_in_other = True
                    if found_in_other:
                        continue
                    else:
                        host_group_to_seq_dd[svd_group].add(seq)    
                else:
                    continue
        print("Sequences that are unique diagnostic of the host group:")
        print(host_group_to_seq_dd)
        
    def investigate_diagnostic_clusters(self):
        """
        Here we will investigate whether here is some form of structure within the ITS2 samples that would allow us
        to correlate to the species grouping.

        First thing to do will be to do will be to plot up the samples accordig to similarity
        and annotate them with the groupings so that we can get an idea of the correlation.
        It will probably be useful to do with only those ITS2 samples that have a group annotation available
        and also with all of the samples.
        """
        # To start with we will do a 1  up one down plot where we have the
        # hierarchical on top of the 99 ITS2 sequences that have the host data associated with them
        # on bottom we will plot an annotation of the host group. We will hope to see clustering.
        # We will need to do this for each of the C and D clades. Start with C as this is the most abundant
        if self.coral == 'Pocillopora':
            fig = plt.figure(figsize=(11, 6))
            # 4 down 1 across
            gs = gridspec.GridSpec(11, 2)
            axes = []
            plot_tyes = ['hier', 'anot']
            hier_ax = plt.subplot(gs[0:4,:])
            seq_bars_ax = plt.subplot(gs[4:6, :])
            seq_leg_ax = plt.subplot(gs[6:7, :])
            anot_ax = plt.subplot(gs[7:8,:])
            anot_leg_ax = plt.subplot(gs[8:9, :])
            island_ax = plt.subplot(gs[9:10, :])
            island_leg_ax = plt.subplot(gs[10:11, :])
        elif self.coral == "Porites":
            fig = plt.figure(figsize=(11, 6))
            # 4 down 1 across
            gs = gridspec.GridSpec(13, 2)
            axes = []
            plot_tyes = ['hier', 'anot']
            hier_ax = plt.subplot(gs[0:4, :])
            seq_bars_ax = plt.subplot(gs[4:6, :])
            seq_leg_ax = plt.subplot(gs[6:7, :])
            anot_ax = plt.subplot(gs[7:8, :])
            anot_leg_ax = plt.subplot(gs[8:9, :])
            anot_sub_ax = plt.subplot(gs[9:10, :])
            anot_sub_leg_ax = plt.subplot(gs[10:11, :])
            island_ax = plt.subplot(gs[11:12, :])
            island_leg_ax = plt.subplot(gs[12:13, :])
        if self.genus == 'Cladocopium':
            dist_df_path = os.path.join(self.input_dir, "2020-05-19_01-11-37.777185.braycurtis_sample_distances_C_sqrt.dist")
        elif self.genus == 'Durusdinium':
            dist_df_path = os.path.join(self.input_dir, "2020-05-19_01-11-37.777185.braycurtis_sample_distances_D_sqrt.dist")

        sph_plot = SPHierarchical(dist_output_path=dist_df_path, no_plotting=True)
        sample_names_in_current_dist = [sph_plot.obj_uid_to_obj_name_dict[_] for _ in sph_plot.dist_df.index.values]
        samples_to_keep = [_ for _ in sample_names_in_current_dist if _ in self.counts_df_with_host.index.values]
        sph_plot = SPHierarchical(
            dist_output_path=dist_df_path, ax=hier_ax,
                                  sample_names_included=samples_to_keep)
        sph_plot.plot()
        hier_ax.spines['right'].set_visible(False)
        hier_ax.spines['top'].set_visible(False)
        hier_ax.set_ylabel('Dissimilarity')
        hier_ax.set_title(f'{self.coral} - {self.genus}')

        spb_plot = SPBars(
            seq_count_table_path=os.path.join(self.input_dir, "98_20200331_DBV_2020-05-19_01-11-37.777185.seqs.absolute.abund_and_meta.txt"),
            profile_count_table_path=os.path.join(self.input_dir, "98_20200331_DBV_2020-05-19_01-11-37.777185.profiles.absolute.abund_and_meta.txt"),
            plot_type='seq_only', legend=True, relative_abundance=True, sample_uids_included=sph_plot.dendrogram_sample_order_uid, bar_ax=seq_bars_ax, seq_leg_ax=seq_leg_ax, limit_genera=[f'{self.genus[0]}']
        )
        spb_plot.plot()
        self._turn_off_spine_and_ticks(seq_bars_ax)
        seq_bars_ax.set_ylabel("ITS2\nseqs")

        # Finally we want to plot up some rectanles that will be the host_group annotations
        # And the island annotations
        # Problem TARA_CO-0000697 anot_sub_ax
        if self.coral == 'Porites':
            self._plot_annotations_and_legends(anot_ax=anot_ax, color_map_name='Dark2', leg_ax=anot_leg_ax,
                                               sample_to_annotation_dict={s: g[0] for s, g in self.sample_to_host_group_dict.items()}, sph_plot=sph_plot)
            anot_ax.set_ylabel("HostGroup")
            self._plot_annotations_and_legends(anot_ax=anot_sub_ax, color_map_name='Set1', leg_ax=anot_sub_leg_ax, sample_to_annotation_dict=self.sample_to_host_group_dict, sph_plot=sph_plot)
            anot_sub_ax.set_ylabel("HostGroup\nsub")
        elif self.coral == 'Pocillopora':
            self._plot_annotations_and_legends(anot_ax=anot_ax, color_map_name='Set1', leg_ax=anot_sub_ax,
                                               sample_to_annotation_dict=self.sample_to_host_group_dict,
                                               sph_plot=sph_plot)
            anot_ax.set_ylabel("Host group")
        self._plot_annotations_and_legends(anot_ax=island_ax, color_map_name='Set3', leg_ax=island_leg_ax,
                                           sample_to_annotation_dict=self.sample_to_island_dict, sph_plot=sph_plot)
        island_ax.set_ylabel("Island")
        plt.savefig(os.path.join(self.figure_dir, f"host_diagnostic_{self.coral}_{self.genus}.png"), dpi=600)
        plt.savefig(
            os.path.join(self.figure_dir, f"host_diagnostic_{self.coral}_{self.genus}.svg"),
            dpi=600)
        foo = 'bar'

    def _plot_annotations_and_legends(self, anot_ax, color_map_name, leg_ax, sample_to_annotation_dict, sph_plot):
        cat_cols = {svd: ctup for svd, ctup in
                    zip(list(set(sample_to_annotation_dict.values())),
                        list(cm.get_cmap(color_map_name).colors[:len(set(sample_to_annotation_dict.values()))]))}
        # Plot the annotations
        x_pos = 0
        patch_list = []
        for sample_name in sph_plot.dendrogram_sample_order_name:
            cat = sample_to_annotation_dict[sample_name]
            patch_list.append(
                Rectangle((x_pos, 0), width=1, height=1, facecolor=cat_cols[cat], edgecolor='black')
            )
            x_pos += 1
        collection = PatchCollection(patch_list, match_original=True)
        anot_ax.add_collection(collection)
        anot_ax.set_xlim(0, len(sph_plot.dendrogram_sample_order_name))
        anot_ax.set_ylim(0, 1)
        self._turn_off_spine_and_ticks(anot_ax)
        # Plot the annotation legends
        x_pos = 0.5
        patch_list = []
        for cat in sorted(list(cat_cols.keys())):
            col = cat_cols[cat]
            patch_list.append(Rectangle((x_pos-0.1, 0), height=0.2, width=0.2, facecolor=col))
            leg_ax.text(x=x_pos, y=0.7, ha='center', va='center', s=cat)
            x_pos += 1
        leg_ax.add_collection(PatchCollection(patch_list, match_original=True))
        leg_ax.set_xlim(0, len(cat_cols))
        leg_ax.set_ylim(0, 1)
        self._turn_off_spine_and_ticks(leg_ax)

    def _turn_off_spine_and_ticks(self, ax):
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_yticks([])
        ax.set_xticks([])


# The look_for_diagnostic_sequences can be run using pre_post='pre' to use the pre_med seqs.
# However this still doesn't produce any viable diagnostic sequences and takes a lot longer due to the number of sequnces
# NB there is no Durusdinium (roughly speaking) in the Porites so it won't plot.

hdz = HostDiagnosticZooxs(pre_post='post', genus='Cladocopium', coral='Porites')
hdz.look_for_diagnostic_sqeuences()
# hdz.investigate_diagnostic_clusters()

"""
Code that I was using to check why there were some ITS2 samples missing with relation to the
corals that Didier has host groupings for
# Read in and process the provenance table
        # NB we checked the provenance table to verify that the 4 samples that we can't match between Didiers hosts indeed
        # don't exist in the ITS2 table.
        # I looked to see if this is because we didn't output info for these samples, but it is not. Rather
        # the sequencing files for these samples never existed so the issue with the samples goes back to before SymPortal
        # processing.
        # As such we will move forwards with only the samples that we have matches for in the Didier table.
        # NB the provenance tables don't actually have the name format that exactly matches Didiers name format
        # So if we want to verfiy the the golden tables then we will need to convert this name format so that
        # it matches on of the fields of provenance.
        # Probably best to match to sampling-design_label. This means going from e.g. I01S01C001POC --> OA001-I00-S00-C000
        # Actually because we don't kow what the OA part is we will try to match to the sample_label e.g. TARA_SAMPLE_20160717T1545Z_OA000-I01-S01-C001_C-COL_CORAL_IMG-CREP_R01_CO-1000004
        # And we are looking for the SEQ-CSL4 versino of this sample_label - this is what was used for the ITS2, 18S and 16S.
        provenance = pd.read_table("/home/humebc/projects/tara/tara_full_dataset_processing/input/TARA-PACIFIC_samples-provenance_20200731d.txt", header=1)
        regex = re.compile(r"^I\d+S\d+Cd+")
        I01-S01-C001

        # Pull out the sampling-design_label for the sample-id
        sampling_design_label_to_sample_id_dict = {k: v for k, v in zip(provenance['sampling-design_label'].values, provenance['sample-id_source'].values) if v in self.counts_df['sample-id'].values}
        # prolem = I01S02C001POC --> TARA_CO-0000065
        self.host_groupings.index = [sampling_design_label_to_sample_id_dict[_] for _ in self.host_groupings['name'].values]
"""