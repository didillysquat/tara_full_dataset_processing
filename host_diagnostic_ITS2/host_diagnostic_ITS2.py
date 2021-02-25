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
import matplotlib.pyplot as plt

class HostDiagnosticZooxs:
    def __init__(self, pre_post='post'):
        self.pre_post = pre_post
        
        # Read in the host group data
        # 'sampling-design_label' 'sample-id_source'
        self.host_groupings = pd.read_table("/home/humebc/projects/tara/tara_full_dataset_processing/host_diagnostic_ITS2/host_groupings.txt", names=["INDV", "clus", "svd_cluster"])
        self.host_groupings.index = self.host_groupings["INDV"]
        self.host_groupings = self.host_groupings["svd_cluster"]
        
        # Read in Guillems golden table
        poc_golden = pd.read_table("/home/humebc/projects/tara/tara_full_dataset_processing/input/data_available_gold_dataset_11_islands_POC.tsv")
        poc_golden.index = poc_golden["INDV"]
        # We want to work with the subset of samples for which there is host grouping info and ITS2 sequences
        poc_golden = poc_golden[~pd.isnull(poc_golden["dataset_ITS2_sample_id"])]
        
        if self.pre_post == 'post':
            self.path_to_zooxs_counts = "/home/humebc/projects/tara/tara_full_dataset_processing/host_diagnostic_ITS2/TARA_PACIFIC_METAB_ITS2_v1/coral/post_med_seqs/TARA_PACIFIC_METAB_ITS2_coral_post_med_seqs_sequences_absolute_abund_and_meta_v1.csv"
            # Get zooxs count data
            self.counts_df = pd.read_csv(self.path_to_zooxs_counts)
            self.counts_df = self.counts_df.iloc[:-1,]
            self.counts_df.index = self.counts_df["sample-id"]
            self.counts_df_with_host = self.counts_df.loc[[_ for _ in self.counts_df.index if _ in poc_golden["dataset_ITS2_sample_id"].values],]
            # For the count table this will likely leave us with a lot of columns where it is only 0s. we will want to get rid of these
            self.counts_df_with_host = self.counts_df_with_host.loc[:, (self.counts_df_with_host != 0).any(axis=0)]
        
        else:
            self.path_to_zooxs_counts = "/home/humebc/projects/tara/tara_full_dataset_processing/host_diagnostic_ITS2/TARA_PACIFIC_METAB_ITS2_v1/coral/pre_med_seqs/TARA_PACIFIC_METAB_ITS2_coral_pre_med_seqs_sequences_absolute_abund_v1.csv"
            with open(self.path_to_zooxs_counts, 'r') as f:
                header = None
                rows = []
                for i, line in enumerate(f):
                    if i == 0:
                        header = line.rstrip().split(',')
                    else:
                        split_line_list = line.rstrip().split(',')
                        if split_line_list[1] in poc_golden["dataset_ITS2_sample_id"].values:
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
        didier_to_sp_name_dict = {k:v for k, v in zip(poc_golden["INDV"].values, poc_golden["dataset_ITS2_sample_id"]) if not pd.isnull(v)}
        self.host_groupings = self.host_groupings.loc[[_ for _ in self.host_groupings.index.values if _ in didier_to_sp_name_dict],]
        self.host_groupings.index = [didier_to_sp_name_dict[_] for _ in self.host_groupings.index.values]
        self.group_to_sample_list_dict = {k: list(self.host_groupings[self.host_groupings==k].index.values) for k in self.host_groupings.unique() if not pd.isnull(k)}

        
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
        foo = 'bar'

         
# The look_for_diagnostic_sequences can be run using pre_post='pre' to use the pre_med seqs.
# However this still doesn't produce any viable diagnostic sequences and takes a lot longer due to the number of sequnces
HostDiagnosticZooxs(pre_post='pre').look_for_diagnostic_sqeuences()

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