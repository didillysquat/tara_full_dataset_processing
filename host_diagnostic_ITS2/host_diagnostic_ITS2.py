"""
Script for investigating whether we can find diagnostic ITS2 sequences that can 'predict' host genotypes as determined by Didier
First approach is simply to look for single diagnostic sequences. This means that for each ITS2 sequences we look to
see if it is present in all samples of a given host grouping but not present in all other host groupings.
Perhaps easiest to produce a multiway venn for this.
"""

import pandas as pd
import re

class HostDiagnosticZooxs:
    def __init__(self):
        self.path_to_zooxs_counts_pre_med = "/home/humebc/projects/tara/tara_full_dataset_processing/host_diagnostic_ITS2/TARA_PACIFIC_METAB_ITS2_v1/coral/pre_med_seqs/TARA_PACIFIC_METAB_ITS2_coral_pre_med_seqs_sequences_absolute_abund_v1.csv"
        self.path_to_zooxs_counts_post_med = "/home/humebc/projects/tara/tara_full_dataset_processing/host_diagnostic_ITS2/TARA_PACIFIC_METAB_ITS2_v1/coral/post_med_seqs/TARA_PACIFIC_METAB_ITS2_coral_post_med_seqs_sequences_absolute_abund_and_meta_v1.csv"
        self.path_to_zooxs_fasta = "/home/humebc/projects/tara/tara_full_dataset_processing/host_diagnostic_ITS2/TARA_PACIFIC_METAB_ITS2_v1/coral/pre_med_seqs/TARA_PACIFIC_METAB_ITS2_coral_pre_med_seqs_sequences_v1.fasta"
        # Get zooxs count data
        self.counts_df = pd.read_csv(self.path_to_zooxs_counts_post_med)
        self.counts_df = self.counts_df.iloc[:-1,]
        self.counts_df.index = self.counts_df["sample-id"]
        
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
        self.counts_df_with_host = self.counts_df.loc[[_ for _ in self.counts_df.index if _ in poc_golden["dataset_ITS2_sample_id"].values],]
        
        # For the count table this will likely leave us with a lot of columns where it is only 0s. we will want to get rid of these
        self.counts_df_with_host = self.counts_df_with_host.loc[:, (self.counts_df_with_host != 0).any(axis=0)]
        
        # Need to convert between didiers naming system and the coral sample-id format
        didier_to_sp_name_dict = {k:v for k, v in zip(poc_golden["INDV"].values, poc_golden["dataset_ITS2_sample_id"]) if not pd.isnull(v)}
        self.host_groupings = self.host_groupings.loc[[_ for _ in self.host_groupings.index.values if _ in didier_to_sp_name_dict],]
        foo = 'bar'
        
    def look_for_diagnostic_sqeuences(self):
        """ For every sequence in the ITS2 count table, look to see if it is found in every individual
            of a given group if it is, then count this sequence. Keep track of which groups a given
            sequence is found in and then we will output this info as a Venn or table. We will obviously
            be particularly interested in those sequences that are only found in one of the groups.
         """

HostDiagnosticZooxs()

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