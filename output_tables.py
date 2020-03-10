import os
import compress_pickle
import sys
from base_18s import EighteenSBase
import pandas as pd
from collections import defaultdict
import numpy as np

class EighteenSOutputTables(EighteenSBase):
    """ Class to produce the three output tables
    1 - Table that is the raw seq abundances (post mothur QC processing) in each of the samples (absolute abundance)
    2 - Table that holds the blast annotation for each of these sequences (Order, Family and Genus)
    3 - Table that is the post consolidation walk sequences from taxonomic origin of one of the
    three target genera. Absolute proportions in the samples.

    4 - Table that will hold meta information for each of the samples that will enable researchers to make a
    decision about whether they want to use given samples are not. The columns will be:

    status - This will be binary as use OR do_not_use

    primary_taxonomic_annotation - This will be either Pocrites, Millepora, Pocillopora or other.
        It will be the cummulatively most abundant of the categories in the sample.
        If this is 'other', it will cause the sample to have a do_not_use status.

    provenance_tax_annotation_is_correct - This will hold a boolean value. This is False if the primary_taxonomic annotation differs
        from the annotation in the tara provenance information file.

    is_inter_coral_contamination - This will be True if the collective abundance of the 'other_coral' sequences is
        greater than 0.01 of the sample. If True this will give a status of do_not_use

    inter_coral_contamination_rel_abund - Float. The actual relative abundance of the 'other_coral' sequences

    is_different_primary_sequence Boolean. Whether the most abundant sequence of the genus in
        question is different to the 'usual' most abundant sequence for the genus in question
        (considering all other samples). This will only apply to samples that have a primary_taxonomic_annotation
        as one of our three target genera. (~50 samples will have this true). Causes do_not_use.

    primary_sequence. The most abundant sequence of the sample that annotates as the primary_taxonomic_annotation.

    low_host_rel_abund. Boolean. If the relative abundance of the host genus in question sequences is
        below 30% as a proportion of all sequence in the sample, this will be true. Will cause do_not_use.

    host_rel_abund
        The ratio of the above.

    is_putative_intra_genus_contamination. Boolean. Whether there is the possibility of intra genus contamination.
        This covers the cases where we have samples that contain 2 abundant sequences from the same genus.
        It could of course be the case that this is due to intragenomic sequence diversity rather than intergenomic
        diversity. But to be on the safe side, I think it is a good idea to flag it up in this table.
        I will set the threshold ratio (ratio of the abundance of the first and second most abundant sequences)
        to be 0.3 for this category to be set to true. (< 10 samples will have this as True). Will cause do_not_use

    putative_intra_genus_contamination_ratio - Ratio that the above is based on


    """
    def __init__(self):
        super().__init__()
        self.info_df = compress_pickle.load(os.path.join(self.cache_dir, 'info_df_with_additional_info.p.bz'))

        # Produce the dictionary for making the coral meta info table
        # This will have sample name as key and a list in order of the df columns given in the comments
        self.primary_seq_dict = self._make_primary_seq_dict()
        self.coral_meta_info_table_dict = self._make_coral_meta_info_table_dict()

        # This will be a dict where full sequence is the key
        # The value will be another dict holding cumulative relative abundance
        # and the taxonomic tup
        self.master_seq_info_dict = {}
        self._populate_master_seq_info_dict()
        # Now get a master list of the sequences in order of cummulative abundance
        self.master_seq_abund_order_list = self._make_master_seq_abund_order_list()

        # Now popualte the dictionary that will be used to create the abundance df
        self.abundance_df_dict = {}
        self._populate_abundance_df_dict()

        self.tax_annotation_df_dict = self._populate_tax_annotation_df_dict()

        # Dict for collecting the sequencing information for making the host only consolidated sequences
        # absolute abundance dataframe
        self.host_only_master_seq_info_dict = self._populate_host_only_master_seq_info_dict()
        self.host_only_master_seq_abund_order_list = self._make_host_only_master_seq_abund_order_list()
        self.consolidated_df_dict = {}
        self._populated_consolidated_df_dict()


    def _make_primary_seq_dict(self):
        # Get the primary sequences for each of the genera
        primary_seq_dict = {}  # genus key to primary seq value
        for genus in ['Porites', 'Millepora', 'Pocillopora']:
            primary_seq_dict[genus] = self.info_df[
                self.info_df['most_abund_coral_genus'] == genus
                ]['most_abund_seq_of_coral_genus'].mode().values.tolist()[0]
        return primary_seq_dict

    def _make_coral_meta_info_table_dict(self):
        if os.path.isfile(os.path.join(self.cache_dir, 'coral_meta_info_table_dict.p.bzX')):
            return compress_pickle.load(os.path.join(self.cache_dir, 'coral_meta_info_table_dict.p.bz'))

        coral_meta_info_table_dict = {}
        print('Populating coral meta info table dict')
        for sample_name in self.info_df.index:
            sys.stdout.write(f'\r{sample_name}')
            # True until False
            use = True

            # Load up the dictionaries that we need
            sample_qc_dir = os.path.join(self.qc_dir, sample_name)
            coral_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'coral_annotation_dict.p.bz'))
            rel_all_seq_abundance_dict = compress_pickle.load(
                os.path.join(sample_qc_dir, 'rel_all_seq_abundance_dict.p.bz'))

            # Get a relative abund for each of the genera in the sample and the 'other' annotation
            # Use the most abundant tax as the primary annotation
            # Use the relative abundances to populate inter_coral_contamination_>_0.01 and
            # inter_coral_contamination_rel_abund.
            coral_tax_rel_count_dd = defaultdict(float)
            for coral_seq_name, tax_designation in coral_annotation_dict.items():
                coral_tax_rel_count_dd[coral_annotation_dict[coral_seq_name]] += rel_all_seq_abundance_dict[coral_seq_name]

            sorted_dict_keys = sorted(coral_tax_rel_count_dd, key=coral_tax_rel_count_dd.get, reverse=True)

            # Get primary taxonomic annotation and check to see if this matches the provenance annotation
            # If the primary annotation was other (i.e. not one of the three target genera) then log
            # False
            primary_taxonomic_annotation = sorted_dict_keys[0]
            provenance_annotation = self.info_df.at[sample_name, 'species']
            if primary_taxonomic_annotation == provenance_annotation:
                is_provenance_tax_annotation_correct = True
            else:
                is_provenance_tax_annotation_correct = False

            # Check to see whether the summed relative abundances of the coral sequences (other than the
            # sequences orginating from the primary_taxonomic_annotation are above 1% of the total coral sample
            inter_coral_contamination_rel_abund = sum([coral_tax_rel_count_dd[k] for k in sorted_dict_keys[1:]])
            if inter_coral_contamination_rel_abund > 0.01:
                is_inter_coral_contamination = True
                use = False
            else:
                is_inter_coral_contamination = False

            # Check to see if the sample has the primary sequences that the majority of samples of its genus have
            # if the sample annotates as other then make this true
            primary_sequence =  self.info_df.at[sample_name, 'most_abund_seq_of_coral_genus']
            try:
                if primary_sequence == self.primary_seq_dict[primary_taxonomic_annotation]:
                    is_different_primary_sequence = False
                else:
                    is_different_primary_sequence = True
                    use=False
            except KeyError:
                # If the primary taxonomic annotation is other
                use = False
                is_different_primary_sequence = True

            # Check to see what proportion the sample sequences the host sequences (of the primary_taxonomic annotation)
            # represent. If less than 0.3, mark low_host_rel_abund as True
            host_rel_abund = coral_tax_rel_count_dd[primary_taxonomic_annotation]
            if host_rel_abund < 0.3:
                is_low_host_rel_abund = True
                use = False
            else:
                is_low_host_rel_abund = False

            # Finally, to check for the putative intragenus contamination
            # To do this we will do a second pass through the coral_annotation_dict and get the two most abundant
            # sequences for the given primary_taxonomic_annotation and work out their ratio
            # If the primary taxaonomic annoatation is other, do the calculation anyway.
            abunds_of_genus = []
            for coral_seq_name, tax_designation in coral_annotation_dict.items():
                if tax_designation == primary_taxonomic_annotation:
                    abunds_of_genus.append(rel_all_seq_abundance_dict[coral_seq_name])
            # sort the abunds to get the top two and calculate a ratio
            sorted_abunds = sorted(abunds_of_genus, reverse=True)
            if len(sorted_abunds) > 1:
                putative_intra_genus_contamination_ratio = sorted_abunds[1]/sorted_abunds[0]
                if putative_intra_genus_contamination_ratio > 0.3:
                    is_putative_intra_genus_contamination = True
                    use = False
                else:
                    is_putative_intra_genus_contamination = False
            else:
                is_putative_intra_genus_contamination = False
                putative_intra_genus_contamination_ratio = 0


            # Add in the island and site information as well because this is useful to have
            island = self.info_df.at[sample_name, 'island']
            site = self.info_df.at[sample_name, 'site']
            # Get the individual
            if sample_name[-2] == '_':
                # then this is a techrep
                sample_base = '_'.join(sample_name.split('_')[:-1])
                sample_base_id = self.sample_provenance_df.at[sample_base, 'C###, F###, MA##, SG##'].replace('C0', '')
                individual = sample_base_id + '_' + sample_name.split('_')[-1]
            else:
                individual = self.sample_provenance_df.at[sample_name, 'C###, F###, MA##, SG##'].replace('C0', '')

            coral_meta_info_table_dict[sample_name] = [
                use, primary_taxonomic_annotation, is_provenance_tax_annotation_correct,
                is_inter_coral_contamination, inter_coral_contamination_rel_abund,
                is_different_primary_sequence, primary_sequence,
                is_low_host_rel_abund, host_rel_abund,
                is_putative_intra_genus_contamination, putative_intra_genus_contamination_ratio,
                island, site, individual
            ]
        compress_pickle.dump(coral_meta_info_table_dict, os.path.join(self.cache_dir, 'coral_meta_info_table_dict.p.bz'))
        return coral_meta_info_table_dict

    def make_and_write_coral_meta_info_output_table(self):
        print('Constructing coral meta info output table')
        df = pd.DataFrame.from_dict(
            self.coral_meta_info_table_dict,
            orient='index',
            columns=['status', 'primary_taxonomic_annotation', 'is_provenance_tax_annotation_correct',
                'is_inter_coral_contamination', 'inter_coral_contamination_rel_abund',
                'is_different_primary_sequence', 'primary_sequence',
                'is_low_host_rel_abund', 'host_rel_abund',
                'is_putative_intra_genus_contamination', 'putative_intra_genus_contamination_ratio',
                'island', 'site', 'individual']
        )
        print('Writing coral meta info output table')
        df.to_csv(os.path.join(self.output_dir, 'coral_meta_info_table.csv.bz'), index=True, compression='bz2')

    def make_and_write_raw_abund_output_table(self):
        # Here we have the self.abundance_df_dict populated and we can now create the dataframe from this dict
        print('Constructing raw abundance table')
        df = pd.DataFrame.from_dict(
            self.abundance_df_dict,
            orient='index',
            columns=self.master_seq_abund_order_list
        )
        print('Writing raw abundance table')
        df.to_csv(os.path.join(self.output_dir, 'raw_seq_abund.csv.bz'), index=True, compression='bz2')

    def make_and_write_tax_output_table(self):
        print('Constructing taxonomy table')
        df = pd.DataFrame.from_dict(
            self.abundance_df_dict,
            orient='index',
            columns=self.master_seq_abund_order_list
        )
        print('Writing taxonomy table')
        df.to_csv(os.path.join(self.output_dir, 'tax_annotation.csv.bz'), index=True, compression='bz2')

    def make_and_write_consolidated_host_output_table(self):
        print('Constructing consolidated host output table')
        df = pd.DataFrame.from_dict(
            self.consolidated_df_dict, orient='index',
            columns=self.host_only_master_seq_abund_order_list
        )
        print('Writing consolidated host output table')
        df.to_csv(os.path.join(self.output_dir, 'consolidated_host.csv.bz'), index=True, compression='bz2')

    def _populate_tax_annotation_df_dict(self):
        # The key should be sequence and the list should be order, family genus, in that order
        tax_annotation_df_dict = {}
        for seq, info_dict in self.master_seq_info_dict.items():
            annotation_tup = info_dict['tax_annotation']
            tax_annotation_df_dict[seq] = [annotation_tup[2], annotation_tup[1], annotation_tup[0]]
        return tax_annotation_df_dict

    def _populate_abundance_df_dict(self):
        if os.path.isfile(os.path.join(self.cache_dir, 'abundance_df_table_output.p.bz')):
            self.abundance_df_dict = compress_pickle.load(os.path.join(self.cache_dir, 'abundance_df_table_output.p.bz'))
        else:
            print('Collecting sequence information (original_seqs): second pass')
            for sample_name in self.info_df.index:
                sys.stdout.write(f'\r{sample_name}')
                # read in the fasta file
                fasta_file = self._read_in_fasta_file(sample_name)

                # read in the name file and make an abundance dictionary
                name_abs_abund_dict = self._make_abs_abund_dict_from_names_path(sample_name)

                # create a seq_to_abs_abund dictionary
                seq_to_abs_abund_dict = {
                    fasta_file[i+1]: name_abs_abund_dict[fasta_file[i].split('\t')[0][1:]] for
                    i in range(0, len(fasta_file), 2)}

                # In the order of the self.master_seq_abund_order_list
                # populate the abundances for the given sample
                temp_abund_list = [seq_to_abs_abund_dict[seq] if seq in seq_to_abs_abund_dict else 0 for seq in self.master_seq_abund_order_list]
                self.abundance_df_dict[sample_name] = temp_abund_list
            compress_pickle.dump(self.abundance_df_dict, os.path.join(self.cache_dir, 'abundance_df_table_output.p.bz'))

    def _populate_host_only_master_seq_info_dict(self):
        # To populate the consolidated count table that is just the
        # sequences found in a sample that are the same genus as the most abundant host sequence
        # We should do similar to before and do two passes
        # First pass we can just use the consolidated_host_seqs_abund_dict that is found in every
        # sample directory. And collect this into the master collection dict
        # For the second pass we need to get an absolute abundance. This is a little trickier.
        # We basically need to create a reverse dictionary. Each consolidated sequence may be representative
        # of several original sequence. We will go through each fasta to see if it has a consolidated sequence
        # representative. If it does or if it doesn't add this to the reverse dict that will be a default dict
        # Once we have the revese dict populated we can go seq by seq in the consolidated_host_seqs_abund_dict
        # and summate the original seqs that the given consolidated seq represents and add this information
        # to the list in the order of the master ordered consolidated sequences.
        if os.path.isfile(os.path.join(self.cache_dir, 'host_only_master_seq_info_dict.p.bz')):
            return compress_pickle.load(os.path.join(self.cache_dir, 'host_only_master_seq_info_dict.p.bz'))
        host_only_master_seq_info_dict = defaultdict(float)
        print('Collecting sequence information (consolidated seqs): first pass')
        for sample_name in self.info_df.index:
            sys.stdout.write(f'\r{sample_name}')

            # Dict that is sequence key to relative abundance in the sample (of only the given genus sequences)
            # I.e. dict adds up to one. We will use this only for the keys
            # to see which seqs we should be concerned with
            consolidated_host_seqs_abund_dict = compress_pickle.load(
                os.path.join(self.qc_dir, sample_name, 'consolidated_host_seqs_abund_dict.p.bz'))

            for seq, rel_abund in consolidated_host_seqs_abund_dict.items():
                host_only_master_seq_info_dict[seq] += rel_abund
        compress_pickle.dump(host_only_master_seq_info_dict, os.path.join(self.cache_dir, 'host_only_master_seq_info_dict.p.bz'))
        return host_only_master_seq_info_dict

    def _populated_consolidated_df_dict(self):
        # For the second pass we need to get an absolute abundance. This is a little trickier.
        # We basically need to create a reverse dictionary. Each consolidated sequence may be representative
        # of several original sequence. We will go through each fasta to see if it has a consolidated sequence
        # representative. If it does or if it doesn't add this to the reverse dict that will be a default dict
        # Once we have the revese dict populated we can go seq by seq in the consolidated_host_seqs_abund_dict
        # and summate the original seqs that the given consolidated seq represents and add this information
        # to the list in the order of the master ordered consolidated sequences.
        if os.path.isfile(os.path.join(self.cache_dir, 'consolidated_df_dict_output_tables.p.bz')):
            self.consolidated_df_dict = compress_pickle.load(os.path.join(self.cache_dir, 'consolidated_df_dict_output_tables.p.bz'))
        else:
            print('Collecting sequence information (consolidated seqs): second pass')
            coral_blasted_seq_to_consolidated_seq_dict = compress_pickle.load(
                os.path.join(self.cache_dir, 'coral_blasted_seq_to_consolidated_seq_dict.p.bz'))
            for sample_name in self.info_df.index:
                sys.stdout.write(f'\r{sample_name}')
                # read in the fasta file
                fasta_file = self._read_in_fasta_file(sample_name)
                fasta_seq_to_name_dict = {fasta_file[i+1]: fasta_file[i].split('\t')[0][1:] for i in range(0, len(fasta_file), 2)}
                # read in the name file and make an abundance dictionary
                name_abs_abund_dict = self._make_abs_abund_dict_from_names_path(sample_name)

                # the dict we are making for each sample that maps consolidated sequence to the
                # original sequences it represents
                consol_seq_to_orig_seq_dict = defaultdict(list)

                for i in range(0, len(fasta_file), 2):
                    seq = fasta_file[i+1]
                    try:
                        # If there is a representative seq for this seq, then log it
                        consol_seq_to_orig_seq_dict[coral_blasted_seq_to_consolidated_seq_dict[seq]].append(seq)
                    except KeyError:
                        # else just map it to its self in the dict
                        consol_seq_to_orig_seq_dict[seq].append(seq)

                consolidated_host_seqs_abund_dict = compress_pickle.load(
                    os.path.join(self.qc_dir, sample_name, 'consolidated_host_seqs_abund_dict.p.bz'))

                temp_abund_list = []
                for master_consolidated_seq in self.host_only_master_seq_abund_order_list:
                    if master_consolidated_seq in consolidated_host_seqs_abund_dict:
                        temp_abund_list.append(sum([name_abs_abund_dict[fasta_seq_to_name_dict[repped_seq]] for repped_seq in consol_seq_to_orig_seq_dict[master_consolidated_seq]]))
                    else:
                        temp_abund_list.append(0)
                self.consolidated_df_dict[sample_name] = temp_abund_list
            compress_pickle.dump(self.consolidated_df_dict, os.path.join(self.cache_dir, 'consolidated_df_dict_output_tables.p.bz'))

    def _make_host_only_master_seq_abund_order_list(self):
        return [tup[0] for tup in sorted([_ for _ in self.host_only_master_seq_info_dict.items()],
                                                                     key=lambda x: x[1],
                                                                     reverse=True)]

    def _make_master_seq_abund_order_list(self):
        return [tup[0] for tup in sorted([_ for _ in self.master_seq_info_dict.items()],
                                                                     key=lambda x: x[1]['cummulative_abund'],
                                                                     reverse=True)]
    def _populate_master_seq_info_dict(self):
        if os.path.isfile(os.path.join(self.cache_dir, 'master_seq_info_dict.p.bz')):
            self.master_seq_info_dict = compress_pickle.load(os.path.join(self.cache_dir, 'master_seq_info_dict.p.bz'))
        else:
            print('Collecting sequence information (original seqs): first pass')
            for sample_name in self.info_df.index:
                sys.stdout.write(f'\r{sample_name}')
                # read in the fasta file
                fasta_file = self._read_in_fasta_file(sample_name)

                # read in the name file and make an abundance dictionary
                name_rel_abund_dict = self._make_rel_abund_dict_from_names_path(sample_name)

                # read in the sample taxonomy dictionary
                sample_annotation_dict = compress_pickle.load(
                    os.path.join(self.qc_dir, sample_name, 'sample_annotation_dict.p.bz'))

                # for each sequence in the fasta file
                # if not already in the dict, init with the rel abund and tax info
                # else simply add to the cumulative abund
                for i in range(0, len(fasta_file), 2):
                    seq_name = fasta_file[i].split('\t')[0][1:]
                    seq = fasta_file[i+1]
                    try:
                        self.master_seq_info_dict[seq]['cummulative_abund'] += name_rel_abund_dict[seq_name]
                    except KeyError:
                        try:
                            tax_tup = sample_annotation_dict[seq_name]
                        except KeyError:
                            tax_tup = ('not_annotated', 'not_annotated', 'not_annotated')
                        self.master_seq_info_dict[seq] = {'cummulative_abund': name_rel_abund_dict[seq_name],
                                                                 'tax_annotation': tax_tup}
            compress_pickle.dump(self.master_seq_info_dict, os.path.join(self.cache_dir, 'master_seq_info_dict.p.bz'))

    def _read_in_fasta_file(self, sample_name):
        with open(
                os.path.join(self.qc_dir, sample_name, 'stability.trim.contigs.good.unique.abund.pcr.unique.fasta'),
                'r') as f:
            return [line.rstrip() for line in f]

    def _make_abs_abund_dict_from_names_path(self, sample_name):
        with open(os.path.join(self.qc_dir, sample_name, 'stability.trim.contigs.good.unique.abund.pcr.names'), 'r') as f:
            name_file = [line.rstrip() for line in f]
        return  {line.split('\t')[0]: len(line.split('\t')[1].split(',')) for line in name_file}

    def _make_rel_abund_dict_from_names_path(self, sample_name):
        abs_abund_dict = self._make_abs_abund_dict_from_names_path(sample_name)
        tot = sum(abs_abund_dict.values())
        return {seq_name: abund/tot for seq_name, abund in abs_abund_dict.items()}

if __name__ == "__main__":
    ot = EighteenSOutputTables()
    ot.make_and_write_coral_meta_info_output_table()
    ot.make_and_write_raw_abund_output_table()
    ot.make_and_write_tax_output_table()
    ot.make_and_write_consolidated_host_output_table()
