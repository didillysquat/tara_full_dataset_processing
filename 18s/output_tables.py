# Then we can continue working on the output_tables.py.
# If we don't do this now. Its just going to be confusing in the future.
import os
import compress_pickle
import sys
from base_18s import EighteenSBase
import pandas as pd
from collections import defaultdict
import numpy as np
import datetime

class EighteenSOutputTables(EighteenSBase):
    """ Class to produce the four output tables
    1 - Table that is the raw seq abundances (post mothur QC processing) in each of the samples (absolute abundance)
    2 - Table that holds the blast annotation for each of these sequences (Order, Family and Genus)
    3 - Table that is the post consolidation walk sequences from taxonomic origin of one of the
    three target genera. Absolute proportions in the samples.
    4 - Table that will hold meta information for each of the samples that will enable researchers to make a
    decision about whether they want to use given samples are not. The columns will be:

    # TODO
    sample-id: This is the barcode_id with TARA_ prepended to it.

    #TODO change this to 'use' as a Bool
    status - This will be binary as use OR do_not_use

    # TODO Change name but also check to see if this is accurate. I think we give the exact genus name
    genetic_18S_genus_taxonomic_annotation - This will be either Pocrites, Millepora, Pocillopora or other.
        It will be the cummulatively most abundant of the categories in the sample.
        If this is 'other', it will cause the sample to have a do_not_use status.

    provenance_tax_annotation_is_correct - This will hold a boolean value. This is False if the primary_taxonomic annotation differs
        from the annotation in the tara provenance information file.

    is_inter_coral_contamination - This will be True if the collective abundance of the 'other_coral' sequences is
        greater than 0.01 of the sample. If True this will give a status of do_not_use

    inter_coral_contamination_rel_abund - Float. The actual relative abundance of the 'other_coral' sequences

    is_different_primary_sequence Boolean. Whether the most abundant sequence of the genus in
        question is different to the 'usual' most abundant sequence for the genus in question
        (considering all other samples). This will only apply to samples that have a genus_18S_taxonomic_annotation
        as one of our three target genera. (~50 samples will have this true). Causes do_not_use.

    primary_sequence. The most abundant sequence of the sample that annotates as the genus_18S_taxonomic_annotation.

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
        # Because there are both coral and non-coral 18s seq pairs in the fastq_info_df, it will be useful
        # to have a list that is just the coral readsets from the table. I'll make this in the base 18s
        # at the same time I make the fastq_info_df.
        
        # For the output tables, on a sample by sample basis we want to have access to what the most abundant
        # coral genus was for a given sample, and what the most abundant sequence was for that genus.
        # For every readset in the fastq_info_df we will gather this information and populate it in a new
        # df called abundance_info_df. This df may later contain addional information such as the QC info.
        # For now we will also add a column called post_qc_seq_depth. This will be the total number of 
        # sequences returned for a given readset.
        self.abundance_info_df = self._make_or_get_abundance_info_df()

        # Produce the dictionary for making the coral meta info table
        # This will have sample name as key and a list in order of the df columns given in the comments
        self.primary_seq_dict = self._make_primary_seq_dict()
        self.coral_meta_info_table_dict = {}
        #TODO We need to incorporate the changes requested by Stephane. These are largely
        # aimed at keeping the same terms that are given in the provenance table
        # in our output table. This is a good point.
        # We had originally hoped to be able to get rid of technical replicates but we will need to include
        # these as there are still technical replicates hosted on the genoscope server.
        self._populate_coral_meta_info_table_dict()

        # TODO move back to being a public method run independently
        self.make_and_write_coral_meta_info_output_table()

        # # This will be a dict where full sequence is the key
        # # The value will be another dict holding cumulative relative abundance
        # # and the taxonomic tup
        # self.master_seq_info_dict = {}
        # self._populate_master_seq_info_dict()
        # # Now get a master list of the sequences in order of cummulative abundance
        # self.master_seq_abund_order_list = self._make_master_seq_abund_order_list()

        # # Now popualte the dictionary that will be used to create the abundance df
        # self.abundance_df_dict = {}
        # self._populate_abundance_df_dict()

        # self.tax_annotation_df_dict = self._populate_tax_annotation_df_dict()

        # # Dict for collecting the sequencing information for making the host only consolidated sequences
        # # absolute abundance dataframe
        # self.host_only_master_seq_info_dict = self._populate_host_only_master_seq_info_dict()
        # self.host_only_master_seq_abund_order_list = self._make_host_only_master_seq_abund_order_list()
        # self.consolidated_df_dict = {}
        # self._populated_consolidated_df_dict()

    def _make_or_get_abundance_info_df(self):
        if os.path.isfile(os.path.join(self.cache_dir, 'abundance_info_df.p.bz')):
            return compress_pickle.load(os.path.join(self.cache_dir, 'abundance_info_df.p.bz'))
        else:
            print('Building abundance_info_df')
            abundance_df_dict = {}
            for readset in self.coral_readsets:
                sys.stdout.write(f'\r{readset}')
                sample_qc_dir = os.path.join(self.qc_dir, readset)
                rel_all_seq_abundance_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'rel_all_seq_abundance_dict.p.bz'))
                coral_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'coral_annotation_dict.p.bz'))
                most_abund_coral_genus = self._identify_most_abund_coral_genus(rel_all_seq_abundance_dict, coral_annotation_dict)
                consolidated_host_seqs_abund_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
                most_abund_seq_of_coral_genus = sorted([_ for _ in consolidated_host_seqs_abund_dict.items()], key=lambda x: x[1], reverse=True)[0][0]
                absolute_seqs_count = sum(self._make_abs_abund_dict_from_names_path(readset).values())
                abundance_df_dict[readset] = [most_abund_coral_genus, most_abund_seq_of_coral_genus, absolute_seqs_count]
            df = pd.DataFrame.from_dict(data=abundance_df_dict, orient='index', columns=['most_abund_coral_genus', 'most_abund_seq_of_coral_genus', 'post_qc_seq_depth'])
            compress_pickle.dump(df, os.path.join(self.cache_dir, 'abundance_info_df.p.bz'))
            return df

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

    def _make_primary_seq_dict(self):
        # Get the primary sequences for each of the genera
        # This is the sequence that is found as the most abundant for the largest number of readsets
        # for a given genus.
        primary_seq_dict = {}  # genus key to primary seq value
        for genus in ['Porites', 'Millepora', 'Pocillopora']:
            primary_seq_dict[genus] = self.abundance_info_df[
                self.abundance_info_df['most_abund_coral_genus'] == genus
                ]['most_abund_seq_of_coral_genus'].mode().values.tolist()[0]
        return primary_seq_dict

    def _populate_coral_meta_info_table_dict(self):
        if os.path.isfile(os.path.join(self.cache_dir, 'coral_18S_meta_info_table_dict.p.bz')):
            return compress_pickle.load(os.path.join(self.cache_dir, 'coral_18S_meta_info_table_dict.p.bz'))

        print('Populating coral meta info table dict')
        tot = len(self.coral_readsets)
        count = 0
        for readset in self.coral_readsets:
            sys.stdout.write(f'\r{readset}: {count}/{tot}')
            tbfour = TableFour(parent=self, readset=readset)
            tbfour.populate_coral_meta_info_table_dict()
            count += 1
        compress_pickle.dump(self.coral_meta_info_table_dict, os.path.join(self.cache_dir, 'coral_18S_meta_info_table_dict.p.bz'))

    def make_and_write_coral_meta_info_output_table(self):
        print('Constructing coral meta info output table')
        column_order = [
                'sample-id', 'use', 'genetic_18S_genus_taxonomic_annotation', 
                'Sample Material label, organismal system level, taxonomic, nominal',
                'Sample Material label, organismal system level, taxonomic, label' ,
                'is_provenance_taxonomic_annotation_correct',
                'inter_coral_contamination_rel_abund', 'is_inter_coral_contamination',
                'primary_sequence', 'is_different_primary_sequence',
                'host_rel_abund', 'is_low_host_rel_abund',
                'putative_intra_genus_contamination_ratio', 'is_putative_intra_genus_contamination',
                'is_representative_for_sample', 'post_qc_seq_depth',
                'fwd_read_name', 'rev_read_name',
                'SAMPLING DESIGN LABEL',
                'ISLAND#', 'SITE#', 'COLONY# (C000) FISH# (F000) MACROALGAE# (MA00)',
                'SAMPLE PROTOCOL LABEL, level 1', 'SAMPLE PROTOCOL LABEL, level 2',
                'Sample Material label - trait'
                ]
        df = pd.DataFrame.from_dict(
            self.coral_meta_info_table_dict,
            orient='index',
            columns=column_order
        )
        # Now we need to move the index to a column
        df.index.name = 'not_readset'
        df['readset'] = list(df.index)
        # swap the readset column with the sample-id column
        column_order.insert(1, 'readset')
        df = df.reindex(columns=column_order)
        print('Writing coral meta info output table')
        dat_string = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f%Z").replace(':', '_')
        df.to_csv(os.path.join(self.output_dir, f'coral_18S_meta_info_table_{dat_string}.csv.bz'), index=False, compression='bz2')

    def make_and_write_raw_abund_output_table(self):
        # Here we have the self.abundance_df_dict populated and we can now create the dataframe from this dict
        print('Constructing raw abundance table')
        df = pd.DataFrame.from_dict(
            self.abundance_df_dict,
            orient='index',
            columns=self.master_seq_abund_order_list
        )
        print('Writing raw abundance table')
        df.to_csv(os.path.join(self.output_dir, 'raw_seq_abund.csv.bz'), index=True, compression='bz2', index_label='sample_id')

    def make_and_write_tax_output_table(self):
        print('Constructing taxonomy table')
        df = pd.DataFrame.from_dict(
            self.tax_annotation_df_dict,
            orient='index',
            columns=['order', 'family', 'genus']
        )
        print('Writing taxonomy table')
        df.to_csv(os.path.join(self.output_dir, 'tax_annotation.csv.bz'), index=True, compression='bz2', index_label='sequence')

    def make_and_write_consolidated_host_output_table(self):
        print('Constructing consolidated host output table')
        df = pd.DataFrame.from_dict(
            self.consolidated_df_dict, orient='index',
            columns=self.host_only_master_seq_abund_order_list
        )
        print('Writing consolidated host output table')
        df.to_csv(os.path.join(self.output_dir, 'consolidated_host.csv.bz'), index=True, compression='bz2', index_label='sample_id')

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
            for sample_id in self.info_df.index:
                sys.stdout.write(f'\r{sample_id}')
                # read in the fasta file
                fasta_file = self._read_in_fasta_file(sample_id)

                # read in the name file and make an abundance dictionary
                name_abs_abund_dict = self._make_abs_abund_dict_from_names_path(sample_id)

                # create a seq_to_abs_abund dictionary
                seq_to_abs_abund_dict = {
                    fasta_file[i+1]: name_abs_abund_dict[fasta_file[i].split('\t')[0][1:]] for
                    i in range(0, len(fasta_file), 2)}

                # In the order of the self.master_seq_abund_order_list
                # populate the abundances for the given sample
                temp_abund_list = [seq_to_abs_abund_dict[seq] if seq in seq_to_abs_abund_dict else 0 for seq in self.master_seq_abund_order_list]
                self.abundance_df_dict[sample_id] = temp_abund_list
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
        for sample_id in self.info_df.index:
            sys.stdout.write(f'\r{sample_id}')

            # Dict that is sequence key to relative abundance in the sample (of only the given genus sequences)
            # I.e. dict adds up to one. We will use this only for the keys
            # to see which seqs we should be concerned with
            consolidated_host_seqs_abund_dict = compress_pickle.load(
                os.path.join(self.qc_dir, sample_id, 'consolidated_host_seqs_abund_dict.p.bz'))

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
            for sample_id in self.info_df.index:
                sys.stdout.write(f'\r{sample_id}')
                # read in the fasta file
                fasta_file = self._read_in_fasta_file(sample_id)
                fasta_seq_to_name_dict = {fasta_file[i+1]: fasta_file[i].split('\t')[0][1:] for i in range(0, len(fasta_file), 2)}
                # read in the name file and make an abundance dictionary
                name_abs_abund_dict = self._make_abs_abund_dict_from_names_path(sample_id)

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
                    os.path.join(self.qc_dir, sample_id, 'consolidated_host_seqs_abund_dict.p.bz'))

                temp_abund_list = []
                for master_consolidated_seq in self.host_only_master_seq_abund_order_list:
                    if master_consolidated_seq in consolidated_host_seqs_abund_dict:
                        temp_abund_list.append(sum([name_abs_abund_dict[fasta_seq_to_name_dict[repped_seq]] for repped_seq in consol_seq_to_orig_seq_dict[master_consolidated_seq]]))
                    else:
                        temp_abund_list.append(0)
                self.consolidated_df_dict[sample_id] = temp_abund_list
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
            for sample_id in self.info_df.index:
                sys.stdout.write(f'\r{sample_id}')
                # read in the fasta file
                fasta_file = self._read_in_fasta_file(sample_id)

                # read in the name file and make an abundance dictionary
                name_rel_abund_dict = self._make_rel_abund_dict_from_names_path(sample_id)

                # read in the sample taxonomy dictionary
                sample_annotation_dict = compress_pickle.load(
                    os.path.join(self.qc_dir, sample_id, 'sample_annotation_dict.p.bz'))

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

    def _read_in_fasta_file(self, sample_id):
        with open(
                os.path.join(self.qc_dir, sample_id, 'stability.trim.contigs.good.unique.abund.pcr.unique.fasta'),
                'r') as f:
            return [line.rstrip() for line in f]

    def _make_abs_abund_dict_from_names_path(self, readset):
        name_path = os.path.join(self.qc_dir, readset, 'stability.trim.contigs.good.unique.abund.pcr.names')
        name_file = EighteenSBase.decompress_read_compress(name_path)
        return  {line.split('\t')[0]: len(line.split('\t')[1].split(',')) for line in name_file}

    def _make_rel_abund_dict_from_names_path(self, readset):
        abs_abund_dict = self._make_abs_abund_dict_from_names_path(readset)
        tot = sum(abs_abund_dict.values())
        return {seq_name: abund/tot for seq_name, abund in abs_abund_dict.items()}

class TableFour():
    def __init__(self, parent, readset):
        # TODO here we will not assign all of the information we want from the provenance table
        # to a variable. We will only assign those things that we need to make our host-related
        # columns from. I think this is just the provenance_annotation.
        self.parent = parent
        self.readset = readset
        self.sample_id = self.parent.fastq_info_df.at[readset, 'sample-id']
        self.status = True
        self.sample_qc_dir = os.path.join(self.parent.qc_dir, readset)
        self.coral_annotation_dict = compress_pickle.load(os.path.join(self.sample_qc_dir, 'coral_annotation_dict.p.bz'))
        self.consolidated_host_seqs_abund_dict = compress_pickle.load(
            os.path.join(self.sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
        self.rel_all_seq_abundance_dict = compress_pickle.load(
            os.path.join(self.sample_qc_dir, 'rel_all_seq_abundance_dict.p.bz'))
        self.coral_tax_rel_count_dd = self._make_coral_tax_rel_count_dd()
        self.sorted_coral_tax_dict_keys = sorted(self.coral_tax_rel_count_dd, key=self.coral_tax_rel_count_dd.get, reverse=True)
        self.genus_18S_taxonomic_annotation = self.sorted_coral_tax_dict_keys[0]
        self.provenance_annotation = self.parent.sample_provenance_df.at[self.sample_id, 'Sample Material label, organismal system level, taxonomic, nominal']

        # The remainder of the variables that we need to populate
        self.is_provenance_tax_annotation_correct = None
        self.inter_coral_contamination_rel_abund = None
        self.is_inter_coral_contamination = None
        self.primary_sequence = None
        self.host_rel_abund = None
        self.putative_intra_genus_contamination_ratio = None
        self.is_putative_intra_genus_contamination = None
        self.is_representative_for_sample = None
        self.post_qc_seq_depth = None

        # Variables that are only associated with processing a Heliopora samples
        self.sample_annotation_dict = None
        self.fasta_dict = None
        self.all_tax_count_dd = None

    def _make_coral_tax_rel_count_dd(self):
        coral_tax_rel_count_dd = defaultdict(float)
        for coral_seq_name, tax_designation in self.coral_annotation_dict.items():
            coral_tax_rel_count_dd[tax_designation] += self.rel_all_seq_abundance_dict[coral_seq_name]
        return coral_tax_rel_count_dd

    def _pop_for_normal_sample(self):
        self.genus_18S_taxonomic_annotation = self.sorted_coral_tax_dict_keys[0]
        self._set_is_provenance_tax_annotation_correct()
        if self.genus_18S_taxonomic_annotation not in ["Porites", "Millepora", "Pocillopora"]:
            self.status = False

        self._set_inter_coral_contamination()

        self._set_primary_sequence()

        self._set_host_rel_abund_normal()

        self._set_intragenus()

        # self._set_island_site_individual()

        self._set_post_qc_seq_depth()

        self._set_is_representative_for_sample()

        self._populate_coral_meta_info_table_dict()

    def _set_post_qc_seq_depth(self):
        self.post_qc_seq_depth = self.parent.abundance_info_df.at[self.readset, 'post_qc_seq_depth']

    def _set_is_representative_for_sample(self):
        """Some of the samples contained multiple pairs of fastq.gz files. We have kept these pairs separate.
        But we want to be able to filter the output count table so that we can only consider one fastq.gz pair
        per coral sample. To facilitate this we will have a is_representative_of_sample field.
        This field will be True False.
        For samples that  contain no tech replicates, the sample will automatically be True. If the sample
        is a sample that contains sample replicates then we will use the fastq.gz pair that produced the greatest
        aboslute number of post-QC sequences (i.e. read depth) as the singular representative for the sample. To
        do this we will set this fastq.gz pair to True in the is_representative_of_sample. For the
        other fastq.gz pairs is_representative_of_sample will be False"""
        if self.parent.fastq_info_df.at[self.readset, 'is_replicate'] == True:
            # Get a list of the sample names that share this base name and make a dict of their abundances
            abs_count_dict = {}
            for readset, ser in self.parent.fastq_info_df.iterrows():
                if ser['sample-id'] == self.sample_id:
                    abs_count_dict[readset] = self.parent.abundance_info_df.at[readset, 'post_qc_seq_depth']
            # Here we have the abs_count_dict populated
            # Now sort it by key and check to see if the current sample_id matches the most abundant
            if self.readset == sorted(abs_count_dict, key=abs_count_dict.get, reverse=True)[0]:
                self.is_representative_for_sample = True
            else:
                self.is_representative_for_sample = False
        elif self.parent.fastq_info_df.at[self.readset, 'is_replicate'] == False:
            self.is_representative_for_sample = True

    def _populate_coral_meta_info_table_dict(self):
        # This is where we do the population of the row
        # We use the host-related values that we've just collected for our columns
        # Then we use columns directly from the output info table and from the provenance table
        # for the other columns.
        # We should have sample_id as the first item to keep people happy
        # but readset that is more useful should be second
        # However, because this is a dictionary we need to use something unique so we will use
        # the readset.
        self.parent.coral_meta_info_table_dict[self.readset] = [
            self.sample_id,
            self.status,
            self.genus_18S_taxonomic_annotation,
            #provenance_annotation
            self.parent.sample_provenance_df.at[self.sample_id, 'Sample Material label, organismal system level, taxonomic, nominal'],
            self.parent.sample_provenance_df.at[self.sample_id, 'Sample Material label, organismal system level, taxonomic, label'],
            self.is_provenance_tax_annotation_correct,
            self.inter_coral_contamination_rel_abund, self.is_inter_coral_contamination,
            self.primary_sequence, self.is_different_primary_sequence,
            self.host_rel_abund, self.is_low_host_rel_abund,
            self.putative_intra_genus_contamination_ratio, self.is_putative_intra_genus_contamination,
            self.is_representative_for_sample, self.post_qc_seq_depth,
            self.parent.fastq_info_df.at[self.readset, 'fwd_read_name'],
            self.parent.fastq_info_df.at[self.readset, 'rev_read_name'],
            self.parent.sample_provenance_df.at[self.sample_id, 'SAMPLING DESIGN LABEL'],
            # Concatenated version of the island, site, individual
            self.parent.sample_provenance_df.at[self.sample_id, 'ISLAND#'],
            self.parent.sample_provenance_df.at[self.sample_id, 'SITE#'],
            self.parent.sample_provenance_df.at[self.sample_id, 'COLONY# (C000) FISH# (F000) MACROALGAE# (MA00)'],
            # Sample protocol labels
            self.parent.sample_provenance_df.at[self.sample_id, 'SAMPLE PROTOCOL LABEL, level 1'],
            self.parent.sample_provenance_df.at[self.sample_id, 'SAMPLE PROTOCOL LABEL, level 2'],
            # Trait
            self.parent.sample_provenance_df.at[self.sample_id, 'Sample Material label - trait']
        ]

    # def _set_island_site_individual(self):
    #     # Add in the island and site information as well because this is useful to have
    #     self.island = self.parent.fastq_info_df.at[self.sample_id, 'island']
    #     self.site = self.parent.info_df.at[self.sample_id, 'site']
    #     # Get the individual
    #     if self.sample_id[-2] == '_':
    #         # then this is a techrep
    #         sample_base = '_'.join(self.sample_id.split('_')[:-1])
    #         sample_base_id = int(
    #             self.parent.sample_provenance_df.at[
    #                 sample_base, 'COLONY# (C000) FISH# (F000) MACROALGAE# (MA00)'].replace('C', ''))
    #         self.individual = str(sample_base_id) + '_' + self.sample_id.split('_')[-1]
    #     else:
    #         self.individual = str(int(
    #             self.parent.sample_provenance_df.at[
    #                 self.sample_id, 'COLONY# (C000) FISH# (F000) MACROALGAE# (MA00)'].replace('C', '')))

    def _set_intragenus(self):
        # Finally, to check for the putative intragenus contamination
        # This needs to be done using the consolidation dictionary
        # Simply sort keys and look for the two most abundant sequences
        # sort the abunds to get the top two and calculate a ratio
        sorted_consolidated_abunds = sorted(self.consolidated_host_seqs_abund_dict.values(), reverse=True)
        if len(sorted_consolidated_abunds) > 1:
            self.putative_intra_genus_contamination_ratio = sorted_consolidated_abunds[1] / sorted_consolidated_abunds[
                0]
            if self.putative_intra_genus_contamination_ratio > 0.3:
                self.is_putative_intra_genus_contamination = True
                self.status = False
            else:
                self.is_putative_intra_genus_contamination = False
        else:
            self.is_putative_intra_genus_contamination = False
            self.putative_intra_genus_contamination_ratio = 0

    def _set_host_rel_abund_normal(self):
        # Check to see what proportion the sample sequences the host sequences (of the primary_taxonomic annotation)
        # represent. If less than 0.3, mark low_host_rel_abund as True
        self.host_rel_abund = self.coral_tax_rel_count_dd[self.genus_18S_taxonomic_annotation]
        if self.host_rel_abund < 0.3:
            self.is_low_host_rel_abund = True
            self.status = False
        else:
            self.is_low_host_rel_abund = False

    def _set_host_rel_abund_heliopora(self):
        # Check to see what proportion the sample sequences the host sequences (of the primary_taxonomic annotation)
        # represent. If less than 0.3, mark low_host_rel_abund as True
        self.host_rel_abund = self.all_tax_count_dd[self.genus_18S_taxonomic_annotation]
        if self.host_rel_abund < 0.3:
            self.is_low_host_rel_abund = True
            self.status = False
        else:
            self.is_low_host_rel_abund = False

    def _set_primary_sequence(self):
        # Check to see if the sample has the primary sequences that the majority of samples of its genus have
        # if the sample annotates as other then make this true
        self.primary_sequence = self.parent.abundance_info_df.at[self.readset, 'most_abund_seq_of_coral_genus']
        try:
            if self.primary_sequence == self.parent.primary_seq_dict[self.genus_18S_taxonomic_annotation]:
                self.is_different_primary_sequence = False
            else:
                self.is_different_primary_sequence = True
                self.status = False
        except KeyError:
            # If the primary taxonomic annotation is other
            self.status = False
            self.is_different_primary_sequence = np.nan

    def _set_inter_coral_contamination(self):
        # Check to see whether the summed relative abundances of the coral sequences (other than the
        # sequences orginating from the genus_18S_taxonomic_annotation are above 1% of the total coral sample
        self.inter_coral_contamination_rel_abund = sum(
            [self.coral_tax_rel_count_dd[k] for k in self.sorted_coral_tax_dict_keys[1:]])
        if self.inter_coral_contamination_rel_abund > 0.01:
            self.is_inter_coral_contamination = True
            self.status = False
        else:
            self.is_inter_coral_contamination = False

    def _pop_for_heliopora_sample(self):
        self.status = False
        self.sample_annotation_dict = compress_pickle.load(
            os.path.join(self.sample_qc_dir, 'sample_annotation_dict.p.bz'))
        self.fasta_dict = self._make_fasta_dict()
        self.all_tax_count_dd = self._make_all_tax_count_dd()
        self.genus_18S_taxonomic_annotation = sorted(self.all_tax_count_dd, key=self.all_tax_count_dd.get, reverse=True)[0]
        self._set_is_provenance_tax_annotation_correct()

        # set inter coral contamination fields to nan
        self.inter_coral_contamination_rel_abund = np.nan
        self.is_inter_coral_contamination = np.nan
        # primary seq should be the most abundant seq in the sample
        self.primary_sequence = self.fasta_dict[sorted(self.rel_all_seq_abundance_dict, key=self.rel_all_seq_abundance_dict.get, reverse=True)[0]]
        # But we set is different to nan
        self.is_different_primary_sequence = np.nan
        # host rel abund will be the abund of heliopora
        self._set_host_rel_abund_heliopora()
        # Set the intragenus to nan
        self.putative_intra_genus_contamination_ratio = np.nan
        self.is_putative_intra_genus_contamination = np.nan
        self._set_post_qc_seq_depth()
        self._set_is_representative_for_sample()
        self._populate_coral_meta_info_table_dict()


    def _set_is_provenance_tax_annotation_correct(self):
        if self.genus_18S_taxonomic_annotation == self.provenance_annotation:
            self.is_provenance_tax_annotation_correct = True
        else:
            self.is_provenance_tax_annotation_correct = False
            self.status = False

    def _make_all_tax_count_dd(self):
        all_tax_count_dd = defaultdict(float)
        for seq_name, tax_designation in self.sample_annotation_dict.items():
            all_tax_count_dd[tax_designation[0]] += self.rel_all_seq_abundance_dict[seq_name]
        return all_tax_count_dd

    def _make_fasta_dict(self):
        fasta_path = os.path.join(self.parent.qc_dir, self.readset,
                               'stability.trim.contigs.good.unique.abund.pcr.unique.fasta')
        fasta_file_as_list = EighteenSBase.decompress_read_compress(fasta_path)
        # with open(os.path.join(self.parent.qc_dir, self.sample_id,
        #                        'stability.trim.contigs.good.unique.abund.pcr.unique.fasta'), 'r') as f:
        #     fasta_file_as_list = [line.rstrip() for line in f]
        fasta_dict = {}
        i = 0
        while i < len(fasta_file_as_list):
            sequence_name = fasta_file_as_list[i][1:].split('\t')[0]
            fasta_dict[sequence_name] = fasta_file_as_list[i + 1]
            i += 2
        return fasta_dict

    def populate_coral_meta_info_table_dict(self):
        if self.provenance_annotation in ["Porites", "Millepora", "Pocillopora"]:
            self._pop_for_normal_sample()
        elif self.provenance_annotation == "Heliopora":
            self._pop_for_heliopora_sample()
        else:
            if self.provenance_annotation.split(' ')[0] in ["Porites", "Millepora", "Pocillopora"]:
                self.provenance_annotation = self.provenance_annotation.split(' ')[0]
                self._pop_for_normal_sample()
            else:
                raise RuntimeError("unexpected provenance_annotation")



if __name__ == "__main__":
    ot = EighteenSOutputTables()
    # ot.make_and_write_coral_meta_info_output_table()
    # ot.make_and_write_tax_output_table()
    # ot.make_and_write_raw_abund_output_table()
    # ot.make_and_write_consolidated_host_output_table()
