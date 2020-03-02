import os
import compress_pickle
import sys
from base_18s import EighteenSBase
import pandas as pd


class EighteenSOutputTables(EighteenSBase):
    """ Class to produce the three output tables
    1 - Table that is the raw seq abundances (post mothur QC processing) in each of the samples (absolute abundance)
    2 - Table that holds the blast annotation for each of these sequences (Order, Family and Genus)
    3 - Table that is the post consolidation walk sequences from taxonomic origin of one of the
    three target genera. Relative proportions in the samples."""
    def __init__(self):
        super().__init__()
        # This will be a dict where full sequence is the key
        # The value will be another dict holding cumulative relative abundance
        # and the taxonomic tup
        self.master_seq_info_dict = self._populate_master_seq_info_dict()
        # Now get a master list of the sequences in order of cummulative abundance
        self.master_seq_abund_order_list = self._make_master_seq_abund_order_list()
        self.raw_seq_abund_path = os.path.join(self.output_dir, 'raw_seq_abund.csv')
        self.tax_df_path = os.path.join(self.output_dir, 'tax_annotation.csv')
        # Now popualte the dictionary that will be used to create the abundance df
        self.abundance_df_dict = self._populate_abundance_df_dict()
        self.abundance_df = None
        self.tax_annotation_df_dict = self._populate_tax_annotation_df_dict()
        self.tax_df = None

    def make_and_write_tax_output_table(self):
        self.tax_df = pd.DataFrame.from_dict(self.abundance_df_dict, orient='index',
                                                   columns=self.master_seq_abund_order_list)
        self.abundance_df.to_csv(self.raw_seq_abund_path, index=True, columns=True)

    def _populate_tax_annotation_df_dict(self):
        # The key should be sequence and the list should be order, family genus, in that order
        tax_annotation_df_dict = {}
        for seq, info_dict in self.master_seq_info_dict.items():
            annotation_tup = info_dict['tax_annotation']
            tax_annotation_df_dict[seq] = [annotation_tup[2], annotation_tup[1], annotation_tup[3]]
        return tax_annotation_df_dict

    def make_and_write_raw_abund_output_table(self):
        # Here we have the self.abundance_df_dict populated and we can now create the dataframe from this dict
        self.abundance_df = pd.DataFrame.from_dict(self.abundance_df_dict, orient='index',
                                                   columns=self.master_seq_abund_order_list)
        self.abundance_df.to_csv(self.raw_seq_abund_path, index=True, columns=True)

    def _populate_abundance_df_dict(self):
        print('Collecting sequence information: second pass')
        for sample_name in self.info_df.index:
            # read in the fasta file
            fasta_file = self._read_in_fasta_file(sample_name)

            # read in the name file and make an abundance dictionary
            name_abs_abund_dict = self._make_abs_abund_dict_from_names_path(sample_name)

            # create a seq_to_rel_abund dictionary
            seq_to_abs_abund_dict = {
                fasta_file[i+1]: name_abs_abund_dict[fasta_file[i].split('\t')[0][1:]] for
                i in range(0, len(fasta_file), 2)}

            # In the order of the self.master_seq_abund_order_list
            # populate the abundances for the given sample
            temp_abund_list = [seq_to_abs_abund_dict[seq] if seq in seq_to_abs_abund_dict else 0 for seq in self.master_seq_abund_order_list]
            self.abundance_df_dict[sample_name] = temp_abund_list

    def _make_master_seq_abund_order_list(self):
        self.master_seq_abund_order_list = [tup[0] for tup in sorted([_ for _ in self.master_seq_info_dict.items()],
                                                                     key=lambda x: x[1]['cummulative_abund'],
                                                                     reverse=True)]
    def _populate_master_seq_info_dict(self):
        print('Collecting sequence information: first pass')
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
                try:
                    self.master_seq_info_dict[fasta_file]['cummulative_abund'] += name_rel_abund_dict[seq_name]
                except KeyError:
                    try:
                        tax_tup = sample_annotation_dict[seq_name]
                    except KeyError:
                        tax_tup = ('not_annotated', 'not_annotated', 'not_annotated')
                    self.master_seq_info_dict[fasta_file] = {'cummulative_abund': name_rel_abund_dict[fasta_file[i]],
                                                             'tax_annotation': tax_tup}

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

