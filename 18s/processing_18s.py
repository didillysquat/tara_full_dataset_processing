"""This script performs the quality control of the 18S sequences and does the taxonomic annotation."""
import os
import sys
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool, current_process
import subprocess
import compress_pickle
from base_18s import EighteenSBase

class EighteenSProcessing(EighteenSBase):
    def __init__(self):
        super().__init__()
        # Perform the QC for the samples
        # Mothur processing that will be fed into the taxonomic screening
        seq_qc = SequenceQC(
            fastq_info_df=self.fastq_info_df, 
            qc_dir=self.qc_dir, 
            cache_dir=self.cache_dir, 
            root_dir=self.root_dir,
            sample_provenance_df=self.sample_provenance_df,
            seq_dir = self.seq_dir
            ).do_mothur_qc()

        # After the mothur QC we will have a set of .name and .fasta files for each sample
        # Now we need to annotate these sequences according to taxonomy
        # The results of this are a set of three pickled out annotation dictionaries
        # One is all sequnces, one is only coral, one is only symbiodiniaceae
        seq_qc.do_bast_qc()

        # Do seq consolidation. This will make and pickle out various abundance
        # dictionaries that we will use in the plotting and in the distance matrix creation
        SeqConsolidator(qc_dir=self.qc_dir, cache_dir=self.cache_dir, fastq_info_df=self.fastq_info_df, coral_readsets=self.coral_readsets).do_consolidation()

        # make the seq_to_total_abund_dict
        if os.path.isfile(os.path.join(self.cache_dir, 'seq_to_total_abund_dict.p.bz')):
            pass
        else:
            seq_to_total_abund_dict = defaultdict(int)
            for readset in self.coral_readsets:
                sample_qc_dir = os.path.join(self.qc_dir, readset)
                consolidated_host_seqs_abund_dict = compress_pickle.load(
                    os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
                for seq_name in consolidated_host_seqs_abund_dict.keys():
                    seq_to_total_abund_dict[seq_name] += 1
            compress_pickle.dump(seq_to_total_abund_dict, os.path.join(self.cache_dir, 'seq_to_total_abund_dict.p.bz'))
        
class SeqConsolidator:
    # In order to do the all_coral_sequence, we are going to need to move
    # the creation of abundance dataframes up here, as we need to keep track of global abundance
    # We also need to have the global abundances to assign the color_dict
    # 1 - Make a consolidation dictionary by walking in order of shortest first
    # 1a - Compress pickle out the abudnace dictionaries for each sample as you go (all seqs)
    # 2 - Go back through all of the sequences collecting abundances and relating them
    # to the corresponding representative consolidated sequences
    # 3 - Compress pickle out the 'translated' abudnace dictionaries for each sample (only coral seqs)
    # 3a - At the same time pickle out the 'translated' abundance dictionaires for each sample that 
    # contain only the 'host sequences': These are the sequencs that are of the genus of the most abundant
    # sequence. We can easily make the minor abundance dictionaries from this.
    # 4 - Make a color dictionary according to the abundance order. Make sure that the yellow
    # red and blue still represent the most abundant sequences in the corresponding corals
    def __init__(self, qc_dir, cache_dir, fastq_info_df, coral_readsets):
        self.qc_dir = qc_dir
        self.cache_dir = cache_dir
        self.fastq_info_df = fastq_info_df
        self.consolidated_host_seqs_rel_abundance_dict = None
        self.coral_blasted_seq_to_consolidated_seq_dict = {}
        self.coral_readsets = coral_readsets

    def do_consolidation(self):
        if not self._check_if_consolidation_already_complete():
            # Firstly get a set of all sequences that are of one of the three genera
            # In the process of doing this, pickle out the all seq abunance dict for
            # the sample
            self.consolidated_host_seqs_rel_abundance_dict = self._make_host_seqs_dict()
            # Then create the consolidation path and
            # consolidate the sequence insitu in the self.host_seqs_dict
            # Also produce a dict that maps blasted_seq to representative consolidated sequence
            # Then pickle out the consolidated_host_seqs_rel_abund_dict
            # and the coral_blasted_seq_to_consolidated_seq_dict
            self._create_consolidation_path_list()
            
            # Revisit the sample directories and make pickle out the:
            # 1 - consolidated_coral_seqs_abund_dict
            # 2 - consolidated_host_seqs_abund_dict
            # The first contains all coral sequences
            # The second contains only coral sequences that are of the genus of the most
            # abundant coral sequence
            # In both cases we will use the consolidated representatives of each sequence
            self._create_and_write_sample_coral_consolidated_rel_abund_dicts()
            self._make_all_coral_sequence_color_dict()
        else:
            return
    
    def _make_all_coral_sequence_color_dict(self):
        """We want to create a color dict where the sequences are coloured in order of abundance
        and when we run out of colors we will use the greys as usual.
        We want to ensure that the most abundant Porites, Millepora and Pocillopora
        sequences are still the same, yellow red and blue, so we'll want to hard code this
        """
        color_dict = {}
        color_list = self._get_colour_list()
        greys = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
        # remove the three colours of interest from the list and add them to the beginning
        # We will need to adjust the order of these three sequences as we don't know
        # which was most abundant. I.e. was porites seq most abundant or millepora etc.
        self._curate_color_list(color_list)
        sorted_seqs_tups = sorted(
            [(seq_name, rel_abund) for seq_name, rel_abund in self.consolidated_host_seqs_rel_abundance_dict.items()], 
            key=lambda x: x[1], 
            reverse=True
            )
        sorted_names = [tup[0] for tup in sorted_seqs_tups]
        for i, seq_name in enumerate(sorted_names):
            if i < len(color_list):
                color_dict[seq_name] = color_list[i]
            else:
                color_dict[seq_name] = greys[i%6]
        
        compress_pickle.dump(color_dict, os.path.join(self.cache_dir, 'all_coral_sequence_color_dict.p.bz'))
        return color_dict

    @staticmethod
    def _curate_color_list(color_list):
        if "#FFFF00" in color_list: color_list.remove("#FFFF00")
        if "#87CEFA" in color_list: color_list.remove("#87CEFA") 
        if "#FF6347" in color_list: color_list.remove("#FF6347")
        if "#00FF00" in color_list: color_list.remove("#00FF00")
        color_list.insert(0, "#FF6347")
        color_list.insert(0, "#87CEFA")
        color_list.insert(0, "#FFFF00")
        # Swap out the 5 and 11 elements as the 5 element is currently a similar color
        # to the 2 element
        five = color_list[5]
        color_list[5] = color_list[11]
        color_list[11] = five
        three = color_list[3]
        color_list[3] = color_list[22]
        color_list[22] = three


    def _check_if_consolidation_already_complete(self):
        if os.path.isfile(os.path.join(self.cache_dir, 'final_consolidated_host_seqs_rel_abundance_dict.p.bz')):
            if os.path.isfile(os.path.join(self.cache_dir, 'coral_blasted_seq_to_consolidated_seq_dict.p.bz')):
                if os.path.isfile(os.path.join(self.cache_dir, 'all_coral_sequence_color_dict.p.bz')):
                    return True
        return False

    def _create_and_write_sample_coral_consolidated_rel_abund_dicts(self):
        print('\nWriting out sample abundance dictionaries\n')
        for readset in self.coral_readsets:
            sys.stdout.write(f'\r{readset}')
            sample_qc_dir = os.path.join(self.qc_dir, readset)
            # Load the already created abundance dictionary
            rel_all_seq_abundance_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'rel_all_seq_abundance_dict.p.bz'))
            # Load the already created taxonomy annotation dictoinaries
            sample_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'sample_annotation_dict.p.bz'))
            coral_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'coral_annotation_dict.p.bz'))
            seq_name_to_seq_dict = self._make_seq_name_to_seq_dict(readset)
            # # Firstly we can write out the consolidated_coral_seqs_abund_dict
            # # There is the possibility that several sequences could be represented by
            # # the same sequences. We will need to check for this and combine these abundances if so
            # self._make_write_consolidated_coral_seqs_abund_dict(
            #     rel_all_seq_abundance_dict, sample_annotation_dict, sample_qc_dir, seq_name_to_seq_dict)

            # Now it is time to do the consolidated_host_seqs_abund_dict
            # Firstly identify the most abundant coral sequence
            most_abundant_coral_genus = self._identify_most_abund_coral_genus(
                rel_all_seq_abundance_dict=rel_all_seq_abundance_dict, 
                coral_annotation_dict=coral_annotation_dict
                )
            
            # Here we have the most abundant genus identified
            self._make_write_consolidated_host_seqs_abund_dict(
                rel_all_seq_abundance_dict, coral_annotation_dict, 
                most_abundant_coral_genus, sample_qc_dir, seq_name_to_seq_dict)
        sys.stdout.write('\n')

    # def _make_write_consolidated_coral_seqs_abund_dict(
    #     self, rel_all_seq_abundance_dict, sample_annotation_dict, sample_qc_dir, seq_name_to_seq_dict):
    #     # Firstly we can write out the consolidated_coral_seqs_abund_dict
    #     # There is the possibility that several sequences could be represented by
    #     # the same sequences. We will need to check for this and combine these abundances if so
    #     consolidated_coral_seqs_abund_dict = {}
    #     for seq_name, rel_abund in rel_all_seq_abundance_dict.items():
    #         try:
    #             if sample_annotation_dict[seq_name] == 'Scleractinia_Anthoathecata':
    #                 try:
    #                     rep_consol_seq = self.coral_blasted_seq_to_consolidated_seq_dict[seq_name_to_seq_dict[seq_name]]
    #                 except KeyError:
    #                     # If key error, then the seq was not consolidated and is itself a representative consolidation
    #                     # sequence and so can be used directly in the consolidated_coral_seqs_abund_dict
    #                     # TODO here we can check that the sequence can be found in the
    #                     if seq_name_to_seq_dict[seq_name] not in self.consolidated_host_seqs_rel_abundance_dict:
    #                         foo = 'bar'
    #                     rep_consol_seq = seq_name_to_seq_dict[seq_name]
    #                 if rep_consol_seq in consolidated_coral_seqs_abund_dict:
    #                     # Then there were multiple seuqence represented by this consol sequence
    #                     # and we need to combine the relative abunances
    #                     current_abund = consolidated_coral_seqs_abund_dict[rep_consol_seq]
    #                     new_abund = current_abund + rel_abund
    #                     consolidated_coral_seqs_abund_dict[rep_consol_seq] = new_abund
    #                 else:
    #                     consolidated_coral_seqs_abund_dict[rep_consol_seq] = rel_abund
    #             else:
    #                 # Then this was not a coral sequence and we are not concerned with it
    #                 pass
    #         except KeyError:
    #             # There was no meaningful annotation available for this sequence in the nt database
    #             pass
    #
    #     # Finally we need to renormalise the realtive abundances for this dictionary as we
    #     # have removed the non_coral samples
    #     tot = sum(consolidated_coral_seqs_abund_dict.values())
    #     consolidated_coral_seqs_abund_dict = {k: v/tot for k, v in consolidated_coral_seqs_abund_dict.items()}
    #     compress_pickle.dump(consolidated_coral_seqs_abund_dict, os.path.join(sample_qc_dir, 'consolidated_coral_seqs_abund_dict.p.bz'))

    def _make_write_consolidated_host_seqs_abund_dict(self, rel_all_seq_abundance_dict, coral_annotation_dict, most_abundant_coral_genus, sample_qc_dir, seq_name_to_seq_dict):
        consolidated_host_seqs_abund_dict = {}
        for seq_name, rel_abund in rel_all_seq_abundance_dict.items():
            # if seq_name == 'GTCGCTACTACCGATTGAATGGTTTAGTGAGGCCTCCTGACTGGCGCCGACACTCTGTCTCGTGCAGAGAGTGGGAGGCCGGGAAGTTGTTCAAACTTGATCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTC':
            #     foo = 'bar'
            try:
                if coral_annotation_dict[seq_name] == most_abundant_coral_genus:
                    # Then this is of the genus that we are interested in
                    # See if there if the sequence has a representative consolidated sequence
                    # If not then use the sequence its self as the key
                    try:
                        rep_consol_seq = self.coral_blasted_seq_to_consolidated_seq_dict[seq_name_to_seq_dict[seq_name]]
                    except KeyError:
                        rep_consol_seq = seq_name_to_seq_dict[seq_name]
                        # TODO test to see that the seq can be found in the abundance dict.
                        if not rep_consol_seq in self.consolidated_host_seqs_rel_abundance_dict:
                            foo = 'bar'
                    if rep_consol_seq in consolidated_host_seqs_abund_dict:
                        # Then there were multiple seuqence represented by this 
                        # consolidated sequence and we need to combine the relative abunances
                        current_abund = consolidated_host_seqs_abund_dict[rep_consol_seq]
                        new_abund = current_abund + rel_abund
                        consolidated_host_seqs_abund_dict[rep_consol_seq] = new_abund
                    else:
                        consolidated_host_seqs_abund_dict[rep_consol_seq] = rel_abund
                else:
                    # Then this is not of the genus we are interested in

                    pass
            except KeyError:
                # Then this was not a coral sequence and we are not concerned with it
                pass
        # Finally we need to renormalise the realtive abundances for this dictionary as we
        # have removed the non_coral samples
        tot = sum(consolidated_host_seqs_abund_dict.values())
        consolidated_host_seqs_abund_dict = {k: v/tot for k, v in consolidated_host_seqs_abund_dict.items()}
        compress_pickle.dump(consolidated_host_seqs_abund_dict, os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))

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
  
    def _write_out_dicts(self):
        compress_pickle.dump(
            self.consolidated_host_seqs_rel_abundance_dict, 
            os.path.join(self.cache_dir, 'final_consolidated_host_seqs_rel_abundance_dict.p.bz')
            )
        compress_pickle.dump(
        self.coral_blasted_seq_to_consolidated_seq_dict, 
        os.path.join(self.cache_dir, 'coral_blasted_seq_to_consolidated_seq_dict.p.bz')
        )

    def _consolidate(self, consolidation_path_list, representative_to_seq_list):
        """In this method we will do two things.
        1 - We will walk the path consolidating the sequences by modifying
        the self.host_seqs_dict insitu.
        2 - We will also create a new dict that is original sequence, to final 
        representative consolidated sequence. To create this dictionary, we will
        use an additional utility dictionary that will keep track of for a given
        sequence, which sequenecs it is the current representative consolidation
        sequence for during the walk of the consolidation path. That way, every time
        we get to a new seq consolidation tuple, for the seq being consolidated, we will
        check to see which seqs it is representative of, and also transfer all of these
        sequences to being represented by the super sequence in question

        NB now that we are doing a bottom up and a top down approach"""
        print('Doing sequence consolidation\n')
        count = 1
        tot = len(consolidation_path_list)
        # search_seq = 'GTCGCTACTACCGATTGAATGGTTTAGTGAGGCCTCCTGACTGGCGCCGACACTCTGTCTCGTGCAGAGAGTGGGAGGCCGGGAAGTTGTTCAAACTTGATCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTC'

        for query_seq, match_seq in consolidation_path_list:
            # if query_seq == search_seq:
            #     foo = 'asdf'
            sys.stdout.write(f'\r{count} out of {tot}')
            count += 1
            query_seq_cummulative_rel_abund = self.consolidated_host_seqs_rel_abundance_dict[query_seq]
            match_seq_cummulative_rel_abund = self.consolidated_host_seqs_rel_abundance_dict[match_seq]
            self.consolidated_host_seqs_rel_abundance_dict[match_seq] = query_seq_cummulative_rel_abund + match_seq_cummulative_rel_abund
            del self.consolidated_host_seqs_rel_abundance_dict[query_seq]
            if query_seq in representative_to_seq_list:
                # then the seq being consolidated was a representative sequence for other sequences
                # and these sequences should also be transfered under the new representative
                for seq_to_transfer in representative_to_seq_list[query_seq]:
                    self.coral_blasted_seq_to_consolidated_seq_dict[seq_to_transfer] = match_seq
                    representative_to_seq_list[match_seq].append(seq_to_transfer)
                del representative_to_seq_list[query_seq]
            self.coral_blasted_seq_to_consolidated_seq_dict[query_seq] = match_seq
            representative_to_seq_list[match_seq].append(query_seq)
        print('\nconsolidation complete\n')

    def _make_host_seqs_dict(self):
        """This is the first step in the process. We need to go through all of the 
        seuences and get a list of all host sequences that are of one of the genera
        we are concerned with. Collect a cummulative abundance for the sequences too so that we can
        use this information when making the consolidation path. We will then create a 
        consolidation path from this."""
        # for every sample, create an complete abundance dict of seq to abundance and pickle it out
        # for those sequences that are of one of the three genera, then add them to the list
        if os.path.isfile(os.path.join(self.cache_dir, 'initial_consolidated_host_seqs_rel_abundance_dict.p.bz')):
            return compress_pickle.load(os.path.join(self.cache_dir, 'initial_consolidated_host_seqs_rel_abundance_dict.p.bz'))
        consolidated_host_seqs_rel_abundance_dict = defaultdict(float)
        for readset in self.coral_readsets:
            print(f'Getting coral sequences from readset {readset}')
            sample_qc_dir = os.path.join(self.qc_dir, readset)
            abs_all_seq_abundance_dict = self._make_abund_dict_from_names_path(readset)
            seq_name_to_seq_dict = self._make_seq_name_to_seq_dict(readset)
            compress_pickle.dump(
                abs_all_seq_abundance_dict, 
                os.path.join(sample_qc_dir, 'abs_all_seq_abundance_dict.p.bz'))
            tot = sum(abs_all_seq_abundance_dict.values())
            rel_all_seq_abundance_dict = {k: v/tot for k, v in abs_all_seq_abundance_dict.items()}
            compress_pickle.dump(
                rel_all_seq_abundance_dict, 
                os.path.join(sample_qc_dir, 'rel_all_seq_abundance_dict.p.bz'))
            sample_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'sample_annotation_dict.p.bz'))
            coral_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'coral_annotation_dict.p.bz'))

            for blasted_seq, annotation in sample_annotation_dict.items():
                if annotation[2] in ['Scleractinia', 'Anthoathecata']:
                    # Then this is a coral seq
                    # If it is of one of the three genera in question then we should
                    # add it to the list of sequences
                    if coral_annotation_dict[blasted_seq] in ['Porites', 'Pocillopora', 'Millepora']:
                        consolidated_host_seqs_rel_abundance_dict[
                            seq_name_to_seq_dict[blasted_seq]
                            ] += rel_all_seq_abundance_dict[blasted_seq]
        compress_pickle.dump(
            consolidated_host_seqs_rel_abundance_dict, 
            os.path.join(self.cache_dir, 'initial_consolidated_host_seqs_rel_abundance_dict.p.bz')
            )
        return consolidated_host_seqs_rel_abundance_dict

    def _make_seq_name_to_seq_dict(self, readset):
        fasta_path = os.path.join(self.qc_dir, readset, 'stability.trim.contigs.good.unique.abund.pcr.unique.fasta')
        fasta_file_as_list = EighteenSBase.decompress_read_compress(fasta_path)
        # with open(os.path.join(self.qc_dir, readset, 'stability.trim.contigs.good.unique.abund.pcr.unique.fasta'), 'r') as f:
        #     fasta_file_as_list = [line.rstrip() for line in f]
        temporary_dictionary = {}
        i = 0
        while i < len(fasta_file_as_list):
            sequence_name = fasta_file_as_list[i][1:].split('\t')[0]
            temporary_dictionary[sequence_name] = fasta_file_as_list[i + 1]
            i += 2
        return temporary_dictionary

    def _create_consolidation_path_list(self):
        """We will do a run from short to long and vice versa. We will only consolidate one sequence
        into another if the sequence to be consolidated has a lower abundance than the putative representative
        sequence."""
        representative_to_seq_list = defaultdict(list)
        for i in range(2):
            # First get a list of the sequences to work with sorted by order of length
            if i == 0:
                # Work short to long
                if os.path.isfile(os.path.join(self.cache_dir, 'consolidation_path_list_s_to_l.p.bz')):
                    consolidation_path_list = compress_pickle.load(
                        os.path.join(self.cache_dir, 'consolidation_path_list_s_to_l.p.bz'))
                    # search_seq = 'GTCGCTACTACCGATTGAATGGTTTAGTGAGGCCTCCTGACTGGCGCCGACACTCTGTCTCGTGCAGAGAGTGGGAGGCCGGGAAGTTGTTCAAACTTGATCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTC'
                    # for tup in consolidation_path_list:
                    #     if search_seq in tup:
                    #         if tup[0] == search_seq:
                    #             print('search seq is tup 0')
                    #         else:
                    #             print('search seq is tup 1')
                else:
                    seq_list = sorted(list(self.consolidated_host_seqs_rel_abundance_dict.keys()), key=len)
                    consolidation_path_list = self._walk_and_test_match_s_to_l(seq_list)
                    compress_pickle.dump(consolidation_path_list, os.path.join(self.cache_dir, 'consolidation_path_list_s_to_l.p.bz'))
            else:
                # work long to short
                seq_list = sorted(list(self.consolidated_host_seqs_rel_abundance_dict.keys()), key=len, reverse=True)
                consolidation_path_list = self._walk_and_test_match_l_to_s(seq_list)
            self._consolidate(
                consolidation_path_list, 
                representative_to_seq_list)
        
        # Pickle out the host_seqs_dict and the consolidated_seq_dict
        self._write_out_dicts()

    def _walk_and_test_match_s_to_l(self, seq_list):
        consolidation_path_list = []
        # Go from start length to end length - 1
        start_n = len(seq_list[0])
        finish_n = len(seq_list[-1])
        print('\nMaking consolidation path short to long')
        for n in range(start_n, finish_n):
            query_seqs = [seq for seq in seq_list if len(seq) == n]
            if not query_seqs:
                continue
            putative_match_seqs = [seq for seq in seq_list if len(seq) > n]
            sub_putative_match_abund_dict = {p_match_seq: self.consolidated_host_seqs_rel_abundance_dict[p_match_seq] for p_match_seq in putative_match_seqs}
            count = 0
            tot_query_seqs = len(query_seqs)
            tot_p_match_seqs = len(putative_match_seqs)
            for q_seq in query_seqs:
                q_seq_abund = self.consolidated_host_seqs_rel_abundance_dict[q_seq]
                matches = [seq for seq in putative_match_seqs if (sub_putative_match_abund_dict[seq] > q_seq_abund) and (q_seq in seq)]
                count += 1
                sys.stdout.write(f'\rseq {count} out of {tot_query_seqs} for level n={n} of {finish_n}. {tot_p_match_seqs} seqs to check against.')
                
                if matches:
                    if len(matches) > 1:
                        # Here we create a list of tuples where the match sequence and the number of
                        # DataSetSamples that contained that sample.
                        # We then sort it according to the cumulative abundance
                        # Finally, we use this sequence as the consolidation representative for the
                        # consolidation path
                        representative_seq = sorted(
                            [(match_seq, self.consolidated_host_seqs_rel_abundance_dict[match_seq]) for match_seq in matches],
                            key=lambda x: x[1], reverse=True)[0][0]
                        consolidation_path_list.append((q_seq, representative_seq))
                    else:
                        consolidation_path_list.append((q_seq, matches[0]))
                else:
                    # If there are no matches then there is no entry required in the consolidation path
                    pass
        return consolidation_path_list

    def _walk_and_test_match_l_to_s(self, seq_list):
        consolidation_path_list = []
        # Go from start length to end length - 1
        start_n = len(seq_list[0])
        finish_n = len(seq_list[-1])
        print('\nMaking consolidation path long to short')
        for n in range(start_n, finish_n, -1):
            query_seqs = [seq for seq in seq_list if len(seq) == n]
            if not query_seqs:
                continue
            putative_match_seqs = [seq for seq in seq_list if len(seq) < n]
            sub_putative_match_abund_dict = {p_match_seq: self.consolidated_host_seqs_rel_abundance_dict[p_match_seq] for p_match_seq in putative_match_seqs}
            count = 0
            tot_query_seqs = len(query_seqs)
            tot_p_match_seqs = len(putative_match_seqs)
            for q_seq in query_seqs:
                q_seq_abund = self.consolidated_host_seqs_rel_abundance_dict[q_seq]
                matches = [seq for seq in putative_match_seqs if (sub_putative_match_abund_dict[seq] > q_seq_abund) and (seq in q_seq)]
                count += 1
                sys.stdout.write(f'\rseq {count} out of {tot_query_seqs} for level n={n} of {finish_n}. {tot_p_match_seqs} seqs to check against.')
                
                if matches:
                    if len(matches) > 1:
                        # Here we create a list of tuples where the match sequence and the number of
                        # DataSetSamples that contained that sample.
                        # We then sort it according to the cumulative abundance
                        # Finally, we use this sequence as the consolidation representative for the
                        # consolidation path
                        representative_seq = sorted(
                            [(match_seq, self.consolidated_host_seqs_rel_abundance_dict[match_seq]) for match_seq in matches],
                            key=lambda x: x[1], reverse=True)[0][0]
                        consolidation_path_list.append((q_seq, representative_seq))
                    else:
                        consolidation_path_list.append((q_seq, matches[0]))
                else:
                    # If there are no matches then there is no entry required in the consolidation path
                    pass
        return consolidation_path_list


    def _make_abund_dict_from_names_path(self, readset):
        name_path = os.path.join(self.qc_dir, readset, 'stability.trim.contigs.good.unique.abund.pcr.names')
        name_file = EighteenSBase.decompress_read_compress(name_path)
        # with open(os.path.join(self.qc_dir, readset, 'stability.trim.contigs.good.unique.abund.pcr.names'), 'r') as f:
        #     name_file = [line.rstrip() for line in f]
        return {line.split('\t')[0]: len(line.split('\t')[1].split(',')) for line in name_file}

    @staticmethod
    def _get_colour_list():
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



class SequenceQC:
    """Quality control of the sequences using Mothur firstly and then running through BLAST
    to remove all non-coral sequences"""
    def __init__(self, fastq_info_df, qc_dir, cache_dir, root_dir, sample_provenance_df, seq_dir):
        self.fastq_info_df = fastq_info_df
        self.qc_dir = qc_dir
        self.root_dir = root_dir
        self.cache_dir = cache_dir
        # These dicts will be used for getting the taxonomic level, and name of a
        # taxonomic match.
        self.taxdump_dir = '/share/databases/nt/taxdump'
        self.node_dict_pickle_path = os.path.join(cache_dir, 'node_dict.p.bz')
        self.name_dict_pickle_path = os.path.join(cache_dir, 'name_dict.p.bz')
        self.node_dict, self.name_dict = self._generate_taxa_and_name_dicts()
        self.sample_provenance_df = sample_provenance_df
        # The directory where the sequencing files are
        self.seq_dir = seq_dir

    def do_mothur_qc(self):
        # For some reason the make.contigs doesn't seem to be using multiple processors and thus
        # this is wasting a lot of time. As such we should try to multiprocess this and use 1
        # process for each. Let's try to implement with Pool
        # We will need a list to pass to Pool for doing the Apply with
        apply_list = []
        for readset, ser in self.fastq_info_df.iterrows():
            # Check to see if this is a coral sample
            if self.sample_provenance_df.at[ser['sample-id'], 'SAMPLE ENVIRONMENT, short'] == 'C-CORAL':
                if self._check_if_qc_already_complete(readset):
                    continue
                else:
                    apply_list.append([readset, os.path.join(self.seq_dir, ser['fwd_read_name']), os.path.join(self.seq_dir, ser['rev_read_name']), self.qc_dir])
        
        if apply_list:
            # NB the Pool does not play well with the debugger that VS code is using.
            # Pool does work OK if run normally (i.e. no debug)
            # By having the Pool operation covered by the if apply_list:
            # we should still be able to debug as we won't have to aceess Pool.
            with Pool(200) as p:
                p.map(self._init_indi_mothur_qc, apply_list)
        return self

    def _init_indi_mothur_qc(self, sample_info):
        """We have to pass the Pool a single function. This is the function we will use to instantiate 
        an IndividualMothurQC class and run the qc."""
        indi = IndividualMothurQC(
            readset=sample_info[0], 
            fwd_path=sample_info[1], 
            rev_path=sample_info[2], 
            qc_base_dir=sample_info[3],
            num_proc=1)
        indi.start_qc()

    def _check_if_qc_already_complete(self, readset):
        return os.path.exists(os.path.join(self.qc_dir, readset, 'qc_complete.txt'))

    def do_bast_qc(self):
        """
        For every set of sequence files run a blast against the nt database.
        For each sequence get 10 blast results.
        Look at the matches and see which taxon was highest represented.
        We will produce three annotation dictionaries.
        In the first, sample_annotation_dict, taxa will be annotated to a genus level
        unless they are Scleractinia, Anthoathecata or Symbiodiniaceae in which case
        they will be annotated thusly.
        
        We will also produce two additional dicts. One for only Symbiodiniaceae and one for only Scleractinia or Anthoathecata.
        In both of these dicts, genus level annotations will be used for the sequences.

        We will run this in serial but use multiple processors for running the blast.

        # TODO given that we are running on zygote and can easily use up the 200 cores, I think that best speed
        # will likely be realised through the use of both multiple processes and using multiple cores for the blast.
        # Perhaps we could try 20 core blasts and use a pool of 10 processes.
        """
        # self._write_out_ncbi_db_config_file()

        # Get an apply list
        apply_list = []
        for readset, ser in self.fastq_info_df.iterrows():
            if self.sample_provenance_df.at[ser['sample-id'], 'SAMPLE ENVIRONMENT, short'] == 'C-CORAL':
                if not os.path.exists(os.path.join(self.qc_dir, readset, 'taxonomy_complete.txt')):
                    apply_list.append(readset)
                
        if apply_list:
            # NB having this set at using a pool of size 12 and using 40 'threads' for the blast
            # works well for doing the blast.
            # There is some bottle neck somewhere for the second annotation part of the 
            # taxonomy processing. You can turn the pool upto 120 but only 6 processes get used.
            # We will not spend the time optimising this now as it is running fast enough.
            # Completes in approx 1 hour. (after having already done the blasts)
            # But if you come back to this in the future you may want to optimise.
            with Pool(120) as p:
                p.map(self._set_tax_running, apply_list)
        
        # for readset in apply_list:
        #     self._set_tax_running(readset)

        foo = 'bar'
    
    def _set_tax_running(self, readset):
        try:
            print(f'{current_process().name}: readset is {readset}')
            sample_id = self.fastq_info_df.loc[readset]['sample-id']
            tax = self.sample_provenance_df.at[sample_id, 'SAMPLE ENVIRONMENT, short']
            indi_tax_an = IndividualTaxAnnotationQC(
                    cache_dir=self.cache_dir, 
                    qc_dir=self.qc_dir, 
                    readset=readset
                    )
            indi_tax_an.do_taxonomy_annotation()
        except Exception as e:
            print(f'{current_process().name}: something went wrong with {readset}')
            print(e)
            raise RuntimeError
        # indi_tax_an.do_taxonomy_annotation()

    def _generate_taxa_and_name_dicts(self):
            if os.path.isfile(self.node_dict_pickle_path):
                node_dict = compress_pickle.load(self.node_dict_pickle_path)
            else:
                # read in the .nodes file. This file tells us which tax level the node is and which node is the parent level
                with open(f'{self.taxdump_dir}/nodes.dmp', 'r') as f:
                    node_file = [line.rstrip() for line in f]
                # now make a dict from this where key is the tax id and the value is a tup where 0 = parent 1 = tax level
                node_dict = {line.split('\t|\t')[0]: (line.split('\t|\t')[1], line.split('\t|\t')[2]) for line in node_file}
                compress_pickle.dump(node_dict, self.node_dict_pickle_path)

            if os.path.isfile(self.name_dict_pickle_path):
                name_dict = compress_pickle.load(self.name_dict_pickle_path)
            else:
                # next read in the names file. This file hold the name of the node.
                with open(f'{self.taxdump_dir}/names.dmp', 'r') as f:
                    name_file = [line.rstrip() for line in f]
                # now make a dict from the names file where the key is the staxid and the value is the name
                name_dict = {line.split('\t|\t')[0]: line.split('\t|\t')[1] for line in name_file if
                            line.split('\t|\t')[3].replace('\t|', '') == 'scientific name'}
                compress_pickle.dump(name_dict, self.name_dict_pickle_path)

            return node_dict, name_dict

class IndividualTaxAnnotationQC:
    def __init__(self, cache_dir, qc_dir, readset):
        self.node_dict_pickle_path = os.path.join(cache_dir, 'node_dict.p.bz')
        self.name_dict_pickle_path = os.path.join(cache_dir, 'name_dict.p.bz')
        self.node_dict, self.name_dict = self._generate_taxa_and_name_dicts()
        self.genus = None
        self.family = None
        self.order = None
        self.readset = readset
        self.qc_dir = qc_dir
        self.fasta_path = os.path.join(qc_dir, readset, 'stability.trim.contigs.good.unique.abund.pcr.unique.fasta')
        # print(f'{current_process().name}: The fasta path is {self.fasta_path}')
        self.blast_out_path = os.path.join(qc_dir, readset, 'blast.out')
        self.output_format = "6 qseqid sseqid staxids evalue pident qcovs staxid stitle ssciname"
        self.blast_out_pickle_path = os.path.join(qc_dir, readset, 'blast.out.p.bz')
        self.blast_output_file = self._get_or_do_blast()
        # This is a dict with seq name as key, and the 10 blast results as a list of lists as value
        self.blast_output_dict = self._make_blast_out_dict()
        # We will have three dictionaries that will hold the results of the taxonomic annotation
        # This dict will hold the annotation for all sequenecs
        self.sample_annotation_dict = {}
        self.coral_annotation_dict = {}
        self.symbiodiniaceae_annotation_dict = {}
        # The paths that we will pickle out the annotation dictionaries to
        self.sample_annotation_dict_pickle_path = os.path.join(qc_dir, readset, 'sample_annotation_dict.p.bz')
        self.coral_annotation_dict_pickle_path = os.path.join(qc_dir, readset, 'coral_annotation_dict.p.bz')
        self.coral_symbiodiniaceae_dict_pickle_path = os.path.join(qc_dir, readset, 'symbiodiniaceae_annotation_dict.p.bz')

    def _generate_taxa_and_name_dicts(self):
            if os.path.isfile(self.node_dict_pickle_path):
                node_dict = compress_pickle.load(self.node_dict_pickle_path)
            else:
                # read in the .nodes file. This file tells us which tax level the node is and which node is the parent level
                with open(f'{self.taxdump_dir}/nodes.dmp', 'r') as f:
                    node_file = [line.rstrip() for line in f]
                # now make a dict from this where key is the tax id and the value is a tup where 0 = parent 1 = tax level
                node_dict = {line.split('\t|\t')[0]: (line.split('\t|\t')[1], line.split('\t|\t')[2]) for line in node_file}
                compress_pickle.dump(node_dict, self.node_dict_pickle_path)

            if os.path.isfile(self.name_dict_pickle_path):
                name_dict = compress_pickle.load(self.name_dict_pickle_path)
            else:
                # next read in the names file. This file hold the name of the node.
                with open(f'{self.taxdump_dir}/names.dmp', 'r') as f:
                    name_file = [line.rstrip() for line in f]
                # now make a dict from the names file where the key is the staxid and the value is the name
                name_dict = {line.split('\t|\t')[0]: line.split('\t|\t')[1] for line in name_file if
                            line.split('\t|\t')[3].replace('\t|', '') == 'scientific name'}
                compress_pickle.dump(name_dict, self.name_dict_pickle_path)
                
            return node_dict, name_dict

    def do_taxonomy_annotation(self):
        sys.stdout.write(f'\r{self.readset}')
        for blasted_seq, result_list in self.blast_output_dict.items():
            self.IndividualTaxAnnotationQCSub(
                blasted_seq_name=blasted_seq, 
                results_list=result_list, 
                sample_dict=self.sample_annotation_dict, 
                symbiodiniaceae_dict=self.symbiodiniaceae_annotation_dict, 
                coral_dict=self.coral_annotation_dict, node_dict=self.node_dict, names_dict=self.name_dict
                ).process_matches()
        
        # Here we have populated the taxonomy dictionaries for the sample in question
        # Now pickle out the dictionaries
        compress_pickle.dump(self.sample_annotation_dict, self.sample_annotation_dict_pickle_path)
        compress_pickle.dump(self.coral_annotation_dict, self.coral_annotation_dict_pickle_path)
        compress_pickle.dump(self.symbiodiniaceae_annotation_dict, self.coral_symbiodiniaceae_dict_pickle_path)
        with open(os.path.join(self.qc_dir, self.readset, 'taxonomy_complete.txt'), 'w') as f:
            f.write(f'{self.readset} complete')

    def _make_blast_out_dict(self):
        # now create a dict that is the list of 10 result items using the seq name as key
        blast_output_dict = defaultdict(list)
        for output_line in self.blast_output_file:
            components = output_line.split('\t')
            blast_output_dict[components[0]].append(components)
        return blast_output_dict

    def _get_or_do_blast(self):
        """
        Check to see if the blast has already been performed and can be read
        in from file. Else perform the blast from scratch.
        """
        # if os.path.isfile(self.blast_out_pickle_path) and os.path.isfile(os.path.join(self.qc_dir, self.readset, 'blast_complete.txt')): # TODO add second verification file
        if os.path.isfile(self.blast_out_pickle_path):
            return compress_pickle.load(self.blast_out_pickle_path)
        else:
            # Run local blast
            print(f'{current_process().name}: Running blast for sample {self.readset}')
            # Decompress
            if os.path.exists(self.fasta_path + '.gz'):
                subprocess.run(['gzip', '-d', self.fasta_path + '.gz'])
            result = subprocess.run(
                ['blastn', '-out', self.blast_out_path, '-outfmt', self.output_format, '-query', self.fasta_path, '-db', '/share/databases/nt/nt',
                    '-max_target_seqs', '10', '-num_threads', '4'])
            if os.path.exists(self.fasta_path):
                subprocess.run(['gzip', self.fasta_path])

            # Read in blast output
            with open(self.blast_out_path, 'r') as f:
                blast_output_file = [line.rstrip() for line in f]

            # pickle out the blast results here to possibly save us time in the future
            compress_pickle.dump(blast_output_file, self.blast_out_pickle_path)
            # Write out a second file for verification purposes, so that we know the
            # write out of the blast out pickle was completed
            with open(os.path.join(self.qc_dir, self.readset, 'blast_complete.txt'), 'w') as f:
                    f.write(f'{self.readset} blast complete')
            # And then delete the blast.out file
            os.remove(self.blast_out_path)
            return blast_output_file
    
    class IndividualTaxAnnotationQCSub:
        """This is a utility class that is used for counting the genus hits for a given sequence"""
        def __init__(self, blasted_seq_name, results_list, sample_dict, symbiodiniaceae_dict, coral_dict, node_dict, names_dict):
            self.node_dict = node_dict
            self.names_dict = names_dict
            # Name of the sequence that was blasted
            self.blasted_seq_name = blasted_seq_name
            # This is the list that will hold 10 lists, 
            # one for each of the top 10 hits for the given sequence
            self.super_results_list = results_list
            # We will use this to hold the result list for the current hit in question
            self.current_result_list = None
            # The annotation dictionaries from the parent IndividualTaxAnnotation class instance
            # we will eventually add the annotation from this given sequence to one of these
            # depending on how sequence is matched in the nt database.
            self.sample_genus_dict = sample_dict
            self.symbiodiniceae_dict = symbiodiniaceae_dict
            self.coral_dict = coral_dict
            # The counters for the three categories of taxa
            self.other_genus_count_dict = defaultdict(int)
            self.symbiodiniaceae_genus_count_dict = defaultdict(int)
            self.coral_genus_count_dict = defaultdict(int)

        def process_matches(self):
            self._count_genus_hits()
            self._populate_annotation_dicts()

        def _populate_annotation_dicts(self):
            # here we have been through each of the 10 hits for the sample and we should populate the
            # sample_tax_dict and the other dicts if the sample is either coral or Symbiodinium
            # see which is the biggest the Sym count or the scler or the other
            
            # If all dicts are empty, then we have had no valid annotations for this seqeunce
            if not any([self.other_genus_count_dict, self.symbiodiniaceae_genus_count_dict, self.coral_genus_count_dict]):
                return

            o_count = sum(self.other_genus_count_dict.values())
            s_count = sum(self.symbiodiniaceae_genus_count_dict.values())
            c_count = sum(self.coral_genus_count_dict.values())
            
            # Find most abundant taxa and populate respective parent dictionary
            if o_count > max(s_count, c_count):
                # log the annotation tup
                self.sample_genus_dict[self.blasted_seq_name] = [a[0] for a in sorted(self.other_genus_count_dict.items(), key=lambda x: x[1], reverse=True)][0]
            
            elif c_count > max(o_count, s_count):
                # log the annotation tup
                self.sample_genus_dict[self.blasted_seq_name] = [a[0] for a in sorted(self.coral_genus_count_dict.items(), key=lambda x: x[1], reverse=True)][0]
                # log the scleractinian genus
                self.coral_dict[self.blasted_seq_name] = [a[0] for a in sorted(self.coral_genus_count_dict.items(), key=lambda x: x[1], reverse=True)][0][0]
            
            elif s_count > max(o_count, c_count):
                # log the annotation tup
                self.sample_genus_dict[self.blasted_seq_name] = [a[0] for a in sorted(self.symbiodiniaceae_genus_count_dict.items(), key=lambda x: x[1], reverse=True)][0]
                # log the symbiodiniaceae genus
                self.symbiodiniceae_dict[self.blasted_seq_name] = [a[0] for a in sorted(self.symbiodiniaceae_genus_count_dict.items(), key=lambda x: x[1], reverse=True)][0][0]

        def _count_genus_hits(self):
            """For each of the 10 results, check to see which genus it was matched to and 
            count to see which genus was most commonly matched to

            log the order, family and genus so that we can output this info inthe
            annotation tables
            """
            for comp_list in self.super_results_list:
                # Get the taxonomic annotations for the sequence
                try:
                    genus_level, family_level, order_level = self._get_taxa_designation_from_staxid(
                        staxid=comp_list[6])
                except:
                    # there are some tax_ids that we don't seem to be able to find
                    # for the time being we will just continue over this seqeunce
                    continue
                if False in [genus_level, family_level, order_level]:
                    # Then we were unable to find one or more of the annotations
                    continue
                
                # Log the match
                annotation_tup = (genus_level, family_level, order_level)
                if family_level == 'Symbiodiniaceae':
                    self.symbiodiniaceae_genus_count_dict[annotation_tup] += 1
                elif order_level == 'Scleractinia' or order_level == 'Anthoathecata':
                    self.coral_genus_count_dict[annotation_tup] += 1
                else:
                    self.other_genus_count_dict[annotation_tup] += 1

        def _get_taxa_designation_from_staxid(self, staxid, tax_level_list=['genus', 'family', 'order']):
            # I have set this up so that you can either feed in the already made node_dict and name_dict or
            # you they can be generated by this method. If you are running this method many times it will make
            # sense to run the above generate_taxa_designation_from_staxid_dicts methods
            # to generate the dicts to feed into this
            # This will take in a list as its argument and find the levels that are listed in the list
            list_to_return = [False for i in tax_level_list]

            # now we can do the searching
            # go through staxid nodes until we get to the tax_level required
            # then grab the name
            while True:
                if staxid == '1':
                    # If we get here then we have not been able to find names for each of the
                    # taxonomic levels for some reason and we should return the list containing a False
                    # for those levels for which an identity could not be found
                    return list_to_return
                current_tax_level = self.node_dict[staxid][1]
                if current_tax_level in tax_level_list:
                    # then this is the taxonomic level we want, grab the name and return
                    list_to_return[tax_level_list.index(current_tax_level)] = self.names_dict[staxid]
                    if False not in list_to_return:
                        # then we have found all of the taxa_levels
                        return list_to_return
                    else:
                        # else we have not and we should continue the procession through the tax ids
                        staxid = self.node_dict[staxid][0]
                else:
                    staxid = self.node_dict[staxid][0]

class IndividualMothurQC:
    def __init__(self, readset, fwd_path, rev_path, qc_base_dir, num_proc=20):
        self.readset = readset
        self.fwd_path = fwd_path
        self.rev_path = rev_path
        self.sample_qc_dir = os.path.join(qc_base_dir, self.readset)
        os.makedirs(self.sample_qc_dir, exist_ok=True)
        self.num_proc=num_proc
        self._write_out_oligo_file()
        self.stability_file_path = self._write_out_file_file()
        self.mothur_batch_file_path = self._write_out_mothur_batch_file()
        
    def start_qc(self):
        result = subprocess.run(['mothur', self.mothur_batch_file_path])
        print(f'Mothur QC complete for {self.readset}')
        # To save on space we will compress all output files
        for file_to_compress in os.listdir(self.sample_qc_dir):
            subprocess.run(['gzip', os.path.join(self.sample_qc_dir, file_to_compress)])
        # Write out a log file to show that qc has been completed
        with open(os.path.join(self.sample_qc_dir, 'qc_complete.txt'), 'w') as f:
            f.write('qc_complete.txt\n')
        # Output file examples:
        # /home/humebc/projects/tara/full_18s_data/seq_qc/TARA_CO-0001684/stability.trim.contigs.good.unique.abund.pcr.names
        # /home/humebc/projects/tara/full_18s_data/seq_qc/TARA_CO-0001684/stability.trim.contigs.good.unique.abund.pcr.unique.fasta

    def _write_out_oligo_file(self):
        oligo_file_path = os.path.join(self.sample_qc_dir, 'primers.oligos')
        oligo_file = ['forward\tTTGTACACACCGCCC', 'reverse\tCCTTCYGCAGGTTCACCTAC']
        with open(oligo_file_path, 'w') as f:
            for line in oligo_file:
                f.write(f'{line}\n')

    def _write_out_file_file(self):
        """This is the file that will be supplied to make the contigs"""
        stability_file = [f'{self.fwd_path}\t{self.rev_path}']
        stability_file_path = os.path.join(self.sample_qc_dir, 'stability.files')

        # write out stability file
        with open(stability_file_path, 'w') as f:
            for line in stability_file:
                f.write(f'{line}\n')

        return stability_file_path

    def _write_out_mothur_batch_file(self):
        base = f'{self.sample_qc_dir}/stability'
        mothur_batch_file = [
            f'set.dir(input={self.sample_qc_dir})',
            f'set.dir(output={self.sample_qc_dir})',
            f'make.contigs(file={self.stability_file_path}, processors={self.num_proc})',
            f'summary.seqs(fasta={base}.trim.contigs.fasta)',
            f'screen.seqs(fasta={base}.trim.contigs.fasta, maxambig=0)',
            f'summary.seqs(fasta={base}.trim.contigs.good.fasta)',
            f'unique.seqs(fasta={base}.trim.contigs.good.fasta)',
            f'summary.seqs(fasta={base}.trim.contigs.good.unique.fasta, '
            f'name={base}.trim.contigs.good.names)',
            f'split.abund(cutoff=2, fasta={base}.trim.contigs.good.unique.fasta, '
            f'name={base}.trim.contigs.good.names)',
            f'summary.seqs(fasta={base}.trim.contigs.good.unique.abund.fasta, '
            f'name={base}.trim.contigs.good.abund.names)',
            f'summary.seqs(fasta={base}.trim.contigs.good.unique.rare.fasta, '
            f'name={base}.trim.contigs.good.rare.names)',
            f'pcr.seqs(fasta={base}.trim.contigs.good.unique.abund.fasta, '
            f'name={base}.trim.contigs.good.abund.names, '
            f'oligos={self.sample_qc_dir}/primers.oligos, pdiffs=2, rdiffs=2, processors={self.num_proc})',
            f'unique.seqs(fasta={base}.trim.contigs.good.unique.abund.pcr.fasta, '
            f'name={base}.trim.contigs.good.abund.pcr.names)',
            f'summary.seqs(fasta={base}.trim.contigs.good.unique.abund.pcr.fasta, '
            f'name={base}.trim.contigs.good.abund.pcr.names)'
        ]

        mothur_batch_file_path = os.path.join(self.sample_qc_dir, 'mothur_batch_file')

        # write out batch file
        with open(mothur_batch_file_path, 'w') as f:
            for line in mothur_batch_file:
                f.write('{}\n'.format(line))

        return mothur_batch_file_path


if __name__ == "__main__":
    esa = EighteenSProcessing()