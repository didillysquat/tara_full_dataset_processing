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
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import time
import numpy as np
import operator
import matplotlib.gridspec as gridspec
from base_18s import EighteenSBase

# TODO later we will be able to write this as a subclass of EighteenSProcessing
# But for the time being we don't want to interfere with any of that code because 
# it is currently running to create the 18S taxa annotations.
class EighteenSAnalysis(EighteenSBase):
    def __init__(self):
        super().__init__()
        # This will hold the additional variables used only in this analysis class
        self.island_site_dict = self._determine_sites_and_island()
        self.islands = sorted(list(self.island_site_dict.keys()))
        self.host_species = ["Porites", "Millepora", "Pocillopora"]

        # Do seq consolidation. This will make and pickle out various abundance
        # dictionaries that we will use in the plotting and in the distance matrix creation
        SeqConsolidator(qc_dir=self.qc_dir, cache_dir = self.cache_dir, info_df=self.info_df).do_consolidation()

    def _determine_sites_and_island(self):
        # For each sample pull out the island and site from the info_df
        # Dict that will have island names as key and a list of sites as value
        island_site_dict = defaultdict(set)
        for sample_index in self.info_df.index:
            island = self.info_df.at[sample_index, 'island']
            site = self.info_df.at[sample_index, 'site']
            island_site_dict[island].add(site)
        return island_site_dict

    def do_stacked_bar_plots(self, plot_type, in_sample_cutoff=None):
        sbp = StackedBarPlotter(
            plot_type=plot_type, islands=self.islands,
            island_site_dict=self.island_site_dict, host_species=self.host_species, 
            fig_output_dir=self.fig_output_dir, qc_dir=self.qc_dir, info_df=self.info_df, cache_dir=self.cache_dir,
            in_sample_cutoff=in_sample_cutoff)
        sbp.plot()


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
    def __init__(self, qc_dir, cache_dir, info_df):
        self.qc_dir = qc_dir
        self.cache_dir = cache_dir
        self.info_df = info_df
        self.consolidated_host_seqs_rel_abundance_dict = None
        self.coral_blasted_seq_to_consolidated_seq_dict = {}

    def do_consolidation(self):
        if not self._check_if_consolidation_already_complete():
            # Firstly get a set of all sequences that are of one of the three genera
            # In the process of doing this, pickle out the all seq abunance dict for
            # the sample
            self.consolidated_host_seqs_rel_abundance_dict = self._make_host_seqs_dict()
            # Then create the consolidation path and
            # Consolidate the sequence insitu in the self.host_seqs_dict
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
        for sample_name in self.info_df.index:
            sys.stdout.write(f'\r{sample_name}')
            sample_qc_dir = os.path.join(self.qc_dir, sample_name)
            # Load the already created abundance dictionary
            rel_all_seq_abundance_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'rel_all_seq_abundance_dict.p.bz'))
            # Load the already created taxonomy annotation dictoinaries
            sample_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'sample_annotation_dict.p.bz'))
            coral_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'coral_annotation_dict.p.bz'))
            seq_name_to_seq_dict = self._make_seq_name_to_seq_dict(sample_name)
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
            if seq_name == 'GTCGCTACTACCGATTGAATGGTTTAGTGAGGCCTCCTGACTGGCGCCGACACTCTGTCTCGTGCAGAGAGTGGGAGGCCGGGAAGTTGTTCAAACTTGATCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTC':
                foo = 'bar'
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
        search_seq = 'GTCGCTACTACCGATTGAATGGTTTAGTGAGGCCTCCTGACTGGCGCCGACACTCTGTCTCGTGCAGAGAGTGGGAGGCCGGGAAGTTGTTCAAACTTGATCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTC'

        for query_seq, match_seq in consolidation_path_list:
            if query_seq == search_seq:
                foo = 'asdf'
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
        for sample_name in self.info_df.index:
            print(f'Getting coral sequences from sample {sample_name}')
            sample_qc_dir = os.path.join(self.qc_dir, sample_name)
            abs_all_seq_abundance_dict = self._make_abund_dict_from_names_path(sample_name)
            seq_name_to_seq_dict = self._make_seq_name_to_seq_dict(sample_name)
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

    def _make_seq_name_to_seq_dict(self, sample_name):
        with open(os.path.join(self.qc_dir, sample_name, 'stability.trim.contigs.good.unique.abund.pcr.unique.fasta'), 'r') as f:
            fasta_file_as_list = [line.rstrip() for line in f]
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
                    search_seq = 'GTCGCTACTACCGATTGAATGGTTTAGTGAGGCCTCCTGACTGGCGCCGACACTCTGTCTCGTGCAGAGAGTGGGAGGCCGGGAAGTTGTTCAAACTTGATCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTC'
                    for tup in consolidation_path_list:
                        if search_seq in tup:
                            if tup[0] == search_seq:
                                print('search seq is tup 0')
                            else:
                                print('search seq is tup 1')
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

    def _make_abund_dict_from_names_path(self, sample_name):
        with open(os.path.join(self.qc_dir, sample_name, 'stability.trim.contigs.good.unique.abund.pcr.names'), 'r') as f:
            name_file = [line.rstrip() for line in f]
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


class StackedBarPlotter:
    """
    This class will be responsible for plotting the stacked bar plots
    that will be arranged in a large 18 x 18 matrice of the islands sites and coral species.
    We will use this plot to get an overview of the sequencing results.
    """
    def __init__(
            self, plot_type, islands, island_site_dict, host_species,
            fig_output_dir, qc_dir, info_df, cache_dir, in_sample_cutoff):
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
        self.info_df = info_df
        self.qc_dir = qc_dir
        self.cache_dir = cache_dir
        self.islands = islands
        self.island_site_dict = island_site_dict
        self.host_species = host_species
        self.fig_output_dir = fig_output_dir
        self.in_sample_cutoff = in_sample_cutoff
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
            for j, lab in enumerate(['POR', 'MIL', 'POC']):
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
            return ['Porites', 'Millepora', 'Pocillopora', 'other_coral',
                    'Symbiodiniaceae', 'other_taxa', 'not_annotated'], col_dict
        elif self.plot_type == 'all_coral_genus':
            col_dict = {'Porites': '#FFFF00', 'Pocillopora': '#87CEFA', 'Millepora': '#FF6347',
                            'other_coral': '#C0C0C0'}
            return ['Porites', 'Millepora', 'Pocillopora', 'other_coral'], col_dict
        elif self.plot_type in ['all_coral_sequence', 'minor_coral_sequence']:
            col_dict = compress_pickle.load(os.path.join(self.cache_dir, 'all_coral_sequence_color_dict.p.bz'))
            return None, col_dict
        else:
            raise NotImplementedError()

    def plot(self):
        # we will go in order of: for island, for site, for species
        for island in self.islands:
            site_list = sorted(list(self.island_site_dict[island]))
            # Make the subgridspec here
            ax_row = int(self.islands.index(island)/6) + 1
            ax_col = self.islands.index(island)%6
            # Number of rows will be the number of sites + 1 for the title
            # Number of columns will be constant and the number of hosts
            sub_gs = self.gs[ax_row, ax_col].subgridspec(len(site_list) + 1, 3)
            # Put the island name in the title plot.
            title_ax = self.fig.add_subplot(sub_gs[0, :])
            self._do_island_title_plot(island=island, ax=title_ax)
            for site in site_list:
                for species in self.host_species:
                    # The last part to the row index is to incorporate the plot that we will do the naming in.
                    # The row will be the site
                    ax_row_index = site_list.index(site) + 1
                    # The col will always be the species
                    ax_col_index = self.host_species.index(species)
                    # In here we can do the actual plotting
                    ax = self.fig.add_subplot(sub_gs[ax_row_index, ax_col_index])
                    # Do the plotting for a given island, site, species set of samples
                    sbip = StackedBarIndiPlot(parent=self, ax=ax, island=island, 
                    site=site, species=species, in_sample_cutoff=self.in_sample_cutoff)
                    if sbip.samples:
                        sbip.do_plotting()
                    else:
                        self._remove_axes_but_allow_labels(ax)

        

        self.fig.suptitle(f'18s {self.plot_type}', fontsize=16)
        if self.in_sample_cutoff:
            svg_path = os.path.join(self.fig_output_dir, f'stacked_bar_18s_{self.plot_type}_in_sample_cutoff_{self.in_sample_cutoff}.svg')
        else:
            svg_path = os.path.join(self.fig_output_dir, f'stacked_bar_18s_{self.plot_type}.svg')
        print(f'\nWriting .svg to {svg_path}')
        plt.savefig(svg_path)
        if self.in_sample_cutoff:
            png_path = os.path.join(self.fig_output_dir, f'stacked_bar_18s_{self.plot_type}_in_sample_cutoff_{self.in_sample_cutoff}.png')
        else:
            png_path = os.path.join(self.fig_output_dir, f'stacked_bar_18s_{self.plot_type}.png')
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


class StackedBarIndiPlot(EighteenSBase):
    def __init__(self, parent, ax, island, site, species, in_sample_cutoff):
        super().__init__()
        self.parent = parent
        self.ax = ax
        self.island = island
        self.site = site
        self.species = species
        self.in_sample_cutoff = in_sample_cutoff
        # As well as the sample names, get the sample ids (i.e. C001).
        # We will plot these very small underneath the plots so that
        # we can see exactly which individual we are working with
        self.samples, self.sample_ids = self._get_sample_name_list()
        self.patches_list = []
        self.ind = 0
        self.color_list = []
        self.num_smp_in_this_subplot = len(self.samples)
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
        
    def _get_sample_name_list(self):
        init_sample_list = self.parent.info_df[
            (self.parent.info_df['island'] == self.island) &
            (self.parent.info_df['site'] == self.site) &
            (self.parent.info_df['species'] == self.species)
        ].index.values.tolist()


        id_list = []
        for sample_name in init_sample_list:
            if sample_name[-2] == '_':
                # then this is a techrep
                sample_base = '_'.join(sample_name.split('_')[:-1])
                sample_base_id = self.sample_provenance_df.at[sample_base, 'C###, F###, MA##, SG##']
                id_list.append(self._make_numeric(sample_base_id + '_' + sample_name.split('_')[-1]))
            else:
                id_list.append(self._make_numeric(self.sample_provenance_df.at[sample_name, 'C###, F###, MA##, SG##']))

        return init_sample_list, sorted(id_list)

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
        if os.path.isfile(os.path.join(self.cache_dir, 'seq_to_total_abund_dict.p.bz')):
            seq_to_total_abund_dict = compress_pickle.load(os.path.join(self.cache_dir, 'seq_to_total_abund_dict.p.bz'))
        else:
            seq_to_total_abund_dict = defaultdict(int)
            for sample_name in self.info_df.index:
                sample_qc_dir = os.path.join(self.qc_dir, sample_name)
                consolidated_host_seqs_abund_dict = compress_pickle.load(
                    os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
                for seq_name in consolidated_host_seqs_abund_dict.keys():
                    seq_to_total_abund_dict[seq_name] += 1
            compress_pickle.dump(seq_to_total_abund_dict, os.path.join(self.cache_dir, 'seq_to_total_abund_dict.p.bz'))
        return seq_to_total_abund_dict

    def _make_abundance_dicts(self):
        # If we are using an in_sample_cutoff (minimum number of samples a sequence must be found in)
        # load up the master abundance dictionary and use this to screen the sample_consolidated_abund_dict
        if self.in_sample_cutoff:
            seq_to_total_abund_dict = self._get_seq_to_total_abund_dict()
            threshold_set = {k for k, v in seq_to_total_abund_dict.items() if v > self.in_sample_cutoff}

        df_dict = {}
        for sample_name in self.samples:
            sample_qc_dir = os.path.join(self.parent.qc_dir, sample_name)
            sample_consolidated_abund_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'consolidated_host_seqs_abund_dict.p.bz'))
            # When we do the minor seqs only then we will want to get rid of the most abundant seq and work with this
            if self.parent.plot_type == 'all_coral_sequence':
                df_dict[sample_name] = sample_consolidated_abund_dict
            elif self.parent.plot_type == 'minor_coral_sequence':
                # Remove the most abundant sequence from the dict
                del sample_consolidated_abund_dict[max(sample_consolidated_abund_dict, key=sample_consolidated_abund_dict.get)]
                if self.in_sample_cutoff:
                    sample_consolidated_abund_dict = {
                        k: v for k, v in sample_consolidated_abund_dict.items() if k in threshold_set}
                tot = sum(sample_consolidated_abund_dict.values())
                sample_consolidated_abund_dict = {k: v/tot for k, v in sample_consolidated_abund_dict.items()}
                df_dict[sample_name] = sample_consolidated_abund_dict
            else:
                raise NotImplementedError
        return df_dict

    def _make_abundance_df(self):
        # Dict that we will populate and then use to make the abundance_df
        df_dict = {}
        for sample_name in self.samples:
            sample_qc_dir = os.path.join(self.parent.qc_dir, sample_name)
            # make a seq_name to abundance dict from the fasta and .names pair
            sample_abund_dict = self._make_abund_dict_from_names_path(sample_name=sample_name)
            # For the all_taxa, we will go sequence by sequence through the fasta file
            with open(
                    os.path.join(sample_qc_dir, 'stability.trim.contigs.good.unique.abund.pcr.unique.fasta'),
                    'r') as f:
                fasta_file_as_list = [line.rstrip() for line in f]
            fasta_names = [line.split('\t')[0][1:] for line in fasta_file_as_list if line[0] == '>']
            # then load the three dictionaries
            sample_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'sample_annotation_dict.p.bz'))
            coral_annotation_dict = compress_pickle.load(os.path.join(sample_qc_dir, 'coral_annotation_dict.p.bz'))

            sample_count_dict = {cat: 0 for cat in self.parent.plotting_categories}
            if self.parent.plot_type == 'all_taxa':
                self._log_abundances_all_taxa(sample_annotation_dict, sample_count_dict, sample_abund_dict, coral_annotation_dict, fasta_names)
            elif self.parent.plot_type == 'all_coral_genus':
                self._log_abundances_all_coral_genus(sample_annotation_dict, sample_count_dict, sample_abund_dict, coral_annotation_dict)
            else:
                raise NotImplementedError

            # Now add the collected abundances to the sample df_dict
            # Making them relative by dividing by the total of the sample_count_dd
            df_dict[sample_name] = [
                sample_count_dict[cat_key]/sum(sample_count_dict.values()) for
                cat_key in self.parent.plotting_categories
            ]

        # Now create the df from the df_dict
        return pd.DataFrame.from_dict(data=df_dict, orient='index', columns=self.parent.plotting_categories)

    def _log_abundances_all_taxa(
            self, sample_annotation_dict, sample_count_dict, sample_abund_dict, coral_annotation_dict, fasta_names):
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
            sample_count_dict[key] += sample_abund_dict[fasta_name]

    def _log_abundances_all_coral_genus(
            self, sample_annotation_dict, sample_count_dict, sample_abund_dict, coral_annotation_dict):
        for blasted_seq, annotation in sample_annotation_dict.items():
            if annotation[2] in ['Scleractinia', 'Anthoathecata']:
                # Then this is a coral seq and we should add the count to either one of the target genera
                # or to an other coral count
                # TODO this is where we can change our logic to according to what type of plot we are doing
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
                sample_count_dict[key] += sample_abund_dict[blasted_seq]

    def _make_abund_dict_from_names_path(self, sample_name):
        with open(os.path.join(
                self.parent.qc_dir, sample_name, 'stability.trim.contigs.good.unique.abund.pcr.names'), 'r') as f:
            name_file = [line.rstrip() for line in f]
        return {line.split('\t')[0]: len(line.split('\t')[1].split(',')) for line in name_file}

    def do_plotting(self):
        for sample_to_plot in self.samples:
            sys.stdout.write(f'\rPlotting sample: {self.island} {self.site} {self.species} {sample_to_plot}')
            if self.parent.plot_type in ['all_taxa', 'all_coral_genus']:
                self._plot_bars_from_df(sample_to_plot)
            else:
                self._plot_bars_from_dicts(sample_to_plot)
            self.ind += 1
        # If <= 10 samples, we want to set the x lim to 10 so that bars are constant width
        # Unless there are more than 10 samples then we want to use this number
        if len(self.samples) > 10:
            self._paint_rect_to_axes(max_num_smpls_in_subplot=len(self.samples))
        else:
            self._paint_rect_to_axes()

    def _plot_bars_from_dicts(self, sample_to_plot):
        bottom_div = 0
        sample_abund_dict = self.abundance_dicts[sample_to_plot]
        order_to_plot = [seq_name for seq_name in self.parent.ordered_seq_name_list if seq_name in sample_abund_dict]
        # In order of the master consolidated seqs
        for seq_name in order_to_plot:
            seq_rel_abund = sample_abund_dict[seq_name] 
            self.patches_list.append(Rectangle((self.ind - 0.5, bottom_div), 1, seq_rel_abund, color=self.parent.color_dict[seq_name]))
            self.color_list.append(self.parent.color_dict[seq_name])
            bottom_div += seq_rel_abund

    def _plot_bars_from_df(self, sample_to_plot):
        bottom_div = 0
        # In the order of the plotting categories
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
        self._remove_axes_but_allow_labels(x_labels=self.sample_ids)
        self.ax.set_ylabel(self.site, fontsize=5, labelpad=0)

    def _remove_axes_but_allow_labels(self, ax=None, x_labels=None):
        if ax is None:
            self.ax.set_frame_on(False)
            self.ax.set_yticks([])
            if x_labels:
                self.ax.set_xticks([])
                self.ax.set_xticks([_-0.25 for _ in range(len(self.sample_ids))])
                self.ax.set_xticklabels(self.sample_ids, fontsize=2, rotation=-90)
                self.ax.xaxis.set_tick_params(which='major',  length=0, pad=0)
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
    for cutoff in [150, 200, 250, 500, 1000 ]:
        EighteenSAnalysis().do_stacked_bar_plots(plot_type='minor_coral_sequence', in_sample_cutoff=cutoff)
    # for plot_type in ['all_taxa', 'all_coral_genus', 'all_coral_sequence', 'minor_coral_sequence']:
    #     if plot_type == 'minor_coral_sequence':
    #         for cutoff in [5,25,50,100]:
    #             EighteenSAnalysis().do_stacked_bar_plots(plot_type=plot_type, in_sample_cutoff=cutoff)
    #     else:
    #         EighteenSAnalysis().do_stacked_bar_plots(plot_type)
    # EighteenSAnalysis().do_stacked_bar_plots('all_taxa')