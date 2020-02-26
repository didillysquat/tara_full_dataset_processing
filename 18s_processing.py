"""Python script responsible for the processing of the 18s seq data to produce distance matrices
of the coral samples"""
import os
import sys
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool
import subprocess
import compress_pickle

class EighteenSAnalysis:
    def __init__(self):
        self.root_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.seq_dir = os.path.join(self.root_dir, 'seq_files')
        # The directory where QC files will be written
        self.qc_dir = os.path.join(self.root_dir, 'seq_qc')
        self.cache_dir = os.path.join(self.root_dir, 'cache')
        os.makedirs(self.cache_dir, exist_ok=True)
        # The directory where the finalised post qc and post taxa screening files will be written
        self.processed_seqs_dir = os.path.join(self.root_dir, 'processed_seqs_dir')
        self.sample_provenance_path = os.path.join(self.root_dir, "tara_samples_provenance.csv")
        self.sample_provenance_df = self._make_sample_provenance_df()
        # The main info df that we will use
        # Sample name as key, fwd and rev path to seq files, coral species
        self.info_df = self._make_info_df()
        
        # Perform the QC for the samples
        # Mothur processing that will be fed into the taxonomic screening
        seq_qc = SequenceQC(info_df=self.info_df, qc_dir=self.qc_dir).do_mothur_qc()

        # After the mothur QC we will have a set of .name and .fasta files for each sample
        # Now we need to annotate these sequences according to taxonomy
        # The results of this are a set of three pickled out annotation dictionaries
        # One is all sequnces, one is only coral, one is only symbiodiniaceae
        seq_qc.do_bast_qc()

        # At this point I think we have the sequences as we want them.
        # We can now work from individual methods.
        # Probably a good idea to run some PCoA work. Just to see how what the spread looks like
        # for the sequnces of the genus in question.

        # Once we've looked at that we can get an idea of whether it may be a good idea to cluster or not
        # Then it will be time to compute the distance matrices per species.
        # We should probably complete a set of unifrac and a set of braycurtis.
        
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

class SequenceQC:
    """Quality control of the sequences using Mothur firstly and then running through BLAST
    to remove all non-coral sequences"""
    def __init__(self, info_df, qc_dir):
        self.info_df = info_df
        self.qc_dir = qc_dir
        # These dicts will be used for getting the taxonomic level, and name of a
        # taxonomic match.
        self.node_dict, self.name_dict = self._generate_taxa_and_name_dicts()

    def do_mothur_qc(self):
        """ 
        for each of the samples in the info df we will make a directory
        in the qc_dir. We will produce all of qc files in this dir
        """
        # For some reason the make.contigs doesn't seem to be using multiple processors and thus
        # this is wasting a lot of time. As such we should try to multiprocess this and use 1
        # process for each. Let's try to implement with Pool
        # We will need a list to pass to Pool for doing the Apply with
        apply_list = []
        for sample_name, ser in self.info_df.iterrows():
            if self._check_if_qc_already_complete(sample_name):
                continue
            else:
                apply_list.append([sample_name, ser['fwd_path'], ser['rev_path'], self.qc_dir])
        
        if apply_list:
            # NB the Pool does not play well with the debugger that VS code is using.
            # Pool does work OK if run normally (i.e. no debug)
            # By having the Pool operation covered by the if apply_list:
            # we should still be able to debug as we won't have to aceess Pool.
            with Pool(24) as p:
                p.map(self._init_indi_mothur_qc, apply_list)

    def _init_indi_mothur_qc(self, sample_info):
        """We have to pass the Pool a single function. This is the function we will use to instantiate 
        an IndividualMothurQC class and run the qc."""
        IndividualMothurQC(
            sample_name=sample_info[0], 
            fwd_path=sample_info[1], 
            rev_path=sample_info[2], 
            qc_base_dir=sample_info[3],
            num_proc=1).start_qc()

    def _check_if_qc_already_complete(self, sample_name):
        return os.path.exists(os.path.join(self.qc_dir, sample_name, 'qc_complete.txt'))

    def do_bast_qc(self):
        """
        For every sample run a blast against the nt database.
        For each sequence get 10 blast results.
        Look at the matches and see which taxon was highest represented.
        We will produce three annotation dictionaries.
        In the first, sample_annotation_dict, taxa will be annotated to a genus level
        unless they are Scleractinia, Anthoathecata or Symbiodiniaceae in which case
        they will be annotated thusly.
        
        We will also produce two additional dicts. One for only Symbiodiniaceae and one for only Scleractinia or Anthoathecata.
        In both of these dicts, genus level annotations will be used for the sequences.

        We will run this in serial but use multiple processors for running the blast.
        """
        for sample_name in self.info_df.index:
            IndividualTaxAnnotationQC(
                node_dict=self.node_dict, 
                names_dict=self.name_dict, 
                qc_dir=self.qc_dir, 
                sample_name=sample_name
                ).do_taxonomy_annotation()

    def _generate_taxa_and_name_dicts(self,
                taxdump_dir='/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload/taxdump'):
            # read in the .nodes file. This file tells us which tax level the node is and which node is the parent level
            with open('{}/nodes.dmp'.format(taxdump_dir), 'r') as f:
                node_file = [line.rstrip() for line in f]
            # now make a dict from this where key is the tax id and the value is a tup where 0 = parent 1 = tax level
            node_dict = {line.split('\t|\t')[0]: (line.split('\t|\t')[1], line.split('\t|\t')[2]) for line in node_file}

            # next read in the names file. This file hold the name of the node.
            with open('{}/names.dmp'.format(taxdump_dir), 'r') as f:
                name_file = [line.rstrip() for line in f]

            # now make a dict from the names file where the key is the staxid and the value is the name
            name_dict = {line.split('\t|\t')[0]: line.split('\t|\t')[1] for line in name_file if
                         line.split('\t|\t')[3].replace('\t|', '') == 'scientific name'}

            return node_dict, name_dict

class IndividualTaxAnnotationQC:
    def __init__(self, node_dict, names_dict, qc_dir, sample_name):
        self.node_dict = node_dict
        self.names_dict = names_dict
        self.genus = None
        self.family = None
        self.order = None
        self.fasta_path = os.path.join(qc_dir, sample_name, 'stability.trim.contigs.good.unique.abund.pcr.unique.fasta')
        self.blast_out_path = os.path.join(qc_dir, sample_name, 'blast.out')
        self.output_format = "6 qseqid sseqid staxids evalue pident qcovs staxid stitle ssciname"
        self.blast_out_pickle_path = os.path.join(qc_dir, sample_name, 'blast.out.p.bz')
        self.blast_output_file = self._get_or_do_blast()
        # This is a dict with seq name as key, and the 10 blast results as a list of lists as value
        self.blast_output_dict = self._make_blast_out_dict()
        # We will have three dictionaries that will hold the results of the taxonomic annotation
        # This dict will hold the annotation for all sequenecs
        self.sample_annotation_dict = {}
        self.coral_annotation_dict = {}
        self.symbiodiniaceae_annotation_dict = {}
        # The paths that we will pickle out the annotation dictionaries to
        self.sample_annotation_dict_pickle_path = os.path.join(qc_dir, sample_name, 'sample_annotation_dict.p.bz')
        self.coral_annotation_dict_pickle_path = os.path.join(qc_dir, sample_name, 'coral_annotation_dict.p.bz')
        self.coral_symbiodiniaceae_dict_pickle_path = os.path.join(qc_dir, sample_name, 'symbiodiniaceae_annotation_dict.p.bz')

    def do_taxonomy_annotation(self):
        for blasted_seq, result_list in self.blast_output_dict.items():
            self.IndividualTaxAnnotationQCSub(
                blasted_seq_name=blasted_seq, 
                results_list=result_list, 
                sample_dict=self.sample_annotation_dict, 
                symbiodiniaceae_dict=self.symbiodiniaceae_annotation_dict, 
                coral_dict=self.coral_annotation_dict
                ).process_matches()
        
        # Here we have populated the taxonomy dictionaries for the sample in question
        # Now pickle out the dictionaries
        compress_pickle.dump(self.sample_annotation_dict, self.sample_annotation_dict_pickle_path)
        compress_pickle.dump(self.coral_annotation_dict, self.coral_annotation_dict_pickle_path)
        compress_pickle.dump(self.symbiodiniaceae_annotation_dict, self.coral_symbiodiniaceae_dict_pickle_path)

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

    def _make_blast_out_dict(self):
        # now create a dict that is the list of 10 result items using the seq name as key
        blast_output_dict = defaultdict(list)
        for output_line in blast_output_file:
            components = output_line.split('\t')
            blast_output_dict[components[0]].append(components)
        return blast_output_dict

    def _get_or_do_blast(self):
        """
        Check to see if the blast has already been performed and can be read
        in from file. Else perform the blast from scratch.
        """
        if os.path.isfile(self.blast_out_pickle_path):
            return compress_pickle.load(self.blast_out_pickle_path)
        else:
            # Run local blast
            subprocess.run(
                ['blastn', '-out', self.blast_out_path, '-outfmt', self.output_format, '-query', self.fasta_path, '-db', 'nt',
                    '-max_target_seqs', '10', '-num_threads', '20'])

            # Read in blast output
            with open(self.blast_out_path, 'r') as f:
                blast_output_file = [line.rstrip() for line in f]

            # pickle out the blast results here to possibly save us time in the future
            compress_pickle.dump(blast_output_file, self.blast_out_pickle_path)
            return blast_output_file
    
    class IndividualTaxAnnotationQCSub:
        """This is a utility class that is used for counting the genus hits for a given sequence"""
        def __init__(self, blasted_seq_name, results_list, sample_dict, symbiodiniaceae_dict, coral_dict):
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
            
            o_count = len(self.other_genus_count_dict.items())
            s_count = len(self.symbiodiniaceae_genus_count_dict.items())
            c_count = len(self.coral_genus_count_dict.items())
            
            # Find most abundant taxa and populate respective parent dictionary
            if o_count > max(s_count, c_count):
                self.sample_genus_dict[self.blasted_seq_name] = [a[0] for a in sorted(self.other_genus_count_dict.items(), key=lambda x: x[1], reverse=True)][0]
            
            elif c_count > max(o_count, s_count):
                self.sample_genus_dict[self.blasted_seq_name] = 'Scleractinia_Anthoathecata'
                # for the scleractinian dictionary we should associate to the most abundant genus
                self.coral_dict[self.blasted_seq_name] = [a[0] for a in sorted(self.coral_genus_count_dict.items(), key=lambda x: x[1], reverse=True)][0]
            
            elif s_count > max(o_count, c_count):
                # Then this is a Symbiodiniaceae and we should store the genus in the specific dictionary
                self.sample_genus_dict[self.blasted_seq_name] = 'Symbiodiniaceae'            
                self.symbiodiniceae_dict[self.blasted_seq_name] = [a[0] for a in sorted(self.symbiodiniaceae_genus_count_dict.items(), key=lambda x: x[1], reverse=True)][0]

        def _count_genus_hits(self):
            """For each of the 10 results, check to see which genus it was matched to and 
            count to see which genus was most commonly matched to"""
            for comp_list in self.super_results_list:
                # Get the taxonomic annotations for the sequence
                try:
                    genus_level, family_level, order_level = self._get_taxa_designation_from_staxid(
                        staxid=comp_list[6], tax_level_list=['genus', 'family', 'order'])
                except:
                    # there are some tax_ids that we don't seem to be able to find
                    # for the time being we will just continue over this seqeunce
                    continue
                if False in [genus_level, family_level, order_level]:
                    # Then we were unable to find one or more of the annotations
                    continue
                
                # Log the match
                if family_level == 'Symbiodiniaceae':
                    self.symbiodiniaceae_genus_count_dict[genus_level] += 1
                elif order_level == 'Scleractinia' or order_level == 'Anthoathecata':
                    self.coral_genus_count_dict[genus_level] += 1
                else:
                    self.other_genus_count_dict[genus_level] += 1

class IndividualMothurQC:
    def __init__(self, sample_name, fwd_path, rev_path, qc_base_dir, num_proc=20):
        self.sample_name = sample_name
        self.fwd_path = fwd_path
        self.rev_path = rev_path
        self.sample_qc_dir = os.path.join(qc_base_dir, self.sample_name)
        os.makedirs(self.sample_qc_dir, exist_ok=True)
        self.num_proc=num_proc
        self._write_out_oligo_file()
        self.stability_file_path = self._write_out_file_file()
        self.mothur_batch_file_path = self._write_out_mothur_batch_file()
        self.mothur_exe_path = "/home/humebc/phylogeneticSoftware/mothur1.43.0/mothur/mothur"
        
    def start_qc(self):
        subprocess.run([self.mothur_exe_path, self.mothur_batch_file_path])
        print(f'Mothur QC complete for {self.sample_name}')
        # We are running very short on space so let's clean up
        for file_to_del in os.listdir(self.sample_qc_dir):
            if file_to_del not in ['stability.trim.contigs.good.unique.abund.pcr.names', 'stability.trim.contigs.good.unique.abund.pcr.unique.fasta']:
                os.remove(os.path.join(self.sample_qc_dir, file_to_del))
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
            f'screen.seqs(fasta={base}.trim.contigs.fasta, maxambig=0, maxhomop=5)',
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

class MakeInfoDF:
    def __init__(self, seq_dir, sample_provenance_df):
        self.paired_files_dict = defaultdict(list)
        self.seq_dir = seq_dir
        self.seq_files = [file for file in os.listdir(self.seq_dir) if 'fastq' in file]
        self.sample_provenance_df = sample_provenance_df
        # Get a list of the samples that we have seq files for
        self.sample_names = self._get_sample_names_list()

    def make_info_df(self):
        # for each of the sample names, find the respective seq files, pair them up and populate the paired_files dict
        self._find_sample_respective_seq_files()

        # here we have the dictionary that tells us which sample names relate to which seqfiles populated
        # Now populate the info dict splitting up those samples that have multiple pairs of seq files
        # Dict that we will create info df from eventually
        return self._return_info_df()

    def _return_info_df(self):
        info_df_dict = {}
        for sample_name in self.paired_files_dict.keys():
            species = self.sample_provenance_df.at[sample_name, 'SAMPLE MATERIAL taxonomy'].split(' ')[1]
            path_lists = self.paired_files_dict[sample_name]
            if len(path_lists) > 1:
                for i in range(len(path_lists)):
                    fwd_path = os.path.join(self.seq_dir, path_lists[i][0])
                    rev_path = os.path.join(self.seq_dir, path_lists[i][1])
                    info_df_dict[f'{sample_name}_{i}'] = [fwd_path, rev_path, species]
            else:
                info_df_dict[sample_name] = [os.path.join(self.seq_dir, path_lists[0][0]), os.path.join(self.seq_dir, path_lists[0][1]), species]
        return pd.DataFrame.from_dict(info_df_dict, orient='index', columns=['fwd_path', 'rev_path', 'species'])

    def _find_sample_respective_seq_files(self):
        for sample_name in self.sample_names:
            matching_seq_files = []
            for seq_file in self.seq_files:
                if sample_name in seq_file:
                    matching_seq_files.append(seq_file)
            # here we have a list of all of the seq files
            # if there are only two, then this is simple
            if not matching_seq_files:
                raise RuntimeError(f'No sequencing files found for sample {sample_name}')
            if len(matching_seq_files) == 2:
                if '_R1.' in matching_seq_files[0]:
                    self.paired_files_dict[sample_name] = [[matching_seq_files[0], matching_seq_files[1]]]
                else:
                    self.paired_files_dict[sample_name] = [[matching_seq_files[1], matching_seq_files[0]]]
            else:
                # if there are more than two then we need to pair them up
                for seq_file in [seq_file for seq_file in matching_seq_files if '_R1.' in seq_file]:
                    self.paired_files_dict[sample_name].append([seq_file, seq_file.replace('_R1.', '_R2.')])

    def _get_sample_names_list(self):
        sample_names = set()
        for seq_file in self.seq_files:
            sample_names.add('_'.join(seq_file.split('_')[:2]))
        sample_names = list(sample_names)
        return sample_names

esa = EighteenSAnalysis()
