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
            if self.sample_provenance_df.at[ser['barcode_id'], 'SAMPLE ENVIRONMENT, short'] == 'C-CORAL':
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

    def _check_if_qc_already_complete(self, sample_name):
        return os.path.exists(os.path.join(self.qc_dir, sample_name, 'qc_complete.txt'))

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
            if self.sample_provenance_df.at[ser['barcode_id'], 'SAMPLE ENVIRONMENT, short'] == 'C-CORAL':
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
            barcode_id = self.fastq_info_df.loc[readset]['barcode_id']
            tax = self.sample_provenance_df.at[barcode_id, 'SAMPLE ENVIRONMENT, short']
            # print(f'{current_process().name}: barcode is {barcode_id}')
            # print(f'{current_process().name}: sample prov tax is {tax}')
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