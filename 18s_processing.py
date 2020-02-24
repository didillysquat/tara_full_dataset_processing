"""Python script responsible for the processing of the 18s seq data to produce distance matrices
of the coral samples"""
import os
import sys
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool
import subprocess

class EighteenSAnalysis:
    def __init__(self):
        self.root_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.seq_dir = os.path.join(self.root_dir, 'seq_files')
        # The directory where QC files will be written
        self.qc_dir = os.path.join(self.root_dir, 'seq_qc')
        # The directory where the finalised post qc and post taxa screening files will be written
        self.processed_seqs_dir = os.path.join(self.root_dir, 'processed_seqs_dir')
        self.sample_provenance_path = os.path.join(self.root_dir, "tara_samples_provenance.csv")
        self.sample_provenance_df = self._make_sample_provenance_df()
        # The main info df that we will use
        # Sample name as key, fwd and rev path to seq files, coral species
        self.info_df = MakeInfoDF(seq_dir=self.seq_dir, sample_provenance_df=self.sample_provenance_df).make_info_df()

        # Perform the QC for the samples
        # Mothur processing that will be fed into the taxonomic screening
        foo = 'bar'

    def _make_sample_provenance_df(self):
        # The SAMPLE ID is in the sequencing file name so we will be able to use this to get the latitute
        df = pd.read_csv(self.sample_provenance_path)
        df.set_index(keys='SAMPLE ID', drop=True, inplace=True)
        df.rename(columns={'EVENT latitude start (North)': 'lat', 'EVENT longitude start (East)': 'lon'}, inplace=True)
        return df

class SequenceQC:
    """Quality control of the sequences using Mothur"""
    def __init__(self, info_df, qc_dir):
        self.info_df = info_df
        self.qc_dir = qc_dir
        self.mothur_exe = "/home/humebc/phylogeneticSoftware/mothur1.40/mothur/mothur"
        # the list that we will map to the mothur qc
        self.apply_list = self._make_apply_list()

    def do_qc(self):
        """ for each of the samples in the info df we will make a directory
        in the qc_dir. We will produce all of qc files in this dir
        """

        # We should aim to use the Pool version of multiprocessing.
        # We will aim to provide a list lists where each list is n=3 or the sample name,
        # the fwd path and rev path.
        # TODO actually to simplify this let's just do the multiprocessing within the mothur
        # i.e. lets just do one sample at a time with.
        with Pool(5) as p:
            p.map(self._mothur_qc, self.apply_list)

    def _make_apply_list(self):
        apply_list = []
        for sample_name, ser in self.info_df.iterrows():
            apply_list.append([sample_name, ser['fwd_path'], ser['rev_path']])
        return apply_list

    def _mothur_qc(self, sample_info):
        sample_name, fwd_path, rev_path = sample_info
        # make a dir in the qc_dir


        # write out an oligo file and a file file
        self._write_out_oligo_file()

class IndividualMothurQC:
    def __init__(self, sample_info, qc_base_dir):
        self.sample_name, self.fwd_path, self.rev_path = sample_info
        self.sample_qc_dir = os.path.join(qc_base_dir, self.sample_name)
        os.makedirs(self.sample_qc_dir, exist_ok=True)
        self._write_out_oligo_file()
        self.stability_file_path = self._write_out_file_file()
        self.mothur_batch_file_path = self._write_out_mothur_batch_file()

    def start_qc(self):
        subprocess.run(['mothur', self.mothur_batch_file_path])
        print(f'Mothur QC complete for {self.sample_name}')

    def _write_out_oligo_file(self):
        oligo_file_path = os.path.join(self.sample_qc_dir, 'primers.oligos')
        oligo_file = ['forward\tTTGTACACACCGCCC', 'reverse\tCCTTCYGCAGGTTCACCTAC']
        with open(oligo_file_path, 'w') as f:
            for line in oligo_file:
                f.write(f'{line}\n')s

    def _write_out_file_file(self):
        """This is the file that will be supplied to make the contigs"""
        stability_file = [f'{self.fwd_path}\t{self.rev_path}']
        stability_file_path = os.join(self.sample_qc_dir, 'stability.files')

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
            f'make.contigs(file={self.stability_file_path}, processors=20)',
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
            f'oligos={self.sample_qc_dir}/primers.oligos, pdiffs=2, rdiffs=2, processors=20)',
            f'unique.seqs(fasta={base}.trim.contigs.good.unique.abund.pcr.fasta, '
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
                    fwd_path = path_lists[i][0]
                    rev_path = path_lists[i][1]
                    info_df_dict[f'{sample_name}_{i}'] = [fwd_path, rev_path, species]
            else:
                info_df_dict[sample_name] = [path_lists[0][0], path_lists[0][1], species]
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
