#!/usr/bin/env nextflow
/* The purpose of this pipeline will be to go from pull down of the tara 18S coral seq data
to producing distance matrices and cluster calls for each of the species*/

// We will start with pulling down the coral fastq.gz sequencing files
// There are many additional sequencing files for the 18S on the genoscope FTP
// that we want to avoid pulling down if possible. For example those that relate to the
// water samples
basedir = "/home/humebc/projects/tara/full_18s_data"
bin_dir = "${basedir}/bin"
// The directory to hold all of the fastq.gz coral 18s seq files
raw_seq_file_out_dir = "/home/humebc/projects/tara/full_18s_data/seq_files"

// Here we will run the pull_down_seq_files.py that will log in to the genoscope FTP server
// and walk the directories in search of the FTP files
process PullDownSeqFilesFromFTP{
    cache 'lenient'
    tag 'ftp pull down'
    conda "envs/18s.yaml"
    publishDir path: raw_seq_file_out_dir
    storeDir raw_seq_file_out_dir
    
    output: 
    tuple file("TARA_CO-*.fastq.gz"), file("pull_down_stats.txt") into ch_raw_seq_files

    script:
    """
    python3 ${bin_dir}/pull_down_seq_files.py
    """
}