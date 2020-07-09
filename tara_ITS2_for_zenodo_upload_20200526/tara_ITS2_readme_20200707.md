# TARA Pacific ITS2 data release

## Table of contents

- [Introduction](#introduction)
- [Data subdirectories](#data-subdirectories)
  - Available for all datasets
    - [pre_med_seqs](#pre_med_seqs)
    - [post_med_seqs](#post_med_seqs)
    - [between_sample_distances](#between_sample_distances)
  - Available only for the coral dataset
    - [between_profile_distances](#between_profile_distances)
    - [its2_type_profiles](#its2_type_profiles)

## Introduction

The TARA Pacific ITS2 data release is organised into three subdirectories representing different types of samples:

- `20200326_tara_corals`: This directory contains results associated with samples taken directly from a coral animal
- `20200326_tara_non_corals`: This directory contains results for samples not taken directly from a coral animal and includes water- and sediment-derived samples.
- `20200326_tara_negatives`: This directory contains results from the negative control sequencing runs.

 For each sample type, datasets were created by submitting the associated sequencing files to the remote instance of the [SymPortal framework](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13004) ([symportal.org](www.symportal.org)). As such, for each sample type, the associated dataset has been processed through the [SymPortal quality control (QC) pipeline](https://github.com/didillysquat/SymPortal_framework/wiki/The-SymPortal-logic#sequence-quality-control).

 In addition to sequence quality control, the `20200326_tara_corals` dataset was included in a SymPortal analysis to predict ITS2 profiles. This dataset therefore contains some additional resources that are not available for the other two datasets. Prediction of ITS2 profiles is not possible for non-host-derived samples. For more details please refer to the [associated manuscript](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13004).

 Where technical replicates were available for a given sample-id (i.e. multiple pairs of fastq.gz files), only the largest replicate (combined size of the two fastq.gz files) was used.

 ***

## Data subdirectories

A directory tree for this data release (`tara_ITS2_dir_tree.txt`) can be found in the root directory of this release.

The following subdirectories exist for each of the datasets:

### pre_med_seqs

This directory contains the following files:

`pre_med_absolute_abundance_df.csv`: The count table detailing the absolute abundance of sequences returned from each sample after processing in the SymPortal QC pipeline, but before undergoing [Minimum Entropy Decomposition](http://merenlab.org/software/med/) that is conducted as standard when loading a dataset into the SymPortal database.

`pre_med_master_seqs.fasta`: Fasta file containing the DNA sequences referenced in the `pre_med_absolute_abundance_df.csv` count table.

### post_med_seqs

This directory contains the following files:

`....seqs.absolute.abund_and_meta.csv`: The count table detailing the absolute abundance of sequences returned from each sample after processing in the SymPortal QC pipeline and after [Minimum Entropy Decomposition](http://merenlab.org/software/med/).

In addition to sequence abundance, the `....seqs.absolute.abund_and_meta.csv` table contains columns detailing QC-related meta information. These columns are as follows:

- _symportal_datasetsample_uid_ - The UID of the DataSetSample object that relates to this sample in the SymPortal remote database
- _sample-id_ - The sample-id as referenced in the TARA Pacific sample provenance tables
- _raw_contigs_ - The number of reads after mothur's make.contigs command is run
- _post_qc_absolute_seqs_ - The absolute number of reads remaining after basic mothur-based quality control has been completed (not including size discrimination)
- _post_qc_unique_seqs_ - The unique (distinct) number of reads remaining after basic mothur-based quality control has been completed (not including size discrimination). I.e. the number of lines in the the mothur generated .names file.
- _post_taxa_id_absolute_symbiodiniaceae_seqs_ - The absolute number of reads remaining after quality control and taxonomic screening (for Symbiodiniaceae sequences) has been completed.
- _post_taxa_id_unique_symbiodiniaceae_seqs_ - The unique (distinct) number of reads remaining after quality control and taxonomic screening (for Symbiodiniaceae sequences) has been completed.
- _size_screening_violation_absolute_ - The absolute number of reads that were thrown out purely due to size violations (sequences were too long or too short).
- _size_screening_violation_unique_ - The unique (distinct) number of reads that were thrown out purely due to size violations (sequences were too long or too short).
- _post_taxa_id_absolute_non_symbiodiniaceae_seqs_ - The absolute number of reads thrown out during taxonomic screening due to being of a non-Symbiodiniaceae origin.
- _post_taxa_id_unique_non_symbiodiniaceae_seqs_ - The unique (distinct) number of reads thrown out during taxonomic screening due to being of a non-Symbiodiniaceae origin.
- _post_med_absolute_ - The absolute number of sequences deposited into the database having gone through all QC, taxonomic screening and MED analysis.
- _post_med_unique_ - The unique (distinct) number of sequences deposited into the database having gone through all QC, taxonomic screening and MED analysis.

Following the meta information columns, the sequence names and their abundances are listed. For sequences with alphanumeric names associated, these names are used (for example: C3, D1a, C3cc etc.). For thsoe sequences without alphanumeric names, an identifier made up of the UID of the related ReferenceSeqeunce object in the SymPortal database followed by the phylogenetic clade to which the sequence belongs is used in the form \<UID_clade\>.

`....seqs.fasta`: Fasta file containing the DNA sequences referenced in the `....seqs.absolute.abund_and_meta.csv` count table.

`....additional_info.txt`: A plain text file detailing meta information relating to the SymPortal submission and analysis.

### between_sample_distances

Data files associated with between-sample similarity distances.

Distance matrices were computed based on the ITS2 sequence assemblages returned from each sample (using post-MED sequence abundances).
The Symbiodiniaceae is divided into several major phylogenetic clades, most of which relate directly to genera within the Symbiodiniaceae (see [LaJeunesse et al. 2018](https://www.sciencedirect.com/science/article/pii/S0960982218309072)).
To prevent signal from within-clade differences being masked by between-clade signal, and given that sequence representatives from some clades cannot be aligned with sequences from others, similairites are calculated in a clade separated manner (see [Hume et al. 2019](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13004) for further details).

This directory is therefore subdivided by the major clades of the Symbiodinaiceae (where present).

For each of these clade directories, pairwise distances were computed using both a Bray-Curtis- and a UniFrac-based approach. For each approach, distance matrices were computed directly from post-MED sequnce abundances and after applying a square root transformation (to increase the weight of lower abundance intragenomic sequences). A such, distance files for each combination of these factors (i.e. Bray-Curtis, UniFrac, with square root transformation, and without) exist (extension `.dist`).

For each distance matrix a corresponding principal coordinate analysis was conducted. The resulting corrdinates are provided in an associated `.csv` file.

For caluculation of UniFrac distances, a Maximum Likelihood sequence tree was generated. This tree was generated using the program [IQ-TREE](http://www.iqtree.org/). The nucleotide evolutionary model used to build the tree and the final tree are provided in the `.iqtree` file.

The sequences used in the generation of the distance matrices are provided both in a multiple sequence alignment format (`...seqs.aligned.fasta`; aligned using [MAFFT](https://mafft.cbrc.jp/alignment/software/)) and a non-aligned format (`...seqs.unaligned.fasta`).

Bray-Curtis and UniFrac computations were implemented in python using [scipy's Bray-Curtis implementation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.braycurtis.html) and [scikit-bio's weighted UniFrac implementation](http://scikit-bio.org/docs/0.4.2/generated/skbio.diversity.beta_diversity.html).

***

The following subdirectories exist only for the `20200326_tara_corals` dataset:

### between_profile_distances

The distances in this directory represent similarity between the ITS2 profiles predicted for the `20200326_tara_corals` dataset.

The data files of this subdirectory are in exactly the same format, as the between_sample_distances and have been calculated in the same way, only related to abundances of defining intragenomic vairances for the profiles in question.

### its2_type_profiles

This directory contains the following files:

`....profiles.absolute.abund_and_meta.csv`: The count table detailing the summed absolute abundances of the profile in question's defining intragnomic variants for a given sample.

For each ITS2 type profile the following attributes are listed:

- _ITS2 type profile UID_ - The uid of the ITS2 type profile object in the SymPortal database that created the output
- _Clade_ - The related phylogenetic clade of the ITS2 type profile object
- _Majority ITS2 sequence_ - A list of the sequences used to define the ITS2 type profile that were identified as the most abundant sequence in any one of the samples that were found to contain the profile.
- _Associated species_ - A list of formal species descriptions that contain ITS2 information that is related to the ITS2 type profile in question. N.B. This does NOT mean that the ITS2 type profile IS representative of these species.
- _ITS2 type abundance local_ - The number of samples in which the ITS2 type profile was found within the set of samples output.
- _ITS2 type abundance DB_ - The number of samples in which the ITS2 type profile was found within the entire SymPortal database from which the analysis was conducted.
- _ITS2 type profile_ - The name of the ITS2 type profile.

In the first column, under the ITS2 type profile meta information is the UID of the DataSetSample object that relates to this sample in the SymPortal remote database.
In the second column is the sample-id as referenced in the TARA Pacific sample provenance tables.

At the bottom of the table is the SymPortal database UID for each of the ReferenceSequence objects associated to the defining intragenomic sequences used to define the profiles.
The average abundance of each of the sequences and the standard deviation is given.

The references relating to the associated Symbiodiniaceae species is listed in the `additional_info.txt` file as is the meta information relating to the analysis.

`....additional_info.txt`: A plain text file detailing meta information relating to the SymPortal submission and analysis.
