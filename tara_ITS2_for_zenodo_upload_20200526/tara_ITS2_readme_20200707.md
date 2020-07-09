# TARA Pacific ITS2 Symbiodiniaceae data release
## List of authors:
Hume, Benjamin C.C.; Poulain, Julie; Pesant, Stéphane; Belser, Caroline; 
Ruscheweyh, Hans-Joachim; Forcioli, Didier; Armstrong, Eric; Clayssen, Quentin; 
Henry, Nicolas; Klinges, Grace; McMinds, Ryan; Paoli, Lucas; Pogoreutz, Claudia; 
Salazar, Guillem; Ziegler, Maren; Moulin, Clementine; Boissin, Emilie; Bourdin, Guillaume; 
Iwankow, Guillaume; Romac, Sarah; Agostini, Sylvain; Banaigs, Bernard; Boss, Emmanuel; 
Bowler, Chris; de Vargas, Colomban; Douville, Eric; Flores, Michel; Furla, Paola; 
Galand, Pierre E.; Gilson, Eric; Lombard, Fabien; Reynaud, Stéphanie; Sullivan, Matthew B.; 
Thomas, Olivier; Troublé, Romain; Vega Thurber, Rebecca; Wincker, Patrick; Zoccola, Didier; 
Planes, Serge; Allemand, Denis; Sunagawa, Shinichi; Voolstra, Christian R.

## List of affiliations:
TBA TBA TBA

## Release 1, 2020707
Points of contact:

Benjamin C C Hume (didillysquat@gmail.com)

Christian R Voolstra (chris.voolstra@gmail.com)

## Table of contents

- [Introduction](#introduction)
- [Data subdirectories](#data-subdirectories)
  - Available for all datasets
    - [pre_med_seqs](#pre_med_seqs)
    - [post_med_seqs](#post_med_seqs)
  - Available for only coral and non_coral
    - [between_sample_distances](#between_sample_distances)
  - Available only for the coral dataset
    - [between_profile_distances](#between_profile_distances)
    - [its2_type_profiles](#its2_type_profiles)
- [Brief usage guide](#brief-usage-guide)
  - [The ITS2 marker in Symbiodiniaceae](#the-its2-marker-in-symbiodiniaceae)
  - [ITS2 type profiles](#its2-type-profiles-as-a-proxy-for-resolving-genotypes-within-the-symbiodiniaceae)
  - [Common use cases](#common-use-cases)
  - [Bespoke analyses using the pre-MED sequences](#bespoke-analyses-using-the-pre-med-sequences)
  

## Introduction

The TARA Pacific ITS2 data release is organised into three subdirectories representing different types of samples:

- `corals`: This directory contains results associated with samples taken directly from a coral animal
- `non_corals`: This directory contains results for samples not taken directly from a coral animal and includes water- and sediment-derived samples.
- `negatives`: This directory contains results from the negative control sequencing runs.

 For each sample type, datasets were created by submitting the associated sequencing 
 files to the remote instance of 
 the [SymPortal framework](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13004) ([symportal.org](www.symportal.org)).
 As such, for each sample type, the associated dataset has been processed through 
 the [SymPortal quality control (QC) pipeline](https://github.com/didillysquat/SymPortal_framework/wiki/The-SymPortal-logic#sequence-quality-control).

 In addition to sequence quality control, the `corals` dataset was included in a SymPortal analysis to predict ITS2 profiles.
 This dataset therefore contains some additional resources that are not available for the other two datasets.
 Prediction of ITS2 profiles is not possible for non-host-derived samples.
 For more details please refer to the [associated manuscript](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13004).
 
 Between sample distances are not provided for the negative samples,
 as only two samples, for one genus(Cladocopium) returned sufficient sequences for these metrics to be calculated
 (according to the cutoffs implemented within SymPortal).

 Where technical replicates were available for a given sample-id (i.e. multiple pairs of fastq.gz files),
 only the largest replicate (combined size of the two fastq.gz files) was used.

 ***

## Data subdirectories

The following subdirectories exist for each of the datasets:

### pre_med_seqs

This directory contains the following files:

`...sequences_absolute_abund_v1.csv`: The count table detailing the absolute abundance of sequences
returned from each sample after processing in the SymPortal QC pipeline,
but before undergoing [Minimum Entropy Decomposition](http://merenlab.org/software/med/)
that is conducted as standard when loading a dataset into the SymPortal database.

`...sequences_v1.fasta`: Fasta file containing the DNA sequences 
referenced in the `...sequences_absolute_abund_v1.csv` count table.

### post_med_seqs

This directory contains the following files:

`...sequences_absolute_abund_and_meta_v1.csv`: The count table detailing the absolute abundance of 
sequences returned from each sample after processing in the SymPortal QC pipeline 
and after [Minimum Entropy Decomposition](http://merenlab.org/software/med/).

In addition to sequence abundance, the `...sequences_absolute_abund_and_meta_v1.csv` table 
contains columns detailing QC-related meta information. These columns are as follows:

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

Following the meta information columns, the sequence names and their abundances are listed.
For sequences with alphanumeric names associated, these names are used 
(for example: C3, D1a, C3cc etc.). For thsoe sequences without alphanumeric names, 
an identifier made up of the UID of the related ReferenceSeqeunce object in the 
SymPortal database followed by the phylogenetic clade to which the sequence belongs is 
used in the form \<UID_clade\>.

`...sequences_v1.fasta`: Fasta file containing the DNA sequences referenced 
in the `...sequences_absolute_abund_and_meta_v1.csv` count table.

`...symportal_output_additional_info_v1.txt`: A plain text file detailing meta 
information relating to the SymPortal submission and analysis.

The between_sample_distances directory does not exist for the `negatives` dataset:

### between_sample_distances

Data files associated with between-sample dissimilarity distances.

Distance matrices were computed based on the ITS2 sequence assemblages returned from each 
sample (using post-MED sequence abundances).
The Symbiodiniaceae is divided into several major phylogenetic clades, 
most of which relate directly to genera within the Symbiodiniaceae 
(see [LaJeunesse et al. 2018](https://www.sciencedirect.com/science/article/pii/S0960982218309072)).
To prevent signal from within-clade differences being masked by between-clade signal, 
and given that sequence representatives from some clades cannot be aligned with sequences 
from others, similairites are calculated in a clade 
separated manner (see [Hume et al. 2019](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13004) 
for further details).

This directory is therefore subdivided by the major clades of the Symbiodinaiceae (where present).

For each of these clade directories, pairwise distances were computed using both a 
Bray-Curtis- and a UniFrac-based approach. For each approach, distance matrices were 
computed directly from post-MED sequence abundances and after applying a square root 
transformation (to increase the weight of lower abundance intragenomic sequences). 
A such, distance files for each combination of these factors (i.e. Bray-Curtis, UniFrac, 
with square root transformation, and without) exist (extension `.dist`).

For each distance matrix a corresponding principal coordinate analysis was conducted. 
The resulting corrdinates are provided in an associated `.csv` file.

For caluculation of UniFrac distances, a Maximum Likelihood sequence tree was generated. 
This tree was generated using the program [IQ-TREE](http://www.iqtree.org/). 
The nucleotide evolutionary model used to build the tree and the final tree are provided in 
the `.iqtree` file.

The sequences used in the generation of the distance matrices are provided both in a multiple 
sequence alignment format (`...sequences_aligned_v1.fasta`; aligned 
using [MAFFT](https://mafft.cbrc.jp/alignment/software/)) and a non-aligned 
format (`...sequences_unaligned_v1.fasta`).

Bray-Curtis and UniFrac computations were implemented in python using 
[scipy's Bray-Curtis implementation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.braycurtis.html) 
and [scikit-bio's weighted UniFrac implementation](http://scikit-bio.org/docs/0.4.2/generated/skbio.diversity.beta_diversity.html).

***

The following subdirectories exist only for the `corals` dataset:

### between_profile_distances

The distances in this directory represent dissimilarity between the ITS2 profiles 
predicted for the `corals` dataset.

The data files of this subdirectory are in exactly the same format, as the between_sample_distances 
and have been calculated in the same way, but related to abundances of defining 
intragenomic vairances for the profiles in question.

### its2_type_profiles

This directory contains the following files:

`...absolute_abund_and_meta_v1.csv`: The count table detailing the summed absolute abundances 
of the profile in question's defining intragnomic variants for a given sample.

For each ITS2 type profile the following attributes are listed:

- _ITS2 type profile UID_ - The uid of the ITS2 type profile object in the SymPortal 
database that created the output
- _Clade_ - The related phylogenetic clade of the ITS2 type profile object
- _Majority ITS2 sequence_ - A list of the sequences used to define the 
ITS2 type profile that were identified as the most abundant sequence in any one of the 
samples that were found to contain the profile.
- _Associated species_ - A list of formal species descriptions that contain ITS2 
information that is related to the ITS2 type profile in question. N.B. This does NOT mean 
that the ITS2 type profile IS representative of these species.
- _ITS2 type abundance local_ - The number of samples in which the ITS2 type profile was 
found within the set of samples output.
- _ITS2 type abundance DB_ - The number of samples in which the ITS2 type profile was found 
within the entire SymPortal database from which the analysis was conducted.
- _ITS2 type profile_ - The name of the ITS2 type profile.

In the first column, under the ITS2 type profile meta information is the UID of the 
DataSetSample object that relates to this sample in the SymPortal remote database.
In the second column is the sample-id as referenced in the TARA Pacific sample provenance tables.

At the bottom of the table is the SymPortal database UID for each of the 
ReferenceSequence objects associated to the defining intragenomic sequences used to define the profiles.
The average abundance of each of the sequences and the standard deviation is given.

The references relating to the associated Symbiodiniaceae species is listed in 
the `...symportal_output_additional_info_v1.txt` 
file as is the meta information relating to the analysis.

`...symportal_output_additional_info_v1.txt`: A plain text file detailing meta information 
relating to the SymPortal submission and analysis.

## Brief usage guide
### The ITS2 marker in Symbiodiniaceae
It is strongly recommended that users unfamiliar with the use of the ITS2 marker in Symbiodiniaceae take the
time to familiarise themselves with it. Certain characteristics of this marker complicate its interpretation.

Briefly, and most importantly: 
1. the rDNA copy number in Symbiodiniaceae is relatively high
(100-1000s of copies per genome).
2. sequence heterogeneity/diversity within the intragenomic copies is also high. 
3. maximum intragenomic distance, is often greater than average inter genomic distances.

As such, the ITS2 marker in Symbiodiniaceae should generally NOT be analysed/inferred from in the same way 
as other widely used barcode gene markers such as the 16S, 18S, COI etc.

An excellent primer
is _Thornhill et al. 2014: Host-specialist lineages dominate the adaptive radiation of reef coral endosymbionts_

### ITS2 type profiles as a proxy for resolving genotypes within the Symbiodiniaceae
Here, genotype is used to refer to genetically disparate populations of Symbiodiniaceae without
reference to a particular taxonomic level.
However, when working with the ITS2 marker, and making use of intragenomic diversity 
(i.e. predicting ITS2 type profiles),
genetic resolutions can generally be considered equivalent to species or sub-species resolutions.

For the TARA ITS2 Symbiodiniaceae dataset, the SymPortal framework has been used to predict ITS2 type profiles
that act as proxies for genotypes.
They can be thought of as the operational taxonomic unit of SymPortal.
In breif, ITS2 type profiles are characterised by sets of sequences 
(referred to as defining intragenomic sequence variants; DIVs) that have been found to co-occur in multiple samples.
Please refer to the [SymPortal manuscript](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13004)
([symportal.org](www.symportal.org)) for further background.

The remaining sections of this readme will deal with use case examples.

### Common use cases

Below are the two most common questions addressed using SymPortal analyses.

#### What is the predominant Symbiodiniaceae genotype in a given coral sample?
Users should refer to the ITS2 type profile count 
table `...absolute_abund_and_meta_v1.csv` in the `ITS2_type_profile` directory.

This table gives the absolute abundance of returned sequences on a per sample, per profile basis.

For example, if a value of 12345 is detailed for a given sample, e.g. TARA_XXX, and a given profile, e.g. A1-A1a-A1ab,
this means that the summed abundances of the A1, A1a and A1ab sequences in the TARA_XXX sample is 12345.

For the total number of seqeunces returned, at various stages of QC, for a specific sample, user should refer to 
the `...sequences_absolute_abund_and_meta_v1.csv` table in the `post_med_seqs`
directory.

Profile predictions have been made using the post-MED sequence
counts that can be found in the `..sequences_absolute_abund_and_meta_v1.csv` files in the `post_med_seqs` directories.

#### To what degree are a number of samples similar/dissimilar with respect to their Symbiodiniaceae community?
How the user goes about answering this question will depend on whether 
the samples are coral-derived or not.

If the samples are coral-derived, the user may check to see if the 
samples are predicted to contain the same or similar ITS2 type profiles
as well use the between sample distance matrices.

ITS2 type profile dissimilarities may be inferred from the distance files in the `between_profile_distances` directory.
Here the user will find both Bray-Curtis and Unifrac-derived between profile distance
matrices and related principal coordinate analysis (PCoA) results organised in a genus/clade separated manner.

If the samples in question are not coral-derived, ITS2 type profile predictions are not available.
Rather, the user must use the between sample dissimilarity matrices, and corresponding PCoA results, that are based on 
the sequence assemblages returned from each sample (housed in the `between_sample_distances` directories).

These distance matrices have been computed using the post-MED sequence
counts found in the `..sequences_absolute_abund_and_meta_v1.csv` files in the `post_med_seqs` directories.

The between sample dissimilarity computations take into account all sequences of a sample on a genus partitioned basis. 
These metrics therefore have, on a per sample basis, the advantage of taking all Symbiodiniaceae taxa
of a given genus/clade into account rather than only the more predomiant genotypes
(as is the case with predicting discrete ITS2 type profiles).
The disadvantage is that no discrete unit of comparison (i.e. ITS2 type profiles) is available.
 
Whether to work with the square root or non-square root transformation-computed distance matrices
will depend on the context of the questions being asked. Count matrices are square root transfored 
to afford lower abundance intragenomic sequence variants a resolving power
closer to that of higher abundance intragenomic sequence variants while maintaining sequence abundance information.
However, some of the lesser abundant sequences within samples may represent lower abundance 
taxa (rather than lower abundance intragenomic sequence variants from the most
abundant taxon). As such, this transformation could also afford a greater significance
to low-level taxa in determining the similarity/dissimilarity between samples.

### Bespoke analyses using the pre-MED sequences
Should the user wish to conduct their own analyses independent of the SymPortal outputs, 
the pre-MED sequence count tables provided in the `pre_med_seqs` directories are likely of interest.
The sequences have been through the same quality control as the post-MED sequences.
The only difference is that they hae not yet undergone MED.
Of note, these sequences have undergone taxonomic screening and are all Symbiodiniaceae in origin.
