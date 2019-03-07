## Scripts used to produce metagene profiles of ChIP-seq histone marks in Arabidopsis Thaliana.

## Retrieval of data
Will download datasets from SRA specified in a file SRRID.txt and unpack the SRA files to fastq.gz.

    download_data.sh

## Script that calculates the rate of that scription from pNET-seq data

    calc_nascent_expression.sh

## Script that can estimate new gene boundries based on TSS-seq, TIF-seq and PAS-seq data in bedgraph format

    estimate_representative_transcript.py

## Pipeline for mapping and processing of pNET-seq data.

    map_NETseq.sh

## Pipeline for mapping ChIP-seq histone Mark data (and RNA-seq with a few alterations)

    map_RNAseq_data.sh

## Script to estimate the most frequent 3'end adapter sequence in a dataset

    pick_adapter.sh

## Script to produce binnned metagene profiles

    mk_binned_metagene_ex_shortlist_50.sh

## Script to plot a metagene profile of a single dataset for different gneesets selected by rate of transcript

    print_profile_ex_50pct.R

## Script to plot many metagene profiles in one plot

    plot_histone_marks.R



CLIP Analysis Pipeline - CLAP, is used to analyse CLIP-seq (specifically PAR-CLIP, HITS-CLIP and iCLIP) data. The family of protocols where cross-linking immunoprecipitation coupled with high-throughput sequencing are referred to as CLIP. They are used identify binding sites of RNA-binding proteins and exist in the different flavours mentioned above. The data produced by these protocols have different characteristics which are accounted for in this pipeline, such that it can be used to process the above three versions and combination of them like iCLIP using 4SU nucleosides (PAR-iCLIP).

## 1. SYSTEM REQUIREMENTS
A relatively powerful computer is needed, absolute minimun 4 GB memory, the mapper will work much faster with more 6+ GB of memory. We currently support Linux and Mac OS X platforms. The pipeline is tested on different versions of Ubuntu, Debian and Mac OS X Sierra. It requires a BASH Shell and some version of awk installed.

## 2. INSTALLATION AND CONFIGURATION
1. Clone repository, find a suitable directory and use git to dowload the repository

        git clone https://github.com/simras/CLAP.git

2. Install Python (https://www.python.org) and Perl (https://www.perl.org), if not alrady installed.

3. Install bwa-pssm (http://bwa-pssm.binf.ku.dk/). bwa-pssm needs to be build from source and requires gdsl.

4. Install bedTools (http://bedtools.readthedocs.io/en/latest). The easy way is to use package managers, brew on Mac, apt-get on debian/Ubuntu or yum on other Linux distributions

5. Install pyicos (https://bitbucket.org/regulatorygenomicsupf/pyicoteo). Needs to be downloaded and installed.

6. Download mapping indexes and other files from (Share-links provided below for 6 species), copy to the folder CLAP, unpack and merge with folder CLAP/resources (it should happen automatically with wget and tar command described below.)

7. Set paths in scripts/CLAP.sh to the mapper bwa-pssm.

        # Absolute path to binary
        bwa=$BASE/../bwa-pssm/bwa

#### Share-links
<b><i>Homo sapiens</i>, original hg19 configuration.</b><BR>
https://sid.erda.dk/share_redirect/a4shAKTCcf

<b>Ensembl version 87</b><BR>
<b><i>Homo sapiens</i></b><BR>
https://sid.erda.dk/share_redirect/e0a4KihL8C
        
<b><i>Mus musculus</i></b><BR>
https://sid.erda.dk/share_redirect/BlzQzJW8jW

<b><i>Drosophila melanogaster</i></b><BR>
https://sid.erda.dk/share_redirect/DxzUSSuL4a

<b><i>Caenorhabditis elegans</i></b><BR>
https://sid.erda.dk/share_redirect/gEs04xFp1S

<b><i>Rattus norvegicus</i></b><BR>
https://sid.erda.dk/share_redirect/BB5p5bQ1lR

<b><i>Saccharomyces cerevisiae</i></b><BR>
https://sid.erda.dk/share_redirect/DThC7AWM7D

The annotation and sequence files are retrieved to the CLAP directory by

        wget <share-link>
        # unpack the compressed file
        tar -xvzf <tar.gz-file>

Ensure that your annotation and sequence file will be used (open scripts/CLAP.sh), by changing lines to fit your Ensembl version and species name.

        # Mapping index location
        # Ensembl version
        ver=87
        species=homo_sapiens

## 3. FURTHER CONFIGURATION
Currently the pipeline is set up with an ENSEMBL version 87 annotation and supports 6 species. If one wished to analyze data from a different species or use a different annotation it has to be integrated by following a number of steps. The scripts we provide parses an Ensembl annotation (GTF) file, if one wishes to use other annotation standards one will have to maunally make it conform to Ensembl formats.

### Updating annotation and mapping indexes:
The script "scripts/make_annotation.sh" contains commands to download and process annotation and sequence files from Ensembl version 87. As Ensembl alters their data formats slightly across versions these script may need to be updated, they have been tested on selected versions back to version 70. To create annotation for other species or other versions, configure the script by changing lines

        #Ensembl Version
        ver=87
        ...
        # Species name
        species="homo_sapiens"
        
and it will be necessary to find the species/assembly name in the Ensembl database, at their ftp site. It is briefly addressed in the following script. Note that some species only have partial assemblies and they will be hard to use in this constext as we assume alle annotation is name according to chromosome and coordinate, rather than scaffolds. You run the script like this
        
        scripts/make_annotation.sh

After creating new annotation it is neccessay to configure the pipeline and create new BWT-indexes for BWA-PSSM, as described above.
        
## 4. TEST-EXAMPLE
To test that everything works, run:

        scripts/testCLAP.sh

It maps, does peak calling and produces a UCSC custom track from reads that map to chr4 in the PAR-CLIP human dataset SRR248532. Since it is a human dataset, it is advice to test that the pipeline is setup correctly before setting it up with annotation for other species.

## 5. USAGE
To get a help menu run:
        
        scripts/CLAP.sh

14 options have to be specified in the sequence presented in the help menu.<BR>
ARGUMENTS:<BR>

    $1: Filename
    $2: Remove adapters?
        "": Don't remove adapters
        Otherwise: Adaptor sequence, ex: ACCTGCA...
    $3: Sequence fixed barcode
    $4: length of random barcode
    $5: Remove duplicates?
        0: No
        1: Yes
    $6: Type of analysis
        1: fixed barcode, random barcodes
        2: no fixed barcode, no random barcodes
        3: only fixed barcode, no random barcodes
        4: no fixed barcode, only random barcodes
    $7: UCSC Custom Tracks (bed tracks)
        0: No UCSC custom tracks
        1: UCSC custom tracks
    $8: Stranded protocol?
        0: Strandless
        1: Stranded
    $9: Index
        1: Genome index
        2: Genome index + exon junction index
    $10: Model
        0: Model T>C conversions (PAR-CLIP), conversion prob 0.125
        1: No model (RNA-Seq, iCLIP or HITS-CLIP)
    $11: Output name
    $12: Quality scores
        0: Phread 64
        1: Phread 33
    $13: Number of threads?
        Input number of threads
    $14: Peak calling?
        0: No
        1: Yes

Example analyses:

    PAR-CLIP (substitution model and no barcodes)
    scripts/CLAP.sh fastq-file TCGTATGCCGTCTTCTGCTTG "" 0 1 2 1 1 2 0 Analysis_name 1 8 1

    HITS-CLIP (no substitution model and no barcodes)
    scripts/CLAP.sh fastq-file TCGTATGCCGTCTTCTGCTTG "" 0 1 2 1 1 2 1 Analysis_name 1 8 1

    iCLIP (with multiplexing and duplication barcodes)
    scripts/CLAP.sh fastq-file TCGTATGCCGTCTTCTGCTTG GGTT 5 1 1 1 1 2 1 Analysis_name 1 8 1

The default substitution model has a T to C conversion rate at 12,5 %. A substitution model with different conversion probability can be created with the script scripts/mk_errorModel.py or the more general script where conversions from and to any nucleotide can be specified (See the repository of BWA-PSSM). <BR>

## 6. HOW TO CITE<BR>
M Plass, SH Rasmussen and A Krogh. Highly accessible AU-rich regions in 3′ untranslated regions are hotspots for binding of proteins and miRNAs. BioRxiv, 2016<BR>
doi: https://doi.org/10.1101/042986 <BR>

## 7. LICENSE<BR>
Copyright (c) 2017, Simon Horskjær Rasmussen and Mireya Plass. The software is open source and released under the MIT license.
