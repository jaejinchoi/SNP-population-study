# SAT_mix. (Still) under renewal.
SNPhylo + Admixture + Treemix (SAT)

## Purpose
The original tool, SNPhylo (https://github.com/thlee/SNPhylo), is designed to build a phylogenetic tree from genome-wide SNPs of soybean. This updated tool added more feature to study population admixture (Admixture) and migration (Treemix), e.g., in humans.

Most modification done in 2014, at Korean Bioinformation Center (KOBIC).  

## Requirements
Updated availability on March 4, 2023  
- R R http://www.r-project.org/index.html
- PYTHON python http://www.python.org/
- PERL perl http://www.perl.org/get.html
- CLUSTALO clustalo https://github.com/GSLBiotech/clustal-omega
- ADMIXTURE admixture http://dalexander.github.io/admixture/
- TREEMIX treemix https://bitbucket.org/nygcresearch/treemix/downloads/
- PLINK plink1 https://zzz.bwh.harvard.edu/plink/

## Setup
Run 'bash setup.sh' (shell script) to check requirements.  
See ![SAT_mix_manual.pdf](SAT_mix_pack/SAT_mix_manual.pdf) to understand overall file structure.  

## ![SAT_mix.sh](SAT_mix_pack/SAT_mix.sh)
### [Arguments]
* -h or ?  
    Show help / options  


### [Input]
Accept SNP data files of VCF, PED, HapMap, simple SNP or GDS.  

### [Output]
See the ![SAT_mix_manual.pdf](SAT_mix_pack/SAT_mix_manual.pdf) #3 'output'

## Arguments  
Below arguments cover all three analyses: SNPylo, Admixture, and Treemix.  

* -h, show_help()  

### 1. Define input format; the default value in parenthesis
Sub arguments available for different file formats.  
* -v [VCF_file]  
    -p [Maximum_PLCS (5)]  
    -c [Minimum_depth_of_coverage (5)]  
* -H [HapMap_file]  
    -p [Maximum_PNSS (0)]  
* -s [simple_SNP_file]  
    -p [Maximum_PNSS (0)]  
* -d [GDS_file]  
    -l [LD_threshold 90.5)]
    -m [MAF_threshold (0.5)]  
* -a [PED (AGCT) file]  

### 2. Basic parameters  
* -l [float]  
    LD: Linkage disequilibrium search length  
* -m [float]  
    MAF: Minor Allele Frequency  
* -M [float]  
    Missing rate (0)  
* -P [output path]  
    prefix_ouput_path (./output)  

### 3. Functional parameters
Parameters with '*' at the end are necessary to TURN-ON its function  

[ADMIXTURE]  
Analyze population size from;
* -k [int], minimum ancestry *  
* -K [int], maximum ancestry *  

[TREEMIX]  
* -t [tree_mix_group_index file] *  
    See group_index_file format explained in the manual  
* -r [str]  
    Specific an outgroup OTU to reroot TREEMIX tree  
* -R [int], maximum migration  

[SNPhylo]  
* -b [int], Number of bootstrap sampling *  
* -o [str], Specific an outgroup OTU to reroot SNPhylo tree  

### [Note]
Applying linkage disequilibrium (LD) acquires random seeding and may output different results per run.

