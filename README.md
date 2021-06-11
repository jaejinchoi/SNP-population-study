# SAT_mix
SNPhylo + Admixture + Treemix.  

Original template script, SNPhylo: https://github.com/thlee/SNPhylo    
* Lee T-H, Guo H, Wang X, Kim C, Paterson AH (2014) SNPhylo: a pipeline to construct a phylogenetic tree from huge SNP data. BMC Genomics 15(1):162. 

The original script modified and expanded function by JaeJin Choi in 2014, at Korean Bioinformation Center (KOBIC).  

## Requirements  
- Perl, for setup and run: http://www.perl.org/get.html
- Python (version 2.7): https://www.python.org/download/releases/2.7/
- MUSCLE: http://www.drive5.com/muscle/
- DNAML dnaml http://evolution.genetics.washington.edu/phylip/
- ADMIXTURE: http://dalexander.github.io/admixture/
- TREEMIX: https://bitbucket.org/nygcresearch/treemix/wiki/Home
- PLINK: http://zzz.bwh.harvard.edu/plink/download.shtml
- Statistical R, for graphical presentation: http://www.r-project.org/index.html
  
## Setup
Run 'bash setup.sh' (shell script) to check and setup requirements.  
See ![SAT_mix_manual.pdf](SAT_mix_pack/SAT_mix_manual.pdf) to understand overall procedure.  

## ![SAT_mix.sh](SAT_mix_pack/SAT_mix.sh)
### [Arguments]
* -h or ?  
    Show help / options  

### 1. Define input format; the default value in parenthesis for each parameter
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
    LD: Linkage Disequilibrium  
* -m [float]  
    MAF: Minor Allele Frequency  
* -M [float]  
    Missing rate (0)  
* -P [output path]  
    prefix_ouput_path (./output)  

### 3. Functional parameters  
[ADMIXTURE]  
Analyze population size from -k [int] to -K [int]  

[TREEMIX]  
* -t [tree_mix_group_index file]  
    See group_index_file format explained in the manual  

[SNPhylo]  
* -b [int], Number of bootstrap sampling  
* -o [str], Outgroup_sample_name for a rooted tree  

### [Note]
Linkage disequilibrium function acquires random seed and may result various output per run.  

### [Input]
SNP data formats of VCF, PED HapMap, simple SNP or GDS.  

### [Output]
See the ![SAT_mix_manual.pdf](SAT_mix_pack/SAT_mix_manual.pdf) #3 'output'

## Limitation

