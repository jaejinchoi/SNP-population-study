#!/bin/bash

VERSION="2023-3; 2021-2; 2014-7"

# Declare functions
print_help_and_exit () {
    echo -e "Determine three population analyses based on SNP data with a VCF, a PED, a HapMap, a Simple SNP or a GDS file" 1>&2
    echo -e "From SNPhylo code additional modication was done for PAPGI study by JaeJin Choi, 2014\n"
    echo -e "Added features (require external programs)"
    echo -e "1. Admixture"
    echo -e "2. TreeMix"
    echo 1>&2
    echo -e "Version: ${VERSION}, customized by JaeJin Choi, KOBIC" 1>&2
    echo 1>&2
    echo -e "Usage:" 1>&2
    echo -e "$(basename $0)
    \n----------Input file type with related arguments
    [-v VCF_file | -p Maximum_PLCS (5) | -c Minimum_depth_of_coverage (5)]
    [-H HapMap_file | -p Maximum_PNSS (0)]
    [-s Simple_SNP_file |-p Maximum_PNSS (0)]
    [-d GDS_file | -l LD_threshold (0.5) | -m MAF_threshold (0.5)]
    [-a PED(ACGT)_file]
    [-T [int], a number of thread (for sequence alignment)]

    \n----------Basic parameters
    [-l [float], LD-linkage disequilibrium]
    [-m [float], MAF-minor allele frequency]
    [-M [float], Missing_rate(0)]
    [-P [absolute_path], prefix_output_folder (default = ./output)]
    [-T [int], threads when multithreading is available (default = 6)]
    
    \n----------Functional parameters
    Parameters with '*' at the end are necessary to TURN-ON its function
    
    1. Admixture [ Population K from -k to -K
    [-k (int > 1), From] *
    [-K (int > -k), To] *

    2. Treemix
    [-t (treemix_group_index), see group_index_file format] *
    [-r [str], Specific an outgroup OTU to reroot TREEMIX tree] 
    [-R [int], maximum migration]

    3. SNPhylo
    [-b (int >= 100), number of bootstrap_sample, normally use 100] *
    [-o [str], Specific an outgroup OTU to reroot SNPhylo tree] 
    
    \n----------ETC
    [-h, help]" 1>&2
    echo
    echo -e "Acronyms:" 1>&2
    echo -e "\tPLCS: The percent of Low Coverage Sample" 1>&2
    echo -e "\tPNSS: The percent of Sample which has no SNP information" 1>&2
    echo -e "\tLD: Linkage Disequilibrium" 1>&2
    echo -e "\tMAF: Minor Allele Frequency" 1>&2
    echo 1>&2
    echo -e "Simple SNP File Format:" 1>&2
    echo -e "\t#Chrom\tPos\tSampleID1\tSampleID2\tSampleID3\t..." 1>&2
    echo -e "\t1\t1000\tA\tA\tT\t..." 1>&2
    echo -e "\t1\t1002\tG\tC\tG\t..." 1>&2
    echo -e "\t..." 1>&2
    echo -e "\t2\t2000\tG\tC\tG\t..." 1>&2
    echo -e "\t2\t2002\tA\tA\tT\t..." 1>&2
    echo -e "\t..." 1>&2
    [ -n "$2" ] && echo -e "\n$2" 1>&2
    exit $1
}

## Initialize variables
# BASE_DIR=$(dirname $(readlink -f "$0")) # readlink in Mac OS X does not support the '-f' option
BASE_DIR=$(cd "$(dirname "$0")" && pwd)
SCRIPTS_DIR="${BASE_DIR}/scripts"

source "${BASE_DIR}/SAT_mix.cfg"

if [ ! -z "${R_LIBS_DIR}" ]
then
    if [ -z "${R_LIBS}" ]
    then
        export R_LIBS="${R_LIBS_DIR}"
    else
        export R_LIBS="${R_LIBS_DIR}:${R_LIBS}"
    fi
fi


# In order to use Pypy instead of Python
PYPY="$(which pypy 2> /dev/null | tail -n1 | tr -d '\t')"
[ ! -z "${PYPY}" ] && PYTHON="${PYPY}"

#supported raw input file
vcf_file=""
hap_file="" # HapMap file path
gds_file=""
smp_file="" # Simple file path
ped_file="" #ACGT ped file path, later convert to hapmap

#filtered file(type can be .hapmap. vcf or .gds)
#stage-1
mnn_filtered_file="" #output of remove_~ script

#fix all parameters, except LD
min_depth_of_coverage=5
max_plcs_pnss=0 #0 is not accepted, x/100

ld_threshold=0.5 #default ld range, 50
maf_threshold=0.5
missing_rate=0 #missing rate = 0(none), however, there are 'N'; gap bases in sequence
prefix_output_folder="${pwd}"
n_thread=6

admix_from=0
admix_to=0

treemix_index_path=""
treemix_root="San" #be cautious
treemix_migration=10

out_sample_id=""
num_bs_sample=0

admix_from=2

# Parse positional parameters
while getopts "v:H:d:a:c:p:l:k:K:t:r:R:m:M:P:o:s:b:hT:" OPT
do
    case "${OPT}" in
        'v')
            vcf_file="${OPTARG}"
            ;;
        'H')
            hap_file="${OPTARG}"
            ;;
        'd')
            gds_file="${OPTARG}"
            ;;
        's')
            smp_file="${OPTARG}"
            ;;
        'a')
            ped_file="${OPTARG}"
            ;;
        'c')
            min_depth_of_coverage="${OPTARG}"
            ;;
        'p')
            max_plcs_pnss="${OPTARG}"
            ;;
        'l')
            ld_threshold="${OPTARG}"
            ;;
        'm')
            maf_threshold="${OPTARG}"
            ;;
        'M')
            missing_rate="${OPTARG}"
            ;;
        'k')
            admix_from="${OPTARG}"
            ;;
        'K')
            admix_to="${OPTARG}"
            ;;
        't')
            treemix_index_path="${OPTARG}"
            ;;
        'r')
            treemix_root="${OPTARG}"
            ;;
        'R')
            treemix_migration="${OPTARG}"
            ;;
        'P')
            prefix_output_folder="${OPTARG}"
            ;;
        'o')
            out_sample_id="${OPTARG}"
            ;;
        'b')
            num_bs_sample="${OPTARG}"
            ;;
        'T')
            n_thread="${OPTARG}"
            ;;        
        'h')
            print_help_and_exit 0
            ;;
        '?')
            print_help_and_exit 1
            ;;
        *)
            print_help_and_exit 1
            ;;
    esac
done


#echo ${prefix_output}
#ret_output=$(cd $(dirname ${prefix_output}) && pwd)
#echo ${ret_output}

#BASE_DIR=$(cd "$(dirname "$0")" && pwd)
#exit 1

#Main processes
[ -z "${vcf_file}${hap_file}${gds_file}${smp_file}${ped_file}" ] && print_help_and_exit 1

#check if any function is TURN-ON


#define path
prefix_output="${prefix_output_folder}/output"

if [ ! -e ${prefix_output_folder} ]
then
    mkdir ${prefix_output_folder}
else
    rm -f ${prefix_output_folder}/*
fi


#single run
if [ ! -z "${vcf_file}" ]
then
    [ ! -e "${vcf_file}" ] && print_help_and_exit 1 "VCF file (${vcf_file}) was not found!"
    [ $(wc -l < "${vcf_file}") -lt 10000 ] && print_help_and_exit 1 "VCF file (${vcf_file}) is too small to run this script!"

    # Remove VCF data which have many low depth of coverage samples
    echo "Start to remove low quality data."
    mnn_filtered_file="${prefix_output}.filtered.vcf"

    #force feed whole, by copying
    echo "Force-feed ${vcf_file} to ${mnn_filtered_file}, ignore -c and -p"
    cp ${vcf_file} ${mnn_filtered_file}
   
    #"${PYTHON}" "${SCRIPTS_DIR}/remove_low_depth_genotype_data.py" "${vcf_file}" ${min_depth_of_coverage} ${max_plcs_pnss} > "${mnn_filtered_file}"
    [ $? != 0 ] && exit 1

    # Determine and show the number of removed SNP data
    echo -e "\n$[$(wc -l < "${vcf_file}") - $(wc -l < "${prefix_output}.filtered.vcf")] low quality lines were removed."
    echo

    # 
    [ $(wc -l < "${mnn_filtered_file}") -lt 10000 ] && print_help_and_exit 1 "Error: There are too small number of SNP data in a file (${mnn_filtered_file})!\nPlease restart this script with different parameter values (-p and/or -c)."

elif [ ! -z "${hap_file}" ]
then
    [ ! -e "${hap_file}" ] && print_help_and_exit 1 "HapMap file (${hap_file}) was not found!"
    [ $(wc -l < "${hap_file}") -lt 10000 ] && print_help_and_exit 1 "HapMap file (${hap_file}) is too small to run this script!"


    # Remove HapMap data which have many no genotype data
    echo "Start to remove low quality data."
    mnn_filtered_file="${prefix_output}.filtered.hapmap"

    "${PYTHON}" "${SCRIPTS_DIR}/remove_no_genotype_data.py" "${hap_file}" ${max_plcs_pnss} > "${mnn_filtered_file}"
    [ $? != 0 ] && exit 1

    # Determine and show the number of removed SNP data
    echo -e "\n$[$(wc -l < "${hap_file}") - $(wc -l < "${mnn_filtered_file}")] low quality lines were removed"
    echo

    # 
    [ $(wc -l < "${mnn_filtered_file}") -lt 10000 ] && print_help_and_exit 1 "Error: There are too small number of SNP data in a file (${mnn_filtered_file})!\nPlease restart this script with different parameter values (-p)."

elif [ ! -z "${smp_file}" ]
then
    [ ! -e "${smp_file}" ] && print_help_and_exit 1 "Simple SNP file (${smp_file}) was not found!"
    [ $(wc -l < "${smp_file}") -lt 10000 ] && print_help_and_exit 1 "Simple data file (${smp_file}) is too small to run this script!"


    # Convert Simple SNP data file to HapMap file
    echo "Start to convert the simple SNP file to a HapMap file."

    "${PYTHON}" "${SCRIPTS_DIR}/convert_simple_to_hapmap.py" "${smp_file}" > "${prefix_output}.hapmap"
    [ $? != 0 ] && print_help_and_exit 1 "Error: The simple SNP file could not be converted to a HapMap file.\nPlease check the file and restart this script."

    # Remove HapMap data which have insufficient genotype data
    echo "Start to remove low quality data."

    mnn_filtered_file="${prefix_output}.filtered.hapmap"
    "${PYTHON}" "${SCRIPTS_DIR}/remove_no_genotype_data.py" "${mnn_filtered_file}" ${max_plcs_pnss} > "${mnn_filtered_file}"
    [ $? != 0 ] && exit 1

    # Determine and show the number of removed SNP data
    echo -e "\n$[$(wc -l < "${prefix_output}.hapmap") - $(wc -l < "${mnn_filtered_file}")] low quality lines were removed"
    echo

    # 
    [ $(wc -l < "${mnn_filtered_file}") -lt 10000 ] && print_help_and_exit 1 "Error: There are too small number of SNP data in a file (${mnn_filtered_file})!\nPlease restart this script with different parameter values (-p)."

elif [ ! -z "${gds_file}" ]
then
    [ ! -e "${gds_file}" ] && print_help_and_exit 1 "GDS file (${gds_file}) was not found!"

elif [ ! -z "${ped_file}" ]
then
    [ ! -e "${ped_file}.ped" ] && print_help_and_exit 1 "PED(ACGT) file (${ped_file}.ped) was not found!"
    [ ! -e "${ped_file}.map" ] && print_help_and_exit 1 "MAP file (${ped_file}.map) was not found!"

    #convert ped to tped
    echo
    echo "Convert ${ped_file} -> ${ped_file}.hapmap"
    ${PLINK} --noweb --file "${ped_file}" --recode --transpose --out "${ped_file}" --missing-genotype 'N' > /dev/null #this output tped, tfam
    [ $? != 0 ] && print_help_and_exit 1 "Try input file name without file extension"

    #convert tped to hapmap
    ${PERL} ${SCRIPTS_DIR}/convert_tped_to_hapmap.pl --tped ${ped_file}.tped --tfam ${ped_file}.tfam --build=papgi
    [ $? != 0 ] && exit 1

   
    #change name
    mv ${ped_file}.tped.hapmap ${ped_file}.hapmap
    hap_file=${ped_file}.hapmap 
    sed -i "s/0__//g" ${hap_file}

    echo "Obtain; ${hap_file}"
    echo

    
    # Remove HapMap data which have many no genotype data
    echo "Start to remove low quality data."
    mnn_filtered_file="${prefix_output}.filtered.hapmap"
    "${PYTHON}" "${SCRIPTS_DIR}/remove_no_genotype_data.py" "${hap_file}" ${max_plcs_pnss} > "${mnn_filtered_file}"
    [ $? != 0 ] && exit 1

    # Determine and show the number of removed SNP data
    echo -e "\n$[$(wc -l < "${hap_file}") - $(wc -l < "${mnn_filtered_file}")] low quality lines were removed"
    echo

    # 
    [ $(wc -l < "${mnn_filtered_file}") -lt 10000 ] && print_help_and_exit 1 "Error: There are too small number of SNP data in a file (${mnn_filtered_file})!\nPlease restart this script with different parameter values (-p)."

else
    print_help_and_exit 1 "SNP data file was not found!"

fi



# Generate sequences from SNP data, and ${prefix_output}.picked.hapmap 
# there are three accepted file formats (vcf, hapmap, and gds), and other formats are converted to one of three before feeded
if [ ! -z "${vcf_file}" ]
then
    "${R}" --slave --vanilla --file="${SCRIPTS_DIR}/generate_snp_sequence.R" --args -v "${mnn_filtered_file}" -l "${ld_threshold}" -m "${maf_threshold}" -M "${missing_rate}" -o "${prefix_output}" > /dev/null
    [ $? != 0 ] && exit 1
   
elif [ ! -z "${hap_file}" ]
then
    "${R}" --slave --vanilla --file="${SCRIPTS_DIR}/generate_snp_sequence.R" --args -H "${mnn_filtered_file}" -l "${ld_threshold}" -m "${maf_threshold}" -M "${missing_rate}" -o "${prefix_output}" > /dev/null
    [ $? != 0 ] && exit 1

elif [ ! -z "${gds_file}" ]
then
    "${R}" --slave --vanilla --file="${SCRIPTS_DIR}/generate_snp_sequence.R" --args -d "${mnn_filtered_file}" -l "${ld_threshold}" -m "${maf_threshold}" -M "${missing_rate}" -o "${prefix_output}" > /dev/null
    [ $? != 0 ] && exit 1

fi

#"${R}" --slave --vanilla --file="${SCRIPTS_DIR}/generate_snp_sequence.R" --args -H "${mnn_filtered_file}" -l "${ld_threshold}" -m "${maf_threshold}" -M "${missing_rate}" -o "${prefix_output}" #hapmap
#"${R}" --slave --vanilla --file="${SCRIPTS_DIR}/generate_snp_sequence.R" --args -v "${mnn_filtered_file}" -l "${ld_threshold}" -m "${maf_threshold}" -M "${missing_rate}" -o "${prefix_output}" #VCF
#[ $? != 0 ] && exit 1


#output .fasta, and .picked.hapmap
picked_file="${prefix_output}.picked"


#print sequence length
#seq_len=$(sed -n "2p" ${prefix_output}.fasta | wc -c)
seq_len=$(head -n 2 "${prefix_output}.fasta" | tail -n 1 | tr -d "\n" | wc -c)
echo -e "!ret(ld:maf:miss)\t${ld_threshold}\t${maf_threshold}\t${missing_rate}\t${seq_len}" > "${prefix_output}.snp_num"
echo "Finally picked; ${seq_len} SNPs"



#SNP size restriction -pseudo function
if [ ${ignore_kill_code} -ne 1 ]; then

    #additional options for sequence length check
    if [ ${seq_len} -gt 20000 ] || [ ${seq_len} -lt 5000 ]; then  #[ ${ld_threshold} -ne 0 ]; then # && [ ${ld_threshold} -ne 0 ]; then

        rm -f -r "${prefix_output_folder}" #remove current dir
        exit 1 #stop further process

    fi
        
fi


#stage-2, throught series of work, faster the first
#1. admixture
#2. treemix
#3. snphylo

if [[ ${admix_from} -ge 1 && ${admix_to} -ge ${admix_from} ]]
then

    admixture_folder="${prefix_output_folder}/admixture"

    if [ ! -e ${admixture_folder} ]
    then
        mkdir ${admixture_folder}
    else
        rm -f ${admixture_folder}/*.*
    fi

    #move to working folder
    cd ${admixture_folder}
    admixture_prefix="${admixture_folder}/out"

    echo
    echo "--Admixture start"
    echo "Prepare Admixture..."
    #convert ACGT format -> 1,2 format by plink
    ${PLINK} --noweb --file "${picked_file}" --recode12 --out "${admixture_prefix}_12" > /dev/null #log file is not necessary
    [ $? != 0 ] && exit 1
    echo "Obtain; ${admixture_prefix}_12.ped(map), --recode12" 

    echo
    echo "Admixture analysis proceed..."
    #population structure analysis using "Admixture", 2014.5.16
    echo "cv_error_trend" > "${admixture_prefix}.cv_error"

    #for k in {${admix_from}..${admix_to}}
    for ((k=${admix_from}; k<=${admix_to}; k++)) #2 to 7
    do
        #"${SCRIPTS_DIR}/admixture" --cv "${prefix_output}_12" $k, admixture require 1,2 based ped
        ${ADMIXTURE} --cv "${admixture_prefix}_12.ped" ${k} > "${admixture_prefix}_12.${k}.log"; [ $? != 0 ] && exit 1

        echo -e "${k}\t$( grep "CV" "${admixture_prefix}_12.${k}.log")" >> "${admixture_prefix}.cv_error"
        echo "1- tree k=${k}" 

        rm -f "${admixture_prefix}_12.${k}.log"

        #draw barplot
        #"${R}" --slave --vanilla --file="${SCRIPTS_DIR}/barplot_multi.R" --args "${admixture_prefix}_12.${k}.Q"
        "${R}" --slave --vanilla --file="${SCRIPTS_DIR}/barplot_multi_index.R" --args -n "${admixture_prefix}_12.nosex" "${admixture_prefix}_12.${k}.Q" #sorted barplot
        echo "2- obtain figure ${admixture_prefix}_12.${k}.Q.png"
        [ $? != 0 ] && exit 1

    done


    cd .. #admixture folder out
    echo "--Admixture done"

else
    echo "--Admixture skip, -K and -k = 0"

fi



#run treemix if treemix_input_path is specificed, 2014.6.23
if [ ! -z "${treemix_index_path}" ]
then
    [ ! -e "${treemix_index_path}" ] && print_help_and_exit 1 "Treemix index file (${treemix_index_path}) was not found!"
    #[ ! -e "${vcf_file}" ] && print_help_and_exit 1 "treemix_index_file (${treemix_index_file}) was not found!"
    #convert hapmap to treemix input format

    echo
    echo "TreeMix analysis proceed..."

    #check if specified root exists in index file
    grep ${treemix_root} ${treemix_index_path}
    [ $? != 0 ] && print_help_and_exit 1 "Root ${treemix_root} was not found in (${treemix_index_path})"

    treemix_folder="${prefix_output_folder}/treemix"

    if [ ! -e ${treemix_folder} ]
    then
        mkdir ${treemix_folder}
    else
        rm -f ${treemix_folder}/*.*
    fi

    #move to working folder
    cd ${treemix_folder}
    treemix_prefix="${treemix_folder}/out"

    #convert pciked.ped to hapmap
    echo
    echo "--Treemix start"
    echo "Prepare TreeMix..."
    #echo "Convert ${prefix_output}.picked.ped(.map) -> ${prefix_output}.picked.hapmap"
    #1. ped -> tped

    echo "Convert ${picked_file}.ped(map) -> ${treemix_prefix}.hapmap"
    ${PLINK} --noweb --file "${picked_file}" --recode --transpose --out "${treemix_prefix}" > /dev/null #log file is not necessary
    [ $? != 0 ] && exit 1

    #2. tped -> hapmap
    ${PERL} ${SCRIPTS_DIR}/convert_tped_to_hapmap.pl --tped "${treemix_prefix}.tped" --tfam "${treemix_prefix}.tfam" --build=papgi
    [ $? != 0 ] && exit 1

    mv ${treemix_prefix}.tped.hapmap ${treemix_prefix}.hapmap

    #sed -i "s/0__//g" ${treemix_prefix}.hapmap #for just in case?

    echo "Obtain; ${treemix_prefix}.hapmap"

    #exit 1 #out

    "${PYTHON}" "${SCRIPTS_DIR}/treemix_input.py" -i "${treemix_index_path}" "${treemix_prefix}.hapmap" > "${treemix_prefix}.treemix_input"
    [ $? != 0 ] && exit 1
    echo "1- convert hapmap -> treemix input format"

    gzip ${treemix_prefix}.treemix_input 
    echo "2- gzip compress ${treemix_prefix}.treemix_input -> ${treemix_prefix}.treemix_input.gz"

    #run treemix(should be installed, root is 'San' as a default
    ${TREEMIX} -i "${treemix_prefix}.treemix_input.gz" -o "${treemix_prefix}" -m "${treemix_migration}" -root "${treemix_root}" > /dev/null #cast process print
    echo "3- run treemix, -m ${treemix_migration} -root ${treemix_root}"

    #visualization, require R lib, "RColorBrewer"
    "${R}" --slave --vanilla --file="${SCRIPTS_DIR}/plotting_funcs_image.R" --args "${treemix_prefix}" > /dev/null
    [ $? != 0 ] && exit 1
    echo "4- obtain figure ${treemix_prefix}.pnd"

    cd .. #treemix out
    echo "--Treemix done"

else
    echo "Proceed without TreeMix analysis(no -t specified or ${tree_index_path} not exists}"

fi


if [ ${num_bs_sample} -ge 100 ] #since Bootstrap minimum is 100
then

    snphylo_folder="${prefix_output_folder}/snphylo"

    if [ ! -e ${snphylo_folder} ]
    then
        mkdir ${snphylo_folder}
    else
        rm -f ${snphylo_folder}/*.*
    fi

    #move to working folder
    cd ${snphylo_folder}
    snphylo_prefix="${snphylo_folder}/out"

    echo
    echo "--SNPhylo start"
    echo "MSA proceed using ${seq_len} SNPs"

    #SNP sequence alignment using MUSCLE, takes much of time, should consider accuracy depend on argument usage
        
    #using clustalo (clustal omega) that supports multithreading for large datasets
    if [ ${seq_len} -ne 0 ]; then #clustal omega alignment, output phylip format --outfmt=phy
        echo "clustalo --threads=${n_thread} -i ${prefix_output}.fasta -o ${snphylo_prefix}.align.txt --outfmt=fasta"
        clustalo --threads=${n_thread} -i ${prefix_output}.fasta -o ${snphylo_prefix}.align.txt --outfmt=fasta; [ $? != 0 ] && exit 1

    else
        echo "Zero length to align using clustalo. Exit"
	exit 1

    fi
    #-maxiters, default=16
    #-diags, dynamic programming rely on k-mer clustering, faster process, but lower accuracy


    echo
    echo "Bootstrap tree draw proceed"

    ## Determine a phylogenetic tree(maximum likelihood)
    for rm_file in outfile outtree
    do
        [ -e "${rm_file}" ] && rm -f "${rm_file}"
    done

    #ln -sf ${snphylo_prefix}.phylip.txt infile #symbolic link

    if [ -z "${out_sample_id}" ]
    then
        #echo "raxmlHPC-PTHREADS-AVX -T 6 -s ./infile -n ./outfile -m GTRCAT -f a" #print what is using
        #raxmlHPC-PTHREADS-AVX -T 6 -s ./infile -n ./outfile -m GTRCAT -f a; [ $? != 0 ] && exit 1 #using raxMLHPC    
        echo "raxml-ng --threads ${n_thread} --msa ${snphylo_prefix}.align.txt  --msa-format FASTA --prefix ${snphylo_prefix} --model GTR+G --all --bs-trees ${num_bs_sample}"
        raxml-ng --threads ${n_thread} --msa ${snphylo_prefix}.align.txt --msa-format FASTA --prefix ${snphylo_prefix} --model GTR+G --all --bs-trees ${num_bs_sample} > /dev/null; [ $? != 0 ] && exit 1

        echo "# run-log saved in ${snphylo_prefix}.raxml.log"

    else
        out_sample_no=$[$(grep -ne "\<${out_sample_id}\>" infile | cut -f1 -d":") - 1]
        
            if [ ${out_sample_no} -eq -1 ]
        then
            echo "Error!!! There is no sample name (${out_sample_id}) to use as a outgroup in a tree input file (${snphylo_prefix}.align.txt)." 1>&2
            echo "Please check the name and restart this script." 1>&2
            rm -f infile
            exit 1
        else
            echo -e "o\n${out_sample_no}\ny\n" | "${DNAML}"; [ $? != 0 ] && exit 1
        fi
    fi

    #mv outfile "${snphylo_prefix}.ml.txt"
    #mv outtree "${snphylo_prefix}.ml.tree"
    #rm -f infile #remove symbolic link

    cd .. #snphylo out
    echo "--SNPhylo done"

fi

echo "!End without notable errors"
                                                                           
