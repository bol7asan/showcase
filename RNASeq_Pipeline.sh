#!/bin/bash

#Author: Mohammed Alhusayan
#Email: mohammed.alhusayan@pennmedicine.upenn.edu

##Parse options from command line
threads=8
function usage {
    echo -e "\n\e[1;31mAvailable Options:\e[1;m"
    echo -e "___________________________ \n"
    echo -e "\e[92m-g       \e[96mDirectory for the reference genome indexed by STAR \e[31m(required)\e[1;m"
    echo -e "\e[92m-f       \e[96mDirectory containing samples of fastq files \e[31m(required)\e[1;m"
    echo -e "\e[92m-t       \e[96mNumber of threads \e[1;m"
    echo -e "\e[92m-o       \e[96mAnalysis directory, defaults to current \e[1;m \n"
}


while getopts "g:f:t:h:o:" opt; do

    case "$opt" in
        g) refgenome=$OPTARG ;; 
        
        f) fastqdir=$OPTARG ;;
        
        t) threads=$OPTARG ;;
        
        h) help=$OPTARG usage ;;

        o) outdir=$OPTARG ;; 

        \?) usage
        
    esac
done
shift $((OPTIND -1))

refgenome=$( echo $refgenome | tr [a-z] [A-Z] )
if [ -z "$fastqdir" ] && [ -z "$refgenome" ] ; then
    echo -e  "\e[35m-f is not supplied\e[1;m    -> \e[94mPlease input the directory containing the fastq files \e[1;m"
    echo -e  "\e[35m-g is not supplied\e[1;m    -> \e[94mPlease input the directory for the reference genome \e[1;m"
    usage
    exit 1
fi

if [ -z "$refgenome" ]; then
    echo -e  "\e[35m-g is not supplied\e[1;m    -> \e[94mPlease input the directory for the reference genome \e[1;m"
    usage
    exit 1
fi

if [[ "$refgenome" == "GRCH38" ]]; then
    refgenome=HG38
fi

if [[ "$refgenome" == "HG38" ]]; then
    gtffile=/media/asangani2/genome/HG38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf
    genomeDir=/media/asangani2/genome/${refgenome}/STAR
fi


if [[ "$refgenome" == "T2T" ]]; then
    # gtffile=/media/asangani2/genome/T2T-CHM13v2.0/CHM13.v2.0.gff3
    gtffile=/media/asangani2/genome/T2T-CHM13v2.0/CHM13_V2.gff3
    genomeDir=/media/asangani2/genome/T2T-CHM13v2.0/STAR
    effsize=3054815472
fi

if [ -z "$fastqdir" ]; then
    echo -e  "\e[35m-f is not supplied\e[1;m    -> \e[94mPlease input the directory containing the fastq files \e[1;m"
    usage
    exit 1
fi

if [ -z "$outdir" ]; then
    
    #Get the parent directory containing the subdirectory that contains the fastq files.
    dirname=$( dirname "$fastqdir" )
    outdir=${dirname}
    
    echo -e  "\e[35m-o is not supplied\e[1;m    -> \e[94mPre-Analysis will be in the parent directory of the fastq files: ${outdir}/pre-analysis \e[1;m"
    
fi
analysisdir=${outdir}/pre-analysis
source activate ngsmo
#Create log file for complete analysis
touch ${analysisdir}/analysis.log


for sample in "$fastqdir"/*/{*PREC*,*LNCAP*,*VCAP*,*DU145*}; do
    
    name=$(basename "$sample")
    sample_dir=${analysisdir}/$name
    fastqFiles=($(ls ${sample}/*gz))

    

    echo -e "$(date)\tRunning pipeline on ${name}." | tee -a "$analysisdir"/analysis.log
        
    #FASTQC##
    Check the quality of the fastqc files
    FASTQC=$sample_dir/FASTQC
    [ -d $FASTQC ] || mkdir -p $FASTQC
    echo -e "$(date)\tStarting FASTQC on ${name}." | tee -a "$analysisdir"/analysis.log
    fastqc "$sample"/*gz -o $FASTQC -t $threads 2>&1 | tee "$FASTQC"/FASTQC.log

    #########STAR Alignment##########
    STARDIR=$sample_dir/STAR
    [ -d $STARDIR ] || mkdir -p $STARDIR
    
    if [ ${#fastqFiles[@]} -eq 2 ]; then

        echo -e "$(date)\tAligning paired end reads." | tee -a "$analysisdir"/analysis.log
        STAR --runThreadN $threads --sjdbGTFfile $gtffile --genomeDir $genomeDir --readFilesIn ${sample}/*gz --readFilesCommand zcat \
        --clip3pAdapterSeq CTGTCTCTTATACACATCT CTGTCTCTTATACACATCT  --clip3pAdapterMMp 0.1 0.1 --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes All --outFileNamePrefix $STARDIR/   --quantMode GeneCounts 2>&1 | tee  "$STARDIR"/STAR.log
    else
        echo -e "$(date)\tAligning single end reads." | tee -a "$analysisdir"/analysis.log
        STAR --runThreadN $threads --sjdbGTFfile $gtffile --genomeDir $genomeDir --readFilesIn ${sample}/*gz --readFilesCommand zcat \
            --clip3pAdapterSeq "CTGTCTCTTATACACATCT" --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes All --outFileNamePrefix $STARDIR/   --quantMode GeneCounts 2>&1 | tee  "$STARDIR"/STAR.log
    fi
    


    #Filter out duplicates and low quality reads using samtools
    echo -e "$(date)\tFiltering out duplicates and low quality reads using samtools." | tee -a "$analysisdir"/analysis.log
    samtools view -F 1024 -q 30 -b  -o ${STARDIR}/Aligned.sortedByCoord_filtered.out.bam ${STARDIR}/Aligned.sortedByCoord.out.bam
    
    echo -e "$(date)\tIndexing Bam File." | tee -a "$analysisdir"/analysis.log
    samtools index ${STARDIR}/Aligned.sortedByCoord_filtered.out.bam 
    
    echo -e "$(date)\tGetting raw read counts using HTSEQ for later use with DeSEQ2." | tee -a "$analysisdir"/analysis.log
    htseq-count -s reverse -m union -i gene_id --additional-attr chr --additional-attr start --additional-attr end --additional-attr gene_name --additional-attr gene_biotype  ${STARDIR}/Aligned.sortedByCoord_filtered.out.bam $gtffile > ${STARDIR}/gene_count.txt


    #######BAM_COVERAGE############

    BIGWIG="$sample_dir"/BIGWIG
    [ -d $BIGWIG ] || mkdir -p $BIGWIG
    echo -e "$(date)\tCreating BigWig file" | tee -a "$sample_dir"/analysis.log
    bamCoverage --bam ${STARDIR}/Aligned.sortedByCoord_filtered.out.bam \
    -o "$BIGWIG"/${name}-CPM.bw  \
    -p $threads \
    --normalizeUsing CPM \
    --effectiveGenomeSize $effsize \
    --ignoreForNormalization chrX chrM
    echo -e "$(date)\tFinished creating BigWig file" | tee -a "$sample_dir"/analysis.log

    ##################################
done
conda deactivate