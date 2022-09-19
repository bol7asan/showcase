#!/bin/bash

#Author: Mohammed Alhusayan
#Email: mohammed.alhusayan@pennmedicine.upenn.edu



##Parse options from command line
function usage {
    echo -e "\n\e[1;31mAvailable Options:\e[1;m"
    echo -e "___________________________ \n"
    echo -e "\e[92m-g       \e[96mDirectory for the reference genome indexed by bowtie \e[31m(required)\e[1;m"
    echo -e "\e[92m-f       \e[96mDirectory containing samples of fastq files \e[31m(required)\e[1;m"
    echo -e "\e[92m-t       \e[96mNumber of threads"
    echo -e "\e[92m-n       \e[96mGenome name"
    echo -e "\e[92m-h       \e[96mPrint help message \e[1;m"
    echo -e "\e[92m-T       \e[96mWhether to trim the fastq files or not \e[1;m"
    echo -e "\e[92m-o       \e[96mAnalysis directory, defaults to current \e[1;m"
    echo -e "\e[92m-q       \e[96mEnter yes to run quality check only, default is no \e[1;m \n"
}

##Set defaults
outdir=.
threads=8
trim="no"
qc="no"
while getopts "g:f:t:n:h:T:o:q:" opt; do

    case "$opt" in
        g) refgenome=$OPTARG ;;

        f) fastqdir=$OPTARG ;;

        t) threads=$OPTARG ;;

        n) genome_name=$OPTARG ;;

        h) help=$OPTARG usage ;;

        T) trim=$OPTARG ;;

        o) outdir=$OPTARG ;;

        q) qc=$OPTARG ;;



        \?) usage

    esac
done
shift $((OPTIND -1))

refgenome=$( echo ${refgenome^^})
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

if [ -z "$fastqdir" ]; then
    echo -e  "\e[35m-f is not supplied\e[1;m    -> \e[94mPlease input the directory containing the fastq files \e[1;m"
    usage
    exit 1
fi

if [[ "$refgenome" == "GRCH38" ]]; then
    refgenome=HG38
fi

if [[ "$refgenome" == "HG38" ]]; then
    gtffile=/media/asangani2/genome/HG38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf
    lncRnaAnnot=/media/asangani2/genome/HG38/gencode.v38.long_noncoding_RNAs.gtf
    genomeDir=/media/asangani2/genome/${refgenome}/STAR
    bowtiedir=/media/asangani2/genome/HG38/BOWTIE/bowtie
    blacklist=/media/asangani2/genome/HG38/hg38-blacklist.v2.bed
    effsize=2747877777
fi

if [[ "$refgenome" == "HG19" ]]; then
    gtffile=/m
    genomeDir=/media/asangani2/genome/${refgenome}/STAR
    bowtiedir=/media/asangani2/genome/HG19/bowtie/hg19
    blacklist=/media/asangani2/genome/HG38/hg38-blacklist.v2.bed
    effsize=2736124973
fi

if [[ "$refgenome" == "MM39" ]]; then
    gtffile=/media/asangani2/genome/GRCm39/gencode.vM28.primary_assembly.annotation.gff3
    lncRnaAnnot=/media/asangani2/genome/HG38/gencode.v38.long_noncoding_RNAs.gtf
    genomeDir=/media/asangani2/genome/GRCm39/STAR
fi


if [[ "$refgenome" == "T2T" ]]; then
    gtffile=/media/asangani2/genome/T2T-CHM13v2.0/CHM13_V2.gff3
    bowtiedir=/media/asangani2/genome/T2T-CHM13v2.0/BOWTIE2/bowtie2
    effsize=3054815472
fi

## Category of marks based on ENCODE
broadMarks=( H3F3A H3K27ME3 H3K36ME3 H3K4ME1 H3K79ME2 H3K9ME1 H3K9ME2 H4K20ME1 MED1 )
narrowMarks=( H2AFZ H3AC H3K27AC H3K4ME2 H3K4ME3 H3K9AC AR )


echo $refgenome
# Loop over directory and get fastq files
for dir in "$fastqdir"/*; do
    
    sampleName=$( basename ${dir} )
    name=$( echo $sampleName | cut -f1,2 -d'_' )
    sampleDir="$outdir/pre-analysis/$name"

    # if [[ "$name" != *"Input"* ]]; then
    # echo -e "$(date)\tAnalyzing $name"
    # Fastq trimming and quality check##
    # if [[ $trim == "yes" ]]; then
        
    #     ######TRIM_GALORE#########
    #     fastqFiles=($(ls "$dir"/{*fastq.gz,*_L*_*.fq.gz} ))
    #     source activate trim_galore
    #     echo -e "$(date)\tTrimming fastq files for adapter and low quality reads" | tee "$sampleDir"/analysis.log
    #     trim_galore --paired --illumina -j "$threads" --basename "$name" -o "$dir"  ${fastqFiles[@]}
    #     echo -e "$(date)\tFinished trimming and quality check" | tee -a "$sampleDir"/analysis.log
    #     conda deactivate

    #     ###########################


    #     ##Capture trimmed files
    #     fastqFiles=($(ls "$dir"/*val*gz ))

    #     ##########FASTQC###########
        
    #     source activate ngsmo
    #     FASTQC="$sampleDir"/FASTQC
    #     [ -d "$FASTQC" ] || { mkdir -p "$FASTQC" && echo -e "$(date)\tRunning FASTQC."; }  | tee -a "$sampleDir"/analysis.log
    #     fastqc ${fastqFiles[@]} -o $FASTQC -t 2
    #     echo -e "$(date)\tFinished FASTQC." >> "$sampleDir"/analysis.log

    #     ##Continue to next iteration if only quality check is needed.
    #     if [ "$qc" == "yes" ]; then
    #         continue
    #     fi
    #     ############################

    #     ##########BOWTIE2###########

    #     BOWTIE2="$sampleDir"/BOWTIE2
    #     [ -d "$BOWTIE2" ] || { mkdir -p "$BOWTIE2" && echo -e "$(date)\tCreating BOWTIE2 directory."; }  | tee -a "$sampleDir"/analysis.log
    #     echo -e "$(date)\tStarting alignment to primary genome using bowtie2 using local alignment" | tee -a "$sampleDir"/analysis.log
    #     # Align to primary genome
    #     bowtie2 -1 ${fastqFiles[0]} -2 ${fastqFiles[1]} --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 \
    #     -p "$threads" -x "$bowtiedir" -S "$BOWTIE2"/${name}.sam &> "$BOWTIE2"/${name}-alignment_stats.txt
    #     echo -e "$(date)\tFinished alignment to genome" | tee -a "$sampleDir"/analysis.log

    #     #############################

    # else
    #     echo -e "$(date)\tNo trimming is done" | tee  "$sampleDir"/analysis.log

    #     ########FASTQC#############

    #     source activate ngsmo
    #     fastqFiles=($(ls "$dir"/{*fastq.gz,*_L*_*.fq.gz} ))
    #     FASTQC="$sampleDir"/FASTQC
    #     [ -d "$FASTQC" ] || { mkdir -p "$FASTQC" && echo -e "$(date)\tRunning FASTQC."; }  | tee -a "$sampleDir"/analysis.log
    #     fastqc ${fastqFiles[@]} -o $FASTQC -t 2
    #     echo -e "$(date)\tFinished FASTQC." >> "$sampleDir"/analysis.log

    #     ##Continue to next iteration if only quality check is needed.
    #     if [ "$qc" == "yes" ]; then
    #         continue
    #     fi
    #     ###########################

    #     ###########BOWTIE2##############

        BOWTIE2="$sampleDir"/BOWTIE2
    #     [ -d "$BOWTIE2" ] || { mkdir -p "$BOWTIE2" && echo -e "$(date)\tCreating BOWTIE2 directory."; }  | tee -a "$sampleDir"/analysis.log
    #     echo -e "$(date)\tStarting alignment to primary genome using bowtie2 using end to end alignment" | tee -a "$sampleDir"/analysis.log
    #     # Align to primary genome
        
    #     bowtie2 -1 ${fastqFiles[0]} -2 ${fastqFiles[1]} --local --very-sensitive --no-mixed --no-discordant --phred33 \
    #         -p "$threads" -x "$bowtiedir" -S "$BOWTIE2"/${name}.sam &> "$BOWTIE2"/${name}-alignment_stats.txt
    #     echo -e "$(date)\tFinished alignment to genome" | tee -a "$sampleDir"/analysis.log

    #     ##################################

    #     fi

    # ###Sort bam file by name
    # echo -e "$(date)\tSorting SAM file by name (for GATK) and converting to BAM using samtools." >> "$sampleDir"/analysis.log
    # source activate ngsmo
    # samtools sort -n -@ $threads -o "$BOWTIE2"/${name}_unfiltered_sortedByName.bam -O BAM "$BOWTIE2"/${name}.sam
    # echo -e "$(date)\tFinished sorting." >> "$sampleDir"/analysis.log


    # ##Mark duplicates
    # echo -e "$(date)\tMarking duplicates in BAM file." >> "$sampleDir"/analysis.log
    # gatk MarkDuplicates --TAGGING_POLICY All -ASO queryname -I "$BOWTIE2"/${name}_unfiltered_sortedByName.bam -O  "$BOWTIE2"/${name}_unfiltered_marked.bam -M "$BOWTIE2"/${name}_markDup.txt
    # echo -e "$(date)\tFinished marking duplicates." >> "$sampleDir"/analysis.log

    # ##Remove low quality reads and duplicates and extract only primary reads based on samflags.
    # echo -e "$(date)\tFiltering out duplicates" | tee -a "$sampleDir"/analysis.log
    # samtools view -@ $threads -b -q 30 -F 3340 -o "$BOWTIE2"/${name}_rmdup.bam "$BOWTIE2"/${name}_unfiltered_marked.bam
    # echo -e "$(date)\tFinished extracting primary reads." >> "$sampleDir"/analysis.log


    # ###Write index and sort by coordinate
    # echo -e "$(date)\tWrite index and sort by coordinates" | tee -a "$sampleDir"/analysis.log
    # samtools sort -@ $threads --write-index -o "$BOWTIE2"/${name}_primarySortedByCoord.bam -O BAM "$BOWTIE2"/${name}_rmdup.bam
    # echo -e "$(date)\tFinished sorting." >> "$sampleDir"/analysis.log
    # gatk MarkDuplicates --TAGGING_POLICY All -ASO queryname -I "$BOWTIE2"/${name}_rmdup.bam -O  "$BOWTIE2"/${name}_filtered_markDup.txt.bam -M "$BOWTIE2"/${name}_filtered_markDup.txt

    
    # ########BAM_COVERAGE############
    source activate ngsmo
    # samtools sort -@ $threads --write-index -o "$BOWTIE2"/${name}_unfiltered_sortedByCoord.bam -O BAM "$BOWTIE2"/${name}_unfiltered_sortedByName.bam
    BIGWIG="$sampleDir"/BIGWIG
    [ -d $BIGWIG ] || mkdir -p $BIGWIG
    echo -e "$(date)\tCreating BigWig file" | tee -a "$sampleDir"/analysis.log
    bamCoverage --bam "$BOWTIE2"/${name}_primarySortedByCoord.bam \
    -o "$BIGWIG"/${name}-BPM.bw  \
    -p $threads \
    --normalizeUsing BPM \
    --effectiveGenomeSize $effsize \
    --ignoreForNormalization chrX chrM -bs 10 --centerReads -e
    echo -e "$(date)\tFinished creating BigWig file" | tee -a "$sampleDir"/analysis.log

    ##################################
    conda deactivate
    ############MACS3#################
    # source activate macs3
    # shopt -s nocasematch
    
    # if  ! [[ "$name" =~ .*"igg".* ]] && ! [[ "$name" =~ .*"input".* ]] ; then
    #     #Set default peak calling mode
        
    #     peak=Both
    #     shopt -s nocaseglob
        
    #     ctrlFiles=$( ls ${outdir}/pre-analysis/{*input*,*igg*}/*BOWTIE2*/*primarySortedByCoord.bam )
    #     # ctrlFiles=/media/asangani1/Reyaz/LCMT1_Project/chipseq/2022-05-05/hg38/LNCAP-LACZ-B56_1/BOWTIE2/LNCAP-LACZ-B56_1_primarySortedByCoord.bam
    #     for mark in "${broadMarks[@]}"; do
    #         if [[ "$name" =~ .*"$mark".* ]];then
                
    #             peak=Broad
    #             echo -e "$(date)\tCalling peaks in broad mode" | tee -a "$sampleDir"/analysis.log 

    #             ##Call peaks with threhold 0.05
    #             MACS3="$sampleDir"/MACS3/broad/0.1
    #             [ -d $MACS3 ] || mkdir -p $MACS3
    #             macs3 callpeak --broad -g hs --broad-cutoff 0.1 -f BAMPE -t "$BOWTIE2"/${name}_primarySortedByCoord.bam --outdir "$MACS3"  -c $ctrlFiles
    #             ##Remove unwanted chromosomes that my skew the results.
    #             sed -i -r '/(chrX|chrM|chrUn)/d' ${MACS3}/NA_peaks.broadPeak
    #             ##Remove blacklisted regions
    #             bedtools intersect -a ${MACS3}/NA_peaks.broadPeak -b $blacklist -v > ${MACS3}/NA_peaks.broadPeakBlFiltered

                
    #             #Call peaks with threshold 0.1
    #             MACS3="$sampleDir"/MACS3/broad/0.05
    #             [ -d $MACS3 ] || mkdir -p $MACS3
    #             macs3 callpeak --broad -g hs --broad-cutoff 0.05 -f BAMPE -t "$BOWTIE2"/${name}_primarySortedByCoord.bam --outdir "$MACS3"  -c $ctrlFiles
    #             sed -i -r '/(chrX|chrM|chrUn)/d' ${MACS3}/NA_peaks.broadPeak
    #             bedtools intersect -a ${MACS3}/NA_peaks.broadPeak -b $blacklist -v > ${MACS3}/NA_peaks.broadPeakBlFiltered

    #             echo -e "$(date)\tFinished peak calling" | tee -a "$sampleDir"/analysis.log     
    #         fi
            
    #     done

        
    #     for mark in "${narrowMarks[@]}"; do

    #         if [[ "$name" =~ .*"$mark".* ]]; then

    #             peak=Narrow
    #             echo -e "$(date)\tCalling peaks in narrow mode" | tee -a "$sampleDir"/analysis.log 
    #             MACS3="$sampleDir"/MACS3/narrow/0.01
    #             [ -d $MACS3 ] || mkdir -p $MACS3
    #             macs3 callpeak -q 0.01 -f  BAMPE -g hs -t "$BOWTIE2"/${name}_primarySortedByCoord.bam --outdir "$MACS3"  -c $ctrlFiles
    #             sed -i -r '/(chrX|chrM|chrUn)/d' ${MACS3}/NA_peaks.narrowPeak
    #             bedtools intersect -a ${MACS3}/NA_peaks.narrowPeak -b $blacklist -v > ${MACS3}/NA_peaks.narrowPeakBlFiltered

    #             MACS3="$sampleDir"/MACS3/narrow/0.05
    #             [ -d $MACS3 ] || mkdir -p $MACS3
    #             macs3 callpeak -q 0.05 -g hs -f BAMPE -t "$BOWTIE2"/${name}_primarySortedByCoord.bam --outdir "$MACS3"  -c $ctrlFiles
    #             sed -i -r '/(chrX|chrM|chrUn)/d' ${MACS3}/NA_peaks.narrowPeak
    #             bedtools intersect -a ${MACS3}/NA_peaks.narrowPeak -b $blacklist -v > ${MACS3}/NA_peaks.narrowPeakBlFiltered
    #             echo -e "$(date)\tFinished peak calling" | tee -a "$sampleDir"/analysis.log 

    #         fi
    #     done

    #     if [[ "$peak" != "Narrow" && "$peak" != "Broad" ]]; then
    #         echo -e "$(date)\tRunning MACS3 in broad and narrow mode since mark is uncategorized." | tee -a "$sampleDir"/analysis.log  

    #         MACS3="$sampleDir"/MACS3/broad/0.05
    #         [ -d $MACS3 ] || mkdir -p $MACS3
    #         macs3 callpeak --broad -g hs --broad-cutoff 0.05 -f BAMPE -t "$BOWTIE2"/${name}_primarySortedByCoord.bam --outdir "$MACS3"  -c $ctrlFiles
    #         ##Remove unwanted chromosomes that my skew the results.
    #         sed -i -r '/(chrX|chrM|chrUn)/d' ${MACS3}/NA_peaks.broadPeak
    #         ##Remove blacklisted regions
    #         bedtools intersect -a ${MACS3}/NA_peaks.broadPeak -b $blacklist -v > ${MACS3}/NA_peaks.broadPeakBlFiltered

    #         ##################

    #         MACS3="$sampleDir"/MACS3/broad/0.1
    #         [ -d $MACS3 ] || mkdir -p $MACS3
    #         macs3 callpeak --broad -g hs --broad-cutoff 0.1 -f BAMPE -t "$BOWTIE2"/${name}_primarySortedByCoord.bam --outdir "$MACS3"  -c $ctrlFiles
    #         sed -i -r '/(chrX|chrM|chrUn)/d' ${MACS3}/NA_peaks.broadPeak
    #         bedtools intersect -a ${MACS3}/NA_peaks.broadPeak -b $blacklist -v > ${MACS3}/NA_peaks.broadPeakBlFiltered

    #         #####################
    #         MACS3="$sampleDir"/MACS3/narrow/0.01
    #         [ -d $MACS3 ] || mkdir -p $MACS3
    #         macs3 callpeak -q 0.01 -f  BAMPE -g hs -t "$BOWTIE2"/${name}_primarySortedByCoord.bam --outdir "$MACS3"  -c $ctrlFiles
    #         sed -i -r '/(chrX|chrM|chrUn)/d' ${MACS3}/NA_peaks.narrowPeak
    #         bedtools intersect -a ${MACS3}/NA_peaks.narrowPeak -b $blacklist -v > ${MACS3}/NA_peaks.narrowPeakBlFiltered

    #         ######################

    #         MACS3="$sampleDir"/MACS3/narrow/0.05
    #         [ -d $MACS3 ] || mkdir -p $MACS3
    #         macs3 callpeak -q 0.05 -g hs -f BAMPE -t "$BOWTIE2"/${name}_primarySortedByCoord.bam --outdir "$MACS3"  -c $ctrlFiles
    #         sed -i -r '/(chrX|chrM|chrUn)/d' ${MACS3}/NA_peaks.narrowPeak
    #         bedtools intersect -a ${MACS3}/NA_peaks.narrowPeak -b $blacklist -v > ${MACS3}/NA_peaks.narrowPeakBlFiltered

    #     fi
    # else
    # echo -e "$(date)\tNot calling peaks on igG." | tee -a "$sampleDir"/analysis.log
    # fi
            
    # conda deactivate
done


# source activate ngsmo
# multiqc "$outdir/pre-analysis"   -o "$outdir/pre-analysis"
# conda deactivate

# if [ "$qc" == "yes" ]; then
#     exit 1
# fi


# set -o noglob

# #Assemble peak stats
# bedFiles="$outdir/pre-analysis/*/MACS3/*/*/*Peak"
# xlsFiles="$outdir/pre-analysis/*/MACS3/*/*/*xls"
# Rscript /media/asangani2/scripts/downstream-analysis/chipseq/assemblePeakStats.R -b $bedFiles -x $xlsFiles -o $outdir/pre-analysis


# #Assemble alignment stats
# alignment_files="$outdir/pre-analysis/"*"/BOWTIE2/"*"alignment_stats.txt"
# dup_files=$outdir/pre-analysis/"*"/BOWTIE2/"*"_"?"_markDup.txt
# rmdup_files=$outdir/pre-analysis/"*"/BOWTIE2/"*"filtered_markDup.txt
# Rscript /media/asangani2/scripts/downstream-analysis/misc/bowtie_alignment_stats.R -a $alignment_files -d $dup_files -o $outdir/pre-analysis -D $rmdup_files


