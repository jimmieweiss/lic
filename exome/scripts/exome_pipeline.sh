#! /bin/bash -l                                                                                                                                                                                                                
set -e

#SBATCH -A b2013064                                                                                                                                                                                       
#SBATCH -o outputs/Halopipeline_140113_01.out                                                                                                                                                            
#SBATCH -e outputs/Halopipeline_140113_01.err                                                                                                                                                            
#SBATCH -J outputs/Halopipeline_140113_01.job                                                                                                                                                            
#SBATCH -p core -n 1                                                                                                                                                                                      
#SBATCH -t 08:00:00                                                                                                                                                                            
#SBATCH --mail-user jimmie.hallman@scilifelab.se                                                                                                                                                          
#SBATCH --mail-type=ALL                                                                                                                                                                                   
                    

source ~/.virtualenvs/virtenv/bin/activate

module unload python/2.6.6
module load python/2.7
module load bioinfo-tools
module load FastQC/0.10.1
module load biopython/1.56
module load htseq/0.5.1
module load bwa/0.7.4
module load samtools/0.1.19

WORKINGDIR=***
INFILE1=$WORKINGDIR'/demultiplexed/lane2_index2_L002_R1_001.fastq'
OUTFILE1=$WORKINGDIR'/demultiplexed/lane2_index2_R1_trimmed.fastq'

INFILE2=$WORKINGDIR'/demultiplexed/lane2_index2_L002_R2_001.fastq'
OUTFILE2=$WORKINGDIR'/demultiplexed/lane2_index2_R2_trimmed.fastq'

TEMP1=$INFILE1'.temp'
TEMP2=$INFILE2'.temp'

ADAPTER_FWD='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
ADAPTER_REV='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'

PREALIGN1=$WORKINGDIR'/demultiplexed/lane2_index2_R1_prealign.fastq'
PREALIGN2=$WORKINGDIR'/demultiplexed/lane2_index2_R2_prealign.fastq'

REF=/bubo/nobackup/uppnex/reference/biodata/genomes/Hsapiens/hg19/bwa/hg19.fa

SAI1=$WORKINGDIR'/aligned/L1_R1_index1.sai'
SAI2=$WORKINGDIR'/aligned/L1_R2_index1.sai'

READGROUP='@RG\tID:Haloplex_index1_lane1\tSM:HaloPlex_index1'
SAMOUT=$WORKINGDIR'/aligned/L1_index1.sam'

SAMHEADER=$WORKINGDIR'/aligned/L1_index2_final.sam'
SAMSORTED=$WORKINGDIR'/aligned/L1_index2.sorted'

MERGEDBAM=$WORKINGDIR'/aligned/index2.sorted.bam'

LANELIST=[LANE1, LAN2, LANE3]

for i in LANELIST:
do
	#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC $INFILE1 -o $OUTFILE1
	#cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT $INFILE2 -o $OUTFILE2
	cutadapt -a $ADAPTER_FWD --minimum-length 50 --paired-output $TEMP2 -o $TEMP1 $INFILE1 $INFILE2
	cutadapt -a $ADAPTER_REV --minimum-length 50 --paired-output $OUTFILE1 -o $OUTFILE2 $TEMP2 $TEMP1
	
	#rm $TEMP1 $TEMP2
	python prealignment.py $OUTFILE1 $OUTFILE2 $PREALIGN1 $PREALIGN2 20

	fastqc $PREALIGN1
	fastqc $PREALIGN2
                                     
	bwa aln -t 8 $REF $PREALIGN1  > $SAI1
	bwa aln -t 8 $REF $PREALIGN2  > $SAI2
	bwa sampe -A -r $READGROUP $REF $SAI1 $SAI2 $PREALIGN1 $PREALIGN2 > $SAMOUT

	python change_sam.py $SAMOUT $SAMHEADER

	samtools view -bS $SAMHEADER | samtools sort - $SAMSORTED

	#rm $OUTFILE1 $OUTFILE2 $SAI1 $SAI2 $SAMOUT
end

samtools merge $MERGEDBAM [LIST OF SORTED BAMS]


rm ../../nobackup/jimmie/aligned/L1_index2_final.sam
rm ../../nobackup/jimmie/aligned/L1_index2.sorted.bam
rm ../../nobackup/jimmie/aligned/L2_index2_final.sam
rm ../../nobackup/jimmie/aligned/L2_index2.sorted.bam

samtools index ../../nobackup/jimmie/aligned/index2.sorted.bam

deactivate







