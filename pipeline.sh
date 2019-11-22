#mkdir creates directories
#we'll make one called RNAseq_analysis
mkdir RNAseq_analysis
#cd is used to move into directories
cd RNAseq_analysis

#Rename no ablation + TBI group
mv V1_R1_001.fastq.gz noAb_TBI_rep1_PE1.fastq.gz
mv V1_R2_001.fastq.gz noAb_TBI_rep1_PE2.fastq.gz
mv V2_R1_001.fastq.gz noAb_TBI_rep2_PE1.fastq.gz
mv V2_R2_001.fastq.gz noAb_TBI_rep2_PE2.fastq.gz
mv V4_R1_001.fastq.gz noAb_TBI_rep3_PE1.fastq.gz
mv V4_R2_001.fastq.gz noAb_TBI_rep3_PE2.fastq.gz
mv V6_R1_001.fastq.gz noAb_TBI_rep4_PE1.fastq.gz
mv V6_R2_001.fastq.gz noAb_TBI_rep4_PE2.fastq.gz

#Rename ablation + TBI group
mv V-A1_R1_001.fastq.gz Ab_TBI_rep1_PE1.fastq.gz
mv V-A1_R2_001.fastq.gz Ab_TBI_rep1_PE2.fastq.gz
mv V-A2_R1_001.fastq.gz Ab_TBI_rep2_PE1.fastq.gz
mv V-A2_R2_001.fastq.gz Ab_TBI_rep2_PE2.fastq.gz
mv V-A5_R1_001.fastq.gz Ab_TBI_rep3_PE1.fastq.gz
mv V-A5_R2_001.fastq.gz Ab_TBI_rep3_PE2.fastq.gz
mv V-A6_R1_001.fastq.gz Ab_TBI_rep4_PE1.fastq.gz
mv V-A6_R2_001.fastq.gz Ab_TBI_rep4_PE2.fastq.gz

#Rename no ablation + sham group
mv V2-noTBI_R1_001.fastq.gz noAb_Sham_rep1_PE1.fastq.gz
mv V2-noTBI_R2_001.fastq.gz noAb_Sham_rep1_PE2.fastq.gz
mv V4-noTBI_R1_001.fastq.gz noAb_Sham_rep2_PE1.fastq.gz
mv V4-noTBI_R2_001.fastq.gz noAb_Sham_rep2_PE2.fastq.gz
mv V5-noTBI_R1_001.fastq.gz noAb_Sham_rep3_PE1.fastq.gz
mv V5-noTBI_R2_001.fastq.gz noAb_Sham_rep3_PE2.fastq.gz
mv V6-noTBI_R1_001.fastq.gz noAb_Sham_rep4_PE1.fastq.gz
mv V6-noTBI_R2_001.fastq.gz noAb_Sham_rep4_PE2.fastq.gz

#Rename ablation + sham group
mv V-A2-noTBI_R1_001.fastq.gz Ab_Sham_rep1_PE1.fastq.gz
mv V-A2-noTBI_R2_001.fastq.gz Ab_Sham_rep1_PE2.fastq.gz
mv V-A4-noTBI_R1_001.fastq.gz Ab_Sham_rep2_PE1.fastq.gz
mv V-A4-noTBI_R2_001.fastq.gz Ab_Sham_rep2_PE2.fastq.gz
mv V-A5-noTBI_R1_001.fastq.gz Ab_Sham_rep3_PE1.fastq.gz
mv V-A5-noTBI_R2_001.fastq.gz Ab_Sham_rep3_PE2.fastq.gz
mv V-A6-noTBI_R1_001.fastq.gz Ab_Sham_rep4_PE1.fastq.gz
mv V-A6-noTBI_R2_001.fastq.gz Ab_Sham_rep4_PE2.fastq.gz

#download the mm10 genome file
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit
#convert it to the correct format
twoBitToFa mm10.2bit mm10.fa
#build genome index
hisat2-build mm10.fa mm10

#download .gtf file
wget ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz

#gtf file when downloaded from ensemble is missing 'chr' prefix for the chromosome column, which will cause an error with htseq-count later on
#Correct gtf file by adding chr prefix at the beginning of every line, except header
awk '{ if($1 !~ /^#/){print "chr"$0} else{print $0} }'  Mus_musculus.GRCm38.98.gtf > Mus_musculus.GRCm38.98.corrected.gtf

#Extract splice sites. This is so hisat2 knows to skip intronic sequences when aligning
python3 ~/pyscripts/hisat2_extract_splice_sites.py Mus_musculus.GRCm38.98.corrected.gtf > Mus_musculus.GRCm38.98.splicesites.txt

#Align reads and convert to bam
for i in *PE1.fastq.gz
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_PE." '{print $1}')
    echo $name
    hisat2 -x mm10 --rna-strandness RF --known-splicesite-infile Mus_musculus.GRCm38.98.splicesites.txt -1 $i -2 ${name}_PE2.fastq.gz | samtools view -bS - | samtools sort -n - -o $name.sorted.bam
done

#Counting alignment before and after junk removal for QC
#Should maybe include chrM depending on biological context...
for i in *sorted.bam 
do
	name=$(echo $i | awk -F"/" '{print $NF}' | awk -F".sorted." '{print $1}')
	echo $name
	echo before junk removal
	samtools view -c $i
	samtools sort $i -o $i
	samtools index $i
	samtools view -bh $i chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > $name.noJunk.bam
	echo after junk removal
	samtools view -c $name.noJunk.bam
done

#Removing PCR Duplicates with fixmate and markdup
for i in *.noJunk.bam
do
	name=$(echo $i | awk -F"/" '{print $NF}' | awk -F".noJunk." '{print $1}')
	echo $name
	echo sorting and removing dups
	samtools sort -n $i -o $name.noJunk.bam
        #fixmate requires the bam be sorted by name
	samtools fixmate -m $name.noJunk.bam $name.fixmate.bam
        #markdup requires the bam be sorted by position (default sorting method)
	samtools sort $name.fixmate.bam -o $name.fixmate.bam
	samtools markdup -rs $name.fixmate.bam $name.noDups.bam		
	rm $name.fixmate.bam
done

#HT-Seq, after PCR duplicates were removed ("-bf 1 $i" is used to output only paired-end alignments. Alignments missing a mate throw an error)
for i in *.noDups.bam
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F".noDups." '{print $1}')
    echo $name
    samtools view -bf 1 $i | htseq-count -r pos -f bam \
        --stranded=reverse - Mus_musculus.GRCm38.98.corrected.gtf > $name.gene.counts.txt
done

sh qc.sh
