rm qc_metrics.txt
echo -e "name\tall.counts\twithout.junk\tunique.counts" >> qc_metrics.txt

for bam in *sorted.bam
do
    name=$(echo $bam | awk -F".sorted.bam" '{print $1}')
    echo $name
    ALL_COUNTS=`samtools view -c $name.sorted.bam`
    WITHOUT_JUNK=`samtools view -c $name.noJunk.bam`
    UNIQUE_COUNTS=`samtools view -c $name.noDups.bam`
    echo -e "${name}\t$ALL_COUNTS\t$WITHOUT_JUNK\t$UNIQUE_COUNTS" >> qc_metrics.txt
done

Rscript qc.R
