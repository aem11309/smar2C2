#!/bin/bash
#SBATCH --job-name=j_v1.ProcessTSR
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --mem=80G
#SBATCH --output=j_ProcessTSR.out
#SBATCH --error=j_ProcessTSR.err
#SBATCH --mail-type=END
#SBATCH --mail-user=aem11309@uga.edu

cd $SLURM_SUBMIT_DIR

module load BEDTools
module load R/3.6.2-foss-2019b

###Sort TSR based on Gene and then counts in peak
awk '$1 !~ /scaf/' MaizeLeaf2hrR1_TSRset.txt > MaizeLeaf2hrR1_TSRset_NoScaf.txt

awk '$4 = /+/ {print $1,$2 - 1,$3,"+",$5,$6,$7,$8,$9,$10,$11,$12}' MaizeLeaf2hrR1_TSRset_NoScaf.txt > MaizeLeaf2hrR1_TSRset.plus.bed
awk '$4 !~ /+/ {print $1,$2 - 1,$3,"-",$5,$6,$7,$8,$9,$10,$11,$12}' MaizeLeaf2hrR1_TSRset_NoScaf.txt > MaizeLeaf2hrR1_TSRset.minus.bed
sed -i '1d' MaizeLeaf2hrR1_TSRset.minus.bed

sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSRset.plus.bed
sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSRset.minus.bed

bedtools sort -i MaizeLeaf2hrR1_TSRset.plus.bed > MaizeLeaf2hrR1_TSRset.plus.sort.bed
bedtools sort -i MaizeLeaf2hrR1_TSRset.minus.bed > MaizeLeaf2hrR1_TSRset.minus.sort.bed

rm MaizeLeaf2hrR1_TSRset.plus.bed
rm MaizeLeaf2hrR1_TSRset.minus.bed

mv MaizeLeaf2hrR1_TSRset.plus.sort.bed MaizeLeaf2hrR1_TSRset.plus.bed
mv MaizeLeaf2hrR1_TSRset.minus.sort.bed MaizeLeaf2hrR1_TSRset.minus.bed

bedtools closest -a MaizeLeaf2hrR1_TSRset.plus.bed -b /scratch/aem11309/Working/Maize_Genome/Maize_v4/Zea_mays.AGPv4.gene.simple.bed -iu -D a > MaizeLeaf2hrR1_TSRPlus.ClosGene.bed
bedtools closest -a MaizeLeaf2hrR1_TSRset.minus.bed -b /scratch/aem11309/Working/Maize_Genome/Maize_v4/Zea_mays.AGPv4.gene.simple.bed -id -D a > MaizeLeaf2hrR1_TSRMinus.ClosGene.bed

awk '$17 < 1000 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$16}' MaizeLeaf2hrR1_TSRPlus.ClosGene.bed > MaizeLeaf2hrR1_TSRPlus.ClosGeneU1K.bed
awk '$17 > -1000 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$16}' MaizeLeaf2hrR1_TSRMinus.ClosGene.bed > MaizeLeaf2hrR1_TSRMinus.ClosGeneU1K.bed

sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSRPlus.ClosGeneU1K.bed
sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSRMinus.ClosGeneU1K.bed

awk '$2 != -1' MaizeLeaf2hrR1_TSRPlus.ClosGeneU1K.bed > MaizeLeaf2hrR1_TSRPlus.ClosGeneU1K.e1.bed
awk '$2 != -1' MaizeLeaf2hrR1_TSRMinus.ClosGeneU1K.bed > MaizeLeaf2hrR1_TSRMinus.ClosGeneU1K.e1.bed

rm MaizeLeaf2hrR1_TSRPlus.ClosGeneU1K.bed
rm MaizeLeaf2hrR1_TSRMinus.ClosGeneU1K.bed

mv MaizeLeaf2hrR1_TSRPlus.ClosGeneU1K.e1.bed MaizeLeaf2hrR1_TSRPlus.ClosGeneU1K.bed
mv MaizeLeaf2hrR1_TSRMinus.ClosGeneU1K.e1.bed MaizeLeaf2hrR1_TSRMinus.ClosGeneU1K.bed

cat MaizeLeaf2hrR1_TSRPlus.ClosGeneU1K.bed MaizeLeaf2hrR1_TSRMinus.ClosGeneU1K.bed > MaizeLeaf2hrR1_TSRAll.ClosGeneU1K.bed
bedtools sort -i MaizeLeaf2hrR1_TSRAll.ClosGeneU1K.bed > MaizeLeaf2hrR1_TSRAll.ClosGeneU1K.sort.bed

###Seperate TSS by strand
awk '$1 !~ /scaf/' MaizeLeaf2hrR1_TSSset.txt > MaizeLeaf2hrR1_TSSset_NoScaf.txt

awk '$3 = /+/ {print $1,$2 -1,$2,$4}' MaizeLeaf2hrR1_TSSset_NoScaf.txt > MaizeLeaf2hrR1_TSSPlus.bed
awk '$3 !~ /+/ {print $1,$2 -1,$2,$4}' MaizeLeaf2hrR1_TSSset_NoScaf.txt > MaizeLeaf2hrR1_TSSMinus.bed
sed -i '1d' MaizeLeaf2hrR1_TSSMinus.bed

sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSSPlus.bed
sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSSMinus.bed

###Find closest gene to TSS
bedtools sort -i MaizeLeaf2hrR1_TSSPlus.bed > MaizeLeaf2hrR1_TSSPlus.sort.bed
bedtools sort -i MaizeLeaf2hrR1_TSSMinus.bed > MaizeLeaf2hrR1_TSSMinus.sort.bed

rm MaizeLeaf2hrR1_TSSPlus.bed
rm MaizeLeaf2hrR1_TSSMinus.bed

mv MaizeLeaf2hrR1_TSSPlus.sort.bed MaizeLeaf2hrR1_TSSPlus.bed
mv MaizeLeaf2hrR1_TSSMinus.sort.bed MaizeLeaf2hrR1_TSSMinus.bed

bedtools closest -a MaizeLeaf2hrR1_TSSPlus.bed -b /scratch/aem11309/Working/Maize_Genome/Maize_v4/Zea_mays.AGPv4.gene.simple.bed -iu -D a > MaizeLeaf2hrR1_TSSPlus.ClosGene.bed
bedtools closest -a MaizeLeaf2hrR1_TSSMinus.bed -b /scratch/aem11309/Working/Maize_Genome/Maize_v4/Zea_mays.AGPv4.gene.simple.bed -id -D a > MaizeLeaf2hrR1_TSSMinus.ClosGene.bed

awk '$9 < 1000 {print $5,$6,$7,$8,$4}' MaizeLeaf2hrR1_TSSPlus.ClosGene.bed > MaizeLeaf2hrR1_TSSPlus.ClosGeneU1K.bed
awk '$9 > -1000 {print $5,$6,$7,$8,$4}' MaizeLeaf2hrR1_TSSMinus.ClosGene.bed > MaizeLeaf2hrR1_TSSMinus.ClosGeneU1K.bed

sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSSPlus.ClosGeneU1K.bed
sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSSMinus.ClosGeneU1K.bed

awk '$2 != -1' MaizeLeaf2hrR1_TSSPlus.ClosGeneU1K.bed > MaizeLeaf2hrR1_TSSPlus.ClosGeneU1K.e1.bed
awk '$2 != -1' MaizeLeaf2hrR1_TSSMinus.ClosGeneU1K.bed > MaizeLeaf2hrR1_TSSMinus.ClosGeneU1K.e1.bed

rm MaizeLeaf2hrR1_TSSPlus.ClosGeneU1K.bed
rm MaizeLeaf2hrR1_TSSMinus.ClosGeneU1K.bed

mv MaizeLeaf2hrR1_TSSPlus.ClosGeneU1K.e1.bed MaizeLeaf2hrR1_TSSPlus.ClosGeneU1K.bed
mv MaizeLeaf2hrR1_TSSMinus.ClosGeneU1K.e1.bed MaizeLeaf2hrR1_TSSMinus.ClosGeneU1K.bed

###Find sum of reads in peak under gene
bedtools sort -i MaizeLeaf2hrR1_TSSPlus.ClosGeneU1K.bed > MaizeLeaf2hrR1_TSSPlus.ClosGeneU1K.sort.bed
bedtools sort -i MaizeLeaf2hrR1_TSSMinus.ClosGeneU1K.bed > MaizeLeaf2hrR1_TSSMinus.ClosGeneU1K.sort.bed
cat MaizeLeaf2hrR1_TSSPlus.ClosGeneU1K.sort.bed MaizeLeaf2hrR1_TSSMinus.ClosGeneU1K.sort.bed > MaizeLeaf2hrR1_TSSAll.ClosGeneU1K.bed
bedtools sort -i MaizeLeaf2hrR1_TSSAll.ClosGeneU1K.bed > MaizeLeaf2hrR1_TSSAll.ClosGeneU1K.sort.bed

bedtools merge -i MaizeLeaf2hrR1_TSSPlus.ClosGeneU1K.sort.bed -c 5,4 -o sum,distinct > MaizeLeaf2hrR1_TSSPlus.GeneSum.bed
bedtools merge -i MaizeLeaf2hrR1_TSSMinus.ClosGeneU1K.sort.bed -c 5,4 -o sum,distinct > MaizeLeaf2hrR1_TSSMinus.GeneSum.bed
bedtools merge -i MaizeLeaf2hrR1_TSSAll.ClosGeneU1K.sort.bed -c 5,4 -o sum,distinct > MaizeLeaf2hrR1_TSSAll.GeneSum.bed

sed -i 's/,/\t/g' MaizeLeaf2hrR1_TSSPlus.GeneSum.bed
sed -i 's/,/\t/g' MaizeLeaf2hrR1_TSSMinus.GeneSum.bed
sed -i 's/,/\t/g' MaizeLeaf2hrR1_TSSAll.GeneSum.bed

awk '{print $1,$2,$3,$4,$5}' MaizeLeaf2hrR1_TSSPlus.GeneSum.bed > MaizeLeaf2hrR1_TSSPlus.GeneSum.e1.bed
awk '{print $1,$2,$3,$4,$5}' MaizeLeaf2hrR1_TSSMinus.GeneSum.bed > MaizeLeaf2hrR1_TSSMinus.GeneSum.e1.bed
awk '{print $1,$2,$3,$4,$5}' MaizeLeaf2hrR1_TSSAll.GeneSum.bed > MaizeLeaf2hrR1_TSSAll.GeneSum.e1.bed

awk '$6 != "" {print $1,$2,$3,$4,$6}' MaizeLeaf2hrR1_TSSPlus.GeneSum.bed > MaizeLeaf2hrR1_TSSPlus.GeneSum.e2.bed
awk '$6 != "" {print $1,$2,$3,$4,$6}' MaizeLeaf2hrR1_TSSMinus.GeneSum.bed > MaizeLeaf2hrR1_TSSMinus.GeneSum.e2.bed
awk '$6 != "" {print $1,$2,$3,$4,$6}' MaizeLeaf2hrR1_TSSAll.GeneSum.bed > MaizeLeaf2hrR1_TSSAll.GeneSum.e2.bed

cat MaizeLeaf2hrR1_TSSPlus.GeneSum.e1.bed MaizeLeaf2hrR1_TSSPlus.GeneSum.e2.bed > MaizeLeaf2hrR1_TSSPlus.GeneSum.FINAL.bed
cat MaizeLeaf2hrR1_TSSMinus.GeneSum.e1.bed MaizeLeaf2hrR1_TSSMinus.GeneSum.e2.bed > MaizeLeaf2hrR1_TSSMinus.GeneSum.FINAL.bed
cat MaizeLeaf2hrR1_TSSAll.GeneSum.e1.bed MaizeLeaf2hrR1_TSSAll.GeneSum.e2.bed > MaizeLeaf2hrR1_TSSAll.GeneSum.FINAL.bed

rm MaizeLeaf2hrR1_TSSPlus.GeneSum.e1.bed
rm MaizeLeaf2hrR1_TSSMinus.GeneSum.e1.bed
rm MaizeLeaf2hrR1_TSSAll.GeneSum.e1.bed
rm MaizeLeaf2hrR1_TSSPlus.GeneSum.e2.bed
rm MaizeLeaf2hrR1_TSSMinus.GeneSum.e2.bed
rm MaizeLeaf2hrR1_TSSAll.GeneSum.e2.bed

sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSSPlus.GeneSum.FINAL.bed
sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSSMinus.GeneSum.FINAL.bed
sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSSAll.GeneSum.FINAL.bed

###Match sum reads in peak under gene to individual peak
R CMD BATCH TPM_Match1.R
R CMD BATCH TPM_Match2.R
R CMD BATCH TPM_Match3.R
R CMD BATCH TPM_Match4.R

sed -i 's/"//g' MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.csv
sed -i 's/"//g' MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.csv
sed -i 's/"//g' MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.all.csv
sed -i 's/"//g' MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.all.csv

sed 's/,/\t/g' MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.csv > MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.tsv
sed 's/,/\t/g' MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.csv > MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.tsv
sed 's/,/\t/g' MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.all.csv > MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.all.tsv
sed 's/,/\t/g' MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.all.csv > MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.all.tsv

sed -i '1d' MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.tsv
sed -i '1d' MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.tsv
sed -i '1d' MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.all.tsv
sed -i '1d' MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.all.tsv

###Find peaks that contain at least 10% of the total same stranded reads in peak under gene
awk '$14 != "NA" {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$7/$14}' MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.tsv > MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.bed
awk '$14 != "NA" {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$7/$14}' MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.tsv > MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.bed

awk '$13 > 0.10' MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.bed > MaizeLeaf2hrR1_TSRGenic.GENESUMg10.plus.bed
awk '$13 > 0.10' MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.bed > MaizeLeaf2hrR1_TSRGenic.GENESUMg10.minus.bed

awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' MaizeLeaf2hrR1_TSRGenic.GENESUMg10.plus.bed > MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.plus.bed
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' MaizeLeaf2hrR1_TSRGenic.GENESUMg10.minus.bed > MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.minus.bed

sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSRGenic.GENESUMg10.plus.bed
sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSRGenic.GENESUMg10.minus.bed
sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.plus.bed
sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.minus.bed

###Find peaks that contain at least 10% of the total reads in peak under gene
awk '$14 != "NA" {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$7/$14}' MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.all.tsv > MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.all.bed
awk '$14 != "NA" {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$7/$14}' MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.all.tsv > MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.all.bed

awk '$13 > 0.10' MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.all.bed > MaizeLeaf2hrR1_TSRGenic.GENESUMg10.plus.all.bed
awk '$13 > 0.10' MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.all.bed > MaizeLeaf2hrR1_TSRGenic.GENESUMg10.minus.all.bed

awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' MaizeLeaf2hrR1_TSRGenic.GENESUMg10.plus.all.bed > MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.plus.all.bed
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' MaizeLeaf2hrR1_TSRGenic.GENESUMg10.minus.all.bed > MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.minus.all.bed

sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSRGenic.GENESUMg10.plus.all.bed
sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSRGenic.GENESUMg10.minus.all.bed
sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.plus.all.bed
sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.minus.all.bed

###Find nongenic peaks
awk '$12 == "NA"' MaizeLeaf2hrR1_TSRset.txt > MaizeLeaf2hrR1_TSRnonGenic.txt

awk '$4 = /+/' MaizeLeaf2hrR1_TSRnonGenic.txt > MaizeLeaf2hrR1_TSRnonGenic.plus.txt
awk '$4 !~ /+/' MaizeLeaf2hrR1_TSRnonGenic.txt > MaizeLeaf2hrR1_TSRnonGenic.minus.txt

sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSRnonGenic.plus.txt
sed -i 's/ //g' MaizeLeaf2hrR1_TSRnonGenic.minus.txt

###Combine genic and nongenic peaks
cat MaizeLeaf2hrR1_TSRnonGenic.plus.txt MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.plus.bed > MaizeLeaf2hrR1_TSRall.plus.bed
cat MaizeLeaf2hrR1_TSRnonGenic.minus.txt MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.minus.bed > MaizeLeaf2hrR1_TSRall.minus.bed

##Find best single nucleotide per peak using same stranded weedout
bedtools intersect -a MaizeLeaf2hrR1_TSSPlus.bed -b MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.plus.bed -wa -wb > MaizeLeaf2hrR1_TSSPlus_TSRInt.bed
bedtools intersect -a MaizeLeaf2hrR1_TSSMinus.bed -b MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.minus.bed -wa -wb > MaizeLeaf2hrR1_TSSMinus_TSRInt.bed

sort -k5 -k6 -k4nr MaizeLeaf2hrR1_TSSPlus_TSRInt.bed > MaizeLeaf2hrR1_TSSPlus_TSRInt.sort.bed
sort -k5 -k6 -k4nr MaizeLeaf2hrR1_TSSMinus_TSRInt.bed > MaizeLeaf2hrR1_TSSMinus_TSRInt.sort.bed

awk '!a[$5 + $6]++ {print $1,$2,$3,$4}' MaizeLeaf2hrR1_TSSPlus_TSRInt.sort.bed > MaizeLeaf2hrR1_TSSPlus_TSRInt.TopSNT.bed
awk '!a[$5 + $6]++ {print $1,$2,$3,$4}' MaizeLeaf2hrR1_TSSMinus_TSRInt.sort.bed > MaizeLeaf2hrR1_TSSMinus_TSRInt.TopSNT.bed

sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSSPlus_TSRInt.TopSNT.bed
sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSSMinus_TSRInt.TopSNT.bed

##Find best single nucleotide per peak using both strand weedout
bedtools intersect -a MaizeLeaf2hrR1_TSSPlus.bed -b MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.plus.all.bed -wa -wb > MaizeLeaf2hrR1_TSSPlus_TSRInt.all.bed
bedtools intersect -a MaizeLeaf2hrR1_TSSMinus.bed -b MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.minus.all.bed -wa -wb > MaizeLeaf2hrR1_TSSMinus_TSRInt.all.bed

sort -k5 -k6 -k4nr MaizeLeaf2hrR1_TSSPlus_TSRInt.all.bed > MaizeLeaf2hrR1_TSSPlus_TSRInt.sort.all.bed
sort -k5 -k6 -k4nr MaizeLeaf2hrR1_TSSMinus_TSRInt.all.bed > MaizeLeaf2hrR1_TSSMinus_TSRInt.sort.all.bed

awk '!a[$5 + $6]++ {print $1,$2,$3,$4}' MaizeLeaf2hrR1_TSSPlus_TSRInt.sort.all.bed > MaizeLeaf2hrR1_TSSPlus_TSRInt.TopSNT.all.bed
awk '!a[$5 + $6]++ {print $1,$2,$3,$4}' MaizeLeaf2hrR1_TSSMinus_TSRInt.sort.all.bed > MaizeLeaf2hrR1_TSSMinus_TSRInt.TopSNT.all.bed

sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSSPlus_TSRInt.TopSNT.all.bed
sed -i 's/ /\t/g' MaizeLeaf2hrR1_TSSMinus_TSRInt.TopSNT.all.bed

###Cleanup
mkdir MaizeLeaf2hrR1_Final_Outputs
mv MaizeLeaf2hrR1_TSSPlus_TSRInt.TopSNT.bed MaizeLeaf2hrR1_Final_Outputs/
mv MaizeLeaf2hrR1_TSSMinus_TSRInt.TopSNT.bed MaizeLeaf2hrR1_Final_Outputs/
mv MaizeLeaf2hrR1_TSSMinus_TSRInt.sort.bed MaizeLeaf2hrR1_Final_Outputs/
mv MaizeLeaf2hrR1_TSSPlus_TSRInt.sort.bed MaizeLeaf2hrR1_Final_Outputs/

rm MaizeLeaf2hrR1_TSSPlus_TSRInt.sort.bed
rm MaizeLeaf2hrR1_TSSMinus_TSRInt.sort.bed
rm MaizeLeaf2hrR1_TSSPlus_TSRInt.bed
rm MaizeLeaf2hrR1_TSSMinus_TSRInt.bed

mv MaizeLeaf2hrR1_TSSPlus.bed MaizeLeaf2hrR1_Final_Outputs/
mv MaizeLeaf2hrR1_TSSMinus.bed MaizeLeaf2hrR1_Final_Outputs/
mv MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.plus.bed MaizeLeaf2hrR1_Final_Outputs/
mv MaizeLeaf2hrR1_TSRGenic.GENESUMg10.OriginalFormat.minus.bed MaizeLeaf2hrR1_Final_Outputs/
mv MaizeLeaf2hrR1_TSRall.plus.bed MaizeLeaf2hrR1_Final_Outputs/
mv MaizeLeaf2hrR1_TSRall.minus.bed MaizeLeaf2hrR1_Final_Outputs/
mv MaizeLeaf2hrR1_TSRnonGenic.plus.txt MaizeLeaf2hrR1_Final_Outputs/
mv MaizeLeaf2hrR1_TSRnonGenic.minus.txt MaizeLeaf2hrR1_Final_Outputs/

rm MaizeLeaf2hrR1_TSRMinus.ClosGene.bed
rm MaizeLeaf2hrR1_TSRPlus.ClosGene.bed
rm MaizeLeaf2hrR1_TSRset_NoScaf.txt
rm MaizeLeaf2hrR1_TSSMinus.GeneSum.bed 
rm MaizeLeaf2hrR1_TSSPlus.GeneSum.bed 

rm MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.bed
mv MaizeLeaf2hrR1_TSRnonGenic.txt MaizeLeaf2hrR1_Final_Outputs
rm MaizeLeaf2hrR1_TSSMinus.ClosGeneU1K.bed
mv MaizeLeaf2hrR1_TSSset.txt MaizeLeaf2hrR1_Final_Outputs
rm MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.csv
mv MaizeLeaf2hrR1_TSRset.minus.bed MaizeLeaf2hrR1_Final_Outputs
mv MaizeLeaf2hrR1_TSSMinus.GeneSum.FINAL.bed MaizeLeaf2hrR1_Final_Outputs
rm MaizeLeaf2hrR1_TSRGenic.GENESUM.minus.tsv
mv MaizeLeaf2hrR1_TSRset.plus.bed MaizeLeaf2hrR1_Final_Outputs
rm MaizeLeaf2hrR1_TSSPlus.ClosGene.bed
rm MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.bed
mv MaizeLeaf2hrR1_TSRset.tab MaizeLeaf2hrR1_Final_Outputs
rm MaizeLeaf2hrR1_TSSPlus.ClosGeneU1K.bed
mv MaizeLeaf2hrR1_TSRGenic.GENESUMg10.minus.bed MaizeLeaf2hrR1_Final_Outputs
rm MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.csv
mv MaizeLeaf2hrR1_TSRset.txt MaizeLeaf2hrR1_Final_Outputs
mv MaizeLeaf2hrR1_TSSPlus.GeneSum.FINAL.bed MaizeLeaf2hrR1_Final_Outputs
mv MaizeLeaf2hrR1_TSRGenic.GENESUMg10.plus.bed MaizeLeaf2hrR1_Final_Outputs
rm MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.tsv
rm MaizeLeaf2hrR1_TSSMinus.ClosGene.bed
mv MaizeLeaf2hrR1_TSSset_NoScaf.txt MaizeLeaf2hrR1_Final_Outputs

