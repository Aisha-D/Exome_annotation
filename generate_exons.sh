## run this script on command line 
# sort the file so that it can be split by chromosome
mysql -h genome-mysql.soe.ucsc.edu -ugenome -A -e "select * from ncbiRefSeq" hg38 > allTranscripts_hg38.bed
sed '1d' allTranscripts_hg38.bed | awk -F "\t" {'print $3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16'} |sort -k 3V,3 -k 5n,5 -o allTranscript_hg38.sorted.bed

mkdir -p split_results
for chr in `cut -f 1 allTranscript_hg38.sorted.bed | sort | uniq`; do
                grep -w $chr allTranscript_hg38.sorted.bed > split_results/$chr.output.bed
done

# in another directiory same level run the script at the the level above the two folders
mkdir -p split_reformatted
cd split_reformatted
for f in ../split_results/*.bed
do
	echo $f
    python3 ../exons.py  -f $f
done

# accidentally saved all header so run this to remove headers
for f in split_reformatted/*.bed
do
    sed -i '1d' $f
done

cat split_reformatted/*bed > allExome_reformatted_hg38.bed

wc -l allExomeReformatted.bed #There are 2,302 exons missing from the bed file generated from ucsc

# clean up the file and order it
# easier way probably exists but for now this is fine
cut -f2 allExome_reformatted_hg38.bed | awk -F "_" {'print $1'} | sed 's/chr//g'  > tmp.bed
paste allExome_reformatted_hg38.bed tmp.bed | awk -F "\t" {'print $8"\t"$5"\t"$6"\t"$4"\t"$1"\t"$7'} | sort -k4 -k5 -k6 > exons_210805.bed