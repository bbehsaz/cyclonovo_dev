allreconstfile=$1

allspectfile=$2

output=$3

cyclonovoDir=$4

err=$5
threshcore=11
summaryfile=$6
cat $allreconstfile | awk '($(NF)>11){print $(NF-3)}' | sort | uniq > "$output"/masses.txt

echo "Calculating P-values and cleaning up ... "
while read line; do grep "$line" "$allreconstfile" | awk '(NF>3 && $(NF-2)>11){n++; if (n<100){print $1}}' | sed 's/,/ /g' > "$output"/"$line"_reconst.txt  ; done < "$output"/masses.txt
while read line; do grep "$line" "$allreconstfile" | awk '(NF>3 && $(NF-2)>11){n++; if (n<100){print $0}}' > "$output"/"$line"_allinfo.txt  ; done < "$output"/masses.txt

cat $allspectfile  | python "$cyclonovoDir"/scripts/get_by_pepmass.py "$output"/masses.txt "$output"/


while read line; do if [ -f "$output"/"$line".mgf ]; then "$cyclonovoDir"/scripts/print_score "$output"/"$line".mgf "$output"/"$line"_reconst.txt --mass_seq_in --no_filter --product_ion_thresh "$err" --no_merge --make_cyclic --multiple_seq_file --concise_output | grep -v "^#" > "$output"/"$line"_pvals.txt; fi; done < "$output"/masses.txt 


while read line; do if [ -f "$output"/"$line".mgf ]; then paste "$output"/"$line"_allinfo.txt "$output"/"$line"_pvals.txt | awk '{print $(NF-1),$0}' | sort -nr | cut -f2- -d' ' |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6+1"\t"$7}'  > "$output"/"$line"_reconst_pvals.txt; fi; done < "$output"/masses.txt 

while read line; do if [ -f "$output"/"$line"_allinfo.txt ]; then rm "$output"/"$line"_allinfo.txt; fi; done < "$output"/masses.txt

while read line; do if [ -f "$output"/"$line"_pvals.txt ]; then rm "$output"/"$line"_pvals.txt; fi; done < "$output"/masses.txt

while read line; do if [ -f  "$output"/"$line".mgf ]; then rm  "$output"/"$line".mgf; fi; done < "$output"/masses.txt

while read line; do if [ -f  "$output"/"$line"_reconst.txt ]; then rm  "$output"/"$line"_reconst.txt; fi; done < "$output"/masses.txt

rm $allreconstfile
while read line; do if [ -f  "$output"/"$line"_reconst_pvals.txt ]; then cat "$output"/"$line"_reconst_pvals.txt >> $allreconstfile ;rm  "$output"/"$line"_reconst_pvals.txt; fi; done < "$output"/masses.txt
while read line; do if [ -f  $allreconstfile ]; then grep -m1 "$line" $allreconstfile | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6}' >> $summaryfile; fi; done < "$output"/masses.txt
rm "$output"/masses.txt