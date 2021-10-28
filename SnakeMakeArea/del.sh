set -euo pipefail
p="/scratch/Shares/layer/workspace/michael_sandbox/layer_lab_chco/results/VariantCalling/SavvySize250/cnv_list_250.bed"
name=$(basename $p)
array=( WES100_S13 WES104_S31 WES88_S29 WES159_S9 WES154_S22 WES98_S3 WES99_S8 WES52_S9 )
for i in "${array[@]}"
do
	echo $i
	cat $p | { grep $i || :; } >> workproband/FindMaxTMP/$name.filtered
done
echo "DONE"
