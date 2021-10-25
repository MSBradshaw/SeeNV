mkdir -p workproband/FindMax/ 2>> workproband/logs/plotting_3.log

string=$(head -1 workproband/SplitParams/3.txt) 2>> workproband/logs/plotting_3.log

region=$(cut -d" " -f1 <<< $string) 2>> workproband/logs/plotting_3.log
svtype=$(cut -d" " -f2 <<< $string) 2>> workproband/logs/plotting_3.log
samplename=$(cut -d" " -f3 <<< $string) 2>> workproband/logs/plotting_3.log
sample=${samplename%.*} 2>> workproband/logs/plotting_3.log
callsfile=$(cut -d" " -f4 <<< $string) 2>> workproband/logs/plotting_3.log
regionclean=$(cut -d" " -f5 <<< $string) 2>> workproband/logs/plotting_3.log
svtype2=$(cut -d" " -f6 <<< $string) 2>> workproband/logs/plotting_3.log

echo "here" 2>> workproband/logs/plotting_3.log

python bin/cnv_plotter.py 			--sample ${samplename%.*} 			--vcf workproband/${samplename%.*}.vcf.gz 			-o workproband/Plots/3.txt 			--scores workproband/MergedAdjZscore/adj_scores.bed.gz 		--exons workproband/Exons/labeled_exons.bed.gz 			--window 100000 			--region ${region} 			--height 7 			--width 5 			--label_exons 			--title "${sample} ${region} ${svtype}" 	--depth workproband/ProbeCoverage/${sample}_probe.cover.mean.stdev.bed 			--max_num_calls workproband/FindMax/all_maxes.txt 	--gnomad_sv workproband/gnomad_sv.bed.gz 			--all_calls workproband/FindMaxTMP/multiple_savvy_calls.txt 
