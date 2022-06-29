python cnviz.py -p -i pos.samples  \
-s ../WES_TargetCoverage_v2.bed \
-c six_pos_calls.txt \
-o ProbandTestOut \
-r SmallNegPanelDB \
-g Homo_sapiens.GRCh37.82.genes.bed.gz \
-a ../gnomad_v2.1_sv.sites.bed.gz \
-t 50
