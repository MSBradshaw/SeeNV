set -euo pipefail
IFS=$'\n\t'
get_rpm_rates.py -m workproband/Mosdepth/WES88_S29.per-base.bed.gz --regions_file workproband/Probes/probes.sorted.bed.gz --num_reads workproband/TotalReads/total_read.txt | bgzip -c > workproband/RPM/WES88_S29.probe.rpm_rate.bed.gz
