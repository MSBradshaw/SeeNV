#!/bin/bash
# read in command line parameters
while getopts ":i:p:c:o:" opt; do
  case $opt in
    i) inputsamples="$OPTARG"
    ;;
    p) probes="$OPTARG"
    ;;
    c) calls="$OPTARG"
    ;;
    o) outputdir="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac

  case $OPTARG in
    -*) echo "Option $opt needs a valid argument"
    exit 1
    ;;
  esac
done

# create the workpanel dir, dump sym linked files in to a location they can be easily accessed by SnakeMake 
# will forably overwrite input files of the same name. All samples named MUST be unique.
#rm -rf workpanel
mkdir -p workpanel
cp $inputsamples workpanel/panel.samples
cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {3} workpanel/{0}.bam"
cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {4} workpanel/{0}.bai"
cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {5} workpanel/{0}.vcf.gz"
cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {6} workpanel/{0}.vcf.gz.tbi"

mkdir -p workpanel/Probes
ln -f -s "$(cd "$(dirname "$probes")"; pwd)/$(basename "$probes")" workpanel/Probes/probes.original.bed
mkdir -p workpanel/Calls
ln -f -s "$(cd "$(dirname "$calls")"; pwd)/$(basename "$calls")" workpanel/Calls/all.calls.original.txt

# put the name of the output dirrectory into a textfile for Snakemake to read in
echo $outputdir > workpanel/outputdir.txt

# start the snakemake pipeline now they input files are in there proper locations
snakemake -c 32 -s build_panel.snake --configfile panel_config.json
