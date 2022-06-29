#!/bin/bash
# in the installation script ensure that #!/bin/bash is the proper path

threads="1"

# read in command line parameters
while getopts ":i:p:c:o:r:g:a:t:" opt; do
  case $opt in
    i) inputsamples="$OPTARG"
    ;;
    p) probes="$OPTARG"
    ;;
    c) calls="$OPTARG"
    ;;
    o) outputdir="$OPTARG"
    ;;
    r) refpanel="$OPTARG"
    ;;
    g) genesfile="$OPTARG"
    ;;
    a) gnomadfile="$OPTARG"
    ;;
    t) threads="$OPTARG"
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

# create the workproband dir, dump sym linked files in to a location they can be easily accessed by SnakeMake 
# will forably overwrite input files of the same name. All samples named MUST be unique.
#rm -rf workproband
mkdir -p workproband
cp $inputsamples workproband/proband.samples
cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {3} workproband/{0}.bam"
cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {4} workproband/{0}.bai"
cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {5} workproband/{0}.vcf.gz"
cat $inputsamples | bin/gargs --sep="\t" "ln -f -s {6} workproband/{0}.vcf.gz.tbi"

mkdir -p workproband/Probes
ln -s "$(cd "$(dirname "$probes")"; pwd)/$(basename "$probes")" workproband/Probes/probes.original.bed
mkdir -p workproband/Calls
ln -s "$(cd "$(dirname "$calls")"; pwd)/$(basename "$calls")" workproband/Calls/all.calls.original.txt
ln -s "$(cd "$(dirname "$genesfile")"; pwd)/$(basename "$genesfile")" workproband/genes.bed.gz
ln -s "$(cd "$(dirname "$gnomadfile")"; pwd)/$(basename "$gnomadfile")" workproband/gnomad_sv.bed.gz
ln -s "$(cd "$(dirname "$gnomadfile")"; pwd)/$(basename "$gnomadfile")".tbi workproband/gnomad_sv.bed.gz.tbi

# put the name of the output dirrectory into a textfile for Snakemake to read in
echo $outputdir > workproband/outputdir.txt
# put the path to the reference panel into a textfie for Snakemake to read in
echo $refpanel > workproband/reference_panel.txt

# start the snakemake pipeline now they input files are in there proper locations
snakemake -c $threads -s run_proband.snake --configfile proband_config.json --printshellcmds --rerun-incomplete
