# Download example samples
The file `download_1KG_samples.sh` contains a series of commands to download 12 samples' .bam and .bai files from 1000 Genomes Project's S3 bucket. The script uses the `aws` command which can be installed [here](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).

It is recommended the sample list files `ref_panel_sample_list.tsv`, and `proband_sample_list.tsv` use absolute file paths for `.bam` and `.bai` files, they include example absolute paths right now and need to be changed to the location on your system.
