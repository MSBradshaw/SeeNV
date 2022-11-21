conda create --name seenv --file Install/seenv.env -y

conda activate seenv

conda install -c bioconda pysam -y
conda install -c bioconda bedtools -y
conda config --add channels bioconda  -y
conda config --add channels conda-forge -y
conda install samtools==1.11 -y
conda install snakemake-minimal>=5.24.1 -y
conda install -c bioconda mosdepth -y

# find where the conda env bin is located
conda_loc=$(which python)
# strip the python on the end of that location
conda_bin=$(dirname $conda_loc)
# add SeeNV's python script to the conda env bin
cp bin/*.py $conda_bin

# install gargs
wget https://github.com/brentp/gargs/releases/download/v0.3.9/gargs_linux
chmod 770 gargs_linux
mv gargs_linux $conda_bin/gargs

cp bin/seenv.py $conda_bin/seenv
cp bin/build_panel.snake $conda_bin/build_panel.snake
cp bin/run_proband.snake $conda_bin/run_proband.snake
cp panel_config.json $conda_bin/panel_config.json
cp proband_config.json $conda_bin/proband_config.json
