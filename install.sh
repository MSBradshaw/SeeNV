conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create --name seenv --file Install/seenv.env -y

conda activate seenv

conda install -c bioconda pysam -y
conda install -c bioconda bedtools -y
conda config --add channels bioconda  -y
conda config --add channels conda-forge -y
conda install samtools==1.11 -y
conda install snakemake-minimal>=5.24.1 -y
conda install -c bioconda mosdepth -y

# add SeeNV's python script to the conda env bin
cp bin/*.py $CONDA_PREFIX/bin

# install gargs
wget https://github.com/brentp/gargs/releases/download/v0.3.9/gargs_linux
chmod 770 gargs_linux
mv gargs_linux $CONDA_PREFIX/bin/gargs

cp bin/seenv.py $CONDA_PREFIX/bin/seenv
cp bin/build_panel.snake $CONDA_PREFIX/bin/build_panel.snake
cp bin/run_proband.snake $CONDA_PREFIX/bin/run_proband.snake
cp panel_config.json $CONDA_PREFIX/bin/panel_config.json
cp proband_config.json $CONDA_PREFIX/bin/proband_config.json
