conda create --name cnviz --file Install/cnviz.env -y

conda activate cnviz

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
# add cnviz's python script to the conda env bin
cp bin/*.py $conda_bin
cp bin/*.snake $conda_bin
cp bin/cnviz $conda_bin

# install gargs
wget https://github.com/brentp/gargs/releases/download/v0.3.9/gargs_linux
chmod 770 gargs_linux
mv gargs_linux $conda_bin/gargs
