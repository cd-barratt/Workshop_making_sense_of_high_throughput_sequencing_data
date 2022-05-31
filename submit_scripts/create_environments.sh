# This is how to make a new virtualenv
# delete previous ones with rm -rf if necessary

#virtualenv
touch ~/.easybuild-yes # then log out and log back in again

module purge
cd /public/barratt #this needs to be the root directory of the user e.g. cd /gpfs0/home/barratt
module load foss/2018b
module load foss/2018b Tkinter/2.7.15-Python-2.7.15
virtualenv virtual_env_RAD_seq
source virtual_env_RAD_seq/bin/activate
pip install numpy
pip install scipy
pip install matplotlib
pip install Ipython
wget https://bitbucket.org/gutenkunstlab/dadi/get/1.7.0.tar.gz
tar xzf 1.7.0.tar.gz
ls
rm -rf 1.7.0.tar.gz
cd gutenkunstlab-dadi-01068e3ad50d/
python setup.py install
cd ../
rm -rf gutenkunstlab-dadi-01068e3ad50d/
pip list
pip freeze > /public/barratt/virtual_env_RAD_seq/virtual_env_RAD_seq_requirements.txt #this needs to be /gpfs0/home/barratt/virtual_env_RAD_seq/virtualenv_requirements.txt
cd /public/barratt/virtual_env_RAD_seq/

#conda
module purge
module load miniconda/3
# add bioconda channel which includes most bioinformatics packages
conda config --add channels bioconda
conda create -n conda_env_RADseq # this will put the environment in your own (hidden) directory e.g. /home/barratt/.conda/envs/conda_env_RADseq
conda install plink
source activate conda_env_RADseq

# You activate these environments from the qsub scripts before running analyses (so they can pick up the bioinformatic tools/packages you need)