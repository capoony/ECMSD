WD=${1}
Output=${2}

echo "*********************"
echo ""
echo "Have a cup of coffee, this may take a while... "
echo '''
   ( (
    ) )
  ........
  |      |]
  \      /
   `----´ '''

echo ""
sleep 2

## activate the conda environment
eval "$(conda shell.bash hook)"

## create a conda environment for the scripts
conda create \
    -p ${WD}/scripts/conda_env \
    -y \
    -c bioconda \
    -c conda-forge \
    mamba 2>${Output}/logs/conda_env.log

## activate the conda environment
conda activate ${WD}/scripts/conda_env

# Force uninstall all OpenJDK versions
${WD}/scripts/conda_env/bin/mamba remove --force openjdk
${WD}/scripts/conda_env/bin/mamba clean --all -y

# Reinstall OpenJDK cleanly
${WD}/scripts/conda_env/bin/mamba install -y -c conda-forge openjdk=17

## install the required packages
${WD}/scripts/conda_env/bin/mamba install \
    -y \
    -c bioconda \
    -c conda-forge \
    minimap2 pigz agbiome::bbtools r r-base r-tidyverse openjdk 2>${Output}/logs/conda_env.log

##
echo "Conda environment created and packages installed successfully."
