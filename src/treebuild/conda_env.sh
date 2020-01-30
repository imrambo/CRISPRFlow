#!/bin/bash
#Conda environment for treebuild SCons builder
# envName=$1
# conda create -n $envName python=3.7 && \
# conda install -n $envName -c conda-forge scons=3.1.1 pillow=6.0.0 matplotlib=3.0.2 && \
# conda install -n $envName -c schrodinger pymol=0.1.0 && \
# conda install -n $envName -c anaconda numpy=1.15.4 pandas=0.23.4
#conda install -n $envName -c bioconda cd-hit cd-hit-auxtools muscle trimal hhsuite=3.2.0


envName=$1
conda create -n $envName python=3.5 && \
conda install -n $envName -c conda-forge scons pillow matplotlib && \
conda install -n $envName -c schrodinger pymol && \
conda install -n $envName -c anaconda numpy pandas
conda install -n $envName -c bioconda cd-hit cd-hit-auxtools muscle trimal hhsuite
