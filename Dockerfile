# add the platform to avoid problems related to apple silicon processors
FROM --platform=linux/amd64 ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive
ARG USER_ID=root
ENV USER_ID $USER_ID

RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*


## Change user to $USER_ID
USER $USER_ID
WORKDIR /home/$USER_ID
ENV HOME /home/$USER_ID

## Git install
RUN apt-get -y update
RUN apt-get -y install git

####### R
RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base 
## enable R package install
RUN chmod a+w /usr/local/lib/R/site-library

RUN apt-get update && apt-get install -y --no-install-recommends \
  libgit2-dev \
  libssl-dev \
  libcurl4-openssl-dev \
  libfontconfig1-dev \
  zlib1g-dev \
  libharfbuzz-dev \
  libfribidi-dev  \
  libfreetype6-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev \
  libxml2-dev \
  libbz2-dev \
  liblapack-dev \
  gfortran \
  cmake


## Install ANACONDA 
ENV ANACONDA_PATH /home/$USER_ID/anaconda3
ENV ANACONDA_INSTALLER Anaconda3-2019.10-Linux-x86_64.sh

# this is to install igraph -> you need the libgfortran.so.4 which is installed in ~/anaconda3/lib
# to find where libgfortran.so.4 is -> find $ANACONDA_PATH -name "libgfortran.so.4"
ENV LD_LIBRARY_PATH $ANACONDA_PATH/lib
ENV PATH $ANACONDA_PATH/bin:$PATH

# RUN echo "export PATH=$ANACONDA_PATH/bin:$PATH" | tee -a /home/.bashrc
# RUN echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ANACONDA_PATH/lib" | tee -a /home/.bashrc

# install 
RUN mkdir /home/root/download && cd /home/root/download && \
  wget https://repo.anaconda.com/archive/$ANACONDA_INSTALLER -nv -q 
RUN /bin/bash ~/download/$ANACONDA_INSTALLER -b

RUN conda update --yes conda && \
  conda update --yes anaconda && \
  conda update --yes --all && \
  conda clean --yes --all && \
  conda init bash
  # conda init bash is to use "conda activate ..." but it requires a restart of the shell

# ## clone the mobster repo
# RUN git clone --branch binomial_noise https://github.com/caravagnalab/mobster.git
RUN git clone https://github.com/Militeee/HMOBSTER.git

## create a conda env
RUN conda create -n mobsterh python=3.10 -y
RUN . $ANACONDA_PATH/bin/activate mobsterh && \
  conda install pip

RUN cd ~/HMOBSTER && \
  $ANACONDA_PATH/envs/mobsterh/bin/python3.10 -m pip --default-timeout=100 install .
RUN rm -R $HOME/HMOBSTER

RUN R -e 'install.packages("remotes"); install.packages("reticulate")'
RUN R -e 'remotes::install_github("caravagnalab/VIBER")'
RUN R -e 'remotes::install_github("caravagnalab/CNAqc")'
RUN R -e 'install.packages("sequenza")'
RUN R -e 'reticulate::use_condaenv("mobsterh"); remotes::install_github("caravagnalab/mobster", ref="binomial_noise")'
RUN R -e 'remotes::install_github("caravagn/evoverse")'
RUN R -e 'remotes::install_github("caravagnalab/CNAqc")'

