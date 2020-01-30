FROM ubuntu:18.04
MAINTAINER Ian Rambo ian.rambo@utexas.edu
<<<<<<< HEAD
=======
#Dockerfile for CRISPRFlow pipeline
>>>>>>> f_crispr_clusters

USER root

RUN apt-get update && apt-get install -y apt-utils \
    default-jdk \
    zlib1g-dev \
    default-jre \
    software-properties-common \
<<<<<<< HEAD
=======
    libdb-dev \
    curl \
>>>>>>> f_crispr_clusters
    python3.7 \
    python3-pip \
    python2.7 \
    python-pip \
    git \
    wget \
    hmmer \
    emboss \
    cd-hit \
    muscle \
    clustalw \
<<<<<<< HEAD
    cd-hit \
=======
>>>>>>> f_crispr_clusters
    ncbi-blast+ \
    prodigal \
    autoconf \
    automake \
    libtool \
    gcc \
    make \
<<<<<<< HEAD
    g++ \
    unzip \
    sudo \
    && apt-get clean

RUN add-apt-repository -y ppa:deadsnakes/ppa

RUN pip3 install biopython pandas scons python-magic
#------------------------------------------------------------------------------
RUN mkdir /build/
=======
    cmake \
    g++ \
    unzip \
    sudo \
    parallel \
    libstdc++6 \
    libgomp1 \
    libpcre3 \
    libpcre3-dev \
    && apt-get clean

RUN pip3 install biopython pandas python-magic
#------------------------------------------------------------------------------
RUN mkdir -p /build
>>>>>>> f_crispr_clusters

#CRISPRDetect_2.4
RUN cd /build && git clone https://github.com/davidchyou/CRISPRDetect_2.4.git \
    && cd CRISPRDetect_2.4 \
    && unzip CRISPRDetect_2.4.zip \
    && rm CRISPRDetect_2.4.zip \
<<<<<<< HEAD
    && cpan Parallel::ForkManager
ENV PATH="/build/CRISPRDetect_2.4/CRISPRDetect_2.4:${PATH}"
=======
    && cpan Parallel::ForkManager \
    && cd CRISPRDetect_2.4 \
    && chmod -R 755 . && chmod 777 tmp

ENV PATH "$PATH:/build/CRISPRDetect_2.4/CRISPRDetect_2.4"
>>>>>>> f_crispr_clusters

#MacSyFinder
RUN cd /build && wget https://bintray.com/gem-pasteur/MacSyFinder/download_file?file_path=macsyfinder-1.0.5.tar.gz \
    && tar -zxf download_file?file_path=macsyfinder-1.0.5.tar.gz \
    && rm download_file?file_path=macsyfinder-1.0.5.tar.gz \
    && cd macsyfinder-1.0.5 \
    && python2.7 setup.py build \
    && python2.7 setup.py test -vv \
    && python2.7 setup.py install
<<<<<<< HEAD
=======

RUN cd /build && git clone -b f_crispr_transposon --single-branch https://github.com/imrambo/CRISPRFlow.git
>>>>>>> f_crispr_clusters
