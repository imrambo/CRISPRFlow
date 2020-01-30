FROM ubuntu:18.04
MAINTAINER Ian Rambo ian.rambo@utexas.edu
#Dockerfile for CRISPRFlow pipeline

USER root

RUN apt-get update && apt-get install -y apt-utils \
    default-jdk \
    zlib1g-dev \
    default-jre \
    software-properties-common \
    libdb-dev \
    curl \
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
    ncbi-blast+ \
    prodigal \
    autoconf \
    automake \
    libtool \
    gcc \
    make \
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

#CRISPRDetect_2.4
RUN cd /build && git clone https://github.com/davidchyou/CRISPRDetect_2.4.git \
    && cd CRISPRDetect_2.4 \
    && unzip CRISPRDetect_2.4.zip \
    && rm CRISPRDetect_2.4.zip \
    && cpan Parallel::ForkManager \
    && cd CRISPRDetect_2.4 \
    && chmod -R 755 . && chmod 777 tmp

ENV PATH "$PATH:/build/CRISPRDetect_2.4/CRISPRDetect_2.4"

#MacSyFinder
RUN cd /build && wget https://bintray.com/gem-pasteur/MacSyFinder/download_file?file_path=macsyfinder-1.0.5.tar.gz \
    && tar -zxf download_file?file_path=macsyfinder-1.0.5.tar.gz \
    && rm download_file?file_path=macsyfinder-1.0.5.tar.gz \
    && cd macsyfinder-1.0.5 \
    && python2.7 setup.py build \
    && python2.7 setup.py test -vv \
    && python2.7 setup.py install

#DIAMOND
#RUN mkdir /build/DIAMOND && cd /build/DIAMOND \
#    && wget http://github.com/bbuchfink/diamond/releases/download/v0.9.26/diamond-linux64.tar.gz \
#    && tar -zxf diamond-linux64.tar.gz \
#    && rm diamond-linux64.tar.gz
#ENV PATH "$PATH:/build/DIAMOND"

#HHSuite
#RUN cd /build \
#    && git clone https://github.com/soedinglab/hh-suite.git \
#    && mkdir -p hh-suite/build \
#    && cd hh-suite/build \
#    && cmake -DCMAKE_INSTALL_PREFIX=. .. \
#    && make -j 4 \
#    && make install

#ENV PATH "$PATH:/build/hh-suite/build/bin"
#ENV PATH "$PATH:/build/hh-suite/build/scripts"


RUN cd /build && git clone -b f_crispr_clusters --single-branch https://github.com/imrambo/CRISPRFlow.git
