FROM ubuntu:18.04
MAINTAINER Ian Rambo ian.rambo@utexas.edu

USER root

RUN apt-get update && apt-get install -y apt-utils \
    default-jdk \
    zlib1g-dev \
    default-jre \
    software-properties-common \
    python3.7 \
    python3-pip \
    python2.7 \
    python-pip \
    scons \
    git \
    wget \
    hmmer \
    emboss \
    cd-hit \
    muscle \
    clustalw \
    cd-hit \
    ncbi-blast+ \
    prodigal \
    autoconf \
    automake \
    libtool \
    gcc \
    make \
    g++ \
    unzip \
    sudo \
    parallel \
    vim \
    && apt-get clean

RUN add-apt-repository -y ppa:deadsnakes/ppa

RUN pip3 install biopython pandas scons python-magic
#------------------------------------------------------------------------------
RUN mkdir -p /build/database

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


#RUN cd /build && git clone -b play_notebook --single-branch https://github.com/imrambo/CRISPRFlow.git


RUN mkdir /build/DIAMOND && cd /build/DIAMOND \
    && wget http://github.com/bbuchfink/diamond/releases/download/v0.9.26/diamond-linux64.tar.gz \
    && tar -zxf diamond-linux64.tar.gz \
    && rm diamond-linux64.tar.gz
ENV PATH "$PATH:/build/DIAMOND"




#CAZymes
#http://bcb.unl.edu/dbCAN2/download/Databases/?C=M;O=D
