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


RUN cd /build && git clone -b play_notebook --single-branch https://github.com/imrambo/CRISPRFlow.git

###---Databases---###
RUN mkdir /build/database/TIGRFAM && cd /build/database/TIGRFAM \
    && wget ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.tar.gz \
    && tar -zxf TIGRFAMs_15.0_HMM.tar.gz \
    && rm TIGRFAMs_15.0_HMM.tar.gz

RUN mkdir /build/database/ANTIFAM && cd /build/database/ANTIFAM \
    && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz \
    && tar -zxf Antifam.tar.gz \
    && rm Antifam.tar.gz

RUN mkdir /build/database/PFAM-A && cd /build/database/PFAM-A \
    && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz
