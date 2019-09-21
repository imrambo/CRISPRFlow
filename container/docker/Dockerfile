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
    && apt-get clean

RUN add-apt-repository -y ppa:deadsnakes/ppa

#RUN cd /build && wget https://www.python.org/ftp/python/2.7.13/Python-2.7.13.tgz \
#    && tar xvzf Python-2.7.13.tgz && cd /build/Python-2.7.13 \
#    && ./configure && make -j4 && make install && make clean \
#    && rm /build/Python-2.7.13.tgz

RUN mkdir /build/

RUN pip3 install biopython \
    && pip3 install pandas \
    && pip3 install scons

RUN useradd -m docker && echo "docker:docker" | chpasswd && adduser docker sudo

USER docker

#CRISPRDetect_2.4
RUN cd /build && git clone https://github.com/davidchyou/CRISPRDetect_2.4.git \
    && cd CRISPRDetect_2.4 \
    && unzip CRISPRDetect_2.4.zip \
    && rm CRISPRDetect_2.4.zip \
    && cpan Parallel::ForkManager
ENV PATH="/build/CRISPRDetect_2.4/CRISPRDetect_2.4:${PATH}"

#CRISPRCasFinder
RUN cd /build && git clone https://github.com/dcouvin/CRISPRCasFinder.git \
    && cd CRISPRCasFinder \
    && sudo bash installer_UBUNTU.sh \
    && ./CRISPRCasFinder.pl -cf CasFinder-2.0.2 -def General -cas -i install_test/sequence.fasta -out Results_test_install –keep \
    && diff Results_test_install/TSV/Cas_REPORT.tsv install_test/Cas_REPORT.tsv \
    && diff Results_test_install/TSV/Crisprs_REPORT.tsv install_test/Crisprs_REPORT.tsv


#------------------------------------------------------------------------------
#Prodigal
#RUN cd /build && git clone https://github.com/hyattpd/Prodigal.git \
#    && cd Prodigal \
#    && make install
#ENV PATH ~/Prodigal:$PATH

#RUN cd /build && wget https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz \
#    && tar -zxf muscle3.8.31_i86linux64.tar.gz \
#    && rm muscle3.8.31_i86linux64.tar.gz
#ENV PATH /build/muscle3.8.31/bin

#RUN cd /build && wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz \
#    && tar -zxf ncbi-blast-2.9.0+-x64-linux.tar.gz \
#    && rm ncbi-blast-2.9.0+-x64-linux.tar.gz
#ENV PATH /build/blast-2.9.0+/bin:$PATH

#RUN cd /build && git clone https://github.com/scapella/trimal.git && cd trimal/source \
#    && make
#ENV PATH /build/trimal/source:$PATH

#ViennaRNA
#RUN cd /build && wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.14.tar.gz \
#    && tar -zxf ViennaRNA-2.4.14.tar.gz \
#    && cd ViennaRNA-2.4.14 \
#    && ./configure && make && make check && make install && make clean \
#    && rm ViennaRNA-2.4.14.tar.gz
