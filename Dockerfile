FROM ubuntu:18.04

LABEL maintainer="ian.rambo@utexas.edu"

#python3.8 \
#python3-pip \

RUN apt-get update && apt-get install -y --no-install-recommends aufs-tools \
    build-essential \
    python2.7 \
    git \
    wget \
    hmmer \
    ncbi-blast+ \
    prodigal \
    cd-hit \
    gcc \
    make \
    cmake \
    g++ \
    unzip \
    libstdc++6 \
    libgomp1 \
    libpcre3 \
    libpcre3-dev \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean


RUN pip install biopython pandas python-magic

RUN mkdir -p /build/bin
WORKDIR /build

#CRISPRDetect_2.4
RUN git clone https://github.com/davidchyou/CRISPRDetect_2.4.git \
    && cd CRISPRDetect_2.4 \
    && unzip CRISPRDetect_2.4.zip \
    && rm CRISPRDetect_2.4.zip \
    && cpan Parallel::ForkManager \
    && cd CRISPRDetect_2.4 \
    && chmod -R 755 . && chmod 777 tmp

ENV PATH "$PATH:/build/CRISPRDetect_2.4/CRISPRDetect_2.4"

WORKDIR /build
#MacSyFinder
RUN wget https://bintray.com/gem-pasteur/MacSyFinder/download_file?file_path=macsyfinder-1.0.5.tar.gz -O macsyfinder-1.0.5.tar.gz \
    && tar -zxvf macsyfinder-1.0.5.tar.gz \
    && cd macsyfinder-1.0.5 \
    && python2.7 setup.py build \
    && python2.7 setup.py test -vv \
    && python2.7 setup.py install \
    && rm -rf ./data \
    && rm -rf ./doc \
    && rm -rf ./macsyview

COPY ./bin /build/bin
