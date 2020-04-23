# CRISPRFlow version 1.0.0.1

Pipeline to search for CRISPR-Cas systems.

Running with Docker
1) Build your docker image

cd container/docker && docker build -t <TAG> .

2) Run your docker container
docker run -it <TAG> "/bin/bash"

#REFERENCES
CRISPRDetect
Biswas, A., Staals, R. H., Morales, S. E., Fineran, P. C. & Brown, C. M. CRISPRDetect: a flexible algorithm to define CRISPR arrays. BMC Genomics 17, 356 (2016)

MacSyFinder
Prodigal
HMMER
CRISPRCasFinder group (DEF and profiles)
Banfield group (Cas14 sequences for HMM)
