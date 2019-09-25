# CRISPRFlow

Pipeline to search for CRISPR-Cas systems.

1) Build your docker container

cd container/docker && docker build -t <TAG> .

2) Run your docker container
docker run -it <TAG> "/bin/bash"

#REFERENCES
CRISPRDetect
Biswas, A., Staals, R. H., Morales, S. E., Fineran, P. C. & Brown, C. M. CRISPRDetect: a flexible algorithm to define CRISPR arrays. BMC Genomics 17, 356 (2016)

MacSyFinder
Prodigal
GNU Parallel
HMMER
CRISPRCasFinder group (DEF and profiles)
Banfield group (Cas14)
