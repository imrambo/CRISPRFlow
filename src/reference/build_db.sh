DBDIR=$1
##---Databases---###
RUN mkdir /build/database/TIGRFAM && cd /build/database/TIGRFAM \
   && wget ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.tar.gz \
   && tar -zxf TIGRFAMs_15.0_HMM.tar.gz \
   && rm TIGRFAMs_15.0_HMM.tar.gz

RUN mkdir /build/database/ANTIFAM && cd /build/database/ANTIFAM \
   && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz \
   && tar -zxf Antifam.tar.gz \
   && rm Antifam.tar.gz

RUN mkdir /build/database/PFAM-A && cd /build/database/PFAM-A \
   && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz \
   && gunzip Pfam-A.hmm.gz
