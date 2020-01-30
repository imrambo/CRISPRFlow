DBDIR=$1
KEYFILE=$2
hmmfetch --index PfamA.gz
hmmfetch -O -f PfamA.hmm $KEYFILE

test -d $DBDIR || mkdir $DBDIR
cd $DBDIR
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0/pfamseq.gz &
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0/uniprot.ssi &

wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz && \
    tar -zxvf Antifam.tar.gz
