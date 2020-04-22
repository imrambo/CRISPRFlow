from Bio import Entrez, SeqIO
import time
import argparse
import sys
import os
#import urllib

"""
Motivation: Batch record download using ENTREZ Efetch.
Author: Ian Rambo
Thirteen... that's a mighty unlucky number... for somebody!
"""
parser = argparse.ArgumentParser()

parser.add_argument('-a', '--acc_file', type=str, dest='acc_file', action='store',
help='File containing a single column of NCBI accessions to fetch')
parser.add_argument('-o', '--output', type=str, dest='output', action='store',
help='output fasta file with fetched sequences')
parser.add_argument('-e', '--email', type=str, dest='email', action='store',
required=True, help='email to use for entrez requests')
parser.add_argument('-b', '--batch', type=int, dest='batch', required=False, nargs='?')
parser.add_argument('-d', '--database', type=str, dest='db', action='store',
default='protein', help='NCBI database to search. Default=protein')
parser.add_argument('-r', '--rettype', type=str, dest='rettype', action='store',
default='fasta', help='Return type for ENTREZ efetch. Default=fasta')
parser.add_argument('-i', '--idtype', type=str, dest='idtype', action='store',
default='acc', help='ENTREZ identifier type.  Default=acc')
parser.add_argument('-w', '--waittime', type=float, dest='waittime', action='store',
default=0.5, help='Time to wait (seconds) between each ENTREZ requests. Default=0.5')
#------------------------------------------------------------------------------
def accession_list(accfile):
    alist = []
    with open(accfile, 'r') as acc:
        for a in acc:
            a = a.rstrip()
            alist.append(a)
    accessions_uniq = list(set(alist))
    accessions_uniq = [d.replace('"', '') for d in accessions_uniq]
    return accessions_uniq
#------------------------------------------------------------------------------
def batch(idlist, size):
    l = len(idlist)
    for start in range(0, l, size):
        yield idlist[start:min(start + size, l)]
#------------------------------------------------------------------------------
def fetch_record(database, id, idtype, rettype, retmax=1, batch=False, waittime=0.5):
    """
    Fetch records with ENTREZ Efetch
    """
    if batch:
        id = ','.join(id)
    else:
        pass

    # try:
    #     handle = Entrez.efetch(db=database,id=id,rettype=rettype,idtype=idtype,retmax=retmax)
    # except urllib.error.HTTPError as err:
    #     if err.code == 502:
    #         print('bad gateway error 502')
    #     else:
    #         raise
    handle = Entrez.efetch(db=database,id=id,rettype=rettype,idtype=idtype,retmax=retmax)
    efetch_out = handle.read()
    handle.close()
    time.sleep(waittime)
    return efetch_out
#------------------------------------------------------------------------------
def main(argv):
    args = parser.parse_args(argv)
    Entrez.email = args.email

    accessions = accession_list(os.path.abspath(args.acc_file))
    print('preparing to fetch entries for %d total accessions...' % len(accessions))
    print('ENTREZ options = rettype:%s, idtype:%s, db:%s' % (args.rettype, args.idtype, args.db))

    with open(os.path.abspath(args.output), 'w') as output:
        if args.batch:
            print('batch efetch mode - %d per request' % args.batch)
            batch_generator = batch(accessions, args.batch)
            for b in batch_generator:
                efetch_out = fetch_record(database=args.db, id=b, idtype=args.idtype, rettype=args.rettype, retmax=args.batch, batch=True, waittime=args.waittime)
                output.write(efetch_out)

        else:
            print('single efetch mode')
            for a in accessions:
                efetch_out = fetch_record(database=args.db, id=a, idtype=args.idtype, rettype=args.rettype, batch=False, waittime=args.waittime)
                output.write(efetch_out)
    print('records written to %s' % os.path.abspath(args.output))

if __name__ == "__main__":
    main(sys.argv[1:])
