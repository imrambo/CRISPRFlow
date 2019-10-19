from Bio import Entrez, SeqIO
import time
import argparse

"""
Motivation: Supply a BLAST output table and fetch matching sequences from ENTREZ.
Author: Ian Rambo
"""
parser = argparse.ArgumentParser()

parser.add_argument('--blasttbl', type=str, dest='blasttbl', action='store',
help='BLAST+ results file')
parser.add_argument('--outfasta', type=str, dest='outfasta', action='store',
help='output fasta file with fetched sequences')
parser.add_argument('--email', type=str, dest='email', action='store',
help='email to use for entrez requests')
parser.add_argument('--dbtype', type=str, dest='dbtype', action='store',
default='protein', help='database to search')

args = parser.parse_args()
#------------------------------------------------------------------------------
def entrez_efetch(db, id, rtp='fasta', idtype='acc', waittime=0.2):
    handle = Entrez.efetch(db=db,
                           id=id,
                           rettype=rtp,
                          idtype=idtype)
    record = SeqIO.read(handle, rtp)
    handle.close()
    #wait so you don't overload the server if doing multiple requests
    time.sleep(waittime)
    return record
#------------------------------------------------------------------------------
Entrez.email = args.email # Always tell NCBI who you are

accessions = []
with open(args.blasttbl, 'r') as dsrb:
    dsrb.readline()
    for line in dsrb:
        line = line.rstrip()
        line_list = line.split('\t')
        blast_hit = line_list[1]
        acc = blast_hit.split(' ')[0]
        accessions.append(acc)

accessions_uniq = list(set(accessions))
accessions_uniq = [d.replace('"', '') for d in accessions_uniq]
print('preparing to fetch fasta sequences for %d accessions...' % len(accessions_uniq))

counter = 0
with open(args.outfasta, 'w') as outfa:
    for i in accessions_uniq:
        seq_record = entrez_efetch(db=args.dbtype,id=i)
        outfa.write(seq_record.format("fasta"))
        print('%d of %d sequences fetched' % (count, accessions_uniq))
        counter += 1
print('%d sequences written to %s' % (counter, args.outfasta))
