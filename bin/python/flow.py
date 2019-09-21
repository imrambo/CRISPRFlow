from Bio import SeqIO
import argparse
import logging
import os
import numpy as np
import pandas as pd
import subprocess
"""
Search for CRISPR-Cas systems using CRISPR arrays as an anchor.
"""
parser = argparse.ArgumentParser()
parser.add_argument('--scaffold', type=str, dest='scaffold_file', action='store',
help='input nucleotide FASTA file. Required.')
parser.add_argument('--out_root', type=str, dest='out_root', action='store',
help='root output directory')
parser.add_argument('--tmp_dir', type=str, dest='tmp_dir', action='store',
default='/tmp',
help='temporary file directory. If not specified, uses /tmp')
parser.add_argument('--crispr_gff', type=str, dest='crispr_gff', action='store', nargs='?',
help='CRISPRDetect GFF file. If supplied, CRISPRDetect will not be run.')
parser.add_argument('--prodigal_amino', type=str, dest='prodigal_amino', action='store', nargs='?',
help='Prodigal amino acid file. Optional. If supplied, Prodigal will not be run.')
parser.add_argument('--window_extent', type=int, dest='window_extent', action='store', default=10000,
help='Number of bp extension to include in a cluster. Default is 10kb.')
parser.add_argument('--logfile_dir', type=str, dest='logfile_dir', action='store',
help='Directory for log files')
parser.add_argument('--crispr_detect_dir', type=str, dest='CRISPRDetectDir', action='store',
help='Directory for CRISPRDetect.pl')
# parser.add_argument('--db_set', type=str, dest='db_set', action='store',
# help='Comma-separated list of HMM databases to search, e.g. KO,PFAM,TIGRFAM. Choose from: {KO, PFAM, TIGRFAM}')
#------------------------------------------------------------------------------
opts = parser.parse_args()


# FORMAT = '%(asctime)-15s %(clientip)s %(user)-8s %(message)s'
# logging.basicConfig(format=FORMAT)
# logging.get_logger()
#==============================================================================
def optstring_join(optdict):
    """
    Join a dictionary of command line options into a single string.
    """
    optstring = ' '.join([str(param) + ' ' + str(val) for param, val in optdict.items()])
    return optstring
#------------------------------------------------------------------------------
def gff_to_pddf(gff, ftype=''):
    """Read in a GFF file as a Pandas data frame. Specify ftype to select
    rows pertaining to a specific feature type."""
    gff_cols = ['source', 'ftype', 'start', 'end', 'score', 'strand',
    'phase', 'attributes']
    if os.path.exists(gff) and os.stat(gff).st_size != 0:
        gff_df = pd.read_csv(crispr_gff, sep='\s+', names=gff_cols, comment='#')

        if ftype:
            gff_df[gff_df['ftype'] == ftype]
            return gff_df
        else:
            return gff_df
#------------------------------------------------------------------------------
def fetch_gene_clusters(gff_df, prodigal_seq_dict, out_fasta, winsize):
    with open(out_fasta, 'w') as fa:
        neighbor_genes = []
        #Loop through records and fetch neighboring genes
        for index, row in gff_df.iterrows():
            seq_objs = [prodigal_seq_dict[key] for key, value in prodigal_seq_dict.items() if row['source'] + '_' in key and prodigal_seq_dict[key]['gene_start'] >= row['start'] - opts.window_extent  and gene_stop <= row['end'] + opts.window_extent]
            neighbor_genes.extend(seq_objs)
        SeqIO.write(neighbor_genes, fa, 'fasta')
#------------------------------------------------------------------------------
def make_seqdict(fasta, prodigal=False):
    """
    Create a SeqIO sequence dictionary from a fasta file. If prodigal is True,
    add gene start and stop coordinates to the dictionary.
    """
    seq_handle = open(fasta, 'r')
    seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle), 'fasta')
    if prodigal:
        for key, value in seq_dict:
            seq_dict[key]['gene_start'] = int(seq_dict[key]['description'].split(' # ')[1])
            seq_dict[key]['gene_stop'] = int(seq_dict[key]['description'].split(' # ')[2])
        return(seq_dict)
    else:
        return(seq_dict)
#==============================================================================
output_paths = {p: os.path.join(opts.out_root, p) for p in ['CRISPRDetect','Prodigal','MacSyFinder','HMMER']}
#for key, value in output_paths:
#    if not os.path.exists(value):
#        os.makedirs(value)
#==============================================================================
#Prodigal options
prodigal_out = os.path.join(output_paths['Prodigal'], os.path.basename(opts.scaffold_file) + '_prodigal.gff')
prodigal_aa = os.path.join(output_paths['Prodigal'], os.path.basename(opts.scaffold_file) + '_prodigal.faa')
prodigal_nt = os.path.join(output_paths['Prodigal'], os.path.basename(opts.scaffold_file) + '_prodigal.fna')

prodigal_opts = {'-i':opts.scaffold_file, '-p':'single', '-o':prodigal_out,
'-a':prodigal_aa, '-d':prodigal_nt, '-f':'gff', '-q':''}
#==============================================================================
###---HMMER suite---###
#Hmmsearch options
hmmsearch_opts = {'--domE':10, '-E':10, '--incE':1e-6, '--incdomE':1e-6, '--seed':42}
hmmsearch_parallel_opts = {'-jobs':5}
hmmsearch_optstring = optstring_join(hmmsearch_opts)
#Jackhmmer options
jackhmmer_opts = {'--mx':'BLOSUM62', '-E':1e-20, '--domE':1e-20, '-N':20,
'--incE':1e-20, '--incdomE':1e-20, '--F1':0.01, '--seed':42, '--cpu':4,
'--noali':''}

###---CRISPRDetect---###
crispr_detect_out = os.path.join(output_paths['CRISPRDetect'], os.path.basename(opts.scaffold_file) + '_CRISPRDetect')
crispr_detect_log = os.path.join(output_paths['CRISPRDetect'], os.path.basename(opts.scaffold_file) + '_CRISPRDetect.log')
#CRISPRDetect options
crispr_detect_opts = {'-f':opts.scaffold_file,
'-o':crispr_detect_out, '-T':2, '-minimum_repeat_length':23,
'-array_quality_score_cutoff':3, '-tmp_dir':opts.tmp_dir,
 '-logfile':crispr_detect_log
 }
crispr_detect_optstring = optstring_join(crispr_detect_opts)

###---CRISPRCasFinder---###
 ccfinder_opts = {'-log':'', '-copyCSS':'', '-repeats':'', '-DBcrispr':'',
'-DIRrepeat':'', '-cas':'', '-ccvRep':'', '-vicinity':10000, '-cluster':20000,
'-minDR':23, '-maxDR':72, '-minSP':8, '-maxSP':64, '-cpuM':1,
'-definition':'SubTyping', '-getSummaryCasfinder':'', '-betterDetectTrunc':'',
'-soFile':ccfinder_sofile, '-keep':''}
ccfinder_optstring = optstring_join(ccfinder_opts)
#==============================================================================
#Read in the nucleotide scaffolds
# scaffold_handle = open(opts.scaffold_file, 'r')
# scaffold_obj = SeqIO.to_dict(SeqIO.parse(scaffold_handle), 'fasta')
#==============================================================================
###---CRISPRDetect---###
if not opts.crispr_gff:
    #Run CRISPRDetect
    crispr_detect_command = os.path.join(opts.CRISPRDetectDir, 'CRISPRDetect.pl %s' % crispr_detect_optstring)
    subprocess.run([crispr_detect_command])
    crispr_gff = crispr_detect_out + '.gff'
else:
    crispr_gff = opts.crispr_gff

crispr_gff_df = gff_to_pddf(gff = crispr_gff, ftype = 'repeat_region')
#==============================================================================
###---Prodigal---###
if not opts.prodigal_aa:
    #Run Prodigal
    prodigal_optstring = optstring_join(prodigal_opts)
    prodigal_command = 'prodigal %s' % prodigal_optstring

    subprocess.run([prodigal_command])

    prodigal_faa_dict = make_seqdict(prodigal_opts['-a'], prodigal=True)

else:
    prodigal_faa_dict = make_seqdict(opts.prodigal_aa, prodigal=True)
#==============================================================================
#Fetch the CRISPR-neighboring genes
neighbor_aa_fasta = os.path.join(output_paths['Prodigal'], os.path.basename(opts.scaffold_file) + '_CRISPR-neighbor-genes.faa')
#neighbor_nt_fasta = os.path.join(prodigal_outdir, os.path.basename(opts.scaffold_file) + '_CRISPR-neighbor-genes.fna')
fetch_gene_clusters(gff_df=crispr_gff_df, prodigal_seq_dict=prodigal_faa_dict, out_fasta=neighbor_aa_fasta, winsize=opts.window_extent)
#==============================================================================
#Subtype the groups of CRISPR-neighboring genes
macsyfinder_opts = {'--sequence_db':neighbor_aa_fasta, '--db_type':'ordered_replicon',
'--e-value-search':1e-6, '--i-evalue-select':1e-4, '--coverage_profile':0.5,
'--def':'../DEF', '--out-dir':output_paths['MacSyFinder'],
'--res-search-suffix':'hmmout', '--res-extract-suffix':'res_hmm_extract',
'--profile-suffix':'search_hmm.out', '--profile-dir':'../profiles/CAS',
'--worker':1, '--verbosity':'-vv'}
if opts.logfile_dir:
    macsyfinder_opts['--log'] = opts.logfile_dir
else:
    macsyfinder_opts['--log'] = output_paths['MacSyFinder']

macsyfinder_optstring = optstring_join(macsyfinder_opts)
macsyfinder_command = 'macsyfinder %s' % macsyfinder_optstring
subprocess.run([macsyfinder_command])
#==============================================================================
# #Run HMMsearches
# hmm_dir = './hmm'
# if not opts.db_set:
#     for database in os.walk(hmm_dir, maxdepth = 1):
#         if parallel:
#             hmmsearch_command = 'find %s -type f -name "*.hmm" | parallel %s hmmsearch --domtblout %s/%s_{/.}.domtbl %s %s %s_{./}' % (database, hmmsearch_parallel_opts, hmmsearch_optstring, neighbor_aa_fasta)
#             subprocess.run([hmmsearch_command])
#         else:
#             hmmsearch_command = 'hmmsearch --domtblout %s/%s_{/.}.domtbl %s %s %s_{./}'
#         #Gather the hits and select the best hits
#         hmmer_best_hits = os.path.join(hmmer_outdir, 'hmmsearch_best_hits.domtbl')
#         hmmsearch_hits_command = 'grep -h ID %s/*.domtbl > %s' % (hmmer_outdir, hmmer_best_hits)
#
# subprocess.run([hmmsearch_hits_command])
#
# dbh_df = domtbl_besthits(hmmer_best_hits)
#=============================================================================

#Create unified table - CRISPR and genes

# def shorten_header_crisprDetect():
#     "if header length+maxsize of crispr > 256"
#Convert to GFF3 format and combine with CRISPR gff df
#==========================================================================
#Parse HMMsearch results

#Nucleotide analyses

#Parse scaffolds - include non-coding regions


#create new features in GFF
