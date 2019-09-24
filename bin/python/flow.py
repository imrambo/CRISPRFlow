"""
Motivation:

Author: Ian Rambo
Contact: ian.rambo@utexas.edu, imrambo@lbl.gov

Thirteen... that's a mighty unlucky number... for somebody!
"""
from Bio import SeqIO
import argparse
import logging
import os
import numpy as np
import pandas as pd
import subprocess
import magic
import gzip
import re
import datetime
import hmmbo
import prodigal

# FORMAT = '%(asctime)-15s %(clientip)s %(user)-8s %(message)s'
# logging.basicConfig(format=FORMAT)
# logging.get_logger()

#REFERENCES
#CRISPRDetect
#MacSyFinder
#Prodigal
#GNU Parallel
#HMMER
#CRISPRCasFinder group (DEF and profiles)
#Banfield group (Cas14)
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
def get_basename(file_path):
    basename = os.path.basename(opts.fasta_file)
    #Remove two extensions, e.g. foo.tar.gz becomes foo
    if re.match(r'^.*?\.[a-z]+\.[a-z]+$', basename):
        basename = re.findall(r'^(.*?)\.[a-z]+\.[a-z]+$', basename)[0]
    else:
        basename = os.path.splitext(basename)[0]
#------------------------------------------------------------------------------
def is_gzipped(file_path):
    """
    Test if a file is gzipped.
    """
    is_gzip = False
    if magic.from_file(file_path).startswith('gzip compressed data'):
        is_gzip = True
        return is_gzip
    else:
        return is_gzip
#------------------------------------------------------------------------------
def make_seqdict(fasta, prodigal=False, gz=False):
    """
    Create a SeqIO sequence dictionary from a fasta file. If prodigal is True,
    add gene start and stop coordinates to the dictionary.
    """
    if gz:
        seq_handle = gzip.open(fasta, 'rb')
    else:
        seq_handle = open(fasta, 'r')

    seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle), 'fasta')
    if prodigal:
        for key in seq_dict.keys():
            seq_dict[key]['gene_start'] = int(seq_dict[key]['description'].split(' # ')[1])
            seq_dict[key]['gene_stop'] = int(seq_dict[key]['description'].split(' # ')[2])
        return seq_dict
    else:
        return seq_dict
#------------------------------------------------------------------------------
def fetch_gene_clusters(gff_anchor, gene_seq_dict, out_fasta, winsize, gff_gene=None, anchor_col='source'):
    """
    Get Prodigal genes within a certain distance from a genomic feature.
    The gff_df must only contain 'anchor' genomic features, i.e. CRISPR array.
    gene_seq_dict is a SeqIO Sequence Dict.
    """
    with open(out_fasta, 'w') as fa:
        neighbor_genes = []
        #Loop through records and fetch neighboring genes
        for index, row in gff_anchor.iterrows():
            #if gff_gene:
            seq_objs = [gene_seq_dict[key] for key in gene_seq_dict.keys() if row[anchor_col] + '_' in key and gene_seq_dict[key]['gene_start'] >= row['start'] - winsize  and gene_seq_dict[key]['gene_stop'] <= row['end'] + winsize]
            neighbor_genes.extend(seq_objs)
        SeqIO.write(neighbor_genes, fa, 'fasta')
#------------------------------------------------------------------------------
# def command_generator(seqdb, profile, prefix=None, optstring, parallel=False, progbar=False, program='hmmsearch', psuffix='.hmm', outdir='/tmp', jobs=1, joblog='joblog'):
#     """
#     Generate a list of bash command strings.
#     Use parallel=True to run parallel jobs using GNU parallel.
#     Pass an ordered dictionary if parameter order matters.
#     """
#     #  BREAK THIS UP INTO FUNCTIONS
#     # if not prefix:
#     #     prefix = get_basename(seqdb)
#
#     commands = []
#
#     if parallel and jobs == 1:
#         parallel = False
#     if not parallel:
#         if os.path.isdir(profile):
#             profile_glob = os.path.join(profile, '*%s' % psuffix)
#             for i in list(glob.glob(profile_glob, recursive = False)):
#                 cmd = '%s %s %s %s' % (program, optstring, i, seqdb)
#                 commands.append(cmd)
#
#         elif os.path.isfile(profile):
#             #Run program for a single profile
#             #optstring = optstring_join(optdict)
#             cmd = '%s %s' % (program, optstring)
#             commands.append(cmd)
#         else:
#             pass
#     elif parallel and jobs > 1:
#         #Build the parallel string
#         parallel_optdict = {'--jobs':jobs}
#         if progbar:
#             parallel_optdict['--bar'] = ''
#         parallel_optdict['--joblog'] = joblog
#         parallel_optstring = opstring_join(parallel_optdict)
#
#         #optstring = optstring_join(optdict)
#         if os.path.isdir(profile):
#             cmd = 'find %s -type f -name "*%s" | parallel %s %s %s {} %s' % (profile, psuffix, parallel_optstring, program, optstring, seqdb)
#             commands.append(cmd)
#         """
#         If a list of paths is provided, break it up into files and directories
#         and generate strings.
#         """
#         if isinstance(profile, list):
#             pfiles = [os.path.isfile(p) for p in profile]
#             pdirs = [os.path.isdir(p) for p in profile]
#             if len(pfiles) + len(pdirs) == len(profile):
#                     if pdirs:
#                         for pdir in pdirs:
#                             cmd = 'find %s -type f -name "*%s" | parallel %s %s %s {} %s' % (pdir, psuffix, parallel_optstring, program, optstring, seqdb)
#                             commands.append(cmd)
#                     if pfiles and len(pfiles) > 1:
#                         cmd = 'parallel %s %s %s {} %s ::: %s' % (parallel_optstring, program, optstring, seqdb, ' '.join(pfiles))
#                         commands.append(cmd)
#     else:
#         pass
#     return commands
#==============================================================================
build_root = '../..'
#==============================================================================
parser = argparse.ArgumentParser()

parser.add_argument('--fasta', type=str, dest='fasta_file', action='store',
help='input nucleotide FASTA file. Required.')
parser.add_argument('--out_root', type=str, dest='out_root', action='store',
help='root output directory. Required.')
parser.add_argument('--tmp_dir', type=str, dest='tmp_dir', action='store',
default='/tmp',
help='temporary file directory.')
parser.add_argument('--threads', type=int, dest='threads', action='store', default=1,
help='number of threads for each command')
parser.add_argument('--jobs', type=int, dest='jobs', action='store', default=1,
help='number of jobs for each step')
parser.add_argument('--crispr_gff', type=str, dest='crispr_gff', action='store', nargs='?',
help='CRISPRDetect GFF file. If supplied, CRISPRDetect will not be run.')
parser.add_argument('--prodigal_amino', type=str, dest='prodigal_amino', action='store', nargs='?',
help='Prodigal amino acid file. Optional. If supplied, Prodigal will not be run.')
parser.add_argument('--window_extent', type=int, dest='window_extent', action='store', default=10000,
help='Number of bp extension to include in a cluster. Default = 10kb.')
parser.add_argument('--joblog_dir', type=str, dest='joblog_dir', action='store',
nargs = '?', help='Directory for log files.')
parser.add_argument('--crispr_detect_dir', type=str, dest='CRISPRDetectDir', action='store',
help='Directory for CRISPRDetect.pl', default='/build/CRISPRDetect_2.4')
parser.add_argument('--database', type=str, dest='database', action='store',
nargs='?', help='Comma-separated list of HMM database paths to use for CRISPR-proximity gene hmmsearch. You can also specify paths to certain profiles, e.g. /path/to/K00001.hmm')
parser.add_argument('--ccs_typing', type=str, dest='ccs_typing', action='store', default='sub-typing',
help='Level of CRISPR-Cas system typing specificity. Choose from: [general, typing, sub-typing]')
parser.add_argument('--profile_suffix', type=str, dest='profile_suffix', action='store', default='.hmm',
help='suffix for HMM gene profiles. Default is ".hmm"')

opts = parser.parse_args()

if opts.jobs < 1:
    logging.warning('--jobs must be >= 1')
#==============================================================================
#Create the output directories
output_paths = {p: os.path.join(opts.out_root, p) for p in ['CRISPRDetect','Prodigal','MacSyFinder','HMMER']}
print(output_paths)
for key, value in output_paths.items():
    if not os.path.exists(value):
        os.makedirs(value)
    else:
        pass

#Create the joblog directory if not specified
if not os.path.exists(opts.joblog_dir):
    os.makedirs(opts.joblog_dir)
else:
    pass
#==============================================================================
#Options for GNU parallel
parallel_optdict = {'--jobs':opts.jobs, '--bar':''}

#Get the file basename to name output files
fasta_basename = get_basename(opts.fasta_file)
#==============================================================================
###---CRISPRDetect---###
crispr_detect_out = os.path.join(output_paths['CRISPRDetect'], fasta_basename + '_CRISPRDetect')
crispr_detect_log = os.path.join(output_paths['CRISPRDetect'], fasta_basename + '_CRISPRDetect.log')
#CRISPRDetect options
crispr_detect_opts = {'-f':opts.fasta_file,
'-o':crispr_detect_out, '-T':opts.threads, '-minimum_repeat_length':23,
'-array_quality_score_cutoff':3, '-tmp_dir':opts.tmp_dir,
 '-logfile':crispr_detect_log}

crispr_detect_optstring = optstring_join(crispr_detect_opts)

###---Run CRISPRDetect---###
if not opts.crispr_gff:
    #Run CRISPRDetect
    crispr_detect_command = os.path.join(opts.CRISPRDetectDir, 'CRISPRDetect.pl %s' % crispr_detect_optstring)
    subprocess.run([crispr_detect_command])
    crispr_gff = crispr_detect_out + '.gff'
else:
    crispr_gff = opts.crispr_gff
#Convert the GFF to a pandas data frame, selecting full CRISPR arrays coords
crispr_gff_df = gff_to_pddf(gff = crispr_gff, ftype = 'repeat_region')
#==============================================================================
###---Prodigal---###
if not opts.prodigal_aa:
    #Generate and run the Prodigal command
    prodigal_command = prodigal_command_generate(fasta=opts.fasta_file, outdir=output_paths['Prodigal'], prefix=fasta_basename)

    subprocess.run([prodigal_command[0]])

    prodigal_faa_dict = make_seqdict(prodigal_command[1]['-a'], prodigal=True)

else:
    prodigal_faa_dict = make_seqdict(opts.prodigal_aa, prodigal=True)
#==============================================================================
#Fetch the CRISPR-neighboring genes
neighbor_aa_fasta = os.path.join(output_paths['Prodigal'], fasta_basename + '_CRISPR-neighbor-genes.faa')
#neighbor_nt_fasta = os.path.join(prodigal_outdir, fasta_basename + '_CRISPR-neighbor-genes.fna')
fetch_gene_clusters(gff_df=crispr_gff_df, gene_seq_dict=prodigal_faa_dict, out_fasta=neighbor_aa_fasta, winsize=opts.window_extent)
#==============================================================================
#Subtype the groups of CRISPR-neighboring genes
macsyfinder_opts = {'--sequence_db':neighbor_aa_fasta, '--db_type':'ordered_replicon',
'--e-value-search':1e-6, '--i-evalue-select':1e-6, '--coverage_profile':0.5,
'--def':'%s/data/definitions/%s' % (build_root, opts.ccs_typing), '--out-dir':output_paths['MacSyFinder'],
'--res-search-suffix':'hmmout', '--res-extract-suffix':'res_hmm_extract',
'--profile-suffix':'hmm', '--profile-dir':'%s/data/profiles/CAS' % build_root,
'--worker':opts.threads, '-vv':''}

if opts.joblog_dir:
    macsyfinder_opts['--log'] = opts.joblog_dir
else:
    macsyfinder_opts['--log'] = output_paths['MacSyFinder']

macsyfinder_command = 'macsyfinder %s %s' % (optstring_join(macsyfinder_opts), 'all')
print(macsyfinder_command)
#subprocess.run([macsyfinder_command])
#==============================================================================
##---HMMSEARCH---###
hmmsearch_opts = {'--domE':10, '-E':10, '--incE':1e-6, '--incdomE':1e-6, '--seed':42, '--cpu':1}

seqdb = neighbor_aa_fasta

#Generate the command strings
if ',' in opts.database:
    database_list = ','.split(opts.database)
    hmmsearch_commands = hmmbo.hmmsearch_command_generator(db_list=database_list, hmmsearch_optdict=hmmsearch_opts, parallel_optdict=parallel_optdict, jobs=1)

else:
    hmmsearch_commands = hmmbo.hmmsearch_command_generator(db_list=[opts.database], hmmsearch_optdict=hmmsearch_opts, parallel_optdict=parallel_optdict, jobs=1)
