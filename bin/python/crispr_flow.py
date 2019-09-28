"""
Motivation:

Author: Ian Rambo
Contact: ian.rambo@utexas.edu, imrambo@lbl.gov

Thirteen... that's a mighty unlucky number... for somebody!
"""
version = 'v0.0.1'

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
from hmmbo import *
from prodigal import *
from shell_tools import *

# FORMAT = '%(asctime)-15s %(clientip)s %(user)-8s %(message)s'
# logging.basicConfig(format=FORMAT)
# logging.get_logger()
#==============================================================================

#------------------------------------------------------------------------------
def get_basename(file_path):
    basename = os.path.basename(file_path)
    #Remove two extensions, e.g. foo.tar.gz becomes foo
    if re.match(r'^.*?\.[a-z]+\.[a-z]+$', basename):
        basename = re.findall(r'^(.*?)\.[a-z]+\.[a-z]+$', basename)[0]
    else:
        basename = os.path.splitext(basename)[0]
    return basename
#------------------------------------------------------------------------------
def make_seqdict(fasta, prodigal=False, gz=False, format='fasta'):
    """
    Create a SeqIO sequence dictionary. If prodigal is True,
    add gene start and stop coordinates from FASTA header to the dictionary.
    """
    if gz:
        seq_handle = gzip.open(fasta, 'rb')
    else:
        seq_handle = open(fasta, 'r')

    seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, format))
    # if prodigal:
    #     for key in seq_dict.keys():
    #         seq_dict[key].gene_start = int(seq_dict[key].description.split('#')[1])
    #         seq_dict[key].gene_stop = int(seq_dict[key].description.split('#')[2])
    #         #seq_dict[key]['gene_start'] = int(seq_dict[key].description.split('#')[1])
    #         #seq_dict[key]['gene_stop'] = int(seq_dict[key].description.split('#')[2])
    #     print(seq_dict)
    #     return seq_dict
    # else:
    return seq_dict
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
help='Directory for CRISPRDetect.pl', default='/build/CRISPRDetect_2.4/CRISPRDetect_2.4')
parser.add_argument('--database', type=str, dest='database', action='store',
nargs='?', help='Comma-separated list of HMM database paths to use for CRISPR-proximity gene hmmsearch. You can also specify paths to certain profiles, e.g. /path/to/K00001.hmm')
parser.add_argument('--ccs_typing', type=str, dest='ccs_typing', action='store', default='sub-typing',
help='Level of CRISPR-Cas system typing specificity. Choose from: [general, typing, sub-typing]')
parser.add_argument('--profile_suffix', type=str, dest='profile_suffix', action='store', default='.hmm',
help='suffix for HMM gene profiles. Default is ".hmm"')
parser.add_argument('--prefix', type=str, dest='prefix', action='store', nargs='?',
help='optional prefix for output files. Uses the input nt fasta basename if not supplied.')
opts = parser.parse_args()

if opts.jobs < 1:
    logging.warning('--jobs must be >= 1')
#==============================================================================
#Create the output directories
output_paths = {p: os.path.join(opts.out_root, p) for p in ['CRISPRDetect','Prodigal','MacSyFinder','HMMER']}

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

nt_fasta = opts.fasta_file
#Get the file basename to name output files
nt_fasta_basename = get_basename(nt_fasta)

prefix = opts.prefix
if not prefix:
    prefix = nt_fasta_basename

#If the nucleotide fasta input is gzipped, gunzip it
if is_gzipped(nt_fasta):
    print('%s is gzip compressed, gunzip file...' % nt_fasta)
    subprocess.run(['gunzip', nt_fasta], shell=False)
    nt_fasta = os.path.splitext(opts.fasta_file)[0]

else:
    pass
#==============================================================================
###---CRISPRDetect---###
crispr_detect_out = os.path.join(output_paths['CRISPRDetect'], nt_fasta_basename + '_CRISPRDetect')
crispr_detect_log = os.path.join(output_paths['CRISPRDetect'], nt_fasta_basename + '_CRISPRDetect.log')
#CRISPRDetect options
crispr_detect_opts = {'-f':nt_fasta,
'-o':crispr_detect_out, '-T':opts.threads, '-minimum_repeat_length':20,
'-array_quality_score_cutoff':3, '-tmp_dir':opts.tmp_dir,
 '-logfile':crispr_detect_log}

###---Run CRISPRDetect---###
if not opts.crispr_gff:
    #Run CRISPRDetect
    crispr_detect_exec = os.path.join(opts.CRISPRDetectDir, 'CRISPRDetect.pl')

    crispr_detect_optlist = exec_cmd_generate(crispr_detect_exec, crispr_detect_opts)
    #subprocess.run([crispr_detect_exec, crispr_detect_optstring], shell=False)
    subprocess.run(crispr_detect_optlist, shell=False)

    crispr_gff = crispr_detect_out + '.gff'
else:
    crispr_gff = opts.crispr_gff
#Convert the GFF to a pandas data frame, selecting full CRISPR arrays coords
crispr_gff_df = gff_to_pddf(gff = crispr_gff, ftype = 'repeat_region')
#==============================================================================
# ###---Prodigal---###
prodigal_outfmt = 'gff'
prodigal_out = os.path.join(output_paths['Prodigal'], prefix + '_prodigal.%s' % prodigal_outfmt)
prodigal_aa = os.path.join(output_paths['Prodigal'], prefix + '_prodigal.faa')
prodigal_nt = os.path.join(output_paths['Prodigal'], prefix + '_prodigal.fna')

prodigal_opts = {'-o':prodigal_out, '-a':prodigal_aa, '-d':prodigal_nt}

if opts.prodigal_gff and opts.prodigal_amino:
    prodigal_gff_df = gff_to_pddf(opts.prodigal_gff)
    prodigal_faa_dict = make_seqdict(opts.prodigal_amino, prodigal=True)
# elif opts.prodigal_gff and not opts.prodigal_amino or not opts.prodigal_gff and opts.prodigal_amino:
#     torun = [k for k in [opts.prodigal_gff, opts.prodigal_amino] if not k]
else:
    #Generate and run the Prodigal command
    prodigal_command = prodigal_command_generate(ntfasta=nt_fasta, optdict=prodigal_opts,
    outfmt=prodigal_outfmt, prodigal='prodigal')

    print('Prodigal will be run in %s mode' % prodigal_command[1]['-p'])
    print(prodigal_command[0])
    subprocess.run(prodigal_command[0], shell=False)

    prodigal_faa_dict = make_seqdict(prodigal_command[1]['-a'], prodigal=True)

else:

# #==============================================================================
#Fetch the CRISPR-neighboring genes
neighbor_aa_fasta = os.path.join(output_paths['Prodigal'], nt_fasta_basename + '_CRISPR-neighbor-genes.faa')
#neighbor_nt_fasta = os.path.join(prodigal_outdir, nt_fasta_basename + '_CRISPR-neighbor-genes.fna')
fetch_gene_clusters(gff_anchor=crispr_gff_df, gene_seq_dict=prodigal_faa_dict, out_fasta=neighbor_aa_fasta, winsize=opts.window_extent, prodigal=True)
#==============================================================================
# #Subtype the groups of CRISPR-neighboring genes
# macsyfinder_opts = {'--sequence_db':neighbor_aa_fasta, '--db_type':'ordered_replicon',
# '--e-value-search':1e-6, '--i-evalue-select':1e-6, '--coverage_profile':0.5,
# '--def':'%s/data/definitions/%s' % (build_root, opts.ccs_typing), '--out-dir':output_paths['MacSyFinder'],
# '--res-search-suffix':'hmmout', '--res-extract-suffix':'res_hmm_extract',
# '--profile-suffix':'hmm', '--profile-dir':'%s/data/profiles/CAS' % build_root,
# '--worker':opts.threads, '-vv':''}
#
# if opts.joblog_dir:
#     macsyfinder_opts['--log'] = opts.joblog_dir
# else:
#     macsyfinder_opts['--log'] = output_paths['MacSyFinder']
#
# macsyfinder_command = 'macsyfinder %s %s' % (optstring_join(macsyfinder_opts), 'all')
# print(macsyfinder_command)
# #subprocess.run([macsyfinder_command])
# #==============================================================================
# ##---HMMSEARCH---###
# hmmsearch_opts = {'--domE':10, '-E':10, '--incE':1e-6, '--incdomE':1e-6, '--seed':42, '--cpu':1}
#
# seqdb = neighbor_aa_fasta
#
# #Generate the command strings
# if ',' in opts.database:
#     database_list = ','.split(opts.database)
#     hmmsearch_commands = hmmbo.hmmsearch_command_generator(db_list=database_list, hmmsearch_optdict=hmmsearch_opts, parallel_optdict=parallel_optdict, jobs=1)
#
# else:
#     hmmsearch_commands = hmmbo.hmmsearch_command_generator(db_list=[opts.database], hmmsearch_optdict=hmmsearch_opts, parallel_optdict=parallel_optdict, jobs=1)
