"""
Motivation: CRISPRCas discovery and subtyping

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
from datetime import datetime
from hmmbo import *
from prodigal import *
from shell_tools import *
from gff3 import gff3_to_pddf
from gene_clusters import *
#==============================================================================
# def scons_command(targets, sources, command, env='env'):
#     """Create a Command builder for SCons"""
#     COMMAND = '%s.Command(["%s"], ["%s"], "%s")\n' % (env, '","'.join(targets), '",'.join(sources), command)
#     return COMMAND
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
# parser.add_argument('--crispr_gff', type=str, dest='crispr_gff', action='store', nargs='?',
# help='CRISPRDetect GFF file. If supplied, CRISPRDetect will not be run.')
# parser.add_argument('--prodigal_amino', type=str, dest='prodigal_amino', action='store', nargs='?',
# help='Prodigal amino acid file. Optional. If supplied, Prodigal will not be run.')
parser.add_argument('--window_extent', type=int, dest='window_extent', action='store', default=10000,
help='Number of bp extension to include in a cluster. Default = 10kb.')
parser.add_argument('--joblog', type=str, dest='joblog', action='store',
nargs = '?', help='Path to logging joblog.')
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

#Set up logger
logging_format = '%(levelname)s %(asctime)s - $(message)s'
if opts.joblog:
    #Create the joblog directory if not specified
    if not os.path.exists(os.path.dirname(opts.joblog)):
        os.makedirs(os.path.dirname(opts.joblog))
    else:
        pass
    logging.basicConfig(filename = opts.joblog, level = logging.DEBUG, format = logging_format)
else:
    logging.basicConfig(filename = os.path.join(opts.tmp_dir, 'CRISPRFlow.%s.joblog' % str(datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))),
    level = logging.DEBUG, format = logging_format)

logger = logging.getLogger()

if opts.jobs < 1:
    logger.warning('--jobs must be >= 1')
#==============================================================================
#Create the output directories
output_paths = {p: os.path.join(opts.out_root, p) for p in ['CRISPRDetect','Prodigal','MacSyFinder','HMMER']}

for key, value in output_paths.items():
    if not os.path.exists(value):
        os.makedirs(value)
    else:
        pass

# #Output SConstruct to rebuild the analysis
# sconstruct_handle = os.path.join(output_paths['SCons'], 'SConstruct')
# SConstruct = open(sconstruct_handle, 'w')
# SConstruct.write('env = Environment()\n')
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
gzip = False
if is_gzipped(nt_fasta):
    gzip = True
    print('%s is gzip compressed, gunzip file...' % nt_fasta)
    subprocess.run(['gunzip', nt_fasta], shell=False)
    nt_fasta = os.path.splitext(opts.fasta_file)[0]

else:
    pass
#==============================================================================
###---CRISPRDetect---###
crispr_detect_out = os.path.join(output_paths['CRISPRDetect'], nt_fasta_basename + '_CRISPRDetect')
#crispr_detect_log = os.path.join(output_paths['CRISPRDetect'], nt_fasta_basename + '_CRISPRDetect.log')
#CRISPRDetect options
crispr_detect_optdict = {'-f':nt_fasta,
'-o':crispr_detect_out, '-T':opts.threads, '-minimum_repeat_length':20,
'-array_quality_score_cutoff':3, '-tmp_dir':opts.tmp_dir}
 #'-logfile':crispr_detect_log}

###---Run CRISPRDetect---###
logger.debug('Run CRISPRDetect')
crispr_detect_exec = os.path.join(opts.CRISPRDetectDir, 'CRISPRDetect.pl')

crispr_detect_optlist = exec_cmd_generate(crispr_detect_exec, crispr_detect_optdict)

# subprocess.run(crispr_detect_optlist, shell=False)
# crispr_detect_gff = crispr_detect_out + '.gff'

if os.path.exists(crispr_detect_gff) and os.stat(crispr_detect_gff).st_size != 0:
    #Convert the GFF to a pandas data frame, selecting full CRISPR arrays coords
    crispr_gff_df = gff3_to_pddf(gff = crispr_detect_gff, ftype = 'repeat_region', index_col=False)
else:
    logger.error('CRISPRDetect GFF file %s not found' % crispr_detect_gff)

# #Write the CRISPRDetect command to SConstruct Command builder
# CRISPR_SOURCES = [nt_fasta]
# CRISPR_TARGETS = [os.path.abspath(os.path.join(root, filename)) for root, dirnames, filenames in os.walk(output_paths['CRISPRDetect']) for filename in filenames]
# #CRISPR_TARGETS.insert(0, CRISPR_TARGETS.pop(CRISPR_TARGETS.index(crispr_detect_log)))
# CRISPR_OPTS = crispr_detect_optdict
# #CRISPR_OPTS['-logfile'] = '${TARGETS}[0]'
# CRISPR_OPTS['-f'] = '$SOURCE'
# CRISPR_COMMAND = ' '.join(exec_cmd_generate(crispr_detect_exec, CRISPR_OPTS))
#
# CRISPR_CMDBLD = scons_command(targets = CRISPR_TARGETS, sources = CRISPR_SOURCES, command = CRISPR_COMMAND)
# SConstruct.write(CRISPR_CMDBLD + '\n')
#==============================================================================
# ###---Prodigal---###
prodigal_outfmt = 'gff'
prodigal_out = os.path.join(output_paths['Prodigal'], prefix + '_prodigal.%s' % prodigal_outfmt)
prodigal_aa = os.path.join(output_paths['Prodigal'], prefix + '_prodigal.faa')
prodigal_nt = os.path.join(output_paths['Prodigal'], prefix + '_prodigal.fna')

prodigal_opts = {'-o':prodigal_out, '-a':prodigal_aa, '-d':prodigal_nt}

#Generate and run the Prodigal command
prodigal_command = prodigal_command_generate(ntfasta=nt_fasta, optdict=prodigal_opts,
outfmt=prodigal_outfmt, prodigal='prodigal')

logger.debug('Prodigal will be run in %s mode' % prodigal_command[1]['-p'])

subprocess.run(prodigal_command[0], shell=False)

if os.path.exists(prodigal_out):
    if os.stat(prodigal_out).st_size != 0:
        prodigal_gff_df = gff3_to_pddf(gff = prodigal_out, ftype = 'CDS', index_col=False)
    else:
        logger.error('Prodigal GFF %s is empty' % prodigal_out)
else:
    logger.error('Prodigal GFF output %s does not exist' % prodigal_out)


if os.path.exists(prodigal_aa):
    if os.stat(prodigal_aa).st_size != 0:
         prodigal_aa_dict = make_seqdict(prodigal_aa, format='fasta')

    else:
        logger.error('Prodigal amino acid fasta %s is empty' % prodigal_aa)
else:
    logger.error('Prodigal amino acid fasta %s does not exist' % prodigal_aa)

# PRODIGAL_SOURCES = [nt_fasta]
# PRODIGAL_TARGETS = [prodigal_out, prodigal_aa, prodigal_nt]
# PRODIGAL_OPTS = prodigal_opts
# PRODIGAL_OPTS['-o'] = '${TARGETS}[0]'
# PRODIGAL_OPTS['-a'] = '${TARGETS}[1]'
# PRODIGAL_OPTS['-d'] = '${TARGETS}[2]'
# PRODIGAL_COMMAND = prodigal_command_generate(ntfasta=PRODIGAL_SOURCES[0], optdict=PRODIGAL_OPTS,
# outfmt=prodigal_outfmt, prodigal='prodigal')
#
# PRODIGAL_CMDBLD = scons_command(targets = PRODIGAL_TARGETS, sources = PRODIGAL_SOURCES, command = PRODIGAL_COMMAND[0])
#
# SConstruct.write(PRODIGAL_CMDBLD + '\n')
#==============================================================================
#Fetch the CRISPR-neighboring genes
neighbor_aa_fasta = os.path.join(output_paths['Prodigal'], nt_fasta_basename + '_CRISPR-neighbor-genes.faa')
# print(crispr_gff_df)
# print(prodigal_gff_df)
# print(prodigal_aa_dict)
#neighbor_nt_fasta = os.path.join(prodigal_outdir, nt_fasta_basename + '_CRISPR-neighbor-genes.fna')
neighbor_gene_clusters = fetch_clusters(anchor_gff_df=crispr_gff_df, gene_gff_df=prodigal_gff_df, gene_seq_dict=prodigal_aa_dict, winsize=opts.window_extent, att_fs=';')
#==============================================================================

#class CRISPR_Cas:

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

#SConstruct.close()
if gzip:
    subprocess.run(['gzip', nt_fasta], shell=False)
else:
    pass
