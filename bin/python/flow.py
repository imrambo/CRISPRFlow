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
            if gff_gene:

            seq_objs = [gene_seq_dict[key] for key in gene_seq_dict.keys() if row[anchor_col] + '_' in key and gene_seq_dict[key]['gene_start'] >= row['start'] - winsize  and gene_seq_dict[key]['gene_stop'] <= row['end'] + winsize]
            neighbor_genes.extend(seq_objs)
        SeqIO.write(neighbor_genes, fa, 'fasta')
#------------------------------------------------------------------------------
def command_generator(seqdb, profile, prefix=None, optdict=None, parallel=False, progbar=False, program='hmmsearch', psuffix='.hmm', outdir='/tmp', jobs=1, joblog='joblog'):
    """
    Generate a list of bash command strings.
    Use parallel=True to run parallel jobs using GNU parallel.
    Pass an ordered dictionary if parameter order matters.
    """
    if not prefix:
        prefix = get_basename(seqdb)

    commands = []

    if parallel and jobs == 1:
        parallel = False
    if not parallel:
        if os.path.isdir(profile):
            profile_glob = os.path.join(profile, '*%s' % psuffix)
            for i in list(glob.glob(profile_glob, recursive = False)):
                cmd = '%s %s %s %s' % (program, optstring, i, seqdb)
                commands.append(cmd)

        elif os.path.isfile(profile):
            #Run program for a single profile
            optstring = optstring_join(optdict)
            cmd = '%s %s' % (program, optstring)
            commands.append(cmd)
        else:
            pass
    elif parallel and jobs > 1:
        #Build the parallel string
        parallel_optdict = {'--jobs':jobs}
        if progbar:
            parallel_optdict['--bar'] = ''
        parallel_optdict['--joblog'] = joblog
        parallel_optstring = opstring_join(parallel_optdict)

        optstring = optstring_join(optdict)
        if os.path.isdir(profile):
            cmd = 'find %s -type f -name "*%s" | parallel %s %s %s {} %s' % (profile, psuffix, parallel_optstring, program, optstring, seqdb)
            commands.append(cmd)
        """
        If a list of paths is provided, break it up into files and directories
        and generate strings.
        """
        if isinstance(profile, list):
            pfiles = [os.path.isfile(p) for p in profile]
            pdirs = [os.path.isdir(p) for p in profile]
            if len(pfiles) + len(pdirs) == len(profile):
                    if pdirs:
                        for pdir in pdirs:
                            cmd = 'find %s -type f -name "*%s" | parallel %s %s %s {} %s' % (pdir, psuffix, parallel_optstring, program, optstring, seqdb)
                            commands.append(cmd)
                    if pfiles and len(pfiles) > 1:
                        cmd = 'parallel %s %s %s {} %s ::: %s' % (parallel_optstring, program, optstring, seqdb, ' '.join(pfiles))
                        commands.append(cmd)
    else:
        pass
    return commands
#==============================================================================
build_root = '../..'
db_ids = ['CAS', 'KO', 'PFAM', 'TIGR']
regex = r'%s\d+' % '|'.join(db_ids)
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
parser.add_argument('--db_set', type=str, dest='db_set', action='store',
nargs='?', help='Comma-separated list of HMM database paths to use for CRISPR-proximity gene hmmsearch. You can also specify paths to certain profiles, e.g. /path/to/K00001.hmm')
parser.add_argument('--ccs_typing', type=str, dest='ccs_typing', action='store', default='sub-typing',
help='Level of CRISPR-Cas system typing specificity. Choose from: [general, typing, sub-typing]')

opts = parser.parse_args()
#==============================================================================
#Create the output directories
output_paths = {p: os.path.join(opts.out_root, p) for p in ['CRISPRDetect','Prodigal','MacSyFinder','HMMER']}
for key, value in output_paths:
    if not os.path.exists(value):
        os.makedirs(value)
    else:
        pass
#==============================================================================
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
'--def':'../data/definitions/%s' % opts.ccs_typing, '--out-dir':output_paths['MacSyFinder'],
'--res-search-suffix':'hmmout', '--res-extract-suffix':'res_hmm_extract',
'--profile-suffix':'hmm', '--profile-dir':'../data/profiles/CAS',
'--worker':opts.threads, '--verbosity':'-vv'}

if opts.joblog_dir:
    macsyfinder_opts['--log'] = opts.joblog_dir
else:
    macsyfinder_opts['--log'] = output_paths['MacSyFinder']

macsyfinder_optstring = optstring_join(macsyfinder_opts)
macsyfinder_command = command_generator(program='macsyfinder', optdict=macsyfinder_opts,
seqdb=prodigal_command[1]['-a'])
seqdb, profile, prefix=None, optdict=None, parallel=False, progbar=False, program='hmmsearch', psuffix='.hmm', outdir='/tmp', jobs=1, joblog='joblog')
#macsyfinder_command = 'macsyfinder %s' % macsyfinder_optstring
#subprocess.run([macsyfinder_command])
#==============================================================================
###---HMMSEARCH---###
# hmmsearch_opts = {'--domE':10, '-E':10, '--incE':1e-6, '--incdomE':1e-6, '--seed':42}
#
# hmmsearch_joblog = os.path.join(output_paths['HMMER'], 'hmmsearch_joblog_' + now.strftime('%D-%M-%Y_%H:%M'))
#
# if ',' in opts.db_set:
#     db_list = ','.split(opts.db_set)
#     pfiles = []
#     for db in db_list:
#         if os.path.isdir(db):
#             if opts.jobs > 1:
#                 command_generate(seqdb=neighbor_aa_fasta, profile=db,
#                 optdict=hmmsearch_opts, jobs=opts.jobs, joblog=hmmsearch_joblog,
#                 outdir=output_paths['HMMER'], psuffix='.hmm', optdict=hmmsearch_opts,
#                 parallel=True, jobs=opts.jobs, joblog=hmmsearch_joblog)
#
#         elif os.path.isfile(db):
#             pfiles.append(db)
#
# hmmer.command_generate(seqdb=neighbor_aa_fasta, profile=, optdict=hmmsearch_opts)
#
# seqdb, profile, outdir, psuffix='.hmm', optdict, prefix*, parallel=False, jobs*, joblog*
#==============================================================================
#transposase_dict = dict()
#==============================================================================
#Read in the nucleotide scaffolds
# scaffold_handle = open(opts.scaffold_file, 'r')
# scaffold_obj = SeqIO.to_dict(SeqIO.parse(scaffold_handle), 'fasta')

# #Jackhmmer options
# jackhmmer_opts = {'--mx':'BLOSUM62', '-E':1e-20, '--domE':1e-20, '-N':20,
# '--incE':1e-20, '--incdomE':1e-20, '--F1':0.01, '--seed':42, '--cpu':4,
# '--noali':''}

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
#==============================================================================
#==============================================================================
# ###---CRISPRCasFinder---###
#  ccfinder_opts = {'-log':'', '-copyCSS':'', '-repeats':'', '-DBcrispr':'',
# '-DIRrepeat':'', '-cas':'', '-ccvRep':'', '-vicinity':10000, '-cluster':20000,
# '-minDR':23, '-maxDR':72, '-minSP':8, '-maxSP':64, '-cpuM':1,
# '-definition':'SubTyping', '-getSummaryCasfinder':'', '-betterDetectTrunc':'',
# '-soFile':ccfinder_sofile, '-keep':''}
# ccfinder_optstring = optstring_join(ccfinder_opts)
