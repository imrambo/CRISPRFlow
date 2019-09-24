#!/usr/bin/python3
"""
Motivation: functions for use with Prodigal.

Author: Ian Rambo
Contact: ian.rambo@utexas.edu, imrambo@lbl.gov

Thirteen... that's a mighty unlucky number... for somebody!
"""
import gzip

def prodigal_mode_select(fasta, version=2):
    """
    Select the correct Prodigal mode based on sequence lengths.
    """
    pgz = is_gzip(fasta)
    seq_dict = make_seqdict(fasta=fasta, gz=pgz)

    if any([len(seq_dict[key]['sequence']) < 20000 for key in seq_dict.keys()]):
        if version == 2:
            mode = 'meta'
        elif version == 3:
            mode = 'anon'
        else:
            pass
    else:
        if version == 2:
            mode = 'single'
        elif version == 3:
            mode = 'normal'
        else:
            pass
    return mode,pgz

#------------------------------------------------------------------------------
def prodigal_command_generate(ntfasta, outdir, prefix, outfmt='gff', version=2):
    """
    Generate a command string to run Prodigal on a fasta file.
    """
    #HOW CAN I EXTRACT THE OUTPUT OF prodigal -v ???
    prodigal_mode = prodigal_mode_select(ntfasta, version=version)

    prodigal_out = os.path.join(outdir, prefix + '_prodigal.gff')
    prodigal_aa = os.path.join(outdir, prefix + '_prodigal.faa')
    prodigal_nt = os.path.join(outdir, prefix + '_prodigal.fna')

    prodigal_opts = {'-p':prodigal_mode[0], '-o':prodigal_out,
    '-a':prodigal_aa, '-d':prodigal_nt, '-f':outfmt}

    #Input file is gzipped, use a pipe and omit input option
    if prodigal_mode[1]:
        if '-i' in prodigal_opts:
            del prodigal_opts['-i']
        prodigal_optstring = optstring_join(prodigal_opts)
        prodigal_command = 'zcat %s | prodigal %s' % (ntfasta, prodigal_optstring)

    else:
        if not prodigal_opts['-i']:
            prodigal_opts['-i'] = ntfasta
        prodigal_optstring = optstring_join(prodigal_opts)
        prodigal_command = 'prodigal %s' % (ntfasta, prodigal_optstring)

    return prodigal_command,prodigal_opts
