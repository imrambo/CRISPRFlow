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
def prodigal_command_generate(ntfasta, outdir, optdict, prefix, outfmt='gff', version=2, prodigal='prodigal'):
    """
    Generate an argument list to run Prodigal using subprocess.run()
    """
    #HOW CAN I EXTRACT THE OUTPUT OF prodigal -v ???
    #GUNZIP PIPE

    if not optdict['-p']:
        prodigal_mode = prodigal_mode_select(ntfasta, version=version)
        optdict['-p'] = prodigal_mode[0]
    else:
        pass

    if not optdict['-i']:
        optdict['-i'] = ntfasta

    prodigal_command = exec_cmd_generate(prodigal, optdict)

    return prodigal_command,optdict

    # #Input file is gzipped, use a pipe and omit input option
    # if prodigal_mode[1]:
    #     if '-i' in optdict:
    #         del optdict['-i']
    #     zcat = subprocess.run(['zcat', ntfasta], stdout=subprocess.PIPE)
    #
    #     #optdicttring = optstring_join(optdict)
    #     #prodigal_command = 'zcat %s | prodigal %s' % (ntfasta, optdicttring)
    #     else:
    #         pass
    #
    # else:
    #     if not optdict['-i']:
    #         optdict['-i'] = ntfasta
    #     prodigal_command = exec_cmd_generate(prodigal, optdict)
    #
    # return prodigal_command,optdict
