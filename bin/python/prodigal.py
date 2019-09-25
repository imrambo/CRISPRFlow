#!/usr/bin/python3
"""
Motivation: functions for use with Prodigal.

Author: Ian Rambo
Contact: ian.rambo@utexas.edu, imrambo@lbl.gov

Thirteen... that's a mighty unlucky number... for somebody!
"""
from shell_tools import *
from Bio import SeqIO
#------------------------------------------------------------------------------
def prodigal_mode_select(fasta, version=2, len_thresh=20000):
    """
    Select the correct Prodigal mode based on input nt fasta sequence lengths.
    """
    pgz = is_gzipped(fasta)
    #seq_dict = make_seqdict(fasta=fasta, gz=pgz)

    if pgz:
        with gzip.open(fasta, 'rb') as fa:
            lines = fa.readlines()
            lines = [lines.decode('ascii') for line in lines]
    else:
        with open(fasta, 'r') as fa:
            lines = fa.readlines()

    #Are there any sequences less than the threshold for 'single/normal' mode?
    #if any([len(seq_dict[key]['sequence']) < len_thresh for key in seq_dict.keys()]):
    #if any([len(line) < len_thresh for line in lines if not line.startswith('>')]):

    #Is the total genome length less than the size threshold for 'single/normal' mode?
    if sum([len(line) for line in lines if not line.startswith('>')]) < len_thresh:
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
def prodigal_command_generate(ntfasta, optdict, outfmt='gff', version=2, prodigal='prodigal'):
    """
    Generate an argument list to run Prodigal using subprocess.run()
    """
    #HOW CAN I EXTRACT THE OUTPUT OF prodigal -v ???
    #GUNZIP PIPE

    if not '-p' in optdict.keys():
        prodigal_mode = prodigal_mode_select(ntfasta, version=version)
        optdict['-p'] = prodigal_mode[0]
    else:
        pass

    if not '-i' in optdict.keys():
        optdict['-i'] = ntfasta
    else:
        pass

    if not '-f' in optdict.keys():
        optdict['-f'] = outfmt
    else:
        pass

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
#------------------------------------------------------------------------------
def fetch_gene_clusters(gff_anchor, gene_seq_dict, out_fasta, winsize, gff_gene=None, prodigal=True):
    """
    Get Prodigal genes within a certain distance from a genomic feature.
    The gff_df must only contain 'anchor' genomic features, i.e. CRISPR array.
    gene_seq_dict is a SeqIO Sequence Dict.

    NOTE: NEED TO CHANGE CODE TO USE GFF FILES FOR COORDINATES
    """
    with open(out_fasta, 'w') as fa:
        neighbor_genes = []
        #Loop through records and fetch neighboring genes
        if prodigal:
            print(gff_anchor.head())
            #for index, row in gff_anchor.iterrows():
            for i in gff_anchor.index:
                start = gff_anchor.get_value(i, 'start')
                end = gff_anchor.get_value(i, 'end')

                print(start)
                print(end)
                #if gff_gene:
                #seq_objs = [gene_seq_dict[key] for key in gene_seq_dict.keys() if row[anchor_col] + '_' in key and int(gene_seq_dict[key].description.split('#')[1] >= row['start']) - winsize  and int(gene_seq_dict[key].description.split('#')[2]) <= row['end'] + winsize]
                #seq_objs = [gene_seq_dict[key] for key in gene_seq_dict.keys() if index + '_' in key and int(gene_seq_dict[key].description.split('#')[1] >= row['start']) - winsize  and int(gene_seq_dict[key].description.split('#')[2]) <= row['end'] + winsize]
                seq_objs = [gene_seq_dict[key] for key in gene_seq_dict.keys() if i + '_' in key and int(gene_seq_dict[key].description.split('#')[1] >= int(start)) - winsize  and int(gene_seq_dict[key].description.split('#')[2]) <= int(end) + winsize]

                neighbor_genes.extend(seq_objs)
            SeqIO.write(neighbor_genes, fa, 'fasta')
        else:
            pass
