
from gff3 import gff3_to_pddf
import gzip
from Bio import SeqIO
from shell_tools import is_gzipped
#------------------------------------------------------------------------------
def make_seqdict(input_file, format='fasta'):
    """
    Create a SeqIO sequence dictionary.
    """
    if is_gzipped(input_file):
        seq_handle = gzip.open(input_file, 'rb')
    else:
        seq_handle = open(input_file, 'r')

    seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, format))

    return seq_dict
#------------------------------------------------------------------------------
def fetch_clusters(anchor_gff_df, gene_gff_df, gene_seq_dict, winsize, att_fs=';'):
    """
    Get genes within a certain distance from a genomic feature.
    The anchor_gff_df must only contain 'anchor' genomic features, i.e. CRISPR array.
    The gene_gff_df contains gene or other feature coordinates you want to fetch.
    The seq_dict is a SeqIO Sequence Dict.
    The winsize will fetch features within x bp of the anchor start and end.
    att_fs is the field separator for the attribute gff column.
    """

    cluster_genes = dict()

    #for i in anchor_gff_df.index:
    for i, row in enumerate(anchor_gff_df.itertuples(), 0):
        anchor_source = anchor_gff_df.at[i, 'source']
        anchor_start = anchor_gff_df.at[i, 'start']
        anchor_end = anchor_gff_df.at[i, 'end']
        anchor_id = anchor_gff_df.at[i, 'attributes'].split(att_fs)[0].split('=')[1].split('_')[0]

        gene_cluster_df = gene_gff_df[(gene_gff_df['source'] == anchor_source) & (gene_gff_df['start'] >= anchor_start - winsize) & (gene_gff_df['end'] <= anchor_end + winsize)]
        gene_cluster_df['gene_id'] = gene_cluster_df['source'].astype(str) + '_' + cluster_df['attributes'].str.split(att_fs).str[0].str.split('=').str[1].str.split('_').str[1]

        #seq_objs = [gene_seq_dict[key] if gid[1] == gene_seq_dict[key].description.split('#')[0] for key in gene_seq_dict.keys() for gid in gene_cluster_df['gene_id'].iteritems()]
        seq_objs = [gene_seq_dict[key] for key in gene_seq_dict.keys() for gid in gene_cluster_df['gene_id'].iteritems() if gid[1] == gene_seq_dict[key].description.split('#')[0]]


        cluster_genes[anchor_id] = seq_objs

    return cluster_genes
