
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
def fetch_clusters(anchor_gff_df, gene_gff_df, gene_seq_dict, winsize):
    """
    Get genes within a certain distance from a genomic feature.
    The anchor_gff_df must only contain 'anchor' genomic features, i.e. CRISPR array.
    The gene_gff_df contains gene or other feature coordinates you want to fetch.
    The seq_dict is a SeqIO Sequence Dict.
    The winsize will fetch features within x bp of the anchor start and end.
    """
    #Do column names match?
    if not list(anchor_gff_df.columns.values) == list(gene_gff_df.columns.values):
        #Reset the column names
        pass

    cluster_genes = []

    for i in anchor_gff_df.index:
        anchor_source = anchor_gff_df.at[i, 'source']
        anchor_start = anchor_gff_df.at[i, 'start']
        anchor_end = anchor_gff_df.at[i, 'end']

        cluster_subset = gene_gff_df[(gene_gff_df.source == anchor_source) & (gene_gff_df.start >= anchor_start - winsize) & (gene_gff_df.end <= anchor_end + winsize)]

        seq_objs = [gene_seq_dict[key] for key in gene_seq_dict.keys() if anchor_source in key and re.match(r'.*?_\d+$', key) and int(gene_seq_dict[key].description.split('#')[1] >= int(start)) - winsize  and int(gene_seq_dict[key].description.split('#')[2]) <= int(end) + winsize]

        cluster_genes.extend(seq_objs)

    return cluster_genes
