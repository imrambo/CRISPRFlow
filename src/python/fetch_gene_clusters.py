
from gff import gff_to_pddf

def fetch_gene_clusters(gff_anchor_df, gene_seq_dict, out_fasta, winsize, gff_gene=None, prodigal=True):
    """
    Get Prodigal genes within a certain distance from a genomic feature.
    The gff_df must only contain 'anchor' genomic features, i.e. CRISPR array.
    gene_seq_dict is a SeqIO Sequence Dict.

    NOTE: NEED TO CHANGE CODE TO USE GFF FILES FOR COORDINATES
    """
    gff_anchor_df['sequence_id'] = gff_anchor_df.index()
    gff_gene_df['sequence_id'] = gff_anchor_df.index().str.extract(r'^(.*?)_\d+$')



    #if all([l in list(gff_gene_df.columns.values) for l in list(gff_anchor_df.columns.values)]) and if len(gff_anchor_df.columns) == len(gff_gene_df.columns):
    # gff_concat_df = pd.concatenate([gff_anchor_df, gff_gene_df])





    df[col].str.extract(r'')

    #for i in gff_anchor_df['sequence_id']:


    def fetch_bin(df, col):
    """
    Extract bin ID.
    """
    df['bin'] = df[col].str.extract(r'(.*?)_scaffold_.*?')
    df['bin1'] = df[col].str.extract(r'^(.*?\.\d+)_\d+$')



    with open(out_fasta, 'w') as fa:
        neighbor_genes = []
        #Loop through records and fetch neighboring genes
        if prodigal:
            print(gff_anchor.head())
            #for index, row in gff_anchor.iterrows():
            for i in gff_anchor.index:
                start = gff_anchor.at[i, 'start']
                end = gff_anchor.at[i, 'end']

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
