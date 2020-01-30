
# coding: utf-8

# In[684]:


"""
Ian Rambo
June 14, 2019
"""


# In[685]:


get_ipython().run_line_magic('matplotlib', 'inline')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re
from functools import reduce


# In[686]:





# In[687]:


def fetch_bin(df, col):
    """
    Extract bin ID.
    """
    df['bin'] = df[col].str.extract(r'(.*?)_scaffold_.*?')
    df['bin1'] = df[col].str.extract(r'^(.*?\.\d+)_\d+$')
    df['bin_id'] = df['bin'].fillna('') + df['bin1'].fillna('')
    df = df.drop(columns=['bin', 'bin1'])
    return df


# In[688]:


def search_bin(df1, df2, idcol, scol1, scol2):
    """
    Using a column of identifiers in df1 (idcol), search for the identifier in
    scol1 and scol2 in df2, which are object columns. If the identifier is found
    in either of the columns, add to the empty dataframe. You can then merge
    this dataframe to your hits by index.
    """
    df_hit = pd.DataFrame()
    for value in df1[idcol].unique():
        bool1 = df2[scol1].str.contains(value)
        bool2 = df2[scol2].str.contains(value)
        booldf = pd.concat([bool1, bool2], axis = 1)
        booldf[idcol] = np.where(booldf[scol1] + booldf[scol2] > 0, value, np.nan)
        hit = booldf.loc[booldf[scol1] | booldf[scol2]]
        df_hit = df_hit.append(hit)
    return df_hit


# In[689]:


#Directory containing CRISPR Cas cluster files
workdir = '/Users/ian/Documents/phd_research/crispr/ccfinder_clusters'
os.chdir(workdir)


# In[690]:


#Load results files for CRISPR Cas clusters
compchrom = pd.read_csv('CRISPR_cas_clusters_compChrom.tsv', sep='\s+')
sfarchaea = pd.read_csv('CRISPR_cas_clusters_sfbay_archaea.tsv', sep='\s+')
sfbacteria = pd.read_csv('CRISPR_cas_clusters_sfbay_bacteria.tsv', sep='\s+')
guayhg = pd.read_csv('CRISPR_cas_clusters_guaymasHG.tsv', sep='\s+')
guay18a = pd.read_csv('CRISPR_cas_clusters_guaymas2018_archaea.tsv', sep='\s+')
guay18b = pd.read_csv('CRISPR_cas_clusters_guaymas2018_bacteria.tsv', sep='\s+')
mbnerr = pd.read_csv('CRISPR_cas_clusters_mb_nerr.tsv', sep='\s+')
ynp = pd.read_csv('CRISPR_cas_clusters_YNP_FS.tsv', sep='\s+')
ncrep = pd.read_csv('CRISPR_cas_clusters_noncomp_repset.tsv', sep='\s+')


# In[691]:


#Add sampling locations
compchrom['Sample_Location'] = 'NCBI complete/chromosome reference'
sfarchaea['Sample_Location'] = 'San Francisco Bay'
sfbacteria['Sample_Location'] = 'San Francisco Bay'
guayhg['Sample_Location'] = 'Guaymas Basin Holy Grail'
ncrep['Sample_Location'] = 'NCBI partial reference'
ynp['Sample_Location'] = 'Yellowstone hot spring'
mbnerr['Sample_Location'] = 'Mesquite Bay'
guay18a['Sample_Location'] = 'Guaymas Basin'
guay18b['Sample_Location'] = 'Guaymas Basin'


# In[692]:


#Extract genome names if applicable
guayhg['BinID'] = guayhg['Sequence'].str.extract(r'(.*?)_scaffold_.*?')
guayhg['BinID'] = guayhg['BinID'].str.replace('Guay', 'G')

guay18a['BinID'] = guay18a['Sequence'].str.extract(r'(.*?)_scaffold_.*?')
guay18b['BinID'] = guay18b['Sequence'].str.extract(r'(.*?)_scaffold_.*?')
guay18ab = pd.concat([guay18a, guay18b], sort=False, axis=0)


# In[693]:


guay18b


# In[694]:



#NCBI reference genome FTP table
#taxa_ncbi = pd.read_csv('taxonomy_ncbi_04-16-2019.csv', sep=',')

#Add extended sample information for guaymas
guay18_archaea_info = pd.read_csv('/Users/ian/Documents/phd_research/guaymas_2019/2019_Guaymas_Basin_Bin_information_archaea.txt', header = 0, sep = '\t')
guay18_bacteria_info = pd.read_csv('/Users/ian/Documents/phd_research/guaymas_2019/2019_Guaymas_Basin_Bin_information_bacteria.txt', header = 0, sep = '\t')
#guay18_archaea_info['BinID'] = guay18_archaea_info['Bin_new_name'].str.replace('_renamed$', '')
#guay18_bacteria_info['BinID'] = guay18_bacteria_info['Bin_new_name'].str.replace('_renamed$', '')

#guay18_bacteria_info


# In[695]:


#Load taxonomy files

#NCBI reference genome FTP table
#taxa_ncbi = pd.read_csv('taxonomy_ncbi_04-16-2019.csv', sep=',')
guay18_archaea_info = pd.read_csv('/Users/ian/Documents/phd_research/guaymas_2019/2019_Guaymas_Basin_Bin_information_archaea.txt', header = 0, sep = '\t')
guay18_bacteria_info = pd.read_csv('/Users/ian/Documents/phd_research/guaymas_2019/2019_Guaymas_Basin_Bin_information_bacteria.txt', header = 0, sep = '\t')
guay18_archaea_info['BinID'] = guay18_archaea_info['Bin_new_name'].str.replace('_renamed$', '')
#guay18_bacteria_info['BinID'] = guay18_bacteria_info['Bin_new_name'].str.replace('_renamed$', '')

#Merge site info for new guaymas
guay18_ab_info = pd.concat([guay18_archaea_info, guay18_bacteria_info], axis=0, sort=False)
guay18_ab_site = guay18_ab_info[['BinID', 'Site']].rename(columns={'Site':'Cluster ID'}).dropna()
guay18ab = guay18ab.merge(guay18_ab_site)

guay18ab['Full_Location'] = guay18ab['Sample_Location'] + '_' + guay18ab['Cluster ID']

guay18a = guay18a.merge(guay18ab).drop(columns=['Cluster ID']).rename(columns={'Full_Location':'Cluster ID'})
guay18b = guay18b.merge(guay18ab).drop(columns=['Cluster ID']).rename(columns={'Full_Location':'Cluster ID'})



# In[696]:


#Guaymas Holy Grail
taxa_guayhg = pd.read_csv('guaymas_hg_genome_information.csv', sep=',')


# In[697]:


#NCBI complete/chromosome
taxa_compchrom = pd.read_csv('compchrom_taxa_hit.csv', sep=',')


# In[698]:


#Create taxonomy mapping dataframe -- SLOW
#compchrom_hit = search_bin(compchrom, taxa_ncbi, 'Sequence', 'Assembly', 'Replicons')
#Write the taxonomy mapping dataframes
#compchrom_hit.to_csv('compchrom_taxa_hit.csv', sep = ',', encoding = 'utf-8', header = True)
guay18ab.head()


# In[699]:


taxa_guay18a = pd.read_csv('2019_Guaymas_Basin_Bin_information_archaea.txt', sep='\s+')
taxa_guay18b = pd.read_csv('2019_Guaymas_Basin_Bin_information_bacteria.txt', sep='\s+')

taxa_guay18a = taxa_guay18a.rename(columns={'Bin_new_name':'BinID', 'Marker':'Cluster ID'})
taxa_guay18a['BinID'] = taxa_guay18a['BinID'].str.replace('_renamed', '')
taxa_guay18a['Cluster ID'] = taxa_guay18a['Cluster ID'].str.replace('^\w__', '')

taxa_guay18b = taxa_guay18b.rename(columns={'Bin_new_name':'BinID', 'Marker':'Cluster ID'})
taxa_guay18b['BinID'] = taxa_guay18b['BinID'].str.replace('_renamed.fna', '')
taxa_guay18b['Cluster ID'] = taxa_guay18b['Cluster ID'].str.replace('^\w__', '')


# In[700]:


#Getting taxonomy information for genbank hits (non-complete references)
#time for i in $(tail -n +2 CRISPR_cas_clusters_noncomp_repset.tsv | cut -f1 | awk 'length($0) > 8' | uniq); do v=$(echo $i | cut -c-6); grep "\"$v\"" prokaryotes_ncbi_06-13-2019.csv | cut -d',' -f2 | sed "s/$/\,${i}/"; done >CRISPR_cas_clusters_noncomp_genbank_taxmap.csv
taxa_noncomp_genbank = pd.read_csv('CRISPR_cas_clusters_noncomp_genbank_taxmap.csv', sep=',', header=None, names=['Organism Groups', 'Sequence'])
taxa_noncomp_genbank['Cluster ID'] = taxa_noncomp_genbank['Organism Groups'].str.split(pat = ';').str[2]
taxa_noncomp_genbank = taxa_noncomp_genbank.rename(columns={'Sequence':'BinID'})
#Getting taxonomy information for non-genbank hits (non-complete references)
#time for i in $(tail -n +2 CRISPR_cas_clusters_noncomp_repset.tsv | cut -f1 | awk 'length($0) <= 8' | uniq); do v=$(echo $i | cut -c-6); grep "\"$v\"" prokaryotes_ncbi_06-13-2019.csv | cut -d',' -f2 | sed "s/$/\,${i}/"; done >CRISPR_cas_clusters_noncomp_taxmap.csv
#taxa_noncomp = pd.read_csv('CRISPR_cas_clusters_noncomp_taxmap.csv')


# In[701]:


#Merge taxonomy information to cluster dfs
guayhg_wtax = pd.merge(guayhg, taxa_guayhg[['BinID', 'Cluster ID']])


# In[702]:


ncrep = ncrep.rename(columns={'Sequence':'BinID'})
#taxa_noncomp_genbank.head()


# In[703]:


ncrep_gbtax = pd.merge(ncrep, taxa_noncomp_genbank[['BinID', 'Cluster ID']])
#ncrep_gbtax = ncrep_gbtax.rename(columns={'Sequence':'BinID'})


# In[704]:


#ncrep_wtax = pd.merge(ncrep, taxa_noncomp[['Sequence', 'Cluster ID']])
#ncrep_wtax = ncrep_tax.rename(columns={'Sequence':'BinID'})
ncrep_gbtax.head()


# In[705]:


ncbi_alltax_gca = pd.read_csv('/Users/ian/Documents/phd_research/crispr/tree/tol/ncbi_taxonomy_mapping_version2.txt', sep='\t', header=None, names=['BinID', 'Organism Groups'])
ncbi_alltax_gca['Cluster ID'] = ncbi_alltax_gca['Organism Groups'].str.split(pat = ';').str[2]
ncbi_alltax_gca.head()


# In[706]:


compchrom_wtax = pd.merge(compchrom, taxa_compchrom[['Sequence', 'Organism Groups']])
compchrom_wtax['Cluster ID'] = compchrom_wtax['Organism Groups'].str.split(pat = ';').str[2]
compchrom_wtax = compchrom_wtax.rename(columns={'Sequence':'BinID'})


# In[707]:


guay18a_wtax = pd.merge(guay18a, taxa_guay18a[['BinID', 'Cluster ID']])
guay18b_wtax = pd.merge(guay18b, taxa_guay18b[['BinID', 'Cluster ID']])


# In[708]:


dataframes = [df[['BinID', 'Cluster ID', 'Nb_CRISPRs', 'Nb_Cas', 'Description']] for df in [guayhg_wtax, ncrep_gbtax, compchrom_wtax, guay18a_wtax, guay18b_wtax]]
dfmaster_wtax = pd.concat(dataframes)


# In[ ]:


dataframes_location = [df[['Sample_Location', 'Nb_CRISPRs', 'Nb_Cas', 'Description']] for df in [guayhg_wtax, ncrep_gbtax, compchrom_wtax, guay18a_wtax, guay18b_wtax, mbnerr, ynp]]

dflocation = pd.concat(dataframes_location)
dflocation_ccs = dflocation[(dflocation['Nb_CRISPRs'] > 0) & (dflocation['Nb_Cas'] > 0)]
#Explode the cluster descriptions into each CCS
s = dflocation_ccs[['Sample_Location', 'Description']]
loc_ccs = pd.DataFrame([[x] + [z] for x, y in zip(s.Sample_Location,s.Description.str.split(',')) for z in y],
                  columns=['Sample_Location', 'CCS Subtype'])

#Filter only CCS subtypes
loc_ccs = loc_ccs[loc_ccs['CCS Subtype'].str.contains('^CAS')]
loc_ccs['CCS num'] = loc_ccs['CCS Subtype'].str.extract(r'^CAS.*?_(n\d+)\[')
loc_ccs['CCS Subtype'] = loc_ccs['CCS Subtype'].str.extract(r'(^CAS.*?)_n\d+\[')

loc_ccs = pd.merge(loc_ccs, dflocation_ccs)
loc_ccs = loc_ccs.drop_duplicates()

# loc_ccs['Cluster ID'] = loc_ccs['Cluster ID'].str.replace(r'unclassified|\(miscellaneous\)', '')
# loc_ccs['Cluster ID'] = loc_ccs['Cluster ID'].str.replace(r'\s+Bacteria\s+', 'Bacteria')
# loc_ccs = loc_ccs[(ccs['Cluster ID'] != 'environmental samples') & (loc_ccs['CCS Subtype'] != 'CAS')]
loc_ccs


# In[728]:


#Select hits for lone Cas genes
dfmaster_wtax_orphan = dfmaster_wtax[(dfmaster_wtax['Nb_CRISPRs'] == 0) & (dfmaster_wtax['Nb_Cas'] > 0)]

o = dfmaster_wtax_orphan[['BinID', 'Description']]
orphan = pd.DataFrame([[x] + [z] for x, y in zip(o.BinID,o.Description.str.split(',')) for z in y],
                  columns=['BinID', 'CCS Subtype'])

orphan = orphan[orphan['CCS Subtype'].str.contains('^CAS')]
orphan['CCS num'] = orphan['CCS Subtype'].str.extract(r'^CAS.*?_(n\d+)\[')
orphan['CCS Subtype'] = orphan['CCS Subtype'].str.extract(r'(^CAS.*?)_n\d+\[')
orphan = pd.merge(orphan, dfmaster_wtax_orphan[['BinID', 'Cluster ID']])


# In[643]:





# In[ ]:


orphan_count = orphan.groupby(['CCS Subtype', 'Cluster ID']).agg('count')


# In[747]:


orphan_entropy_df = orphan.sample(n=100, replace=True, random_state=1)
for i in range(1,101):
    orphan_subset = orphan.sample(n=50, replace=True, random_state=i)
    orphan_entropy = relative_entropy(df = orphan_subset, groupvar = 'Cluster ID', countvar1 = 'BinID', countvar2 = 'CCS Subtype')
    #orphan_entropy_mean = orphan_entropy.fillna(1e-99).groupby('Cluster ID').agg('mean')
    orphan_entropy_df = orphan_entropy_df.append(orphan_entropy, sort=True, ignore_index=False)
orphan_entropy_df = orphan_entropy_df.reset_index()


# In[749]:



orphan_entropy = relative_entropy(df = orphan, groupvar = 'Cluster ID', countvar1 = 'BinID', countvar2 = 'CCS Subtype')


# In[751]:


#orphan_entropy = relative_entropy(df = orphan, groupvar = 'Cluster ID', countvar1 = 'BinID', countvar2 = 'CCS Subtype')
orphan_entropy


# In[93]:


orphan_entropy = orphan_entropy[(orphan_entropy['Cluster ID'] != 'environmental samples') & (orphan_entropy['CCS Subtype'] != 'CAS')]

orphan_entropy_pivot = orphan_entropy.pivot(index='Cluster ID', columns='CCS Subtype', values='H')
orphan_entropy_pivot = orphan_entropy_pivot.fillna(1e-99)

orphan_entropy_pivot_norm = (orphan_entropy_pivot - orphan_entropy_pivot.mean()) / (orphan_entropy_pivot.max() - orphan_entropy_pivot.min())
orphan_entropy_norm_stack = orphan_entropy_pivot_norm.stack().reset_index().rename(columns={0:'H'})
#orphan_entropy_norm_stack['Cluster ID'].unique()


# In[760]:


#oens_candidatus = orphan_entropy_norm_stack[orphan_entropy_norm_stack['Cluster ID'].str.contains('^Candidatus|^Asgard|^Bathyarchaeota|^Crenarchaeota|^Korarchaeota|^Nanohaloarchaeota|^Euryarchaeota|^Altiarchaeales|^Thermotogae')]
oens_candidatus = orphan_entropy_norm_stack[orphan_entropy_norm_stack['Cluster ID'].str.contains('^Candidatus')]
#oens_candidatus = oens_candidatus[(oens_candidatus['CCS Subtype'] != 'CAS-TypeIE')|(oens_candidatus['CCS Subtype'] != 'CAS-TypeIIA')]
plt.figure(figsize=(18,15))
sns.set_style('whitegrid')
plt.xticks(rotation='vertical', fontsize=18)
plt.yticks(fontsize=18)
#sns.color_palette("cubehelix", 20)
violin = sns.violinplot(x='CCS Subtype', y='H', data=oens_candidatus, color = 'k')
#v.set_xticklabels('CRISPR-Cas Subtype', rotation=30)
# v.set_xlabel("CRISPR-Cas type-subtype", fontsize=10)
# v.set_ylabel("Taxonomy",fontsize=10)
# v.tick_params(labelsize=10)
markers = ['o','v','s','p','*']
swarm = sns.swarmplot(x='CCS Subtype', y='H', data=oens_candidatus,
                      hue='Cluster ID', palette="Set2", size=12)
plt.legend(bbox_to_anchor=(1, 1), loc=2, fontsize=32)
swarm.set_ylabel("Relative entropy H'",fontsize=28)
swarm.set_xlabel("CRISPR-Cas subtype association", fontsize=28)

#plt.setp(swarm.get_legend().get_texts(), fontsize='22')
#v.savefig("violin_swarm_orphan_cas_entropy.png")


# In[ ]:





# In[409]:


#CRISPR-Cas system hits
dfmaster_wtax_ccs = dfmaster_wtax[(dfmaster_wtax['Nb_CRISPRs'] > 0) & (dfmaster_wtax['Nb_Cas'] > 0)]


# In[410]:


#Explode the cluster descriptions into each CCS
s = dfmaster_wtax_ccs[['BinID', 'Description']]
ccs = pd.DataFrame([[x] + [z] for x, y in zip(s.BinID,s.Description.str.split(',')) for z in y],
                  columns=['BinID', 'CCS Subtype'])

#Filter only CCS subtypes
ccs = ccs[ccs['CCS Subtype'].str.contains('^CAS')]
ccs['CCS num'] = ccs['CCS Subtype'].str.extract(r'^CAS.*?_(n\d+)\[')
ccs['CCS Subtype'] = ccs['CCS Subtype'].str.extract(r'(^CAS.*?)_n\d+\[')

ccs = pd.merge(ccs, dfmaster_wtax_ccs[['BinID', 'Cluster ID']])
ccs = ccs.drop_duplicates()

ccs['Cluster ID'] = ccs['Cluster ID'].str.replace(r'unclassified|\(miscellaneous\)', '')
ccs['Cluster ID'] = ccs['Cluster ID'].str.replace(r'\s+Bacteria\s+', 'Bacteria')
ccs = ccs[(ccs['Cluster ID'] != 'environmental samples') & (ccs['CCS Subtype'] != 'CAS')]


# In[411]:


ccs_counts = ccs[['BinID', 'CCS Subtype', 'Cluster ID']].groupby(['CCS Subtype', 'Cluster ID']).agg('count')


# In[610]:


#Create a pivot table of CCS per genome

# ccs_pivot = pd.pivot_table(ccs, values='BinID', index=['Cluster ID'],
#                         columns=['CCS Subtype'], aggfunc='count')

ccs_pivot = pd.pivot_table(ccs, values='CCS num', index=['BinID'],
                         columns=['CCS Subtype'], aggfunc='count')

# ccs_pivot = pd.pivot_table(ccs, values='BinID', index=['CCS Subtype'],
#                         columns=['Cluster ID'], aggfunc='count')

ccs_pivot = ccs_pivot.fillna(0)

#ccs_pivot.to_csv('ccs_subtype_pivot.tsv', sep = '\t', encoding = 'utf-8', header = True)


# In[611]:


#Accessions without hits
# ncrep_noccs = 'noncomp_repset_accessions_noccs.txt'
# compchrom_noccs = 'compchrom_accessions_noccs.txt'
# noccs = []
# for dataframe in [ncrep_noccs, compchrom_noccs]:
#     df = pd.read_csv(dataframe, header=None, names=['BinID'] + list(ccs_pivot.columns))
#     df = df.drop_duplicates()
#     print(df.size)
#     noccs.append(df)
# noccs_df = pd.concat(noccs, axis=0, ignore_index=False).fillna(0)
# noccs_df = noccs_df.set_index('BinID')

noccs = pd.read_csv('./noncomp_compchrom_accessions_noccs.txt', header=None, names=['BinID'] + list(ccs_pivot.columns))

noccs = noccs.fillna(0).set_index('BinID')


# In[612]:


ccs_pivot_combined = pd.concat([ccs_pivot, noccs], axis=0, ignore_index=False, sort=False).reset_index(drop=False)


# In[613]:


missing_taxa = pd.read_csv('missing_taxids.csv', sep='\s+', header=None, names=['BinID', 'Organism Groups'])
missing_taxa['Cluster ID'] = missing_taxa['Organism Groups'].str.split(pat = ';').str[2]


# In[622]:


pmerge_dfs = [ncbi_alltax_gca[['BinID','Cluster ID']], taxa_noncomp_genbank[['BinID','Cluster ID']], taxa_guayhg[['BinID','Cluster ID']], missing_taxa[['BinID', 'Cluster ID']]]
combotax = pd.concat(pmerge_dfs, axis=0)
ccs_pivot_taxonomy = pd.merge(ccs_pivot_combined, combotax, on='BinID', how='left')
ccs_pivot_taxonomy['Cluster ID'] = ccs_pivot_taxonomy['Cluster ID'].str.replace(' ', '_')
ccs_pivot_taxonomy['Cluster ID'] = ccs_pivot_taxonomy['Cluster ID'].str.replace('/', '_')
ccs_pivot_taxonomy['Cluster ID'] = ccs_pivot_taxonomy['Cluster ID'].str.replace('\(|\)', '')

ccs_pivot_taxonomy.to_csv('ccs_pivot_taxonomy.csv', sep = ',', encoding = 'utf-8', header = True)
#ccs_pivot_taxonomy = reduce(lambda left,right: pd.merge(left,right,on='BinID'), pmerge_dfs)


# In[ ]:





# In[623]:


#missing_taxa = ccs_pivot_taxonomy[ccs_pivot_taxonomy.isna().any(axis=1)]
#missing_taxa['BinID'].to_csv('missing_taxa.csv', sep = ',', encoding = 'utf-8', header = False)


# In[ ]:





# In[624]:


#ccs_entropy = relative_entropy(df = ccs, groupvar = 'Cluster ID', countvar1 = 'BinID', countvar2 = 'CCS Subtype')


# In[625]:



#ccs_pivot_taxonomy.loc[:, df.columns != 'b']
#ccs_pivot_taxonomy_notnull = ccs_pivot_taxonomy[(ccs_pivot_taxonomy_notnull.columns != 'BinID') | (ccs_pivot_taxonomy['Cluster ID'].notnull())]


#Only the reference genomes and MAGs with taxonomic assignments
ccs_pivot_taxonomy_notnull = ccs_pivot_taxonomy[ccs_pivot_taxonomy['Cluster ID'].notnull()]


# In[641]:


list(ccs_pivot_taxonomy['Cluster ID'].unique())


# In[627]:


def norm_df(df):
    df_norm = (df - df.mean()) / (df.max() - df.min())
    return df_norm


# In[628]:


def relative_entropy_wide(df, genome_col, group_col):
    #Total number of other genomes
    totalG = df[genome_col].unique().size

    #Genomes per group
    #Q = ccs_pivot_taxonomy_notnull[ccs_pivot_taxonomy_notnull['Cluster ID']!= i].sum()
    #P = ccs_pivot_taxonomy_notnull[ccs_pivot_taxonomy_notnull['Cluster ID']== i].sum()

    #dfQ = pd.DataFrame(Q[1:-1,], columns=['Count']).reset_index()
    #dfP = pd.DataFrame(P[1:-1,], columns=['Count'])

    genome_per_group = df.loc[:, df.columns != genome_col].groupby(group_col).agg('count')

    item_per_group = df.loc[:, df.columns != genome_col].groupby(group_col).agg('sum')

    Qi = item_per_group/(totalG - genome_per_group)
    Pi = item_per_group/genome_per_group

    foo = Pi/Qi
    foo = foo.fillna(1e-99)
    entropy = Pi*np.log2(foo)
    entropy_norm = norm_df(entropy)
    return entropy_norm

#Qi = Qi.reset_index()
#Pi = Pi.reset_index()

#Qi_norm = (Qi - Qi.mean()) / (Qi.max() - Qi.min())
#Pi_norm = (Pi - Pi.mean()) / (Pi.max() - Pi.min())

#foo = Pi_norm/Qi_norm
#Qi = Qi.replace(0.0, 1e-99)
#Qi

#ccs_pivot_taxonomy_notnull['CAS-TypeIA'].sum()
# random_entropy(df, sampn, rown)


# #ccs_pivot_taxonomy.shape
# df = ccs_pivot_taxonomy_notnull.loc[:, ccs_pivot_taxonomy_notnull.columns != 'BinID']
# for i in range(0, sampn):
#     bork = df.sample(n = rown, replace=True)



# In[629]:



pd.set_option('display.max_rows', 500)


# In[630]:


# df_sample = ccs_pivot_taxonomy_notnull.sample(n=10000, replace=True, random_state=4)
# df_sample_entropy = relative_entropy_wide(df_sample, genome_col = 'BinID', group_col = 'Cluster ID')
# df_sample_entropy.fillna(1e-99)

#For each bootstrap replicate, take the mean entropy. Create a dataframe of mean entropies per
#group for each sampling.
entropy_bootstrap = None
entropy_bootstrap = relative_entropy_wide(ccs_pivot_taxonomy_notnull.sample(n=10000, replace=True, random_state=0), genome_col = 'BinID', group_col = 'Cluster ID')
entropy_bootstrap = entropy_bootstrap.fillna(1e-99).groupby('Cluster ID').agg('mean')
for i in range(1, 1001):
    df_sample = ccs_pivot_taxonomy_notnull.sample(n=10000, replace=True, random_state=i)
    df_sample_entropy = relative_entropy_wide(df_sample, genome_col = 'BinID', group_col = 'Cluster ID')
    df_entropy_mean = df_sample_entropy.fillna(1e-99).groupby('Cluster ID').agg('mean')
    entropy_bootstrap = entropy_bootstrap.append(df_sample_entropy, sort=True, ignore_index=False)
entropy_bootstrap = entropy_bootstrap.reset_index()


# In[637]:


def violin_swarm(df, plotid, ymin, ymax):
    #Create a violin plot with swarm plot overlay
    plt.figure(figsize=(20,22))
    sns.set_style('whitegrid')
    plt.xticks(rotation='vertical', fontsize=40)
    plt.yticks(fontsize=35)
    #sns.color_palette("cubehelix", 20)
    sns.violinplot(data=df, color = 'k')
    #swarm = sns.swarmplot(data=df, palette="Set2", size=6)
    sns.swarmplot(data=df, palette="Set2", size=6)
    #swarm.set_ylabel("Mean relative entropy H'",fontsize=22)
    #swarm.set_xlabel('CRISPR-Cas subtype', fontsize=22)
    #swarm.set_title(plotid, fontsize=20)

    plt.ylabel("Mean normalized relative entropy (H')", fontsize=35)
    plt.title(plotid, fontsize=45)
    axes = plt.gca()
    axes.set_ylim([ymin,ymax])
    plotname = 'swarm_violin_%s.png' % plotid
    plt.savefig(plotname, bbox_inches='tight')
    plt.close("all")


# In[638]:


#list(entropy_bootstrap['Cluster ID'].unique())[12:]
[i for i,x in enumerate(list(entropy_bootstrap['Cluster ID'].unique())) if x == 'Crenarchaeota']


# In[640]:


emin = -0.25
emax = 1.25

# for tax in list(entropy_bootstrap['Cluster ID'].unique()):
# #for tax in ['Asgard']:
#     tax_entropy = entropy_bootstrap[entropy_bootstrap['Cluster ID'] == tax]
#     tax_entropy = tax_entropy.loc[:, tax_entropy.columns != 'Cluster ID']
#     violin_swarm(df = tax_entropy, plotid = tax, ymin=emin, ymax=emax)

for tax in list(entropy_bootstrap['Cluster ID'].unique())[34:]:
#for tax in ['Asgard']:
    tax_entropy = entropy_bootstrap[entropy_bootstrap['Cluster ID'] == tax]
    tax_entropy = tax_entropy.loc[:, tax_entropy.columns != 'Cluster ID']
    violin_swarm(df = tax_entropy, plotid = tax, ymin=emin, ymax=emax)



# In[24]:



cas14_guaymas19.loc[cas14_guaymas19['eval_dom'].idxmin()]


# In[95]:


cas14_guaymas19[cas14_guaymas19['eval_dom'] <= 1e-20]


# In[44]:


#cas14_guaymas19[cas14_guaymas19['scaffold'] == 'Meg22_810_Bin_139']

cas14_guaymas19[cas14_guaymas19['scaffold'].str.contains('Meg22_810_Bin_139_')]
