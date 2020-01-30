#N Tropy
def relative_entropy(df, groupvar, countvar1, countvar2):
    #Number of z
    ngenome = df[countvar1].unique().size
    #Number of z per group
    n_per_group = df[[countvar1, groupvar]].drop_duplicates().groupby(groupvar).agg('count').rename(columns={countvar1:'n_per_group'})
    #Counts of x per group
    x_counts = df[[countvar1,groupvar,countvar2]].drop_duplicates().groupby([groupvar, countvar2]).agg('count').rename(columns={countvar1:'n_count'})

    counts = x_counts.join(n_per_group, how = 'outer')
    counts['Pi'] = counts['n_count']/counts['n_per_group'] #observed probability
    counts['Qi'] = counts['n_count']/ngenome #expected probability
    counts['H'] = counts['Pi']*np.log2(counts['Pi']/counts['Qi']) #relative entropy
    counts_reset = counts.reset_index()
    return counts_reset
