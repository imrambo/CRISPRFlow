#!/usr/bin/python
def domtbl_besthits(source):
    """
    Get the best hit for each target AA sequence from HMMER3 domain
    table output. The best hit is based on:
    1. min independent e-value
    2. max bitscore
    3. max alignment length
    4. max mean posterior probability of aligned residues
    in the maximum expected accuracy alignment, i.e. reliability of overall
    alignment
    """
    domtbl = pd.read_csv(source, comment='#', header=None,
    names = ['query_alignment_name','target_name','target_accession',
    'target_len','query_name','accession','query_length','evalue_sequence',
    'bitscore_sequence','bias','domain_number','n_domains','evalue_conditional',
    'evalue_independ','bitscore_domain','bias_domain','hmm_from','hmm_to',
    'ali_from','ali_to','env_from','env_to','acc','target_desc'], sep = '\s+')

    domtbl['alignment_length'] = (domtbl['ali_to'] - domtbl['ali_from']).astype(int)

    aggregations = {'query_alignment_name':'first', 'evalue_independ':min, 'acc':max, 'bitscore_domain':max, 'alignment_length':max}
    dom_agg = domtbl.groupby(['target_name'], as_index = False).agg(aggregations)

    #dom_agg.to_csv(target, sep = '\t', encoding='utf-8', header = True, index = False)

    return dom_agg
#------------------------------------------------------------------------------
def hmm_builder(seqs, name, muscle_iter, suffix):
    #Hmmbuild options
    hmmbuild_opts = {'--cpu':2, '--seed':42, '-n':name}
    hmmbuild_optstring = optstring.join(hmmbuild_opts)
    #MUSCLE alignment
    muscle_command = 'muscle -maxiters %d -in %s -out %s_muscle.afa' % (muscle_iter, seqs, name)
    #Use trimal to remove columns with gaps in 95% or more:
    trimal_command = 'trimal -in %s_muscle.afa -out %s_muscle_trimal.afa -gt 0.95' % (name, name)
    #Build the HMM
    hmmbuild_command = 'hmmbuild %s %s.%s %s_muscle_trimal.afa' % (hmmbuild_optstring, name, suffix, name)

    joined_command = ' & '.join([muscle_command, trimal_command, hmmbuild_command])
    subprocess.run([joined_command])
    return
