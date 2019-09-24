#!/usr/bin/python3
"""
Motivation: functions to use with hmmsearch, jackhmmer, etc.
"""
import os

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
    import pandas as pd
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
    """
    Build a profile HMM.
    """
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
#------------------------------------------------------------------------------
def hmmsearch_command_generator(db_list, hmmsearch_optdict, parallel_optdict=None, jobs=1):
    """
    Generate command strings for hmmsearch. Set jobs > 1 for GNU Parallel mode.
    Returns a list of command strings.
    """
    commands = []
    dbfiles = [d for d in db_list if os.path.isfile(d)]
    dbdirs = [d for d in db_list if os.path.isdir(d)]

    if jobs > 1 and not parallel_optdict:
        parallel_optdict = {'--jobs':jobs}
    #Loop through the directory paths, and run hmmsearch for all profiles in
    #the path
    if dbdirs:
        for db in dbdirs:
            db_basename = get_basename(db)
            hmmsearch_joblog = os.path.join(output_paths['HMMER'], 'hmmsearch_joblog_%s_%s' % (db_basename, now.strftime('%D-%M-%Y_%H:%M')))
            out_dir = os.path.mkdirs(os.path.join(output_paths['HMMER'], db_basename))

            if jobs > 1:
                parallel_optdict['--joblog'] = hmmsearch_joblog
                parallel_optstring = optstring_join(parallel_optdict)
                parallel_command = 'parallel %s' % parallel_optstring

                hmmsearch_optdict['--domtbl'] = os.path.join(out_dir, '%s_{/.}.domtbl' % get_basename(seqdb))
                optdicttring = optstring.join(optdict)
                hmmsearch_command = 'hmmsearch %s {} %s' % (optdicttring, seqdb)

                parahmm_command = 'find %s -name "*%s" | %s %s' % (db, opts.profile_suffix, parallel_command, hmmsearch_command)
                commands.append(parahmm_command)

            else:
                #Run each hmmsearch in a for loop
                profile_glob = os.path.join(db, '*%s' % opts.profile_suffix)
                for profile in list(glob.glob(profile_glob, recursive = False)):
                    hmmsearch_optdict['--domtbl'] = os.path.join(out_dir, '%s_%s.domtbl' % (get_basename(seqdb), get_basename(profile)))
                    hmmsearch_command = 'hmmsearch %s %s %s' % (optstring_join(optdict), profile, seqdb)
                    commands.append(hmmsearch_command)
    #Run hmmsearch for individual files
    elif dbfiles:
        if len(dbfiles) > 1:
            if jobs > 1:
                #temporary file containing paths for input HMM profiles
                path_file = os.path.join(opts.tmp_dir, 'profile_paths_%s' % now.strftime('%D-%M-%Y_%H:%M'))
                try:
                    with open(path_file, 'w') as pf:
                        for d in dbfiles:
                            pf.write('%s\n' % d)
                except IOError:
                    'Unable to write hmmsearch input profile file paths to %s' % path_file

                hmmsearch_joblog = os.path.join(output_paths['HMMER'], 'hmmsearch_joblog_%s_%s' % ('individual', now.strftime('%D-%M-%Y_%H:%M')))
                out_dir = os.path.mkdirs(os.path.join(output_paths['HMMER'], 'individual'))

                parallel_optdict['--joblog'] = hmmsearch_joblog
                parallel_optstring = optstring_join(parallel_optdict)
                parallel_command = 'parallel %s' % parallel_optstring

                hmmsearch_optdict['--domtbl'] = os.path.join(out_dir, '%s_{/.}.domtbl' % get_basename(seqdb))
                hmmsearch_optstring = optstring.join(optdict)
                hmmsearch_command = 'hmmsearch %s {} %s' % (hmmsearch_optstring, seqdb)

                parahmm_command = '%s %s :::: %s' % (parallel_command, hmmsearch_command, path_file)
                commands.append(parahmm_command)
            else:
                for profile in dbfiles:
                    hmmsearch_optdict['--domtbl'] = os.path.join(out_dir, '%s_%s.domtbl' % (get_basename(seqdb), get_basename(profile)))
                    hmmsearch_command = 'hmmsearch %s' % optstring_join(hmmsearch_optdict)
                    commands.append(hmmsearch_command)

        #Search one profile at a time
    elif len(dbfiles) == 1:
        hmmsearch_optdict['--domtbl'] = os.path.join(out_dir, '%s_%s.domtbl' % (get_basename(seqdb), get_basename(dbfiles[0])))
        hmmsearch_command = 'hmmsearch %s %s %s' % (optstring_join(hmmsearch_optdict), dbfiles[0], seqdb)
        commands.append(hmmsearch_command)

    else:
        pass
    return commands
#------------------------------------------------------------------------------
#BONEYARD
# def hmmsearch_command_generate(seqdb, profile, outdir, psuffix='.hmm', optdict, prefix*, parallel=False, jobs*, joblog*, progbar=False):
#     """
#     Generate a command string to run hmmsearch.
#     Use parallel=True to run parallel jobs. Useful for large numbers of profiles.
#     """
#     if not prefix:
#         prefix = get_basename(seqdb)
#     commands = []
#
#     if parallel and jobs == 1:
#         parallel = False
#     if not parallel:
#         if os.path.isdir(profile):
#             profile_glob = os.path.join(profile, '*%s' % (id_string, psuffix))
#             for i in list(glob.glob(profile_glob, recursive = False)):
#                 hmmsearch_command = 'hmmsearch %s %s %s' % (hmmsearch_optstring, i, seqdb)
#         elif os.path.isfile(profile):
#             #Run hmmsearch for a single profile
#             hmmsearch_optstring = optstring_join(optdict)
#             if not optstring['--domtbl']:
#                 domtbl = os.path.join(outdir, prefix + '.domtbl')
#             hmmsearch_command = 'hmmsearch %s %s %s' % (hmmsearch_optstring, profile, seqdb)
#     elif parallel and jobs > 1:
#         #Build the parallel string
#         parallel_optdict = {'--jobs':jobs}
#         if progbar:
#             parallel_optdict['--bar'] = ''
#         if not joblog:
#             from pathlib import Path
#             home = str(Path.home())
#             joblog = os.path.join(home, 'hmmsearch_joblog_' + now.strftime('%D-%M-%Y_%H:%M'))
#         parallel_optdict['--joblog'] = joblog
#         parallel_optstring = opstring_join(parallel_optdict)
#
#         if optdict['--domtbl']:
#             del optdict['--domtbl']
#         hmmsearch_optstring = optstring_join(optdict)
#         """
#         If a list of profile paths is provided
#         """
#         if isinstance(profile, list):
#             pfiles = [os.path.isfile(p) for p in profile]
#             pdirs = [os.path.isdir(p) for p in profile]
#             if len(pfiles) + len(pdirs) == len(profile):
#                     if pdirs:
#                         for pdir in pdirs:
#                             hmmsearch_command = 'find %s -type f -name "*%s" | parallel %s hmmsearch --domtbl %s/{/.}_%s.domtbl %s {} %s' % (pdir, psuffix, parallel_optstring, outdir, prefix, hmmsearch_optstring, seqdb)
#                             commands.append(hmmsearch_command)
#                     if pfiles and len(pfiles) > 1:
#                         hmmsearch_command = 'parallel %s hmmsearch --domtbl %s/{/.}_%s.domtbl %s {} %s ::: %s' % (pdir, psuffix, parallel_optstring, outdir, prefix, hmmsearch_optstring, seqdb, ' '.join(pfiles))
#                         commands.append(hmmsearch_command)

#------------------------------------------------------------------------------
