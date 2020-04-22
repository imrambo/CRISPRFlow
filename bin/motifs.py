#!/usr/bin/python
methylation_motifs = {'GATC':r'GATC', 'GANTC':r'GA[ACTG]TC', 'CCGG':r'CCGG', 'GGCC':r'GGCC'}

def sliding_window(sequence, window_extent):
    """
    Brute-force sliding window for calculating pattern frequency.
    """
    raw_data = {motif:[] for motif, pattern in methylation_motifs}
    raw_data['position'] = []
    raw_data['GC'] = []

    position = window_extent
    while position <= len(sequence) - window_extent:
        raw_data['position'].append(position)
        window = sequence[position-window_extent:position+window_extent]

        nucleotides = nucleotide_calculator(window)
        gc = float(((nucleotides['G'] + nucleotides['C'])/len(window))*100)
        raw_data['GC'].append(gc)

        for motif, pattern in methylation_motifs:
            motif_position = re.search(pattern, window).start()
            motif_score = 0
            while motif_position > 0:
                motif_score += float((1/(1 + abs((len(window)/2) - motif_position))))
        position += 1


def nucleotide_calculator(sequence):
    nuc_counts = {}
    nuc_counts['A'] = sequence.count('A')
    nuc_counts['C'] = sequence.count('C')
    nuc_counts['T'] = sequence.count('T')
    nuc_counts['G'] = sequence.count('G')

    return nuc_counts
