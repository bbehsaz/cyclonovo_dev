# def generateKmersIteratively():


def findkthelement(inputlist, N):
    if len(inputlist) == 0:
        return 0
    if len(inputlist) == 1:
        return inputlist[0]
    return sorted(inputlist, reverse=True)[min(len(inputlist) - 1, N)]


def findProperMAXKmer(
        candidate_kmers):  # It returns one of the highest scored kmers as representative. The one with most different amino acids is chosen randonmly
    num_uniq_aas = {}
    kmers_with_highest_numuniqaas = [0, []]
    for kmer in candidate_kmers:
        kmer_seq = kmer.split("-")
        num_uniq_aas = len(set(kmer_seq))
        if num_uniq_aas > kmers_with_highest_numuniqaas[0]:
            kmers_with_highest_numuniqaas[0] = num_uniq_aas
            kmers_with_highest_numuniqaas[1] = [kmer_seq]
        if num_uniq_aas == kmers_with_highest_numuniqaas[0]:
            kmers_with_highest_numuniqaas[1].append(kmer_seq)
    lowest_repeat_element = [20, ""]

    for kmer_seq in kmers_with_highest_numuniqaas[1]:
        repeat_element = max([kmer_seq.count(s) for s in set(kmer_seq)])
        if repeat_element < lowest_repeat_element[0]:
            lowest_repeat_element[0] = repeat_element
            lowest_repeat_element[1] = kmer_seq
    return "-".join(lowest_repeat_element[1])


def generate_All_Kmers(combination, kmerSize, pepMass, spectrumVector, kmer_score_thresh, e):
    # this function generates all the kmers using the amino acids in combination and calculaets their linear
    # score using an iterative method
    aaFreq = {}
    aminoFrequencies = {}
    protonMass = round(1.00728, 3)
    list_combination = list(combination)
    aminos_set = set(list_combination)

    sorted_aas = sorted(combination)
    aa2int = {sorted_aas[i]:i for i in range(len(sorted_aas))}


    def make_binary_addE(spectrumVector):
        binary_spectrumVector_e = [0] * 2000000
        for p in range(len(spectrumVector)):
            if spectrumVector[p] > 0:
                for x in range(int(p - (round(protonMass + e, 3) * 1000)), int(p - (round(protonMass - e, 3) * 1000))):
                    binary_spectrumVector_e[x] = 1
        return binary_spectrumVector_e

    binary_spectrumVector_e = make_binary_addE(spectrumVector)

    def linear_score_fragment(fragmentMass, spectrumVector, e):
        score = 0
        for x in range(int(round(fragmentMass + protonMass - e, 3) * 1000),
                       int(round(fragmentMass + protonMass + e, 3) * 1000)):
            if spectrumVector[x] > 0:
                return 1
        return 0

    def extend_to_k(k, kmerScores, kmerFrequencies, kmerSequences,
                    N):  # all-kmers with their linear scores using the k-1-mers in the kmers_scores
        kmerScores[k] = {}
        kmerFrequencies[k] = {}  # Frequency of all possible aa's in kmers with size k
        kmerSequences[k] = {}
        kmerscorethresh = findkthelement(kmerScores[k - 1].values(), N)
        for aa in aminos_set:
            # for k1mer in [x for x in kmerScores[k - 1] if kmerScores[k - 1][x] >= kmerscorethresh]:
            for k1mer in kmerScores[k - 1]:
                k1merSeq = list(kmerSequences[k - 1][k1mer])
                newKmer_sequence = k1merSeq + [aa]
                if len(newKmer_sequence) > 0:
                    if sum(newKmer_sequence) > pepMass:
                        continue
                if k1mer == '':
                    newKmer = str(aa)
                else:
                    newKmer = k1mer + "-" + str(aa)
                newScore = kmerScores[k - 1][k1mer]
                newScore += sum(
                    [binary_spectrumVector_e[int(round(sum(newKmer_sequence[i:]), 3) * 1000)] for i in range(k)])
                kmerScores[k][newKmer] = newScore
                kmerSequences[k][newKmer] = tuple(newKmer_sequence)
                kmerFrequencies[k][newKmer] = kmerFrequencies[k - 1][k1mer].copy()
                kmerFrequencies[k][newKmer][aa] +=1 
        newk = k + 1

        if newk == kmerSize + 1:
            return kmerScores, kmerFrequencies, kmerSequences
        return extend_to_k(newk, kmerScores, kmerFrequencies, kmerSequences, N)

    kmerScores = {}
    kmer_aminos = set()
    kmerFrequencies = {}
    kmerSequences = {}
    kmerScores[0] = {'': 0}
    kmerFrequencies[0] = {'': dict((aa, 0) for aa in aminos_set)}
    kmerSequences[0] = {'': []}
    N = 50
    all_kmerScores, all_kmerFrequencies, all_kmerSequences = extend_to_k(
        1, kmerScores, kmerFrequencies, kmerSequences, N)
    sorted_scores = sorted(all_kmerScores[kmerSize].values(), reverse=True)

    candidate_kmers = []
    if len(sorted_scores)>0:
        if sorted_scores[0] < kmer_score_thresh:
            return "NA", sorted_scores[0],{},{},{}
    else:
        return "NA", 0,{},{},{}
    #Choose a best-scoring k-mer to report. K-mer with the most number of distinct aa's and the highest score is selected
    for i in range(min(2, len(set(sorted_scores)))): 
        candidate_kmers += [t for t in all_kmerScores[kmerSize] if all_kmerScores[kmerSize][t] == sorted_scores[i]]
    selectedkmer = findProperMAXKmer(candidate_kmers)

    return all_kmerScores[kmerSize][selectedkmer], selectedkmer, all_kmerScores[kmerSize], all_kmerFrequencies[kmerSize],all_kmerFrequencies[kmerSize-1]


def kmerScoreTest(final_clusters_peaks, kmerSize, realPepMass, spectrumVector, kmerScoreIdent, e):
    building_blocks = final_clusters_peaks.keys()
    kmerSize = 5
    maxkmer_score, maxKmerSeq, kmers_scores, kmers_freqs, prefixKmer_freqs = generate_All_Kmers(building_blocks, kmerSize, realPepMass, spectrumVector, kmerScoreIdent, e)
    final_compound_basedon_kmerintesnity = "unclassified"
    if maxkmer_score == "NA":
        final_compound_basedon_kmerintesnity = "unclassified"
    elif maxkmer_score > kmerScoreIdent:
        final_compound_basedon_kmerintesnity = "cyclopeptide"
    return final_compound_basedon_kmerintesnity, maxkmer_score, maxKmerSeq,  kmers_scores, kmers_freqs,prefixKmer_freqs
