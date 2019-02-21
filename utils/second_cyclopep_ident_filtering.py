def calculateCyclopepIntensities(final_clusters_peaks, intensities, polymer_repeat_units):
    polymer_peaks = []
    aa_peaks = []
    twoDigit_intensities = {}
    for p in intensities:
        if round(p, 2) not in twoDigit_intensities:
            twoDigit_intensities[round(p, 2)] = intensities[p]
        else:
            twoDigit_intensities[round(p, 2)] += intensities[p]
    allintensity = sum(intensities.values())
    final_cyclopeptidic_monomers = []

    def calculate_internsity(peaks):
        intensity_sum = 0
        for peak in peaks:
            intensity_sum += twoDigit_intensities[round(peak, 2)]
        return intensity_sum

    for mass in final_clusters_peaks:
        peak_pairs = []
        aa_pairs = []
        for group in final_clusters_peaks[mass]:
            peak_pairs += final_clusters_peaks[mass][group]
        if mass in polymer_repeat_units:
            for pair in peak_pairs:
                polymer_peaks += list(pair)
        else:
            final_cyclopeptidic_monomers.append(mass)
            for pair in peak_pairs:
                aa_peaks += list(pair)

    polymer_intensity = calculate_internsity(list(set(polymer_peaks)))
    aa_intensity = calculate_internsity(list(set(aa_peaks)))

    percentProteinoPeaksIntens = aa_intensity / float(allintensity)

    percentPolymerPeaksIntens = polymer_intensity / float(allintensity)
    return percentProteinoPeaksIntens, percentPolymerPeaksIntens, final_cyclopeptidic_monomers


def cycloIntensityTest(
        final_clusters_peaks, intensities, polymer_repeat_units, cycloinsity_thresh):
    percentProteinoPeaksIntens, percentPolymerPeaksIntens, final_cyclopeptidic_monomers = calculateCyclopepIntensities(
        final_clusters_peaks, intensities, polymer_repeat_units)
    if percentProteinoPeaksIntens > 0.6:
        final_compound = "cyclopeptide"
        if percentProteinoPeaksIntens < 0.7 and percentPolymerPeaksIntens > 0.35:
            percentProteinoPeaksIntens = min(percentProteinoPeaksIntens - 0.1, 0)
            final_compound = "polymer"
    return final_compound, percentProteinoPeaksIntens, percentPolymerPeaksIntens, final_cyclopeptidic_monomers
# if percentProteinoPeaksIntens > percentPolymerPeaksIntens:
# 	final_compound = 'polymer'
# else:
# 	percentProteinoPeaksIntens = min(percentPolymerPeaks,percentProteinoPeaksIntens)

# def generateKmersIteratively():


# def generate_All_Kmers(building_blocks,kmerSize,realPepMass,peaks,e): 
# #this function generates all the kmers using the amino acids in combination and calculaets their linear 
# #score using an iterative method 
# 	aaFreq = {}
# 	aminoFrequencies = {}

# 	list_combination = list(combination)
# 	aminos_set = set(list_combination)
# 	for x in aminos_set:
# 		aminoFrequencies[x] = list_combination.count(x)
# 	def extend_to_k(k,kmerScores, kmerFrequencies, kmerSequences): #all-kmers with their linear scores using the k-1-mers in the kmers_scores
# 		kmerScores[k] = {}
# 		kmerFrequencies[k] = {} #Frequency of all possible aa's in kmers with size k
# 		kmerSequences[k] = {}

# 		for aa in aminos_set:
# 			for k1mer in kmerScores[k-1]:
# 				k1merSeq = kmerSequences[k-1][k1mer][:]
# 				newKmer_sequence = k1merSeq + [aa]
# 				if len(newKmer_sequence)>0:
# 					if sum(newKmer_sequence) >realPepMass:
# 						continue
# 				if k1mer == '':
# 					newKmer = str(aa)
# 				else:
# 					newKmer = k1mer + "-" + str(aa)
# 				newScore =kmerScores[k-1][k1mer][:]
# 				newScore += linear_score_fragments(
# 					[newKmer_sequence[i:] for i in range(0,k)], peaks, e)
# 				kmerScores[k][newKmer] = newScore
# 				kmerFrequencies[k][newKmer] = {}
# 				for x in kmerFrequencies[k-1][k1mer]:
# 					kmerFrequencies[k][newKmer][x] = kmerFrequencies[k-1][k1mer][x]
# 				kmerFrequencies[k][newKmer][aa] = kmerFrequencies[k-1][k1mer][aa] + 1
# 				kmerSequences[k][newKmer] = newKmer_sequence
# 		newk = k+1

# 		if newk== kmerSize + 1:
# 			return kmerScores, kmerFrequencies , kmerSequences
# 		return extend_to_k(newk, kmerScores, kmerFrequencies, kmerSequences)

# 	kmerScores = {}
# 	kmerFrequencies = {}
# 	kmerSequences = {}
# 	kmerScores[0] = {'':[]}
# 	kmerFrequencies[0] = {'': dict((aa, 0) for aa in aminos_set)}
# 	kmerSequences[0] = {'':[]}
# 	all_kmerScores, all_kmerFrequencies , all_kmerSequences  = extend_to_k(
# 		1, kmerScores, kmerFrequencies, kmerSequences)
# 	final_allKmerScores = {}
# 	kmerMatches = all_kmerScores[kmerSize].copy()

# 	for size in all_kmerScores:
# 		final_allKmerScores[size] = {}
# 		for kmer in all_kmerScores[size]:
# 			final_allKmerScores[size][kmer] = len(( [x for x in all_kmerScores[size][kmer]] ))
# 	return final_allKmerScores, all_kmerFrequencies , all_kmerSequences, kmerMatches


# def kmerScoreTest(final_clusters_peaks, allSpectraVector_dic):
# 	building_blocks = final_clusters_peaks.keys()
# 	# building_blocks = [key for key, value in standardAutconvCleaned.iteritems() if value>1]
# 	kmerSize = 5
# 	all_kmerScores, all_kmerFrequencies , all_kmerSequences, all_kmerMatches = generate_All_Kmers(
# 		building_blocks,kmerSize,realPepMass,intensities.keys(),e)
# 	highetScoredKmer = ("",0)
# 	sorted_kmers = sorted(all_kmerScores[5].items(), key=itemgetter(1),reverse=True)
# 	if len(all_kmerScores[5])>0:
# 		highetScoredKmer = sorted_kmers[0]
# 	topkmer = ""
# 	if highetScoredKmer[1]>3 and percentProteinoPeaks>0.6:
# 		final_compound = "cyclopeptide"
# 		topkmer = highetScoredKmer[0]
# 		if percentProteinoPeaks <0.7 and percentPolymerPeaks>0.35:
# 						final_compound="polymer"
# 	if topkmer != "":
# 		roundedKmer = "-".join([str(round(float(km),2)) for km in topkmer.split("-")])
# 	else:
# 		roundedKmer = "-"
