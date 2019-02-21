# def generate_convolutions(
# 	standardAminoMasses,peaks,e,thresh,charge,representative, precursorMass):
# 	constituentMonoMers = []
# 	finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge = generateAllAutconv(peaks,e,thresh,representative,standardAminoMasses,precursorMass,charge)
# 	return finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge, standardAminoMasses, e

def median(lst):
    quotient, remainder = divmod(len(lst), 2)
    if remainder:
        return sorted(lst)[quotient]
    return sum(sorted(lst)[quotient - 1:quotient + 1]) / 2.


def generate_aa_convolutions_vectorized(
        standardAminoMasses, peaks, e, realPepMass):
    sattelites_removed_linked_cyclopep_clusters = {}
    for offset in sorted(standardAminoMasses):
        sattelites_removed_linked_cyclopep_clusters[offset] = convolution(peaks, peaks, offset, e, realPepMass)
    finalconvols = {}
    clustermass_sattelites_removed_linked_cyclopep_clusters = {}
    for aa in sattelites_removed_linked_cyclopep_clusters:
        sumofalldists = 0
        totalnumalldists = 0
        alldistances = []
        for group in sattelites_removed_linked_cyclopep_clusters[aa]:
            if len(sattelites_removed_linked_cyclopep_clusters[aa]) > 0:
                distances = [
                    min([abs(pair[1] - pair[0]) for pair in sattelites_removed_linked_cyclopep_clusters[aa][group]],
                        key=lambda x: abs(aa - x))]
                # distances = [abs(pair[1]-pair[0]) for pair in sattelites_removed_linked_cyclopep_clusters[aa][group]]
                alldistances += distances
                sumofalldists += sum(distances)
                totalnumalldists += len(distances)
        if totalnumalldists != 0:
            # clustermass = sumofalldists/(len(alldistances)*1.0)
            clustermass = median(alldistances)
        else:
            clustermass = 0
        if abs(clustermass - aa) < e:
            clustermass_sattelites_removed_linked_cyclopep_clusters[aa] = sattelites_removed_linked_cyclopep_clusters[
                aa]
        else:
            clustermass_sattelites_removed_linked_cyclopep_clusters[aa] = {}

    return clustermass_sattelites_removed_linked_cyclopep_clusters


def convolution(spectrum1, spectrum2, offset, e, realPepMass):
    convolutionPairPeaks_list = []

    def make_binary_addE(spectrumVector):
        binary_spectrumVector_e = [0] * 2000000
        for p in range(len(spectrumVector)):
            if spectrumVector[p] > 0:
                for x in range(int(p - (5 * round(e, 3) * 1000)), int(p + (5 * round(e, 3) * 1000))):
                    binary_spectrumVector_e[x] = 1
        return binary_spectrumVector_e

    for peak1 in range(min(int(realPepMass * 1000), len(spectrum1))):
        if spectrum1[peak1] > 0:
            for peak2 in range(max(int(peak1 + (offset - 5 * e) * 1000), 0),
                               min(2000000, int(peak1 + (offset + 5 * e) * 1000) + 1)):
                if spectrum2[peak2] > 0:
                    convolutionPairPeaks_list.append(
                        (round(float(peak1) / 1000.0, 3), round(float(peak2) / 1000.0, 3)))  # peakpair
                    break

    # for peak1 in [peak for peak in range(min(int(realPepMass*1000), len(spectrum1))) if spectrum1[peak]>0]:
    # 	peak2 = peak1+int(round(offset,3)*1000)
    # 	if spectrum2_vector_5e[peak2]> 0:
    # 		convolutionPairPeaks_list.append(( round(float(peak1)/1000.0,3),round(float(peak2)/1000.0,3) ) ) #peakpair

    redundancy_removed = {}

    current_group = 0
    convolutionPairPeaks_noRedundancy = {}
    pair_groups = {}
    for pair in convolutionPairPeaks_list:
        new_group = True
        if current_group == 0:
            current_group += 1
            minmax_x = [pair[0], pair[0]]
            minmax_y = [pair[1], pair[1]]
            convolutionPairPeaks_noRedundancy[current_group] = [pair]
            pair_groups[pair] = current_group
            continue
        adduct = 28
        if minmax_x[0] - adduct < pair[0] < minmax_x[1] + adduct:
            if minmax_y[0] - adduct < pair[1] < minmax_y[1] + adduct:
                new_group = False
        if not new_group:
            if pair[0] < minmax_x[0]:
                minmax_x[0] = pair[0]
            if pair[0] > minmax_x[1]:
                minmax_x[1] = pair[0]

            if pair[1] < minmax_y[0]:
                minmax_y[0] = pair[1]
            if pair[1] > minmax_y[1]:
                minmax_y[1] = pair[1]
            convolutionPairPeaks_noRedundancy[current_group].append(pair)
            pair_groups[pair] = current_group
        else:
            current_group += 1
            minmax_x = [pair[0], pair[0]]
            minmax_y = [pair[1], pair[1]]
            convolutionPairPeaks_noRedundancy[current_group] = [pair]
            pair_groups[pair] = current_group
            continue
    min_element = offset
    max_element = offset
    alldistances = {}
    final_group = {}
    for pair in convolutionPairPeaks_list:
        if round(abs(pair[1] - pair[0]), 3) not in alldistances:
            alldistances[round(abs(pair[1] - pair[0]), 3)] = [pair]
        else:
            alldistances[round(abs(pair[1] - pair[0]), 3)].append(pair)
    # if  abs(round(abs(pair[1]-pair[0]),3) - offset) < e:
    # 	final_group.append(round(abs(pair[1]-pair[0]),3))
    sorted_alldistances = sorted(alldistances)
    last_distance = 0
    aa_group = False
    last_group_status = False
    final_group_distanes = []
    for i in range(len(sorted_alldistances)):
        if i == len(sorted_alldistances) - 1:
            if last_group_status:
                if abs(distance - offset) > e:
                    break
                else:
                    final_group_distanes.append(distance)
                    break
        distance = sorted_alldistances[i]
        if abs(distance - offset) < e:
            aa_group = True
        if abs(distance - last_distance) > e:
            if last_group_status:
                if abs(distance - offset) > e:
                    break
            final_group_distanes = [distance]

        else:
            final_group_distanes.append(distance)
        last_distance = distance
        last_group_status = aa_group
    final_groups = {}

    for distance in final_group_distanes:

        for pair in alldistances[distance]:
            final_groups[pair_groups[pair]] = convolutionPairPeaks_noRedundancy[pair_groups[pair]]

    return final_groups


def find_repeat_chains_longerThenx(offset_pairs, x):
    all_xcords = [pair[0] for pair in offset_pairs]
    all_ycords = [pair[1] for pair in offset_pairs]

    def all_indices(value, qlist):
        indices = []
        idx = -1
        while True:
            try:
                idx = qlist.index(value, idx + 1)
                indices.append(idx)
            except ValueError:
                break
        return indices

    def intersect(a, b):
        """ return the intersection of two lists """
        return list(set(a) & set(b))

    def extendChainOneStep(start_points, all_shared_values, all_xcords, all_ycords):
        next_stepX = []
        for ind in start_points:
            yvalue = all_ycords[ind]
            if yvalue in all_shared_values:
                next_stepX.append(yvalue)
        return list(set(next_stepX))

    original_shared_values = intersect(all_xcords, all_ycords)
    maxchain = 2
    last_shared_values = original_shared_values[:]
    while maxchain < x and len(last_shared_values) > 0:
        start_points = []
        new_shared_values = []
        for value in last_shared_values:
            start_points += all_indices(value, all_xcords)

        new_shared_values = extendChainOneStep(start_points, original_shared_values, all_xcords, all_ycords)
        last_shared_values = new_shared_values[:]
        if len(new_shared_values) > 0:
            maxchain += 1
        else:
            break
    if maxchain == x:
        return 1
    else:
        return 0


def output_cyclopeptide_polymers(final_clusters_convolutions, final_clusters_peaks, polymerMasses, charge,
                                 thresholdValue):
    totalCyclopeptide = 0
    compound = "nonpeptidic"
    polymer = False
    # frequent_cyclopep_clusters = {k:final_clusters_convolutions[k] for k in final_clusters_convolutions if final_clusters_convolutions[k]>=thresholdValue}
    from operator import itemgetter
    # if frequent_cyclopep_clusters >-1:
    sorted_convolutions = sorted(final_clusters_convolutions.items(), key=itemgetter(1), reverse=True)
    final_STNDconvolutions = []  # this list will hold all the final convolutions that will contribute to the final convolutions with filtered convoltuions
    checkMasses = []
    numPolymers = []

    for mass in [x[0] for x in sorted_convolutions if x[1] > thresholdValue]:
        if mass in polymerMasses:
            if mass > 26:
                numPolymers.append(mass)
        elif int(mass) not in set(checkMasses):
            if int(mass) == 87:
                continue
            if charge == 2 and (2 * int(mass) in set(checkMasses) or int(mass) / 2 in set(checkMasses)):
                continue
            final_STNDconvolutions.append(mass)
            # checkMasses += [int(mass)-1,int(mass),int(mass)+1, int(mass)+2, int(mass)-2]
            checkMasses += [int(mass) - 1, int(mass), int(mass) + 1]
    if len(numPolymers) > 1:
        if len(final_STNDconvolutions) < 3 * len(numPolymers) + 1:
            polymer = True
    if len(final_STNDconvolutions) < 2:
        compound = "polymer"
    if len(final_STNDconvolutions) < 4 and (not polymer):  # checking if this is polymer of a single amino acid
        found_chains = 0
        chains = []
        for mass in final_STNDconvolutions:
            peak_pairs = []
            for group in final_clusters_peaks[mass]:
                peak_pairs += final_clusters_peaks[mass][group]
            # max_chain =  find_repeat_chain([peak_pairs)
            # chains.append(max_chain)
            found_chains += find_repeat_chains_longerThenx(peak_pairs, 4)
        # if max_chain>3:
        # 	found_chains +=1
        if found_chains > 1:
            compound = "polymer"
        else:
            compound = "cyclopeptide"

    if len(final_STNDconvolutions) > 3 and not polymer:
        totalCyclopeptide += 1
        compound = "cyclopeptide"
    elif len(final_STNDconvolutions) > 3 and polymer:
        found_chains = 0
        n = 0
        chains = []
        for mass in set(numPolymers):
            n += 1
            peak_pairs = []
            for group in final_clusters_peaks[mass]:
                peak_pairs += final_clusters_peaks[mass][group]
            found_chains += find_repeat_chains_longerThenx(peak_pairs, 4)
        # if max_chain>3:
        # 	found_chains +=1
        # find_repeat_chains_longerThenx(mass ,x)

        if found_chains > 1:
            compound = "polymer"
        else:
            compound = "cyclopeptide"

    elif len(numPolymers) > 1:
        compound = "polymer"

    else:
        compoud = "unclassified"

    return compound
