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
            if round(peak, 2) in twoDigit_intensities:
                intensity_sum += twoDigit_intensities[round(peak, 2)]
            elif round(peak, 2) - 0.01 in twoDigit_intensities:
                intensity_sum += twoDigit_intensities[round(peak - 0.01, 2)]
            elif round(peak + 0.01, 2) in twoDigit_intensities:
                intensity_sum += twoDigit_intensities[round(peak + 0.01, 2)]
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
    final_compound = "unclassified"
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
