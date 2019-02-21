def classifySpectrum(
    spectrum_file, readStringMGF,checkConditions,generate_aa_convolutions_vectorized,output_cyclopeptide_polymers,verboseprint,polymer_repeat_units,
    cycloIntensityTest,kmerScoreTest,writeOriginalSpectra,
    building_blocks_main,params,protonMass,file_mapping,outputprefix,get_original_fpath, spectrum_scan,mgfFile,cyclopeptide_spectra_file
    ):
    import math
    from operator import itemgetter

    fixed_spectrum_file = spectrum_file + "END IONS\n"

    peptide = str(spectrum_scan)
    peaksnIntensity, pepMasses, charges, retentions, fileLines, allSpectraVector_dic = readStringMGF(
        fixed_spectrum_file, spectrum_scan)
    if len(peaksnIntensity) == 0:
        return -1
    verboseprint("=======================================================")
    verboseprint("Spectrum mass, rt, charge: {0},{1},{2}".format(pepMasses[peptide], retentions[peptide],
                                                                 charges[peptide]))

    ### Initiation
    
    thisbuildingblock = building_blocks_main.copy()
    standardAutconvCleaned = []
    precursorMass = pepMasses[peptide]
    # num_raw_spectra += 1
    realPepMass = pepMasses[peptide] * charges[peptide] - protonMass
    betta_values = [float(params["beta"])]
    alpha_values = [float(params["alpha"])]
    e = float(params['fragment_ion_thresh'])

    ### Remove spectrum with low num peaks and pepmass out of range
    if not checkConditions(realPepMass, len(peaksnIntensity[peptide])):
        verboseprint("Peptide mass is out of range or number peaks is too low!")
        return -1
        

    # else:
    #     num_spectra_analyzed += 1

    ##################################################################
    ### Analyze spectrum to identify cyclospectrum

    ### Generate convolutions for spectrum
    final_cyclopep_clusters = generate_aa_convolutions_vectorized(thisbuildingblock,
                                                                  allSpectraVector_dic[peptide],
                                                                  float(params["fragment_ion_thresh"]),
                                                                  realPepMass)
    final_cyclopep_convolutions_clusters = dict((aa, len(final_cyclopep_clusters[aa])) for aa in
                                                final_cyclopep_clusters)
    # final_monomers_thresh1 = dict((aa, len(final_cyclopep_clusters[aa])) for aa in
    #                                             final_cyclopep_clusters if len(final_cyclopep_clusters[aa])>int(params['kmer']))
    ### First Round of Cyclospectra Identification
    compound_type = "unclassified"
    for x in alpha_values:
        alpha = round(float(x), 3)
        for betta in betta_values:
            suffix = "_" + str(x) + "_" + str(betta)
            thresholdValue = round(math.ceil(alpha * realPepMass + betta))
            frequent_cyclopep_clusters = dict((k, final_cyclopep_convolutions_clusters[k]) for k in
                                              final_cyclopep_convolutions_clusters if
                                              final_cyclopep_convolutions_clusters[k] >= thresholdValue)
            if len(frequent_cyclopep_clusters) >= int(float(params['num_frequent_clusters'])):
                compound_type = output_cyclopeptide_polymers(
                    final_cyclopep_convolutions_clusters, final_cyclopep_clusters, polymer_repeat_units,
                    charges[peptide], thresholdValue)
            else:
                verboseprint("not enough frequent clusters")

    ### Second Round of Cyclospectra Identification

    # Testing if the cycloIntensity matches the cyclointensityThreshold
    compound_after_cyclointensity = "unclassified"

    if compound_type == "cyclopeptide":
        intensities = peaksnIntensity[peptide]
        convolThresh = int(params['aminoThresh'])
        final_clusters_peaks = dict((m, final_cyclopep_clusters[m]) for m in final_cyclopep_clusters if
                                    len(final_cyclopep_clusters[m]) > convolThresh)
        compound_after_cyclointensity, cycloIntensity_perc, polymerIntensity_perc, final_cyclopeptidic_monomers = cycloIntensityTest(
            final_clusters_peaks, intensities, polymer_repeat_units, convolThresh)
        final_cyclopep_monomers_convol = dict((m, final_cyclopep_convolutions_clusters[m]) for m in
                                              final_clusters_peaks if (m not in polymer_repeat_units and len(final_cyclopep_clusters[m])>int(params['aminoThresh'])))
        for x in final_cyclopep_monomers_convol:
            verboseprint("{0}\t{1}".format(x, final_cyclopep_monomers_convol[x]))
    else:

        cycloIntensity_perc, polymerIntensity_perc = 0, 0
    # CONTINUE FROM HERE
    if compound_after_cyclointensity == "cyclopeptide":
        final_compound_kmerScore, maxKmerScore, maxKmerSeq, dict_kmer_scores,kmers_freqs, prefixKmer_freqs = kmerScoreTest(final_cyclopep_monomers_convol, 5,
                                                                           realPepMass,
                                                                           allSpectraVector_dic[peptide],
                                                                           int(params['kmerScoreIdent']), e)
    else:
        final_compound_kmerScore, maxKmerScore, maxKmerSeq, dict_kmer_scores, kmers_freqs,prefixKmer_freqs = "unclassified", 0, "NA", dict(),dict(),dict()

    verboseprint(
        "cyclointensity,kmerscore,maxKmerSequence:\t{0}\t{1}\t{2}".format(cycloIntensity_perc, maxKmerScore,
                                                                          maxKmerSeq))
    final_compound = final_compound_kmerScore
    ### Write the selected Cyclospectra into the file
    if final_compound == "cyclopeptide":
        oneCyclopepFound = True
        
        # cyclopeptide_spectra_file = open(outputprefix + "cyclospectra.mgf", "a")
        writeOriginalSpectra(fileLines[peptide], cyclopeptide_spectra_file)
    if final_compound == 'polymer':
        num_polymer_spectra_dic[(alpha, betta)] += 1

    if final_compound != "cyclopeptide":
        final_compound = "NC"
    if final_compound == "cyclopeptide":
        building_blocks_multiplicity = ",".join(
            ["(" + str(round(key, 2)) + "," + str(value) + ")" for key, value in
             sorted(final_cyclopep_monomers_convol.items(), key=itemgetter(1), reverse=True)])
    else:
        building_blocks_multiplicity = []
        final_cyclopep_monomers_convol = {}
    
    info = [get_original_fpath(mgfFile, file_mapping), precursorMass, retentions[peptide], charges[peptide],
            final_compound,
            cycloIntensity_perc, maxKmerScore, maxKmerSeq,
            building_blocks_multiplicity]  # reports classification measure for each spectrum
    return info, realPepMass, final_cyclopep_monomers_convol,dict_kmer_scores,kmers_freqs,prefixKmer_freqs