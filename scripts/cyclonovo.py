#!/usr/bin/env python
#################################################
# CycloNovo: algorithm for high-throughput de novo cyclopeptide sequencing using tandem mass spectrometry
# Pevzner lab 
# University of California San Diego, La Jolla, CA, USA
#################################################

import errno
import os
import math as math
from os.path import basename
from utils.filter_apr1 import *
# from utils.denovo_benchmark_apr1 import *
from utils.denovoseq_cyclospectrum import *
from utils.filter_round2_kmers import *

# from scripts.read_string_mgf_vectored import readStringMGF
from scripts.read_string_split_mgf_vectored import readStringMGF

from utils.sharedFunctions import writeOriginalSpectra
from utils.sharedFunctions import getSpecializedBuildingBlocks
from utils.sharedFunctions import checkConditions

from utils.cyclointensity_test import *
from utils.spectral_convolutions_analysis import *
from utils.kmerscore_test import *
from utils.classify_spectrum import *

from operator import itemgetter
from subprocess import Popen, PIPE

from utils.spectra_utils import get_spectra_fpaths
from utils.file_utils import verify_file
from utils.common import parse_params_xml
from utils.file_utils import removeDir
from utils.file_utils import mkdir_p

def is_valid_file(parser, arg):
    """
    Check if arg is a valid file that already exists on the file system.
    Parameters
    """
    arg = os.path.abspath(arg)
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


total = 0


def get_parser():
    """Parse arguments and check if the spectrum file exists. """
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(prog="cyclonovo.py",
                            description="Algorith for analyzing cyclopeptides using Tandem Mass Spectra ...",
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument("-s", "--spectrum", dest="spectra", action='append',
                        help="Input spectra files and dirs with spectra files", metavar="FILE", required=True)
    parser.add_argument("-o", "--output", dest="output", help="Output dir", required=True)
    parser.add_argument("-p", "--params", dest="params", help="params.xml with all settings (GNPS running mode)")

    parser.add_argument("--pname", dest="pname",
                        help="Prefix for generated files", default="")

    parser.add_argument("--preprocess", dest="preprocess", action="store_true",
                        help='Filter and Preprocess the spectra (otherwise it\'s assumed it\'s preprocessed\')')
    parser.add_argument("-d", "--precursor_ion_thresh", dest='precursor_ion_thresh',
                        help='Precursor Ion Thresh (in Da)', default='0.02')
    parser.add_argument("-e", "--fragment_ion_thresh", dest='fragment_ion_thresh',
                        help='Fragment Ion Thresh (in Da)', default='0.02')
    parser.add_argument("--verbosity", dest="verbose", action="store_true")
    parser.add_argument("--monomers", dest="monomers", action="store",
                        help="Monomers to be considered ('standard','top25' or'ripps')",
                        choices=['standard', 'top25', 'ripps'], default='standard')

    # ARGUMENTS FOR CYCLOPEPTIDE IDENTIFICATION
    parser.add_argument("--alpha", dest='alpha',
                        help='Alpha value for calculating multiplicity thresholds',
                        default=float(0.0067), type=float)
    parser.add_argument("--beta", dest='beta',
                        help='Beta value for calculating multiplicity threshold',
                        default=float(-1), type=float)
    parser.add_argument("--kmer_score", dest='kmerScoreIdent',
                        help='k-mer score threshold for identifying cyclopeptidic spectra',
                        default=int(4), type=int)
    parser.add_argument("--cyclointensity", dest='cyclointensityThresh',
                        help='cycloIntensity threshold for identifying cyclopeptidic spectra',
                        default=int(60), type=int)
    parser.add_argument("--num_frequent_clusters", dest='num_frequent_clusters',
                        help='cycloIntesity threshold for identifying cyclopeptidic spectra',
                        default=int(2), type=int)
    parser.add_argument("--aminoacid_multiplicity", dest='aminoThresh',
                        help='Amino acid multiplicity threshold for de novo sequencing',
                        default=int(1), type=int)

    # ARGUMENTS FOR DE NOVO SEQUENCING
    parser.add_argument('--denovo', dest='denovo', help="De novo sequence cyclopeptide spectra", action='store_true')
    parser.add_argument("--kmer_size", dest='kmerSize',
                        help='k-mer size for building de Bruijn graph', default=int(5), type=int)
    parser.add_argument("--kmer_threshold", dest='kmerThreshold',
                        help='Threshold value for selecting k-mers for graph',
                        default=int(2), type=int)

    return parser


def get_params_and_mapping(opts):
    if opts.params is not None:
        verify_file(opts.params, description="params file")
        params, file_mapping = parse_params_xml(opts.params)
    else:
        params = opts.__dict__  # TODO: take a subset of all options
        file_mapping = None
    return params, file_mapping

def get_original_fpath(fpath, file_mapping=None):
    return file_mapping[basename(fpath)] if file_mapping and basename(fpath) in file_mapping else fpath


if __name__ == "__main__":
    args = get_parser().parse_args()
    cyclonovoPath = os.path.dirname(os.path.realpath(__file__))

    if args.verbose:
        def verboseprint(stringtoprint):
            print(stringtoprint)
    else:
        verboseprint = lambda *a, **k: None  # do-nothing function

    spectra_fpaths = get_spectra_fpaths(args.spectra)
    if not spectra_fpaths:
        verboseprint('No input spectra were found!')
        sys.exit(1)
    params, file_mapping = get_params_and_mapping(args)
    verboseprint("++ Arguments: " + str(args))
    verboseprint("++ Running CycloNovo from directory: " + cyclonovoPath)
    protonMass = 1.00728  # This is simply mass of hydrogen --- do not remove!
    digits = 3
    output = args.output
    outputprefix = os.path.join(output, args.pname)
    ### Monomers to use
    building_blocks_main, polymer_repeat_units = getSpecializedBuildingBlocks(params["monomers"], cyclonovoPath)
    removeDir(output)
    mkdir_p(output)
    summary_file = open(outputprefix + "summary.tsv", "w")
    summary_file.write("SpectrumFile\t#RawSpectra\t#Cyclospectra\t#SpectraSequenced\n")
    # info = [get_original_fpath(mgfFile, file_mapping), precursorMass, retentions[peptide], charges[peptide], peptide, final_compound,
    #   cycloIntensity_perc, maxKmerScore, maxKmerSeq, building_blocks_multiplicity]  # reports classification measure for each spectrum
    with open(outputprefix + "classification_scores.tsv", "w") as compound_type_reports_file:
        compound_type_reports_file.write(
            "SpectrumFile\tSpecMass\tRT\tCharge\tCyclopeptide?\tCyclointensity\tKmerScore\tAKmerWMAXScore\tCyclopeptidicClusters\n")
    for mgfFile in spectra_fpaths:
        verboseprint("++ Reading input spectrum: " + mgfFile)
        ### Read spectra with print_spectrum
        if params["preprocess"]:
            Process = Popen(
                ['print_spectrum ' + str(
                    mgfFile) + ' --max_charge 1 --scan_num -1 --print_spectrum ' +
                 '--configs_dir ' + cyclonovoPath],
                shell=True, stdin=PIPE, stdout=PIPE)
        else:
            Process = Popen(['print_spectrum ' + str(mgfFile) +
                             ' --scan_num -1 --print_spectrum --no_filter --no_merge --max_charge 1 ' + '--configs_dir ' + cyclonovoPath],
                            shell=True, stdin=PIPE, stdout=PIPE)
        peaksfile = Process.communicate()[0]

        ### Initiations
        peptides_considered = 0
        num_spectra_analyzed = 0
        num_cyclopeptide_spectra_dic = {}
        betta_values = [float(params["beta"])]
        alpha_values = [float(params["alpha"])]
        alpha = round(float(params["beta"]),3)
        beta = float(params["alpha"])
        e = float(params["fragment_ion_thresh"])
        num_cyclopeptide_spectra_dic[(alpha, beta)] = 0
        kValues = [int(args.kmerSize)]
        kthresholdValues = [int(args.kmerThreshold)]
        allidsconsidered = []
        allsequenced = 0
        representative = 'median'
        oneCyclopepFound = False

        ### starting the analysis one spectrum at a time
        split_input_spectra = peaksfile.split("END IONS")
        thisbuildingblock = building_blocks_main.copy()
        num_raw_spectra = len(split_input_spectra) - 1
        verboseprint("# of spectra in the file:\t{0}".format(num_raw_spectra))
        spectrum_scan = -1
        cyclopeptide_spectra_file = open(outputprefix + "cyclospectra.mgf", "a")
        for spectrum_file in split_input_spectra[:len(split_input_spectra) - 1]:
            if not spectrum_file:
                continue
            spectrum_scan += 1
            total += 1
            fixed_spectrum_file = spectrum_file + "END IONS\n"
            # ##################################################################
            # ### Analyze spectrum to identify cyclospectrum
            # print spectrum_file, readStringMGF,checkConditions,generate_aa_convolutions_vectorized,output_cyclopeptide_polymers,verboseprint,polymer_repeat_units,cycloIntensityTest,kmerScoreTest,writeOriginalSpectra,building_blocks_main,params,protonMass,file_mapping,outputprefix,get_original_fpath,spectrum_scan,mgfFile,cyclopeptide_spectra_file
            classify_spec_output = classifySpectrum(fixed_spectrum_file, readStringMGF,checkConditions,generate_aa_convolutions_vectorized,output_cyclopeptide_polymers,verboseprint,polymer_repeat_units,cycloIntensityTest,kmerScoreTest,writeOriginalSpectra,building_blocks_main,params,protonMass,file_mapping,outputprefix,get_original_fpath,spectrum_scan,mgfFile,cyclopeptide_spectra_file)
                # spectrum_scan,mgfFile,cyclopeptide_spectra_file)
            if classify_spec_output != -1:
                num_spectra_analyzed+=1

                info, realPepMass, dict_clusters_thresh1, dict_kmer_scores,kmers_freqs,prefixKmer_freqs = classify_spec_output
                if info[4] == "cyclopeptide":
                    num_cyclopeptide_spectra_dic[(alpha,beta)] +=1
                    oneCyclopepFound = True
            else:
                continue
            final_compound = info[4]

            if final_compound == "cyclopeptide":
                with open(outputprefix + "classification_scores.tsv", "a") as compound_type_reports_file:
                    compound_type_reports_file.write("\t".join(str(z) for z in info) + "\n")

            ##################################################################
            ### De novo sequence the cyclospectrum if the user specified so
                if final_compound == "cyclopeptide":
                    verboseprint("classified: " + str("cyclospectrum"))
                    if not args.denovo:
                        continue
                    benchmark_file = {}
                    reconstructions_file = {}
                    # for kmerSize in kValues:
                    #     for kmerThreshold in kthresholdValues:
                    kmerSize = int(params['kmerSize'])    
                    kmerThreshold= int(params['kmerThreshold'])   
                    benchmark_file[(kmerSize, kmerThreshold)] = open(
                    outputprefix + "sequencing_k" + str(kmerSize) + "_t" + str(kmerThreshold) + "_summary.txt", "a")
                    reconstructions_file[(kmerSize, kmerThreshold)] = open(
                    outputprefix + "sequencing_k" + str(kmerSize) + "_t" + str(kmerThreshold) + "_reconstructions.txt", "a")

                    selected_kmer = dict((kmer, kmers_freqs[kmer]) for kmer in dict_kmer_scores if
                                            dict_kmer_scores[kmer] > params['kmerThreshold'])
                    
                    building_blocks = []
                    for p in sorted(dict_clusters_thresh1.iteritems(), key=itemgetter(1), reverse=True):
                        if p[1] > args.aminoThresh:
                            # verboseprint("{0}\t{1}".format(round(p[0], 2), p[1]))
                            building_blocks.append(p[0])
                    if len(building_blocks) < 3:
                        verboseprint("Number of predicted cyclopeptidic amino acids is low")
                        break
                    candidateSequences = denovo_sequence(kmerSize, selected_kmer, realPepMass, building_blocks, e,e,kmerThreshold,verboseprint,prefixKmer_freqs)
                    output_denovo_results(candidateSequences, info[1], info[2], info[3],benchmark_file[(kmerSize, kmerThreshold)], 
                        kmerSize,kmerThreshold,reconstructions_file[(kmerSize, kmerThreshold)])

        nameofrecontfile = reconstructions_file[(kmerSize, kmerThreshold)].name
        nameofsummaryfile = benchmark_file[(kmerSize, kmerThreshold)].name
        # nameofcyclospecfile = output+"/cyclospectra.mgf"
        nameofcyclospecfile = cyclopeptide_spectra_file.name
        reconstructions_file[(kValues[0], kthresholdValues[0])].close()
        verboseprint("calculating P-values and cleaning up ...")
        if args.denovo:
            Popen(["sh " + cyclonovoPath + '/scripts/generate_pvalus_recontfile_mgf.sh ' + nameofrecontfile +
                   " " + nameofcyclospecfile + " " + str(output) + " " + str(cyclonovoPath) + " " + str(e) + " " + nameofsummaryfile],
                  shell=True)
        # if oneCyclopepFound:
        cyclopeptide_spectra_file.close()
        # report the number of cyclopeptides identified in output_summary.txt
        argstoprint = [get_original_fpath(mgfFile, file_mapping),
                       str(num_raw_spectra), num_cyclopeptide_spectra_dic[(alpha, beta)], allsequenced]
        for (alpha, betta) in num_cyclopeptide_spectra_dic:
            summary_file.write("\t".join([str(ap) for ap in argstoprint]) + "\n")
    summary_file.close()
    verboseprint("CycloNovo succesfully finished. Find the results in " + str(output))
