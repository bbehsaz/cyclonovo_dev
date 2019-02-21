
def generate_feasible_combinations(building_blocks, pepLengths, pepMass, delta):
    # Generates all feasible combinations for lengths that
    failedLenths = []
    allFeasibleCombinations = {}
    for l in pepLengths:
        allFeasibleCombinations[l] = set()
        allFeasibleCombinations[l] = find_feasible_combintions(
            building_blocks, l, pepMass, delta)
        if len(allFeasibleCombinations) == 0:
            failedLenths.append(l)
    return allFeasibleCombinations, failedLenths

def __combinations_with_replacement(iterable, r):
    # combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
    pool = tuple(iterable)
    n = len(pool)
    if not n and r:
        return
    indices = [0] * r
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != n - 1:
                break
        else:
            return
        indices[i:] = [indices[i] + 1] * (r - i)
        yield tuple(pool[i] for i in indices)

def find_feasible_combintions(building_blocks, l, peptideMass, delta):
    # Among all possible combinations of size l of the selected autconvolutions (building_blocks) returns the ones that
    # their sum matches the precursor mass allowing for delta error
    listfeasibleCombinations = set()
    from operator import itemgetter
    import itertools
    protonMass = 1.00728
    for combination in __combinations_with_replacement(building_blocks, l):
        if abs(sum(combination) - (peptideMass)) < delta:
            listfeasibleCombinations.add(tuple(sorted([float(x) for x in combination])))


    return listfeasibleCombinations

def select_kmers_per_composition(composition_freq,kmers):
    kmer_set= set()   
    for  kmer in kmers:
        include = True
        for aa in kmers[kmer]:
            if aa not in composition_freq:
                if kmers[kmer][aa] >0:
                    include = False
                    break
            else:
                if kmers[kmer][aa] > composition_freq[aa]:
                    include= False
                    break
        if include:
            kmer_set.add(kmer)
    return kmer_set


def select_kmers_compositions(all_compositions,kmers):
    all_selected_kmers= dict()
    compositions_freq = dict()
    for composition in all_compositions:
        composition_freq = {aa_c:composition.count(aa_c) for aa_c in composition}
        all_selected_kmers[tuple(composition)] = select_kmers_per_composition(composition_freq,kmers)    
        compositions_freq[tuple(composition)] = composition_freq
    return all_selected_kmers, compositions_freq


def find_feasible_composition(building_blocks, pepLengths, pepMass, delta):
# Generates all feasible combinations for found lengths
    failedLenths = []
    allFeasibleCombinations = {}
    for l in pepLengths:
        allFeasibleCombinations[l] = set()
        allFeasibleCombinations[l] = find_feasible_combintions(
                                        building_blocks, l, pepMass, delta)
        if len(allFeasibleCombinations) == 0:
            failedLenths.append(l)
    return allFeasibleCombinations, failedLenths



def constructDeBruijn(kmers, kmerSize,prefixKmer_freqs):
    # Constructs de Bruijn graph from k-mers with k-1
    adjacencies = {}
    nodes2aminos = {}
    n = 0
    edges = []
    nodes = {}
    def suffix(string):
        return string[1:]

    def prefix(string):
        return string[:len(string) - 1]
    # To get the deburijn graph that actually nodes are the k-mer size ... only for plotting purposes!
    for kmer in kmers:
        # add prefix and suffix nodes
        kmer_tuple = tuple([float(a) for a in kmer.split("-")])

        #adding prefix node
        prefixK1merSeq = prefix(kmer_tuple)
        prefixK1merName = "-".join(str(aa) for aa in prefix(kmer_tuple))
        prefixK1merFrequencies = prefixKmer_freqs[prefixK1merName]
        # prefixK1merFrequencies = dict((a, prefixK1merSeq.count(a)) for a in set(prefixK1merSeq))
        nodes[prefixK1merName] = (prefixK1merSeq, prefixK1merFrequencies)
        #adding suffix node
        suffixK1merSeq = suffix(kmer_tuple)
        suffixK1merName = "-".join(str(aa) for aa in suffixK1merSeq)
        suffixK1merFrequencies = prefixKmer_freqs[suffixK1merName]
        # suffixK1merFrequencies = dict((a, suffixK1merSeq.count(a)) for a in set(suffixK1merSeq))
        nodes[suffixK1merName] = (suffixK1merSeq, suffixK1merFrequencies)

        if prefixK1merName in adjacencies:
            adjacencies[prefixK1merName].append(suffixK1merName)
            edges.append((prefixK1merName, suffixK1merName))
        else:
            adjacencies[prefixK1merName] = [suffixK1merName]
            edges.append((prefixK1merName, suffixK1merName))

    return adjacencies, nodes, edges

def denovo_sequence(kmerSize, kmers, realPepMass, building_blocks, e, delta, kmerThreshold, verboseprint,prefixKmer_freqs):
    # De novo sequence spectrum with mass pepMass using amino acids in building blocks

    lengths = choose_length(realPepMass)
    allcompos = []

    print "WWHATTGTTTTASDFSDFASDFAAT"
    # exit()
    protonMass = 1.00728  
    sorted_monomers = sorted(building_blocks)
    mapAmino2Int = {str(sorted_monomers[i]):i for i in range(len(building_blocks))}

    allfeasible_compositions, failedLenths = find_feasible_composition(building_blocks, lengths, realPepMass, delta)
    for l in allfeasible_compositions:
        if l < kmerSize + 1:
            continue
        allcompos += list([sorted(x) for x in allfeasible_compositions[l]])
    kmers_per_composition, compositions_freq = select_kmers_compositions(allcompos,kmers)
    # print kmers_per_composition
    checkedSeqs = set()
    for composition in kmers_per_composition:
        verboseprint("Checking amino acid combination:\t{0}".format(composition))

        verboseprint("# of Putative Kmers: {}".format(len(kmers_per_composition[composition])))
        l = len(composition)
        adjacencies, nodes, edges = constructDeBruijn(
                                                    kmers_per_composition[composition],kmerSize,prefixKmer_freqs
                                                    )  # Build the graph

        adjacencies, nodes, edges = pruneGraph(adjacencies, nodes, edges)
        listAllPaths = []
        for node in nodes:
            if (realPepMass + protonMass) - sum(nodes[node][1]) > -e:
                listAllPaths.extend(findWalksLengthK(composition, adjacencies, node, l, nodes, building_blocks))
        allPaths = set(listAllPaths)
        if len(allPaths) < 1:
            verboseprint("No feasible cycle found!")
            continue
    
        for path in allPaths:
            for x in composition:
                if  compositions_freq[composition][x]!= path.count(x):
                    continue
            # if not pathIsClosed(list(path), chosenKmers, kmerSize):
            #     continue
            checked = False

            for permutation in relatedPermutations(list(path)):
                if permutation in checkedSeqs:
                    checked = True
                    break
            if checked:
                continue
            else:
                checkedSeqs.add(path)
        if len(checkedSeqs) > 0:
            verboseprint("Number of feasible cycles: {0}".format(len(checkedSeqs)))
        else:
            verboseprint("No feasible cycle found!")
    return checkedSeqs



def output_denovo_results(candidatePeptides, precursorMass, retention, charge, benchmark_file, kmerSize,
                          kmerThreshold, reconstructions_file):
    # This function just the final reconstructions in the format wanted! This one is for benchamrking CYCLOLIBRARY
    from operator import itemgetter

    # reconstruction_sorted = sorted(candidatePeptides.items(), key=itemgetter(1), reverse=True)
    scores = set()

    peptidesUptoCorrectScore = []
    correctFound = False
    # for reconstruction in candidatePeptides:
    score= 20 

    for reconstruction in candidatePeptides:
        # scores.add(score)
        # if len(scores) > 4:
        #     break
        # elif score > 12:
        #     peptidesUptoCorrectScore.append(reconstruction)
        #     scores.add(score)
        if score > 11:
            reconstructions_file.write(
                "{0}\t{1}\t{2}\t{3}\t{4}\n".format(",".join([str(x) for x in reconstruction]), precursorMass, retention,
                                                   charge, score))


    # allMatches = {}
    # for k in range(1, kmerSize):
    #     currentMatches = matchAllKmers(building_blocks, k, peaks, e)
    #     allMatches.update(currentMatches)
    # # } CHanging back
    # for l in allFeasibleCombinations:
    #     if l < kmerSize + 1:
    #         continue
    #     allcompos += list([sorted(x) for x in allFeasibleCombinations[l]])
    # verboseprint("Checking feasible combinations")
    # kmerScores, kmerFrequencies, kmerSequences = {}, {}, {}

    # candidatePeptides = {}  # tuples (peptide, score)

    # for l in [x for x in allFeasibleCombinations if len(allFeasibleCombinations[x]) > 0]:
    #     for combination in allFeasibleCombinations[l]:
    #         verboseprint("Checking amino acid combination:\t{0}".format(combination))
    #         all_kmerFrequencies = {}
    #         kmers = {}
    #         all_kmerSequences = {}
    #         kmers, all_kmerFrequencies[kmerSize], all_kmerSequences[kmerSize] = generateAllKmerLinearScore(allMatches,
    #                                                                                                        list(
    #                                                                                                            combination),
    #                                                                                                        kmerSize,
    #                                                                                                        pepMass)
    #         chosenKmers = dict((key, value) for key, value in kmers.items() if value > kmerThreshold)

    #         # verboseprint("#Putative Kmers: {}".format(len(chosenKmers)))
    #         adjacencies, nodes, edges = constructDeBruijn(chosenKmers, all_kmerFrequencies, all_kmerSequences,
    #                                                       kmerSize)  # Build the graph
    #         adjacencies, nodes, edges = pruneGraph(adjacencies, nodes, edges)
    #         listAllPaths = []
    #         for node in nodes:
    #             if (pepMass + protonMass) - sum(nodes[node][1]) > -e:
    #                 listAllPaths.extend(findWalksLengthK(combination, adjacencies, node, l, nodes, building_blocks))
    #         allPaths = set(listAllPaths)

    #         if len(allPaths) < 1:
    #             verboseprint("No feasible cycle found!")
    #             continue
    #         checkedSeqs = set()

    #         for path in allPaths:
    #             for x in combination:
    #                 if combination.count(x) != path.count(x):
    #                     continue
    #             if not pathIsClosed(list(path), chosenKmers, kmerSize):
    #                 continue
    #             checked = False

    #             for permutation in relatedPermutations(list(path)):
    #                 if permutation in checkedSeqs:
    #                     checked = True
    #                     break
    #             if checked:
    #                 continue
    #             else:
    #                 checkedSeqs.add(path)
    #         numfeasibleCycles = 0
    #         for path in checkedSeqs:
    #             matchedPeaksFragments = findCyclicScore(path, kmerSize, peaks, e)
    #             score = len(set([x[0] for x in matchedPeaksFragments]))
    #             candidatePeptides[tuple(path)] = score
    #             if score > 11:
    #                 # verboseprint("{}\t{}".format(path, score))
    #                 numfeasibleCycles += 1
    #         if numfeasibleCycles > 0:
    #             verboseprint("Number of feasible cycles: {0}".format(numfeasibleCycles))
    #         else:
    #             verboseprint("No feasible cycle found!")
    # return candidatePeptides


def pruneGraph(adjacencies, nodes, edges):
    backwardEdges = {}
    leaves = {}
    newleaves = []
    singletons = []
    for node in nodes:
        backwardEdges[node] = []
    for node in nodes:
        if node not in adjacencies:
            adjacencies[node] = []
            continue
        for tail in adjacencies[node]:
            backwardEdges[tail].append(node)

    for node in nodes:
        if len(adjacencies[node]) == 0:
            if node not in backwardEdges:
                singletons.append(node)
    [nodes.remove(x) for x in singletons]
    for node in nodes:
        if node not in adjacencies:
            adjacencies[node] = []
    # singletons = set(nodes) - set(set(nodes) - set(adjacencies.keys()) )

    newleaves = findLeaves(adjacencies)
    while len(newleaves) > 0:

        leaves = newleaves[:]
        newleaves = []

        for leaf in leaves:

            for parent in backwardEdges[leaf]:
                # parent = backwardEdges[leaf][0]
                # if len(adjacencies[parent]) ==1:
                #   newleaves.append(parent)
                adjacencies[parent].remove(leaf)
                edges.remove((parent, leaf))
            # backwardEdges[leaf].remove(parent)
            del backwardEdges[leaf]
            del adjacencies[leaf]
            del nodes[leaf]
        newleaves = findLeaves(adjacencies)

    return adjacencies, nodes, edges
def choose_length(pepMass):
    if pepMass > 1100:
        pepLengths = [9, 10, 11]
    if pepMass > 950:
        pepLengths = [8, 9, 10]
    elif pepMass > 800:
        pepLengths = [7, 8, 9]
    elif pepMass > 700:
        pepLengths = [6, 7, 8, 9]
    elif pepMass > 600:
        pepLengths = [5, 6, 7]
    elif pepMass > 500:
        pepLengths = [4, 5, 6, 7]
    elif pepMass > 400:
        pepLengths = [3, 4, 5]
    return pepLengths    
def pathIsClosed(path, chosenKmers, kmerSize):
    closed = True
    for i in range(0, kmerSize):
        kmerSeq = path[-(kmerSize - i):] + path[:i]
        kmerName = "-".join([str(x) for x in kmerSeq])
        kmerName = tuple(kmerSeq)
        if kmerName in chosenKmers:
            continue
        else:
            closed = False
            return closed
    return closed


def relatedPermutations(permutation):
    cyclicPermutations = [permutation[i:] + permutation[:i] for i in range(len(permutation))]
    permutation.reverse()
    reverseCyclicPermutations = [permutation[i:] + permutation[:i] for i in range(len(permutation))]
    allRelated = cyclicPermutations + reverseCyclicPermutations
    allRelatedNames = []
    for x in allRelated:
        name = tuple([y for y in x])
        # name = '-'.join([str(y) for y in x])
        allRelatedNames.append(name)
        x.reverse()
        name = tuple([y for y in x])
        # name = '-'.join([str(y) for y in x])
        allRelatedNames.append(name)

    return allRelatedNames


protonMass = 1.00728


def findWalksLengthK(combination, graph, start, l, nodes, building_blocks):
    # Finds all Walks of length k from the node start in a directed graph. Uses modified BFS for traversal.
    combinationFreq = dict((aa, combination.count(aa)) for aa in set(building_blocks))
    allPaths = []
    pathAAfrequencies = {}
    for aa in building_blocks:
        if aa in nodes[start][1]:
            pathAAfrequencies[aa] = nodes[start][1][aa]
        else:
            pathAAfrequencies[aa] = 0
    startSequence = nodes[start][0]

    queue = [(start, pathAAfrequencies, list(startSequence))    ]
    while True:
        if not queue:
            break
        (vertex, frequencies, path) = queue.pop(0)
        if vertex in graph:
            for nextName in set(graph[vertex]):
                next = nodes[nextName][0][:]
                nextAA = next[-1]
                newFrequencies = frequencies.copy()
                if nextAA in frequencies:
                    newFrequencies[nextAA] = frequencies[nextAA] + 1
                else:
                    newFrequencies[nextAA] = 1
                if newFrequencies[nextAA] > combinationFreq[nextAA]:
                    continue
                newPath = path + [nextAA]
                if len(newPath) == l:
                    feasiblePath = True
                    for aa in building_blocks:
                        if newFrequencies[aa] != combinationFreq[aa]:
                            feasiblePath = False
                            break
                    if feasiblePath:
                        allPaths.append(newPath)
                if len(newPath) < l:
                    queue.append((nextName, newFrequencies, newPath))
                if len(newPath) > l + 1:
                    break
    return [tuple(path1) for path1 in allPaths]


def spellPath(path):
    spelledSeq = list(path[0])
    for i in range(1, len(path)):
        # pointer = len(kmers[0])-1
        spelledSeq = spelledSeq + [path[i][-1]]
    return spelledSeq


def findLeaves(adjacencies1):
    newleaves = []
    for node in adjacencies1:
        if len(adjacencies1[node]) == 0:
            newleaves.append(node)
    return newleaves
