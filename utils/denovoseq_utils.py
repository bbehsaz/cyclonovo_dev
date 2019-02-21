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


def generate_feasible_combinations(foundAutoConv, pepLengths, pepMass, delta):
    # Generates all feasible combinations for lengths that
    failedLenths = []
    allFeasibleCombinations = {}
    for l in pepLengths:
        allFeasibleCombinations[l] = set()
        allFeasibleCombinations[l] = find_feasible_combintions(
            foundAutoConv, l, pepMass, delta)
        if len(allFeasibleCombinations) == 0:
            failedLenths.append(l)
    return allFeasibleCombinations, failedLenths


