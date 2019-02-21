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


def find_repeat_chains_longerThenx(offset, offset_pairs, x):
    all_xcords = [pair[0] for pair in offset_pairs]
    all_ycords = [pair[1] for pair in offset_pairs]
    original_shared_values = intersect(all_xcords, all_ycords)
    maxchain = 2
    last_shared_values = original_shared_values[:]
    while maxchain < x and len(last_shared_values) > 0:
        start_points = []
        new_shared_values = []
        for value in last_shared_values:
            start_points += all_indices(value, all_xcords)

        new_shared_values = extendChainOneStep(start_points, original_shared_values, all_xcords, all_ycords)
        last_shared_values = new_shared_values
        if len(new_shared_values) > 0:
            maxchain += 1
        else:
            break
    if maxchain == x:
        return 1
    else:
        return 0