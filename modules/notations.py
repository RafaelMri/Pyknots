# Notation conversion functions.

__all__ = ['crossing_pattern', 'dowker_thistlethwaite']

def crossing_pattern(braid):
    """ Returns the crossing pattern defined by tracing the knot from
        a starting point and continuing until the starting point is
        reached again.

        At each crossing we create a list representing
        the under and over arcs. Let 'o' be the overarc, 'u' the
        under arc, then the representation is [u, o, u].

        This representation is used in computing the braids DT code,
        as well as the Gauss code and potentially others.
    """
    strands = [i for i in range(braid.braid_index)]
    start_strands = [i for i in strands]
    crossings = []
    valid = True
    while valid:
        for i in braid.braid_notation:
            s1, s2 = strands[abs(i)-1], strands[abs(i)]
            if i > 0:
                strands[i-1], strands[i] = s2, s1
                crossings.append([s1, s2, s1])
            elif i < 0:
                strands[abs(i)-1], strands[abs(i)] = s2, s1
                crossings.append([s2, s1, s2])
        if strands == start_strands:
            valid = False
    return crossings

def _normalize_DT(DT):
    """ Helper for dowker_thistlethwaite function.
        DT notation can be normalized by multiplying -1 to all
        elements. This is because mirror images can not be
        distinguished by DT notation.
    """
    pos = neg = 0
    for i in DT:
        if i > 0:
            pos += 1
        else:
            neg += 1
    if pos >= neg:
        return DT
    else:
        return [-i for i in DT]

def dowker_thistlethwaite(braid):
    """ Returns Dowker-Thistlethwaite notation of the braid.

        Dowker Thistlewaite notation, ie DT notation, is a sequence of
        even integers computed by walking the length of the knot from
        a starting point (until the starting point is reached again)
        counting up and labelling each crossing by the current value
        of the count. If the value is even and we are going over a
        crossing, then we take the negative of the value.

        Each crossing will appear exactly twice, and so each is assigned
        a tuple of 1 even and 1 odd integer. The tuples are ordered by
        the odd integers, then the resulting sequence of even is
        the DT notation.
    """
    DT, DT_temp = [], []
    step = 1
    DT_dict = {}
    for i in range(braid.braid_length):
        DT_dict[i] = []
    for i, j in enumerate(crossing_pattern(braid)):
        crossing = i % braid.braid_length
        if 0 == j[0]:
            DT_dict[crossing].append(step)
            step += 1
        elif 0 == j[1]:
            if step % 2 == 0:
                DT_dict[crossing].append(-step)
            else:
               DT_dict[crossing].append(step)
            step += 1

    for i in DT_dict:
        temp = []
        for j in DT_dict[i]:
            temp.append(j)
        DT_temp.append(temp)

    for i in range(1, (2*braid.braid_length)+1):
        for j, k in enumerate(DT_temp):
            if i % 2 != 0:
                if i in k:
                    DT.append([m for m in k if m != i])
                if -i in k:
                    DT.append([m for m in k if m != -i])

    DT = [i for j in DT for i in j]
    DT = _normalize_DT(DT)
    return DT

def gauss(braid):
    """ Returns Gauss notation of the braid.
    """
    gauss = []
    for i, j in enumerate(crossing_pattern(braid)):
        crossing = i % braid.braid_length
        crossing += 1
        if 0 == j[0]:
            gauss.append(-crossing)
        elif 0 == j[1]:
            gauss.append(crossing)
    return gauss
