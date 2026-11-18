

def moments(start, step, lower, upper):
    """Count, sum, and square sum in start+step*k within [lower, upper]."""
    first = start + max(0, (lower-start+step-1)//step)*step
    if first > upper:
        return 0, 0, 0
    count = (upper-first)//step+1
    s1 = count*(count-1)//2
    s2 = count*(count-1)*(2*count-1)//6
    return count, count*first+step*s1, count*first*first+2*first*step*s1+step*step*s2


def v2(value):
    if value <= 0:
        raise ValueError('2-adic valuation requires a positive integer')
    return (value & -value).bit_length()-1
