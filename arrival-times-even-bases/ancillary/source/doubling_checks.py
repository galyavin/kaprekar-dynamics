

def multiplicity(base, pair):
    """Number of nonconstant four-slot strings with this difference pair."""
    u, v = pair
    if v == 0:
        return (base-u)*(12*u-4)
    if u == v:
        return 6*(base-u)
    return 24*(base-u)*(u-v)
