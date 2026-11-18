"""Seven-slot global clock potential and proved bounds.

The nine-step universal descent is audited in seven_clock_audit.py.
These are bounds, not a capture-time evaluator or a new phase sampler.
"""

def potential(base,q):
    u,v,w=q
    return sum(abs(4*g-base) for g in (base-u,u-v,v-w,w))


def positive_bound(base):
    """Pair first-capture bound if no zero difference is ever visited."""
    if not isinstance(base,int) or base<400:
        raise ValueError('Integer base at least 400 required.')
    return 54*base+85**3


def global_bound(base):
    """Pair first-capture bound, including all zero-face visits."""
    return positive_bound(base)+6
