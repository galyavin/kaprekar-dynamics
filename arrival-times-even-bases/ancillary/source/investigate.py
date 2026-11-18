

def output_digits(base, pair):
    u, v = pair
    if v:
        return (u, v - 1, base - 1 - v, base - u)
    return (u - 1, base - 1, base - 1, base - u)


def transition(base, pair):
    a, b, c, d = sorted(output_digits(base, pair), reverse=True)
    return a - d, b - c


def published_maximum_distance(base):
    """Devlin--Zeng (2021), Theorem 1, for bases divisible by five.

    https://math.colgate.edu/~integers/v97/v97.pdf
    Maximum finite hitting time; excludes trajectories entering other cycles.
    """
    if base < 5 or base % 5:
        raise ValueError('This function requires a positive base divisible by five.')
    odd_part, exponent = base // 5, 0
    while odd_part % 2 == 0:
        odd_part //= 2
        exponent += 1
    if odd_part > 1:
        return exponent + 2
    if exponent < 3:
        return (4, 7, 10)[exponent]
    if exponent % 4 == 0:
        return 4 * exponent + 6
    if exponent % 4 == 3:
        return 5 * exponent + 6
    return 3 * exponent + 5
