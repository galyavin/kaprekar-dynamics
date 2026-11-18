from collections import Counter
from cycle_model import number_from_pair


def digit_step(base, number):
    digits=[]
    for _ in range(4):
        number, digit=divmod(number,base)
        digits.append(digit)
    digits.sort()
    low=high=0
    for digit in digits:
        low=base*low+digit
    for digit in reversed(digits):
        high=base*high+digit
    return high-low


def direct_input_histogram(base, target):
    destinations={number_from_pair(base,p) for p in target}
    distances={x:0 for x in destinations}
    counts=Counter()
    for number in range(base**4):
        if number%(base**3+base**2+base+1)==0:
            continue
        path, seen, x=[],set(),number
        while x not in distances and x not in seen:
            path.append(x); seen.add(x)
            x=digit_step(base,x)
        value=distances.get(x)
        for x in reversed(path):
            if value is not None:
                value+=1
            distances[x]=value
        if distances[number] is not None:
            counts[distances[number]]+=1
    return dict(sorted(counts.items()))
