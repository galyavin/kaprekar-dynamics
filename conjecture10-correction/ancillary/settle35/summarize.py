"""Recompute the frozen finite probe tallies from the integer TSV records."""
from collections import Counter
from pathlib import Path

from comparison_geometry import weight

HERE = Path(__file__).resolve().parent


def rows(name):
    with (HERE/name).open() as stream:
        return [tuple(map(int,line.split())) for line in stream]


def roots(data):
    n = len(data)
    censored = sum(a[2] < 0 for a in data)
    exported = sum(a[7] > 0 for a in data)
    good = sum(a[9] for a in data)
    bounds = sum(a[2]>=0 and not a[7] and not a[9] for a in data)
    return n,censored,exported,good,bounds


def physical(data,n,header,base_column):
    m = (n+1)//2
    counter = Counter()
    for row in data:
        base, p = row[base_column], row[header:]
        detected = p[0]>=0
        counter['resolved'] += detected
        counter['censored'] += not detected
        counter['all_C'] += p[4]
        counter['all_Cstar'] += p[5]
        counter['above_C_1/1000'] += p[4] and p[3]*1000>m*base
        counter['above_C_1/200'] += p[4] and p[3]*200>m*base
        counter['uninterrupted_above_C_1/1000'] += p[4] and p[3]*1000>m*base and p[9]<0
        for j, name in enumerate([
            'height_fail','pair_fail','other_row_fail','leave_Cstar','return_Cstar',
            'admitted_row_breach','return_C_after_row_breach','leave_C','return_C',
        ],6):
            counter[name] += p[j]>=0
        counter['entries_outside_C'] += p[13]==0
        counter['literal_steps_compared'] += p[15]+1
    return dict(counter)


def present(location, counter):
    print(location, ' '.join(f'{k}={v}' for k,v in counter.items()))


def main():
    print('Root cube L / N / censored / positive recurrent export / budgeted / zero export outside budgets')
    cubes = []
    for L in (64,128,256,512,1024,2048,4096):
        result = roots(rows(f'cube-{L}.tsv'))
        cubes.append(result)
        print(L,*result)
    low = [sum(a[i] for a in cubes[:3]) for i in range(5)]
    high = [sum(a[i] for a in cubes[-3:]) for i in range(5)]
    decline = (low[3]*high[0]-(high[3]+high[1])*low[0])*10>=low[0]*high[0]
    print('Falling criterion',decline,'low share floor',f'{low[3]}/{low[0]}',
          'high share ceiling',f'{high[3]+high[1]}/{high[0]}')
    for B in (200000,2000003,20000000):
        data = rows(f'box35-{B}.tsv')
        print('Box35',B,'autonomous root',*roots(data))
        summary = physical(data,35,22,0)
        present('  literal',summary)
        complement = {5,7,10,14,15}
        central = [B//20 + (B//50 if i in complement else 0) for i in range(18)]
        central[0] += B-sum(central)
        print('  entire-box original weight interval',
              *(v['weight'] for v in weight(35,B,B//10000,central)))
        for a in data:
            p = a[22:]
            if p[4] and p[3]*1000>18*B:
                print('  hit',a[1],'entry',p[0],'period',p[1],'exact peak numerator',p[3],
                      'root h2',a[6],'label exports at 5,7,10,14,15',*p[26:31])
    for length in (97,119):
        data = rows(f'comparison-{length}.tsv')
        assert data[0][2]==-1
        present(f'Length {length} centre control',physical(data[:1],length,3,1))
        present(f'Length {length} random box',physical(data[1:],length,3,1))


if __name__ == '__main__':
    main()
