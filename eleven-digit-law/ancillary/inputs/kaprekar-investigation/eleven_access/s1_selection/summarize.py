"""Reduce private indexed panel shards, keeping covariance on the original draws."""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
import json
from math import sqrt
from pathlib import Path
from statistics import NormalDist
import sys

Z95 = NormalDist().inv_cdf(.975)
ZP = NormalDist().inv_cdf(1-.05/(2*7))
ZL = NormalDist().inv_cdf(1-.05/(2*6))


def estimate(values, weight=None):
    n = len(values)
    if weight is None:
        weight = [1]*n
    total = sum(weight)
    mean = sum(x*w for x,w in zip(values,weight))/total
    return mean, [w*(x-mean)*n/total for x,w in zip(values,weight)]


def linear(*terms):
    return sum(k*value for k,(value,_) in terms), [sum(k*v for (k,_),v in zip(terms,row))
                                                  for row in zip(*(a for _,(_,a) in terms))]


def interval(value, crit=Z95):
    estimate, influence = value
    n = len(influence)
    half = crit*sqrt(sum(x*x for x in influence)/(n*(n-1)))
    return dict(estimate=estimate, halftime=half, CI=[estimate-half,estimate+half])


def wilson(m,n):
    p = m/n
    c = (p+Z95**2/(2*n))/(1+Z95**2/n)
    h = Z95*sqrt(p*(1-p)/n+Z95**2/(4*n*n))/(1+Z95**2/n)
    return dict(count=m, denominator=n, Wilson95=[max(0,c-h),min(1,c+h)])


def read(work, kind, count, parse):
    rows, hashes = [], {}
    for path in sorted(work.glob(f'{kind}-*.txt')):
        digest = sha256()
        with path.open('rb') as handle:
            for line in handle:
                digest.update(line)
                row = line.split()
                assert int(row[0])==len(rows)
                rows.append(parse(row[1:]))
        hashes[path.name] = digest.hexdigest()
    assert len(rows)==count, (kind,len(rows),count)
    return rows, hashes


def parse_physical(row):
    entered, *fields = row
    assert len(fields)%10==0
    return int(entered), [tuple(map(int,fields[j:j+9]))+(float(fields[j+9]),)
                          for j in range(0,len(fields),10)]


def observed(columns, j, b):
    keys = ('clock','strict','strict_loss','five','five_loss','skipped','fallback','status','model_loss','model')
    out = {key: list(col) for key,col in zip(keys,zip(*(v[j] for _,v in columns)))}
    out['value'] = [min(1,max(0,t/b)) for t in out['clock']]
    out['captured'] = [t>=0 for t in out['clock']]
    return out


def compare(l,r,values,weights=None):
    return linear((1,estimate(l[values],weights)),(-1,estimate(r[values],weights)))


def count_exclusions(a,leading,base):
    chart = Counter()
    for i,(c,s,l,f,m,model) in enumerate(zip(a['captured'],a['strict'],a['strict_loss'],
                                                a['five'],a['five_loss'],a['status'])):
        flags = (not leading[i],not c,s<0, f<0, 50*s>base, 50*f>base, bool(l&1),bool(l&2),
                 bool(m&1),bool(m&2),model not in (0,5,6))
        names = ('leading_cap','uncaptured','strict_never','five_never','strict_late','five_late',
                 'strict_interlevel','strict_zero','five_interlevel','five_zero','other_model_stop')
        chart[','.join(k for k,v in zip(names,flags) if v) or 'none'] += 1
    n = len(leading)
    return dict(flags={name: wilson(sum(v for k,v in chart.items() if name in k.split(',')),n)
                       for name in names}, overlaps={k: wilson(v,n) for k,v in sorted(chart.items())},
                accelerated=sum(x>0 for x in a['skipped']), skipped_updates=sum(a['skipped']),
                fallback=sum(a['fallback']), model_status=dict(Counter(a['status'])))


def pair_coverage(a,b,flags):
    chart = Counter()
    for i,d in enumerate(flags):
        def status(x):
            return (not x['captured'][i], x['strict'][i]<0, x['five'][i]<0,
                    x['strict_loss'][i],x['five_loss'][i],x['status'][i])
        chart[str((d,status(a),status(b)))] += 1
    return {tag: wilson(n,len(flags)) for tag,n in sorted(chart.items())}


def premise(work):
    columns, hashes = read(work,'premise',60_000,parse_physical)
    rs = (3,5,15,17)
    arms = {r: observed(columns,j,192000+r) for j,r in enumerate(rs)}
    report = {'rows':len(columns),'hashes':hashes,'readers':{}}
    hist = json.loads((Path(__file__).resolve().parent.parent/'r2_routes/diagnostic-results.json').read_text())
    for reader in ('strict','five'):
        pairs = []
        info = {}
        original = hist['readers'][reader]['scales']['192000']
        for lhs,rhs in ((3,5),(15,17)):
            l,r = arms[lhs],arms[rhs]
            chosen = [int(a>=0 and b>=0 and c and d and e==f==0) for a,b,c,d,e,f in
                      zip(l[reader],r[reader],l['captured'],r['captured'],l['status'],r['status'])]
            k = f'pair_{lhs}_{rhs}'
            graph = compare(l,r,'model',chosen)
            actual = compare(l,r,'value',chosen)
            assert abs(graph[0]-original[k]['entry']['model']['estimate'])<1e-13, (reader,k,graph[0])
            assert abs(actual[0]-original[k]['entry']['physical']['estimate'])<1e-13, (reader,k,actual[0])
            assert sum(chosen)==original[k]['retained_entry']['count']
            info[k] = dict(model=interval(graph),physical=interval(actual),retained=sum(chosen))
            pairs.append(graph)
        info['pooled_model'] = interval(linear((.5,pairs[0]),(.5,pairs[1])))
        report['readers'][reader] = info
    return report


def explore_parse(row):
    fields = iter(row)
    entered = int(next(fields))
    if entered<0:
        return entered, []
    arms = []
    for _ in range(24):
        status, t, n = int(next(fields)),float(next(fields)),int(next(fields))
        stages = []
        for _ in range(n):
            duration = float(next(fields))
            ranks = tuple(int(next(fields)) for _ in range(2))
            length, coeff, dec, index, residue, loss = (int(next(fields)) for _ in range(6))
            offset = tuple(int(next(fields)) for _ in range(10))
            after = tuple(int(next(fields)) for _ in range(6))
            stages.append((duration,ranks,length,coeff,dec,index,residue,loss,offset,after))
        arms.append((status,t,stages))
    assert next(fields,None) is None
    return entered,arms


def explore(work):
    data, hashes = read(work,'explore',4000,explore_parse)
    eligible = [(i,a) for i,(e,a) in enumerate(data) if e>=0 and all(not k[0] for k in a)]
    results = {}
    n = len(data)
    for r,s in combinations(range(24),2):
        counts = Counter()
        total, squares, changed = 0.,0.,0
        increments = [0.]*15
        profile = Counter()
        examples = {}
        for row,arms in eligible:
            l,h = arms[r],arms[s]
            delta = l[1]-h[1]
            total += delta
            squares += delta**2
            changed += delta!=0
            first_any = first_macro = 0
            for i in range(max(len(l[2]),len(h[2]))):
                a = l[2][i] if i<len(l[2]) else None
                b = h[2][i] if i<len(h[2]) else None
                difference = (a[0] if a else 0)-(b[0] if b else 0)
                increments[i] += difference
                if not first_any and a!=b:
                    first_any = i+1
                if not first_macro and (difference!=0 or ((a[-1] if a else None)!=(b[-1] if b else None))):
                    first_macro = i+1
                    route = (first_macro,a[1][0] if a else -1,b[1][0] if b else -1)
                    profile[route] += 1
                    if route not in examples:
                        examples[route] = dict(row=row,left=a,right=b)
            counts[first_any,first_macro] += 1
        se = sqrt(max(0,squares-total*total/n)/(n*(n-1)))
        top = profile.most_common(3)
        results[f'{r}-{s}'] = dict(mean=total/n, SE=se, Wilson95_changed=wilson(changed,n),
                                   overlay={str(k):v for k,v in sorted(counts.items())},
                                   stage_clock_mean_differences=[v/n for v in increments if v!=0],
                                   per_stage={k+1:v/n for k,v in enumerate(increments) if v},
                                   top_macro_ranks=[dict(stage=a, incoming_ranks=(b,c), count=v,
                                                         example=examples[a,b,c]) for (a,b,c),v in top])
    return dict(rows=n, complete_all_24=len(eligible), hashes=hashes, by_pair=results)


def local_parse(row):
    assert len(row)%4==0
    return [(int(row[j]),float(row[j+1]),int(row[j+2]),int(row[j+3])) for j in range(0,len(row),4)]


def local(work):
    count = 4096
    data, hashes = read(work,'local',count,local_parse)
    residues = (3,5,9,11,15,17,21,23)
    estimates, arms = {}, {}
    for r,col in zip(residues,zip(*data)):
        steps, weights, skipped, fallback = zip(*col)
        score = [max(0,min(1,((t/(30720000+r))-.010383)/.000003)) for t in steps]
        est = estimate(score,weights)
        estimates[r] = est
        arms[r] = dict(ramp=interval(est), accelerated=sum(x>0 for x in skipped), fallback=sum(fallback),
                       one=sum(s==1 for s in score), partial=sum(0<s<1 for s in score))
    choices = ((3,5),(15,17),(3,9),(5,11),(3,15),(5,17))
    test = {f'{r}-{s}':interval(linear((1,estimates[r]),(-1,estimates[s])),ZL) for r,s in choices}
    separation = all(.09<=test[k]['CI'][0]<=test[k]['CI'][1]<=.19 for k in ('3-5','15-17'))
    equality = all(max(abs(v) for v in test[k]['CI'])<=.035 for k in ('3-9','5-11','3-15','5-17'))
    return dict(rows=count, hashes=hashes, arms=arms, contrasts=test,
                separation=separation,equality=equality)


def holdout(work):
    n = 36_000
    data, hashes = read(work,'holdout',n,parse_physical)
    leading = [t>=0 for t,_ in data]
    anchor = 3072000
    residues = (3,5,9,11,15,17,0,1)
    a = {r: observed(data,j,anchor+r) for j,r in enumerate(residues)}
    metadata = {}
    for r in residues:
        item = a[r]
        for k,flag in [('leading',leading),('strict',[e>=0 for e in item['strict']]),
                       ('five',[e>=0 for e in item['five']])]:
            item[k+'_score'] = [x if y and z else 0 for x,y,z in zip(item['value'],flag,item['captured'])]
        item['model_score'] = [x if y and status==0 else 0 for x,y,status in zip(item['model'],leading,item['status'])]
        metadata[r] = count_exclusions(item,leading,anchor+r)
    choices = ((3,5),(15,17),(3,9),(5,11),(3,15),(5,17),(0,1))
    result = {}
    observed_pairs = {}
    for r,s in choices:
        l,h = a[r],a[s]
        key = f'{r}-{s}'
        result[key] = {m:interval(compare(l,h,m),ZP if m=='leading_score' else Z95)
                       for m in ('leading_score','strict_score','five_score','value')}
        result[key]['cohort_overlaps'] = pair_coverage(l,h,leading)
        observed_pairs[key] = compare(l,h,'leading_score')
        if r not in (9,11) and s not in (9,11):
            model_n = 12_000
            graph = compare({k:v[:model_n] for k,v in l.items()},
                            {k:v[:model_n] for k,v in h.items()},'model_score')
            physical = compare({k:v[:model_n] for k,v in l.items()},
                               {k:v[:model_n] for k,v in h.items()},'leading_score')
            result[key]['model_subpanel'] = interval(graph)
            result[key]['physical_minus_model_subpanel'] = interval(linear((1,physical),(-1,graph)))
    # The seven simultaneous claims use the two subresidue panels separately.
    pool = interval(linear((.5,observed_pairs['3-5']),(.5,observed_pairs['15-17'])))
    separated = all(.00010<=result[k]['leading_score']['CI'][0]<=result[k]['leading_score']['CI'][1]<=.00036
                    and result[k]['leading_score']['halftime']<=.00010 for k in ('3-5','15-17'))
    eq = all(max(abs(x) for x in result[k]['leading_score']['CI'])<=.00012
             for k in ('3-9','5-11','3-15','5-17'))
    parity = result['0-1']['leading_score']['CI'][0]>.00015
    return dict(rows=n, model_subpanel_rows=12_000, hashes=hashes, leading_entry_by_cap=wilson(sum(leading),n), arm_reports=metadata,
                pairs=result, pooled_3_minus_5=pool, forecast=dict(separated=separated,equality=eq,parity=parity))


def main():
    kind,path = sys.argv[1:]
    output = globals()[kind](Path(path))
    print(json.dumps(output,indent=2))


if __name__=='__main__':
    main()
