"""Exact reached release inventories, polytope geometry and independent paths."""
import argparse,copy,json
from fractions import Fraction as Q
from itertools import combinations
from math import ceil
from pathlib import Path
import nine_a_facet as af
import nine_b_facet as bf
import nine_ab as ab
import nine_bb as bb
from nine_entry import evaluate,compose,EXIT_FACETS
from nine_b_facet_checks import require,from_orders
from nine_phase_checks import replay,raw_step
from nine_releases import release,records
from odd_affine import word,progression

ROOT=Path(__file__).resolve().parent
A0=(1,0,0,1,-1);A1=(1,1,-1,0,-1)
L=((1,0,0,0,0),(1,0,0,-1,0),(0,1,0,-1,0),(0,1,-1,0,0),(0,0,0,0,1))
# BA0: (U,V,p,e,B,1). BBA1: (U,p,r,e,B,1).
CHARTS={
 'BA0':((1,0,0,0,0,0),(0,1,0,0,0,0),(0,-1,1,0,1,0),(-1,0,0,1,1,0)),
 'BBA1':((1,0,0,0,0,0),(Q(-1,2),Q(1,2),0,Q(1,2),1,0),
          (Q(1,2),Q(1,2),0,Q(-1,2),0,0),(Q(1,4),Q(1,4),Q(-1,2),Q(1,4),0,0)),
}
ACTIVE={'BA0':(bf.B0,A0),'BBA1':(bf.B0,bb.B1,A1)}


def solve(rows):
    a=[list(map(Q,(*r[:4],-r[4]))) for r in rows]
    for i in range(4):
        j=next((j for j in range(i,4) if a[j][i]),None)
        if j is None:return None
        a[i],a[j]=a[j],a[i];p=a[i][i];a[i]=[x/p for x in a[i]]
        for j in range(4):
            if j!=i:
                p=a[j][i];a[j]=[x-p*y for x,y in zip(a[j],a[i])]
    return tuple(row[-1] for row in a)


def rank(rows):
    rows=[list(map(Q,row)) for row in rows];i=0
    for col in range(4):
        j=next((j for j in range(i,len(rows)) if rows[j][col]),None)
        if j is None:continue
        rows[i],rows[j]=rows[j],rows[i];p=rows[i][col];rows[i]=[x/p for x in rows[i]]
        for j in range(len(rows)):
            if j!=i:
                p=rows[j][col];rows[j]=[x-p*y for x,y in zip(rows[j],rows[i])]
        i+=1
        if i==len(rows):break
    return i


def face_vertices(active):
    others=[r for r in bf.FACETS if r not in active];vs=set()
    for extra in combinations(others,4-len(active)):
        v=solve((*active,*extra))
        if v is not None and all(evaluate(r,(*v,1))>=0 for r in bf.FACETS):vs.add(v)
    return tuple(sorted(vs))


def exit_geometry():
    require({compose(row,L) for row in bf.FACETS}==set(bf.FACETS),
            'Leading rotation permutes the complete core-facet inventory')
    names={r:n for n,r,s in EXIT_FACETS};result={}
    for kind,active,delta,expected in (
        ('A',(af.A2,),af.DELTA_A,{'A0','A1','B1'}),
        ('B',(bf.B0,),bf.DELTA_B,{'A0','A1','B1','B2'}),
        ('BB',(bf.B0,bb.B1),bb.DELTA_BB,{'A1','B2'})):
        out=[];full=[]
        for row in bf.FACETS:
            if evaluate(row[:4],delta)>=0:continue
            vs=face_vertices((*active,row))
            dim=rank([tuple(a-b for a,b in zip(v,vs[0])) for v in vs[1:]]) if vs else -1
            name=names.get(row,str(row));out.append({'face':name,'dimension':dim,'vertices':vs})
            if dim==3-len(active):full.append(name)
        require(set(full)==expected,'Complete full-dimensional outgoing face inventory')
        result[kind]=out
    return result


def extended_chart(kind):return (*CHARTS[kind],(0,0,0,0,1,0),(0,0,0,0,0,1))


def specialize(kind,phase):
    if kind=='BA0':
        p,e=phase
        return tuple((r[0],r[1],r[4],r[2]*p+r[3]*e+r[5]) for r in extended_chart(kind))
    p,r,e=phase
    return tuple((g[0],0,g[4],g[1]*p+g[2]*r+g[3]*e+g[5]) for g in extended_chart(kind))


def incoming(kind):
    phases=bf.PHASES if kind=='BA0' else bb.PHASES;out=[];thresholds=[];count=0
    for phase in phases:
        orders=bf.six_orders(phase) if kind=='BA0' else bb.twelve_orders(phase)
        _,guards,_=replay(from_orders(orders));ts=[]
        for g in guards:
            row=compose(g,extended_chart(kind))
            macro=(row[0],row[1],row[4]) if kind=='BA0' else (row[0],0,row[4])
            if any(macro):continue
            c=row[2]*phase+row[5] if kind=='BA0' else row[1]*phase[0]+row[2]*phase[1]+row[5]
            slope=row[3];count+=1
            if slope>0:ts.append(ceil(-c/slope))
            else:require(slope==0 and c>=0,'No omitted transverse release guard')
        m=max(ts);width=8 if kind=='BA0' else 16
        prefix=(phase,) if kind=='BA0' else phase
        thresholds.append((*prefix,m));out.extend((*prefix,e) for e in range(m-width,m))
    require(len(out)==(32 if kind=='BA0' else 192),'Complete release input inventory')
    return out,thresholds,count


def numeric_source(kind,B,phase):
    U=(21*B//25 if kind=='BA0' else 83*B//100);V=61*B//100
    if kind=='BBA1':
        p,r,e=phase;U+=(-U-p+2*r-e)%4
    source=specialize(kind,phase)
    q=tuple(evaluate(row,(U,V,B,1)) for row in source[:4])
    require(all(Q(x).denominator==1 for x in q),'Integral release source')
    return tuple(map(int,q))


def generate():
    result=[];B=1000000000
    for kind in CHARTS:
        for phase in incoming(kind)[0]:
            q=numeric_source(kind,B,phase);prefix=[];rotation=1 if kind=='BA0' else 2
            rec=word(B,q,rotation);prefix.extend(rec['orders']);q=tuple(evaluate(r,(*q,B,1)) for r in rec['matrix'])
            for _ in range(1000):
                p=af.phase(B,q) if kind=='BA0' else ab.decode(B,q)[2]
                valid=p in (af.PHASES if kind=='BA0' else ab.PHASES)
                if valid:
                    orders=af.twelve_orders(p) if kind=='BA0' else ab.return_orders(p);rec=word(B,q,len(orders))
                    if rec['orders']==orders:
                        n,delta=progression(B,q,rec);expected=af.DELTA_A if kind=='BA0' else ab.DELTA_AB
                        if delta==expected and n>10000:
                            result.append({'kind':kind,'phase':phase,'orders':prefix});break
                rec=word(B,q,3);prefix.extend(rec['orders']);q=tuple(evaluate(r,(*q,B,1)) for r in rec['matrix'])
            else:raise ArithmeticError(('Unresolved finite release',kind,phase))
    data={'records':result}
    (ROOT/'nine_releases_certificates.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    records.cache_clear();return data


def guard_check(row,kind,vertices):
    a,b,c,k=row
    if a or b or c:
        values=[a*v[0]+b*v[1]+c for v in vertices]
        require(min(values)>=0 and sum(values)>0,'Leading release guard on the entire relative face interior')
    else:require(k>=0,'Exact constant release guard')


def verify(data):
    expected=[(kind,phase) for kind in CHARTS for phase in incoming(kind)[0]]
    require([(r['kind'],tuple(r['phase'])) for r in data['records']]==expected,'Complete two-release inventory')
    geometry=exit_geometry();count=0;maximum={kind:0 for kind in CHARTS}
    vertices={kind:face_vertices(ACTIVE[kind]) for kind in CHARTS}
    for record,(kind,phase) in zip(data['records'],expected):
        source=specialize(kind,phase);M,guards,offsets=replay(from_orders(tuple(map(tuple,record['orders']))))
        for g in guards:guard_check(compose(g,source),kind,vertices[kind]);count+=1
        target=tuple(compose(row,source) for row in M)
        leading=tuple(r[:3] for r in source[:4])+((0,0,1),)
        for _ in range(1 if kind=='BA0' else 2):leading=tuple(compose(r,leading) for r in L)
        require(tuple(r[:3] for r in target)==leading[:4],'Exact leading release rotation')
        p=2*target[0][-1]-target[1][-1]
        if kind=='BA0':
            require(p in af.PHASES,'Recurrent outgoing A phase');orders=af.twelve_orders(p);delta=af.DELTA_A
        else:
            r=target[0][-1]+target[1][-1]-2*target[3][-1];p=(p,r)
            require(p in ab.PHASES,'Recurrent outgoing AB phase');orders=ab.return_orders(p);delta=ab.DELTA_AB
        ext=(*target,(0,0,1,0),(0,0,0,1));R,gs,_=replay(from_orders(orders))
        require(tuple(compose(row,ext) for row in R)==tuple((*row[:3],row[-1]+d) for row,d in zip(target,delta)),
                'Outgoing exact recurrent displacement')
        for g in gs:
            row=compose(g,ext);guard_check(row,kind,vertices[kind]);count+=1
            if not any(row[:3]):require(evaluate(g[:4],delta)>=0,'Release has no bounded-repeat obstruction')
        maximum[kind]=max(maximum[kind],len(record['orders']))
    return {'geometry':geometry,'rational_guard_checks':count,'max_prefix':maximum,
            'inventories':{kind:{'count':len(incoming(kind)[0]),'thresholds':incoming(kind)[1]} for kind in CHARTS}}


def actual_checks():
    steps=0;count=0
    for kind in CHARTS:
        for phase in incoming(kind)[0]:
            for B in range(1000000,1000005):
                q=numeric_source(kind,B,phase);out=release(B,q,kind)
                for _ in range(out['elapsed']):q=raw_step(B,q)
                require(q==out['target'],'Independent full boundary release');steps+=out['elapsed'];count+=1
    return {'independent_releases':count,'independent_steps':steps}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--generate',action='store_true');args=parser.parse_args()
    data=generate() if args.generate else json.loads((ROOT/'nine_releases_certificates.json').read_text())
    result=verify(data);controls=0
    for kind in ('missing','word','phase'):
        bad=copy.deepcopy(data)
        if kind=='missing':bad['records'].pop()
        elif kind=='word':bad['records'][0]['orders'].pop()
        else:
            bad['records'][0]['phase']=list(bad['records'][0]['phase']);bad['records'][0]['phase'][0]+=1
        try:verify(bad)
        except (ArithmeticError,AssertionError,ValueError):controls+=1
        else:raise ArithmeticError('Corrupt release certificate accepted')
    result.update(actual_checks());result.update({'status':'passed','corruptions_rejected':controls,
          'scope':'Complete reached BA0 and BBA1 release inventories and leading exit geometry; original-input phase weights remain separate.'})
    (ROOT/'nine_releases_results.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
    report=(f'PASS: {result["rational_guard_checks"]} exact guard checks; both complete release inventories (32 and 192 states).\n'
            f'Maximum prefix lengths: {result["max_prefix"]}.\n'
            f'{result["independent_releases"]} independent complete releases over {result["independent_steps"]} steps; '
            f'{controls} corruptions rejected.\n')
    (ROOT/'nine_releases_verification.txt').write_text(report);print(report,end='')


if __name__=='__main__':main()
