"""Exact boundary releases from B/A0 and BB/A1 into the solved A routes."""
import json
from functools import lru_cache
from pathlib import Path
import nine_a_facet as af
import nine_ab as ab
from nine_entry import evaluate
from odd_affine import word


@lru_cache(None)
def records():
    data=json.loads(Path(__file__).with_name('nine_releases_certificates.json').read_text())
    return {(r['kind'],tuple(r['phase'])):r for r in data['records']}


def phase(B,q,kind):
    u,v,w,z=q
    if kind=='BA0':return v+w-B,u+z-B
    if kind=='BBA1':return v+w-B,u+v-2*z-B,u+v-w-B
    raise ValueError('Unknown release.')


def release(B,q,kind):
    q=tuple(q);record=records().get((kind,phase(B,q,kind)))
    if record is None:raise ValueError('Outside the complete incoming release inventory.')
    expected=tuple(map(tuple,record['orders']));rec=word(B,q,len(expected))
    if rec['orders']!=expected:raise ValueError('Boundary-release prefix guards fail.')
    target=tuple(evaluate(row,(*q,B,1)) for row in rec['matrix'])
    orders=af.twelve_orders(af.phase(B,target)) if kind=='BA0' else ab.return_orders(ab.decode(B,target)[2])
    if word(B,target,len(orders))['orders']!=orders:raise ValueError('Outgoing release return guards fail.')
    return {'elapsed':len(expected),'target':target,'outgoing':'A' if kind=='BA0' else 'AB'}
