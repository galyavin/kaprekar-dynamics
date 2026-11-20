"""Certify maximum clock from exact anchor drift on the accepted guarded cover."""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
import argparse
import fcntl
import json
from pathlib import Path
import sys
import tempfile

from replay import ACCESS, ROOT
sys.path.insert(0,str(ACCESS/"z1b_rank_one"))
import verify
import cover


def one(executable,continuation,root,args):
    exact_reply = cover.reply
    anchors = {}
    reads = {}

    def checked_reply(*a,**kw):
        result = exact_reply(*a,**kw)
        kind, details, exits, _ = result
        if kind == "stage":
            _,_,period,_,rank = details
            assert exits
            drift = {descriptor[1] for _,descriptor,_,_ in exits}
            assert len(drift)==1
            slope = Q(-drift.pop(),period)
            anchors[rank] = min(anchors.get(rank,slope),slope)
            reads[rank] = reads.get(rank,0)+1
        return result

    cover.reply = checked_reply
    try:
        status = verify.run_root(executable,continuation,root,args)
    finally:
        cover.reply = exact_reply
    if status:
        raise RuntimeError(f"rank-one closure failed for {root=}: {status}")
    assert set(anchors) == set(range(1,6))
    data = dict(root=root,minimum_drift=anchors,reads=reads,minimum=min(anchors.values()))
    print("anchor-clock",json.dumps(data,default=str),flush=True)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--roots",type=int,nargs="+",choices=range(6),default=[1,3,5])
    args = parser.parse_args()
    options = argparse.Namespace(shorts=500,ingress_budget=50000,line_budget=100000,
                                 state_dir=None,resume=False,ingress_only=False)
    with (ROOT/".task-monitor/compute.lock").open() as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        with tempfile.TemporaryDirectory(prefix="task146-anchor-bound-") as area:
            work = Path(area)
            executable = verify.build(work/"ingress",False)
            continuation = verify.build(work/"continuation",True)
            with ProcessPoolExecutor(max_workers=len(args.roots)) as pool:
                calls = [pool.submit(one,executable,continuation,r,options) for r in args.roots]
                result = [call.result() for call in calls]
    minima = [r["minimum"] for r in result]
    print("anchor-clock-summary",json.dumps(dict(roots=result,minimum=min(minima)),default=str))


if __name__=="__main__":
    main()
