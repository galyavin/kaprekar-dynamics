"""Build file-producing checks in a private snapshot and serialize substantial runs."""
from concurrent.futures import ThreadPoolExecutor
import fcntl
from pathlib import Path
import shutil
import subprocess
import sys


def main():
    source = Path(__file__).resolve().parent
    repo = source.parents[2]
    mode, path = sys.argv[1:]
    work = Path(path).resolve()
    assert work.is_dir() and not work.is_relative_to(repo)
    if mode == 'holdout':
        subprocess.run(['shasum','-a','256','-c','FREEZE.sha256'],cwd=source,check=True)
    files = [source/'model.cpp',source/'physical.cpp',source.parent/'r2_routes/diagnostic.cpp',
             source.parent/'stage_graph/graph.cpp']
    if mode=='holdout':
        files.append(source/'local.cpp')
    for file in files:
        shutil.copy2(file,work/file.name)
    for key in ('model','physical',*(['local'] if mode=='holdout' else [])):
        subprocess.run(['clang++','-std=c++17','-O3','-Wall','-Wextra','-Wno-return-type',
                        f'{key}.cpp','-o',key],cwd=work,check=True)

    def phase(kind, executable, label, n, chunk):
        def shard(start):
            with (work/f'{kind}-{start:06}.txt').open('wb') as output:
                subprocess.run([f'./{executable}',label,str(start),str(min(chunk,n-start))],
                               cwd=work,stdout=output,check=True)
            print(f'{kind} {start} to {min(n,start+chunk)-1}', flush=True)
        with ThreadPoolExecutor(max_workers=8) as pool:
            list(pool.map(shard,range(0,n,chunk)))

    with (repo/'.task-monitor/compute.lock').open('rb') as lock:
        print(f'{mode} waiting for the exclusive compute window',flush=True)
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX)
        print(f'{mode} acquired the exclusive compute window',flush=True)
        if mode=='derivation':
            phase('premise','physical','premise',60_000,7500)
            phase('explore','model','explore',4000,500)
        elif mode=='holdout':
            subprocess.run(['./physical','verify'],cwd=work,check=True)
            phase('holdout','physical','3072000',36_000,4500)
            phase('local','local','30720000',4096,512)
        else:
            raise ValueError(mode)


if __name__ == '__main__':
    main()
