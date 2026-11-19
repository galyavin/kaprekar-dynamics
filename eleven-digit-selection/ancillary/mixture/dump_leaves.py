"""Print the aggregate record and write optional leaf data to a gzip file."""
import argparse
import fcntl
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('first', type=int)
    parser.add_argument('last', type=int)
    parser.add_argument('output', type=Path)
    parser.add_argument('--all', action='store_true', help='dump every leaf in the selected residues')
    args = parser.parse_args()
    if not 0 <= args.first <= args.last <= 5:
        parser.error('expected a residue interval in 0..5')
    source = Path(__file__).resolve().parent
    env = os.environ.copy()
    if args.all:
        env['KAPREKAR_DUMP_ALL'] = '1'
    else:
        env.pop('KAPREKAR_DUMP_ALL', None)
    with tempfile.TemporaryDirectory(prefix='p7-leaves-') as tmp:
        scratch = Path(tmp)
        lock_path = Path(env.get('KAPREKAR_COMPUTE_LOCK') or scratch / 'compute.lock')
        with lock_path.open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            with (scratch / 'verify.cpp').open('wb') as patched:
                subprocess.run(['patch', '-s', '-o', '-'], cwd=source, stdout=patched,
                               input=(source / 'leaf-dump.patch').read_bytes(), check=True)
            executable = scratch / 'verify'
            subprocess.run([env.get('CXX', 'clang++'), '-O3', '-std=c++17', '-Wall', '-Wextra',
                            scratch / 'verify.cpp', '-o', executable], check=True)
            command = [executable, str(args.first), str(args.last)]
            with gzip.open(args.output, 'xb') as output:
                with subprocess.Popen(command, env=env, stderr=subprocess.PIPE) as process:
                    assert process.stderr is not None
                    try:
                        shutil.copyfileobj(process.stderr, output)
                    except BaseException:
                        process.kill()
                        raise
                    if status := process.wait():
                        raise subprocess.CalledProcessError(status, command)


if __name__ == '__main__':
    main()
