#!/bin/sh
set -e
cd "$(dirname "$0")"
python3 -B even_check.py ./ww threshold 5 271137 1000000
python3 -B even_check.py ./ww exhaustive 5 1201440 2000160 2000004 2998560 2999994
python3 -B even_check.py ./ww sample 5 10002720:4000:1 100000800:500:1 100000002:300:0 1000000800:100:1 1000000002:40:0 10000000800:20:1 100000000800:5:0
