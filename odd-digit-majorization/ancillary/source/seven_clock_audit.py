"""Yices subprocess adapter."""
import shutil, subprocess, time


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def yices_status(smt,timeout=180):
    executable=shutil.which('yices-smt2')
    require(executable is not None,'Yices 2 command yices-smt2 is required')
    start=time.monotonic()
    try:
        r=subprocess.run([executable,f'--timeout={timeout}'],input='(set-logic QF_LIA)\n'+smt,text=True,capture_output=True,timeout=timeout+5)
    except subprocess.TimeoutExpired:return 'unknown (process timeout)',time.monotonic()-start
    status=r.stdout.strip()
    require(r.returncode==0 and status in ('sat','unsat','unknown'),'Unexpected Yices result: '+r.stdout+r.stderr)
    return status,time.monotonic()-start
