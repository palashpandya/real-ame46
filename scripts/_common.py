import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
FAILED = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)
    if not ok:
        FAILED.append(name)


def done():
    print(f"\n{'ALL CHECKS PASSED' if not FAILED else f'{len(FAILED)} CHECK(S) FAILED: {FAILED}'}")
    sys.exit(1 if FAILED else 0)
