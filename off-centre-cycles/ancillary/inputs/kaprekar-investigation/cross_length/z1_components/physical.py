"""Validate a positive original-word lift of an order-8 productive entrance."""

from pathlib import Path
import sys

from weights import ORDER, REFERENCES, step


def main(accepted_verifier):
    sys.path.insert(0, str(accepted_verifier.resolve()))
    from verify_root import exact_lift, verify

    for reference in REFERENCES:
        first = {}
        state = reference["source"]
        while state not in first:
            first[state] = len(first)
            state = step(state)
            assert len(first) <= 500_000
        capture = first[state]
        stop = len(first)
        print(reference["name"], "reference enters cycle at", capture, "return at", stop)
        result = verify({"order": ORDER, "source": reference["source"],
                         "insertions": reference["intervals"], "cycle": (capture, stop)})
        exact_lift(51, result)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
