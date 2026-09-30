"""Prepare a bounded second oblique search from first-round observable scores.

The explicit proposals were fixed after reading first-round O/X and Doppler
scores, before opening any of the three full truth density grids.
"""

from __future__ import annotations

import json
import os
import shutil

from oblique_model_family_cases import RUN, candidate


PROPOSALS = {
    "nequick": {
        "q1p04_hp10": dict(dq1=.04, dq2=0.0, dh=10.0),
        "q1p08": dict(dq1=.08, dq2=0.0, dh=0.0),
    },
    "chapman": {
        "q1p04_hp10": dict(dq1=.04, dq2=0.0, dh=10.0),
        "hp20": dict(dq1=0.0, dq2=0.0, dh=20.0),
    },
    "iri2016": {
        "q2m06_hp10": dict(dq1=0.0, dq2=-.06, dh=10.0),
        "q2m12": dict(dq1=0.0, dq2=-.12, dh=0.0),
        "q2m06_gq2p08": dict(dq1=0.0, dq2=-.06, dh=0.0,
                              q2_gradient=.08),
        "q2m06_gq2m08": dict(dq1=0.0, dq2=-.06, dh=0.0,
                              q2_gradient=-.08),
    },
}


def main() -> None:
    os.umask(0o077)
    tasks = []
    for case, proposals in PROPOSALS.items():
        first = RUN / case / "first_round_scores.json"
        if not first.exists():
            shutil.copyfile(RUN / case / "scores.json", first)
            first.chmod(0o600)
        scores = json.loads(first.read_text())
        if len(scores["candidates"]) != 7:
            raise ValueError(f"Full first-round scores missing for {case}")
        for name, params in proposals.items():
            candidate(case, name, **params)
            tasks.append({"case": case, "candidate": name})
    path = RUN / "second_round_tasks.json"
    path.write_text(json.dumps(tasks, indent=2) + "\n")
    path.chmod(0o600)
    print(f"Prepared {len(tasks)} bounded second-round oblique candidates")


if __name__ == "__main__":
    main()
