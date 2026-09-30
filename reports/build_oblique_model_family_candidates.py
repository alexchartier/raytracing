"""Construct a model-agnostic first-round oblique topside coordinate search."""

from __future__ import annotations

import json
import os
from pathlib import Path

from oblique_model_family_cases import RUN, CASES, candidate

PROPOSALS = {
    "q1p04": {"dq1": .04, "dq2": 0.0, "dh": 0.0},
    "q1m04": {"dq1": -.04, "dq2": 0.0, "dh": 0.0},
    "q2p06": {"dq1": 0.0, "dq2": .06, "dh": 0.0},
    "q2m06": {"dq1": 0.0, "dq2": -.06, "dh": 0.0},
    "hp10": {"dq1": 0.0, "dq2": 0.0, "dh": 10.0},
    "hm10": {"dq1": 0.0, "dq2": 0.0, "dh": -10.0},
}


def main() -> None:
    os.umask(0o077)
    tasks = []
    for case in CASES:
        for name, changes in PROPOSALS.items():
            candidate(case, name, **changes)
            tasks.append({"case": case, "candidate": name})
    path = RUN / "first_round_tasks.json"
    path.write_text(json.dumps(tasks, indent=2) + "\n")
    path.chmod(0o600)
    print(f"Created {len(tasks)} candidate grids from vertical ionograms and in-situ samples")


if __name__ == "__main__":
    main()
