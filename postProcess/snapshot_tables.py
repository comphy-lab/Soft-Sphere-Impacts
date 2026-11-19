"""Extract the force and footprint tables from a run's saved snapshots."""

import math
from pathlib import Path
import subprocess
import tempfile


def weber_label(value):
    weber = float(value)
    if not math.isfinite(weber) or weber <= 0:
        raise ValueError("We must be finite and positive")
    return format(weber, ".0f") if weber.is_integer() else repr(weber)


def force_summary_path(stem, case="bounce"):
    if case not in {"bounce", "bounce-ve", "bounce-ve-wi", "bounce-kv", "bounce-bingham"}:
        raise ValueError(f"Unknown case {case!r}")
    return Path("..") / f"{stem}_{case}.dat"


def extract_snapshot_table(helper, output):
    snapshots = [Path(f"intermediate/snapshot-{0.01 * index:5.4f}")
                 for index in range(100)]
    snapshots = [snapshot for snapshot in snapshots if snapshot.is_file()]
    if not snapshots:
        raise SystemExit("No snapshots found in intermediate/")

    # Put the temporary and final tables on the same filesystem so a failed
    # extraction leaves the previous table intact.
    with tempfile.TemporaryDirectory(prefix=".extract-", dir=".") as work:
        table = Path(work) / "table.dat"
        for snapshot in snapshots:
            subprocess.run([f"./{helper}", snapshot, table], check=True)
        table.replace(output)
