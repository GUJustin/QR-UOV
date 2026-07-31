#!/usr/bin/env python3
"""Deterministic multi-instance regression for the reduced QR-UOV attack.

Each row generates a fresh reduced QR-UOV key from a distinct seed, produces
and verifies an honest signature, runs the local-separator forgery, verifies
the forgery under the unchanged public map, serializes the public transcript,
and invokes the public-only replay checker on that transcript.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path

import qruov_local_separator_forgery as impl


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--count", type=int, default=8)
    ap.add_argument("--seed-start", type=int, default=20260731)
    ap.add_argument("--output", type=Path, default=Path("qruov_local_separator_regression.json"))
    ap.add_argument("--transcript-dir", type=Path, default=Path("regression_transcripts"))
    ap.add_argument("--q", type=int, default=13)
    ap.add_argument("--V", type=int, default=3)
    ap.add_argument("--O", type=int, default=2)
    ap.add_argument("--max-salts", type=int, default=4)
    ap.add_argument("--max-outer-attempts", type=int, default=5)
    ap.add_argument("--max-chart-tries", type=int, default=3)
    ap.add_argument("--max-separator-tries", type=int, default=20)
    args = ap.parse_args()

    args.transcript_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    all_ok = True
    suite_start = time.perf_counter()
    verifier = Path(__file__).resolve().with_name("verify_qruov_local_separator_run.py")

    for offset in range(args.count):
        seed = args.seed_start + offset
        message = f"local-separator regression seed {seed}".encode()
        t0 = time.perf_counter()
        result = impl.run_forgery(
            q=args.q,
            V=args.V,
            O=args.O,
            seed=seed,
            message=message,
            max_salts=args.max_salts,
            max_outer_attempts=args.max_outer_attempts,
            max_chart_tries=args.max_chart_tries,
            max_separator_tries=args.max_separator_tries,
        )
        elapsed = time.perf_counter() - t0
        transcript = args.transcript_dir / f"seed_{seed}.json"
        payload = asdict(result)
        payload["elapsed_seconds"] = elapsed
        transcript.write_text(json.dumps(payload, indent=2), encoding="utf-8")

        replay_ok = False
        replay_stdout = ""
        replay_stderr = ""
        if result.success:
            proc = subprocess.run(
                [sys.executable, str(verifier), str(transcript)],
                cwd=Path(__file__).resolve().parent,
                text=True,
                capture_output=True,
                check=False,
            )
            replay_ok = proc.returncode == 0
            replay_stdout = proc.stdout.strip()
            replay_stderr = proc.stderr.strip()

        successful_record = result.attempt_records[-1] if result.success and result.attempt_records else None
        trace = successful_record["trace"] if successful_record else {}
        row = {
            "seed": seed,
            "success": result.success,
            "honest_signature_verified": result.honest_signature_verified,
            "forgery_verified": result.forgery_verified,
            "public_replay_succeeded": replay_ok,
            "salt_attempts": result.attempts,
            "attack_records": len(result.attempt_records),
            "elapsed_seconds": elapsed,
            "branch_lifts": trace.get("branch_lifts", 0),
            "rational_reconstructions": trace.get("rational_reconstructions", 0),
            "target_eliminant_degree": successful_record.get("target_eliminant_degree", 0) if successful_record else 0,
            "exhaustive_base_roots": trace.get("exhaustive_base_roots", 0),
            "exhaustive_extension_roots": trace.get("exhaustive_extension_roots", 0),
            "transcript": str(transcript),
            "replay_stdout": replay_stdout,
            "replay_stderr": replay_stderr,
        }
        ok = result.success and result.honest_signature_verified and result.forgery_verified and replay_ok
        all_ok = all_ok and ok
        rows.append(row)
        print(
            f"seed={seed} success={result.success} verify={result.forgery_verified} "
            f"replay={replay_ok} salts={result.attempts} records={len(result.attempt_records)} "
            f"elapsed={elapsed:.3f}s"
        )

    summary = {
        "parameters": {"q": args.q, "ell": 3, "V": args.V, "O": args.O},
        "count": args.count,
        "all_succeeded": all_ok,
        "successful_forgery_count": sum(bool(r["success"] and r["forgery_verified"]) for r in rows),
        "public_replay_count": sum(bool(r["public_replay_succeeded"]) for r in rows),
        "total_elapsed_seconds": time.perf_counter() - suite_start,
        "maximum_single_run_seconds": max((float(r["elapsed_seconds"]) for r in rows), default=0.0),
        "maximum_salt_attempts": max((int(r["salt_attempts"]) for r in rows), default=0),
        "rows": rows,
    }
    args.output.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "rows"}, indent=2))
    print(f"wrote {args.output}")
    if not all_ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
