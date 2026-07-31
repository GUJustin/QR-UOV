#!/usr/bin/env python3
"""Independent public replay of a reduced QR-UOV local-separator forgery run.

The verifier uses only the serialized public key, message, forged signature,
and attack certificate.  It reconstructs the affine core, reruns the full
split-before-lift/local-separator solver with the recorded separator, and
checks the final signature under the unchanged public map.
"""
from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Tuple

import qruov_local_separator_forgery as impl


def quad_from_json(d: Dict[str, Any]) -> impl.Quad:
    return impl.Quad(A=d["A"], b=d["b"], c=d["c"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("run", type=Path, nargs="?", default=Path("qruov_local_separator_forgery_run.json"))
    args = ap.parse_args()
    data = json.loads(args.run.read_text(encoding="utf-8"))
    if not data.get("success"):
        raise SystemExit("run does not contain a successful forgery")

    pub = data["public_instance"]
    q = int(pub["q"])
    pk = impl.ReducedPublicKey(
        q=q,
        V=int(pub["V"]),
        O=int(pub["O"]),
        public_matrices=pub["public_matrices"],
        hash_prefix=bytes.fromhex(pub["hash_prefix_hex"]),
    )
    message = bytes.fromhex(pub["message_hex"])
    salt = bytes.fromhex(data["salt_hex"])
    signature = data["signature"]
    target = impl.hash_to_field(pk.hash_prefix, message, salt, q, 3 * pk.O)
    assert target == data["target"]
    assert impl.verify(pk, message, (salt, signature))

    cert = data["successful_certificate"]
    assert cert is not None
    T = cert["target_transform"]
    mixed = impl.h.output_mix(pk.public_matrices, T, q)
    ty = impl.h.matvec(T, target, q)
    assert ty == [0] * (len(target) - 1) + [1]

    U = cert["isotropic_subspace_U"]
    r = len(target) - 3
    assert impl.h.rank(U, q) == r
    for A in mixed[:3]:
        assert all(x == 0 for row in impl.h.restrict_form(A, U, q) for x in row)

    z0 = cert["chart_z0"]
    B = cert["chart_kernel_B"]
    restricted = [impl.h.restrict_form(A, U, q) for A in mixed]
    f = [impl.quad_from_chart(A, z0, B, q) for A in restricted[3:-1]]
    H0 = impl.quad_from_chart(restricted[-1], z0, B, q)
    assert [asdict(x) for x in f] == cert["affine_system"]
    assert asdict(H0) == cert["scale_form_H0"]

    local = cert["local_separator"]
    lam: List[Tuple[int, int]] = [tuple(x) for x in local["lambda"]]  # type: ignore[list-item]
    c: Tuple[int, int] = tuple(local["c"])  # type: ignore[assignment]
    coords, reason, septries, degree, trace, replay = impl.local_separator_solve(
        f,
        H0,
        q,
        random.Random(0),
        max_separator_tries=1,
        forced_lambda=lam,
        forced_c=c,
    )
    assert reason == "success" and septries == 1 and replay is not None
    assert coords == local["affine_coordinates"]
    for key in (
        "C",
        "C0",
        "precision",
        "vandermonde_nodes",
        "start_roots",
        "lambda",
        "c",
        "q_valuation_at_target",
        "qbar",
        "b1bar",
        "psibar",
        "tau",
        "qprime_at_tau",
        "b1_at_tau",
        "psi_at_tau",
        "H0_at_point",
        "affine_coordinates",
        "crosscheck",
    ):
        assert replay[key] == local[key], key
    assert degree == len(local["qbar"]) - 1

    z = [
        (z0[i] + sum(B[i][j] * coords[j] for j in range(len(coords)))) % q
        for i in range(r)
    ]
    assert z == cert["projective_chart_point"]
    scale = int(cert["affine_scale"])
    zscaled = [(scale * x) % q for x in z]
    replay_signature = impl.h.matvec(U, zscaled, q)
    assert replay_signature == signature == cert["forged_signature_vector"]
    assert impl.h.eval_map(mixed, replay_signature, q) == ty
    assert impl.h.eval_map(pk.public_matrices, replay_signature, q) == target

    print("Public replay succeeded")
    print(f"q={q}, V={pk.V}, O={pk.O}, a={len(f)}, eliminant_degree={degree}")
    print(f"branch_lifts={trace.branch_lifts}, rational_reconstructions={trace.rational_reconstructions}")
    print(f"signature verifies={impl.verify(pk, message, (salt, replay_signature))}")


if __name__ == "__main__":
    main()
