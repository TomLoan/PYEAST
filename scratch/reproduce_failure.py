"""Isolate the `in_vivo_assembly` failure mode.

Finding
-------
The failure is NOT caused by the number of fragments. The full 20-part
violacein plasmid assembles fine in ~1.3 s (see build_plasmid.py).

It is caused by REPEATED HOMOLOGY BLOCKS. In vivo / homologous-recombination
assembly needs every junction to be a *unique* stretch of homology. As soon as
a part is reused (the same terminator or promoter driving several genes -- very
common in multi-gene yeast pathways), that duplicated ~250 bp body is a valid
recombination site in more than one place. pydna's assembly graph
(`assembly2.Assembly`) then contains many parallel cycles; it enumerates simple
cycles up to a hard cap of 10,000 (assembly2.py:1789) and raises:

    ValueError: Too many possible paths (more than 10000)

This script runs two otherwise-identical 11-part assemblies at limit=20
(the practical minimum homology for yeast IVA) to show the difference.

    uv run python reproduce_failure.py
"""

from __future__ import annotations

import os
import time

from pydna.readers import read
from pydna.design import primer_design, assembly_fragments
from pydna.amplify import pcr
from pydna.assembly2 import in_vivo_assembly, Assembly

HERE = os.path.dirname(os.path.abspath(__file__))
LIMIT = 20  # minimum homology used for yeast in vivo assembly


def load(name: str):
    rec = read(os.path.join(HERE, f"{name}.fasta"))
    rec.name = name
    return rec


def run_case(label: str, names: list[str]) -> None:
    parts = [load(n) for n in names]
    annealing = [primer_design(p) for p in parts]
    tailed = assembly_fragments(annealing, overlap=35, circular=True)
    amplicons = [pcr(a.forward_primer, a.reverse_primer, t)
                 for a, t in zip(tailed, parts)]

    edges = Assembly(amplicons, limit=LIMIT,
                     use_all_fragments=True).G.number_of_edges()

    print(f"\n{label}")
    print(f"  parts ({len(names)}): {', '.join(names)}")
    print(f"  overlap-graph edges at limit={LIMIT}: {edges}")
    t0 = time.perf_counter()
    try:
        products = in_vivo_assembly(amplicons, limit=LIMIT, circular_only=False)
        dt = time.perf_counter() - t0
        circular = [p for p in products if p.circular]
        print(f"  --> {len(products)} product(s) "
              f"({len(circular)} circular) in {dt:.2f}s  [OK]")
    except Exception as exc:  # noqa: BLE001
        dt = time.perf_counter() - t0
        print(f"  --> {type(exc).__name__}: {exc}  (after {dt:.2f}s)  [FAILED]")


if __name__ == "__main__":
    # Control: three cassettes, each with its own unique terminator.
    run_case(
        "CONTROL  (unique terminators)",
        ["ampR_ColE1", "2micron",
         "pPDC1", "VioA", "tCYC1",
         "pPGI1", "VioB", "tFBA1",
         "pPGK1", "VioC", "tGPM1"],
    )

    # Failure: identical layout, but tCYC1 is reused as the terminator for all
    # three cassettes. Same part count, same limit -- only the repeat differs.
    run_case(
        "FAILURE  (tCYC1 reused x3)",
        ["ampR_ColE1", "2micron",
         "pPDC1", "VioA", "tCYC1",
         "pPGI1", "VioB", "tCYC1",
         "pPGK1", "VioC", "tCYC1"],
    )

    print("\nConclusion: the blow-up tracks the number of repeated homology")
    print("junctions, not the number of fragments. Reused promoters/terminators")
    print("are the trigger; giving every junction unique homology avoids it.")
