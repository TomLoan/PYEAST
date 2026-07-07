"""Sandbox: round-trip PYEAST-designed TAR primers through pydna.

Goal: prove the integration path before touching PYEAST itself.

  PYEAST parts (SeqRecord)  --design_circular_primers-->  tailed primers
        |                                                       |
        +--- resolve freezer template (real, else dummy) -------+
                              |
                          pydna.pcr  --> Amplicon (tails + primer features + template link)
                              |
                     in_vivo_assembly --> circular plasmid (+ provenance source)

Uses PYEAST's OWN functions (design_circular_primers, get_templates,
rationalize_templates) so we validate the real design, not a re-derivation.
Templates model "what's in the freezer": most promoters/terminators are
amplified from the BY4741 genome; synthetic parts (Vio genes) have no template
and get a dummy named <part>_missing_template so history is never empty.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# --- wire in PYEAST (real functions) without installing the whole CLI ---
PYEAST_DATA = Path("C:/Users/u5390772/Desktop/Code/PYEAST_data_TomL")
PYEAST_SRC = Path("C:/Users/u5390772/Desktop/Code/PYEAST Local/src")
os.environ["PYEAST_DATA_DIR"] = str(PYEAST_DATA)
sys.path.insert(0, str(PYEAST_SRC))

# pyeast.utils.__init__ eagerly imports the plotting stack (matplotlib,
# dna_features_viewer) we don't need here. Stub it so importing the primer/
# sequence utils doesn't require the whole viz dependency tree.
import types  # noqa: E402
sys.modules["pyeast.utils.visualisation"] = types.ModuleType("pyeast.utils.visualisation")

from Bio import SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402

from pyeast.utils.primer_utils import design_circular_primers  # noqa: E402
from pyeast.utils.sequence_utils import (  # noqa: E402
    get_templates,
    load_sequences,
    rationalize_templates,
)

from pydna.dseqrecord import Dseqrecord  # noqa: E402
from pydna.primer import Primer  # noqa: E402
from pydna.amplify import pcr  # noqa: E402
from pydna.assembly2 import Assembly, in_vivo_assembly  # noqa: E402

COMPONENT_DIR = PYEAST_DATA / "component_libraries" / "Saccharomyces_cerevisiae"
TEMPLATE_DIR = PYEAST_DATA / "templates"

# The violacein-pathway plasmid, in a circular order.
ORDER = [
    "AmpR_ColE1", "2Micron", "Leu2Mx",
    "pPDC1", "VioA", "tCYC1",
    "pPGI1", "VioB", "tFBA1",
    "pPGK1", "VioC", "tGPM1",
    "pTDH1", "VioD", "tPDC1",
    "pTDH3", "VioE", "tTDH2",
    "tTDH3", "tTef1",
]

OVERHANG = 25    # bp homology per primer -> 50 bp junction overlap (PYEAST default)
LIMIT = 25       # min overlap the assembler will accept
# NB on the PCR anneal limit: pydna.pcr defaults to a 13 bp seed, which is shorter
# than the real designed footprints and produces FALSE "non-specific"/"no product"
# calls on large templates and repeat-containing parts. The right value is the
# actual designed footprint length (primer minus its overhang), computed per part.


def anneal_limit(fp_seq, rp_seq):
    """Designed 3' footprint = primer minus the homology overhang; use the shorter."""
    return max(12, min(len(str(fp_seq)), len(str(rp_seq))) - OVERHANG)


def load_template_index(directory: Path) -> dict[str, SeqRecord]:
    """Index every template record by name (mirrors get_templates' file scan)."""
    index = {}
    for f in directory.iterdir():
        if f.suffix.lower() in (".fasta", ".fa", ".fsa"):
            fmt = "fasta"
        elif f.suffix.lower() in (".gb", ".gbk"):
            fmt = "genbank"
        else:
            continue
        for rec in SeqIO.parse(f, fmt):
            index[rec.name] = rec
    return index


def _construct_amplicon(part, fp_seq, rp_seq, fp, rp):
    """Deterministically build the intended amplicon when PCR is non-specific.

    Product = fwd overhang + full part + rev-comp of rev overhang. This is what
    the tailed primers would yield, and mirrors PYEAST's current concatenation.
    """
    fwd_oh = str(fp_seq)[:OVERHANG]
    rev_oh_rc = str(Seq(str(rp_seq)[:OVERHANG]).reverse_complement())
    amp = Dseqrecord(fwd_oh + str(part.seq) + rev_oh_rc)
    amp.name = part.id
    amp.template = Dseqrecord(str(part.seq), name=f"{part.id}_constructed")
    amp._pyeast_primers = [fp, rp]
    return amp


def build_amplicon(part, fp_name, fp_seq, rp_name, rp_seq, template_id, tindex):
    """Make a pydna Amplicon for one part: real template -> dummy -> constructed."""
    fp = Primer(str(fp_seq), name=fp_name)
    rp = Primer(str(rp_seq), name=rp_name)
    lim = anneal_limit(fp_seq, rp_seq)
    has_real = template_id not in (None, "Not found") and template_id in tindex

    # 1) real freezer template
    if has_real:
        template = Dseqrecord(str(tindex[template_id].seq))
        template.name = template_id
        try:
            amp = pcr(fp, rp, template, limit=lim)
            amp.name = part.id
            return amp, template_id, "PCR (real template)", None
        except Exception as exc:
            real_note = f"{template_id}: {str(exc).splitlines()[0]}"
    else:
        real_note = None

    # 2) dummy template (part sequence only)
    dummy = Dseqrecord(str(part.seq))
    dummy.name = f"{part.id}_missing_template"
    try:
        amp = pcr(fp, rp, dummy, limit=lim)
        amp.name = part.id
        status = "PCR (dummy)" if not has_real else "PCR (dummy; real non-specific)"
        return amp, dummy.name, status, real_note
    except Exception as exc:
        # 3) genuinely non-specific (repeat inside the part) -> construct it
        amp = _construct_amplicon(part, fp_seq, rp_seq, fp, rp)
        note = real_note or f"{part.id}: {str(exc).splitlines()[0]}"
        return amp, amp.template.name, "constructed (non-specific)", note


def main() -> int:
    print(f"Loading parts from {COMPONENT_DIR.name} ...")
    seqs = load_sequences(str(COMPONENT_DIR))
    missing = [n for n in ORDER if n not in seqs]
    if missing:
        print(f"  MISSING parts (not in library): {missing}")
        return 1
    parts = [seqs[n] for n in ORDER]

    print("Designing TAR primers (PYEAST design_circular_primers)...")
    primers = design_circular_primers(parts, target_tm=50, overhang_length=OVERHANG)
    primer_items = list(primers.items())  # ordered: part0 F, part0 R, part1 F, ...

    print("Resolving freezer templates (PYEAST get_templates)...")
    template_dict = get_templates(parts, str(TEMPLATE_DIR))
    chosen = rationalize_templates(template_dict)
    tindex = load_template_index(TEMPLATE_DIR)

    print("\nBuilding amplicons via pydna.pcr:\n")
    print(f"{'part':<12}{'amp':>6}  {'template':<26}{'status'}")
    print("-" * 80)
    amplicons = []
    notes = []
    for i, part in enumerate(parts):
        fp_name, fp_seq = primer_items[2 * i]
        rp_name, rp_seq = primer_items[2 * i + 1]
        template_id = chosen.get(part.name)
        amp, source_note, status, note = build_amplicon(
            part, fp_name, fp_seq, rp_name, rp_seq, template_id, tindex
        )
        amplicons.append(amp)
        print(f"{part.id:<12}{len(amp):>6}  {source_note:<26}{status}")
        if note:
            notes.append((part.id, note))

    if notes:
        print("\nPCR-specificity findings (pydna flagged; PYEAST concatenation would miss these):")
        for pid, note in notes:
            print(f"  - {note[:100]}")

    # --- cheap ambiguity diagnostic (the PYEAST design-time check) ---
    print("\nAmbiguity diagnostic (overlap partners per fragment):")
    G = Assembly(amplicons, limit=LIMIT, use_all_fragments=True).G
    from collections import defaultdict
    partners = defaultdict(set)
    for u, v, _ in G.edges(keys=True):
        if u > 0 and v > 0:
            partners[abs(u)].add(abs(v))
    ambiguous = [(i, ORDER[i - 1], sorted(partners.get(i, [])))
                 for i in range(1, len(parts) + 1) if len(partners.get(i, [])) > 2]
    if ambiguous:
        for _, name, p in ambiguous:
            print(f"  AMBIGUOUS: {name} overlaps fragments {p}")
    else:
        print("  clean - every fragment overlaps exactly its 2 neighbours")

    # --- assemble ---
    print(f"\nRunning in_vivo_assembly (limit={LIMIT}, circular_only)...")
    products = in_vivo_assembly(amplicons, limit=LIMIT, circular_only=True)
    circular = [p for p in products if p.circular]
    print(f"  products: {len(products)} ({len(circular)} circular)")
    for p in products:
        print(f"    {'circular' if p.circular else 'linear'}: {len(p)} bp | "
              f"source: {type(getattr(p, 'source', None)).__name__}")

    if circular:
        best = max(circular, key=len)
        print("\nProvenance chain of the assembled plasmid:")
        print(f"  plasmid: {len(best)} bp circular, source={type(best.source).__name__}")
        print(f"  built from {len(best.source.input)} amplicons; each carries "
              f"template + primer features:")
        for amp in amplicons[:5]:
            prims = amp.primers() if hasattr(amp, "primers") else getattr(amp, "_pyeast_primers", [])
            prset = [pr.name for pr in prims]
            feats = [f.qualifiers.get('label', ['?'])[0]
                     for f in amp.features if f.type == 'primer_bind']
            tname = amp.template.name if hasattr(amp, "template") else "?"
            print(f"    {amp.name:<10} <- template '{tname}' "
                  f"| primers {prset} | primer_bind feats {feats}")
        print("    ... (showing first 5)")

        out = Path("sandbox_assembly.gb")
        best.write(str(out))
        print(f"\n  wrote {out} ({len(best)} bp)")

        # one worked example of a real-template PCR figure (annealing context)
        from pydna.amplicon import Amplicon
        example = next((a for a in amplicons
                        if isinstance(a, Amplicon) and hasattr(a, "template")
                        and "missing_template" not in a.template.name), None)
        if example is not None:
            print(f"\nExample real-template PCR ({example.name} <- {example.template.name}):")
            print(example.figure()[:600])
    return 0


if __name__ == "__main__":
    sys.exit(main())
