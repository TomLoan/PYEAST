"""Practice-informed PCR specificity analysis for the TAR amplicons.

Rather than judging a binding site by how long its exact match is, we judge it
by whether it would actually PRIME at the bench: keep pydna's permissive 13 bp
seed to FIND every candidate site, then score each site's exact 3' footprint by
melting temperature (Bio Tm_NN, matching PYEAST's design Tm scale). A secondary
PCR product is only worth warning about if its *weaker* primer site clears the
off-target Tm cutoff.

Calibration (from real parts):
  - tFBA1 spurious site  = TTTGTAAAAATAT  Tm 19.0 C  -> noise, ignore
  - 2Micron off-target   = ATCTGTGCTTCAT  Tm 32.6 C  -> a 1269 bp side-product
Default cutoff 40 C (tunable); at 40 C the 2Micron 32.6 C site does NOT fire.
"""

from __future__ import annotations

# Reuse the sandbox's PYEAST wiring + template resolution (its module-level code
# sets sys.path / PYEAST_DATA_DIR and stubs the viz import).
from sandbox_pyeast_pydna import (
    ORDER,
    OVERHANG,
    TEMPLATE_DIR,
    anneal_limit,
    load_sequences,
    COMPONENT_DIR,
    design_circular_primers,
    get_templates,
    rationalize_templates,
    load_template_index,
    Dseqrecord,
    Primer,
)

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from pydna.amplify import Anneal

SEED = 13            # permissive: find every candidate site
OFFTARGET_TM = 40.0  # cutoff for warning about a priming off-target (tunable)


def footprint_tm(seq) -> float:
    s = str(seq)
    if len(s) < 2:
        return float("nan")
    try:
        return float(mt.Tm_NN(Seq(s)))
    except Exception:
        return float("nan")


def analyse(fp, rp, template, expected_size, cutoff=OFFTARGET_TM):
    """Return (intended, offtargets, others) products as (size, min_tm, f_tm, r_tm).

    The intended amplicon is the product closest to the expected size (part +
    both overhangs); every other product is a candidate off-target, flagged when
    its weaker primer site clears the Tm cutoff.
    """
    an = Anneal([fp, rp], template, limit=SEED)
    infos = []
    for prod in an.products:
        f_tm = footprint_tm(prod.forward_primer.footprint)
        r_tm = footprint_tm(prod.reverse_primer.footprint)
        infos.append((len(prod), min(f_tm, r_tm), f_tm, r_tm))
    if not infos:
        return None, [], []
    intended = min(infos, key=lambda x: abs(x[0] - expected_size))
    others = [i for i in infos if i is not intended]
    offtargets = [i for i in others if i[1] >= cutoff]
    return intended, offtargets, others


def resolve_template(part, tid, tindex):
    if tid not in (None, "Not found") and tid in tindex:
        t = Dseqrecord(str(tindex[tid].seq)); t.name = tid
        return t, tid
    t = Dseqrecord(str(part.seq)); t.name = f"{part.id}_missing_template"
    return t, t.name


def main():
    seqs = load_sequences(str(COMPONENT_DIR))
    parts = [seqs[n] for n in ORDER]
    primers = list(design_circular_primers(parts, target_tm=50, overhang_length=OVERHANG).items())
    chosen = rationalize_templates(get_templates(parts, str(TEMPLATE_DIR)))
    tindex = load_template_index(TEMPLATE_DIR)

    print(f"PCR specificity (seed={SEED}, off-target Tm cutoff={OFFTARGET_TM} C)\n")
    print(f"{'part':<12}{'template':<24}{'intended':>10}{'2ndary sites (size@minTm)'}")
    print("-" * 88)
    warned = []
    for i, part in enumerate(parts):
        fp = Primer(str(primers[2 * i][1]), name=primers[2 * i][0])
        rp = Primer(str(primers[2 * i + 1][1]), name=primers[2 * i + 1][0])
        template, tname = resolve_template(part, chosen.get(part.name), tindex)
        expected = len(part.seq) + 2 * OVERHANG
        try:
            intended, offtargets, others = analyse(fp, rp, template, expected)
        except Exception as exc:
            print(f"{part.id:<12}{tname:<24}  analyse error: {str(exc).splitlines()[0][:40]}")
            continue
        # show any secondary products, marking those above cutoff
        sec = ", ".join(f"{sz}bp@{mtm:.0f}{'*' if mtm >= OFFTARGET_TM else ''}"
                        for sz, mtm, _, _ in others) or "-"
        flag = "  <== OFF-TARGET" if offtargets else ""
        isz, itm, _, _ = intended
        print(f"{part.id:<12}{tname:<24}{f'{isz}bp@{itm:.0f}':>10}  {sec}{flag}")
        if offtargets:
            warned.append((part.id, offtargets))

    print("\n* = secondary product whose weaker primer site clears the "
          f"{OFFTARGET_TM:.0f} C cutoff (likely visible band).")
    if warned:
        print("\nWarnings to surface to the user:")
        for pid, offs in warned:
            for sz, mtm, ftm, rtm in offs:
                print(f"  - {pid}: possible off-target ~{sz} bp "
                      f"(secondary site Tm {mtm:.1f} C)")
    else:
        print("\nNo off-target products above cutoff across the assembly.")


if __name__ == "__main__":
    main()
