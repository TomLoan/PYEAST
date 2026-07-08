# Copyright CSIRO 2025. Thomas Loan
# See LICENSE for full GpLv2 license.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful;,
# but WITHOUT ANY WARRENTY; without even the implied warranty of
# MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# ===========================================================================

# src/pyeast/utils/pydna_utils.py
"""
pydna adapter for PYEAST.

Bridges PYEAST's own design functions (design_circular_primers / design_linear_primers,
get_templates / rationalize_templates) to pydna's PCR and in-vivo-assembly machinery so
the "assembly" step is a real homologous-recombination simulation rather than a blind
concatenation. Provides:

- amplicon construction with a real-template -> dummy -> constructed fallback ladder,
- a cheap ambiguity diagnostic (overlap-graph partner count) that runs before enumeration,
- a practice-informed PCR-specificity analysis (Tm of the exact 3' footprint, not match
  length),
- a wrapped in_vivo_assembly that degrades gracefully on the "too many paths" explosion,
- re-annotation of the assembled product with PYEAST_component features, handling the
  coordinate-origin difference between pydna and PYEAST.

All pydna imports live here so the rest of the package need not import pydna eagerly.

There are two independent "limit" parameters and they must never be conflated:
  * PCR anneal limit  - per-primer seed for building/finding amplicons. For building we use
    the designed footprint length (primer minus its overhang); for specificity *finding* we
    use pydna's permissive 13 bp seed to catch every candidate site.
  * Assembly overlap limit - the minimum homology the recombination will accept
    (~20-25 bp; 20 is the practical yeast HR minimum). Do NOT raise it to ease computation.
"""

# ===========================================================================

from __future__ import annotations

import os
from collections import defaultdict
from datetime import datetime
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt
from pydna.amplify import Anneal
from pydna.assembly2 import Assembly, in_vivo_assembly, pcr_assembly
from pydna.dseqrecord import Dseqrecord
from pydna.opencloning_models import CloningStrategy
from pydna.primer import Primer

from .path_utils import get_private_equivalent
from .sequence_utils import (
    INTEGRATION_CONSTRUCT_DESCRIPTION,
    PYEAST_COMPONENT_NOTE,
    TAR_CONSTRUCT_DESCRIPTION,
    is_pyeast_component,
)

# The permissive seed used only to *find* candidate priming sites for the specificity
# analysis (short enough to surface weak/partial sites, which are then judged by Tm).
SPECIFICITY_SEED = 13

# Default off-target Tm cutoff (deg C). A secondary product is only worth warning about
# when its weaker primer site clears this melting temperature. Tunable per call.
OFFTARGET_TM = 40.0

# Templates longer than this (bp) are treated as genome/chromosome-scale fallbacks; any
# specificity warning sourced from them is flagged low-confidence rather than trusted.
GENOME_SCALE_BP = 100_000

# Templates larger than this (bp) are trimmed to the amplified region (+ flank) before PCR,
# so the recorded cloning history embeds a local window rather than a whole chromosome.
# The amplicon sequence is unchanged (the part sits within the window), only the recorded
# template does. Genomic coordinates are not preserved.
TRIM_TEMPLATE_BP = 50_000
TRIM_FLANK_BP = 200


def anneal_limit(fp_seq, rp_seq, overhang: int) -> int:
    """Per-primer PCR seed = the designed 3' footprint (primer minus its overhang).

    Using the true footprint length (rather than pydna's fixed 13 bp default) avoids false
    "non-specific"/"no product" calls on large genome templates and repeat-containing parts.
    A floor of 12 keeps short GC-rich footprints working.
    """
    return max(12, min(len(str(fp_seq)), len(str(rp_seq))) - overhang)


def load_template_index(directory: str) -> dict[str, SeqRecord]:
    """Index every template record by record name across public and private dirs.

    get_templates returns only matching template *names*; pydna's pcr needs the actual
    sequence, so this mirrors get_templates' file scan and returns the SeqRecords keyed by
    ``record.name`` (later files win on a name clash, matching load order).

    Args:
        directory: Public template directory. Its private equivalent (if any) is also scanned.

    Returns:
        Dict mapping template record name to its SeqRecord.
    """
    index: dict[str, SeqRecord] = {}

    public_dir = Path(directory)
    try:
        private_dir = get_private_equivalent(public_dir)
    except ValueError:
        private_dir = None

    for search_dir in [public_dir] + ([private_dir] if private_dir else []):
        if not search_dir.exists():
            continue
        for filename in os.listdir(search_dir):
            file_path = search_dir / filename
            lower = str(file_path).lower()
            if lower.endswith((".fasta", ".fa", ".fsa")):
                fmt = "fasta"
            elif lower.endswith((".gb", ".gbk")):
                fmt = "genbank"
            else:
                continue
            for record in SeqIO.parse(file_path, fmt):
                index[record.name] = record

    return index


def _primer_pairs(primers: dict) -> list[tuple[str, object, str, object]]:
    """Split an ordered primer dict (2 per part: part0 F, part0 R, ...) into pairs."""
    items = list(primers.items())
    pairs = []
    for i in range(0, len(items) - 1, 2):
        (fp_name, fp_seq), (rp_name, rp_seq) = items[i], items[i + 1]
        pairs.append((fp_name, fp_seq, rp_name, rp_seq))
    return pairs


def _template_for_pcr(record: SeqRecord, part: SeqRecord, name: str) -> Dseqrecord:
    """Build the PCR template Dseqrecord, trimming genome-scale records to the amplified region.

    For a large template (e.g. a whole chromosome) the part is located within it and only a
    local window (part + TRIM_FLANK_BP either side) is kept, so the recorded history embeds a
    small region rather than the entire chromosome. The window keeps the template's name and
    any features falling inside it. Falls back to the full record when the part is not found.
    """
    def _named(seqrec) -> Dseqrecord:
        dseq = Dseqrecord(seqrec, circular=False)
        dseq.name = name
        return dseq

    if len(record.seq) <= TRIM_TEMPLATE_BP:
        return _named(record)

    tseq = str(record.seq).upper()
    length = len(part.seq)
    idx = tseq.find(str(part.seq).upper())
    if idx == -1:
        idx = tseq.find(str(part.seq.reverse_complement()).upper())
    if idx == -1:
        return _named(record)  # part not locatable; keep full template

    start = max(0, idx - TRIM_FLANK_BP)
    end = min(len(record.seq), idx + length + TRIM_FLANK_BP)
    return _named(record[start:end])


def _construct_amplicon(part: SeqRecord, fp_seq, rp_seq, fp, rp, overhang: int) -> Dseqrecord:
    """Deterministically build the intended amplicon when PCR is non-specific.

    Product = forward overhang + full part + reverse-complement of the reverse overhang,
    i.e. exactly what the tailed primers would yield. Mirrors PYEAST's historical
    concatenation so the assembly still closes.
    """
    fwd_oh = str(fp_seq)[:overhang]
    rev_oh_rc = str(Seq(str(rp_seq)[:overhang]).reverse_complement())
    amp = Dseqrecord(fwd_oh + str(part.seq) + rev_oh_rc)
    amp.name = part.id
    amp.template = Dseqrecord(str(part.seq), name=f"{part.id}_constructed")
    return amp


def _tag_component(amp, part: SeqRecord, index: int) -> None:
    """Tag the amplicon's part body with a PYEAST component feature.

    A standard misc_feature carrying a /note=PYEAST_component marker (so the tag survives
    round-trips through GenBank tools that coerce non-standard feature keys) and a /label
    of f"{part.id}_{index}". It propagates through in_vivo_assembly onto the assembled
    product - into both the .gb and the OpenCloning history JSON - already carrying its
    final label. The part sits forward in the amplicon (amplicons are built
    forward_overhang + part + reverse_overhang_rc, and PCR products are oriented
    forward-primer-first), so a forward search always locates it; skips silently if not.
    """
    idx = str(amp.seq).upper().find(str(part.seq).upper())
    if idx == -1:
        return
    amp.features.append(SeqFeature(
        FeatureLocation(idx, idx + len(part.seq)),
        type="misc_feature",
        qualifiers={"label": [f"{part.id}_{index}"], "note": [PYEAST_COMPONENT_NOTE]},
    ))


def build_amplicons(
    parts: list[SeqRecord],
    primers: dict,
    rationalized_templates: dict | None,
    template_index: dict | None,
    overhang: int,
) -> tuple[list, list[dict]]:
    """Make one pydna Amplicon per part via a real -> dummy -> constructed fallback ladder.

    Primers are paired to parts by dict order (2 per part: forward then reverse). Each part
    is amplified from its rationalized freezer template when available; otherwise from a
    dummy template made of the part sequence; otherwise deterministically constructed.
    History is therefore never empty.

    Args:
        parts: Parts in assembly order.
        primers: Ordered primer dict from design_circular/linear_primers.
        rationalized_templates: Optional {part.name: template_name}. Empty/None -> dummy path.
        template_index: Optional {template_name: SeqRecord} from load_template_index.
        overhang: Homology/overhang length used in the primer design.

    Returns:
        (amplicons, notes) where notes is a list of dicts describing any fallback taken.
    """
    rationalized_templates = rationalized_templates or {}
    template_index = template_index or {}
    pairs = _primer_pairs(primers)

    amplicons: list = []
    notes: list[dict] = []

    for i, (part, (fp_name, fp_seq, rp_name, rp_seq)) in enumerate(zip(parts, pairs)):
        fp = Primer(str(fp_seq), name=fp_name)
        rp = Primer(str(rp_seq), name=rp_name)
        lim = anneal_limit(fp_seq, rp_seq, overhang)

        template_id = rationalized_templates.get(part.name)
        has_real = (
            template_id not in (None, "Not found") and template_id in template_index
        )

        def _amplify(template):
            """Return the single expected amplicon, or None if non-specific / no product.

            Uses pcr_assembly (not pcr) so the amplicon carries a PCRSource - primers +
            template + locations - giving the assembled product a full OpenCloning history.
            add_primer_features keeps the primer_bind features PYEAST annotates with.
            """
            try:
                products = pcr_assembly(template, fp, rp, add_primer_features=True, limit=lim)
            except Exception:
                return None
            return products[0] if len(products) == 1 else None

        # 1) real freezer template. Build from the SeqRecord (not just its sequence string)
        # so the template's own features carry through PCR -> assembly into the final .gb.
        # Genome-scale templates are trimmed to the amplified region to keep history small.
        if has_real:
            template = _template_for_pcr(template_index[template_id], part, template_id)
            amp = _amplify(template)
            if amp is not None:
                amp.name = part.id
                _tag_component(amp, part, i)
                amplicons.append(amp)
                continue
            notes.append({
                "part": part.id,
                "kind": "real_nonspecific",
                "message": f"{part.id}: PCR on {template_id} was non-specific; "
                           f"fell back to the part sequence",
            })

        # 2) dummy template (part sequence only)
        dummy = Dseqrecord(str(part.seq))
        # Cap at GenBank's 16-char LOCUS name limit: OpenCloning serialises this template
        # to a GenBank string when building the cloning history, and an over-long name warns.
        # The part identity is carried by the amplicon (amp.name = part.id) below.
        dummy.name = f"{part.id}_dummy"[:16]
        amp = _amplify(dummy)
        if amp is not None:
            amp.name = part.id
            if not has_real:
                notes.append({
                    "part": part.id,
                    "kind": "dummy_template",
                    "message": f"{part.id}: no freezer template; amplified from part sequence",
                })
            _tag_component(amp, part, i)
            amplicons.append(amp)
            continue

        # 3) genuinely non-specific (repeat inside the part) -> construct it deterministically
        amp = _construct_amplicon(part, fp_seq, rp_seq, fp, rp, overhang)
        notes.append({
            "part": part.id,
            "kind": "constructed",
            "message": f"{part.id}: PCR non-specific; amplicon constructed deterministically",
        })
        _tag_component(amp, part, i)
        amplicons.append(amp)

    return amplicons, notes


def cloning_history_json(product) -> str | None:
    """Serialise the product's cloning history as an OpenCloning CloningStrategy JSON.

    Walks the pydna assembly source DAG (product -> InVivoAssemblySource -> PCRSource per
    amplicon -> primers + templates) via CloningStrategy.from_dseqrecords. Returns the JSON
    string, or None if the product carries no usable history (e.g. a constructed fallback).
    """
    if product is None:
        return None
    try:
        return CloningStrategy.from_dseqrecords([product]).model_dump_json()
    except Exception:
        return None


def _location_kind(loc, fragment_len: int, overhang: int) -> str:
    """Classify an overlap location on a fragment as terminal or internal."""
    margin = max(2, overhang)
    start = int(loc.start)
    end = int(loc.end)
    if start <= margin or end >= fragment_len - margin:
        return "terminal"
    return "internal"


def diagnose_assembly(amplicons: list, limit: int, overhang: int = 25) -> list[dict]:
    """Cheap pre-enumeration ambiguity diagnostic (Feature 1).

    Builds the overlap graph and flags any forward fragment that overlaps MORE than its two
    neighbours (>2 partners), which signals accidental part reuse or shared/internal homology
    that makes recombination ambiguous. Building the graph is fast even when full cycle
    enumeration would explode.

    Args:
        amplicons: pydna Amplicons in assembly order.
        limit: Assembly overlap limit (min homology accepted, ~20-25).
        overhang: Designed overhang length, used to classify terminal vs internal overlaps.

    Returns:
        A list of findings, one per ambiguous fragment:
        {part, fragment, partners: [names], overlaps: [{partner, location, kind}]}.
        Empty list means every fragment overlaps exactly its two neighbours (clean).
    """
    graph = Assembly(amplicons, limit=limit, use_all_fragments=True).G

    def frag_name(idx: int) -> str:
        i = abs(idx) - 1
        return amplicons[i].name if 0 <= i < len(amplicons) else str(idx)

    partners: dict[int, set[int]] = defaultdict(set)
    overlap_info: dict[int, list[dict]] = defaultdict(list)

    for u, v, data in graph.edges(data=True):
        if u > 0 and v > 0 and u != v:
            partners[u].add(v)
            locs = data.get("locations") or []
            if locs:
                frag_len = len(amplicons[u - 1])
                overlap_info[u].append({
                    "partner": frag_name(v),
                    "location": f"{int(locs[0].start)}..{int(locs[0].end)}",
                    "kind": _location_kind(locs[0], frag_len, overhang),
                })

    findings: list[dict] = []
    for u in sorted(partners):
        if len(partners[u]) > 2:
            findings.append({
                "part": frag_name(u),
                "fragment": u,
                "partners": sorted(frag_name(v) for v in partners[u]),
                "overlaps": overlap_info.get(u, []),
            })
    return findings


def _footprint_tm(seq) -> float:
    """Tm of an exact 3' footprint via Bio Tm_NN (matches PYEAST's design Tm scale).

    Returns NaN on the sequences Bio's tables choke on, so callers can skip rather than crash.
    """
    s = str(seq)
    if len(s) < 2:
        return float("nan")
    try:
        return float(mt.Tm_NN(Seq(s)))
    except Exception:
        return float("nan")


def pcr_specificity(
    parts: list[SeqRecord],
    primers: dict,
    rationalized_templates: dict | None,
    template_index: dict | None,
    overhang: int,
    cutoff: float = OFFTARGET_TM,
) -> list[dict]:
    """Practice-informed PCR-specificity analysis (Feature 2).

    Keeps pydna's permissive 13 bp seed to *find* every candidate priming site, then judges
    each site by whether it would actually prime at the bench = melting temperature of the
    exact 3' footprint. The intended product is the one closest to the expected size
    (part + both overhangs); any *other* product whose weaker footprint Tm clears ``cutoff``
    is surfaced as a soft off-target note. Warnings sourced from genome-scale fallback
    templates are flagged low-confidence.

    Returns:
        A list of warnings: {part, template, size, tm, low_confidence, message}.
    """
    rationalized_templates = rationalized_templates or {}
    template_index = template_index or {}
    pairs = _primer_pairs(primers)

    warnings: list[dict] = []

    for part, (fp_name, fp_seq, rp_name, rp_seq) in zip(parts, pairs):
        template_id = rationalized_templates.get(part.name)
        has_real = (
            template_id not in (None, "Not found") and template_id in template_index
        )
        # Specificity is only meaningful against a real freezer template (where the part sits
        # in flanking context). With no real template the "template" would be the bare part,
        # for which off-target scoring is not informative - skip it.
        if not has_real:
            continue

        fp = Primer(str(fp_seq), name=fp_name)
        rp = Primer(str(rp_seq), name=rp_name)
        template = Dseqrecord(str(template_index[template_id].seq))
        template.name = template_id
        tname = template_id

        low_confidence = len(template) > GENOME_SCALE_BP

        try:
            products = Anneal([fp, rp], template, limit=SPECIFICITY_SEED).products
        except Exception:
            continue
        if not products:
            continue

        infos = []
        for prod in products:
            f_tm = _footprint_tm(prod.forward_primer.footprint)
            r_tm = _footprint_tm(prod.reverse_primer.footprint)
            infos.append((len(prod), min(f_tm, r_tm)))

        expected = len(part.seq) + 2 * overhang
        intended = min(infos, key=lambda x: abs(x[0] - expected))
        for info in infos:
            if info is intended:
                continue
            size, min_tm = info
            if min_tm >= cutoff:
                conf = " (low-confidence: genome-scale template)" if low_confidence else ""
                warnings.append({
                    "part": part.id,
                    "template": tname,
                    "size": size,
                    "tm": min_tm,
                    "low_confidence": low_confidence,
                    "message": f"{part.id}: possible off-target ~{size} bp on {tname} "
                               f"(secondary site Tm {min_tm:.1f} C){conf}",
                })

    return warnings


class AssemblyExplosionError(Exception):
    """Raised internally when in_vivo_assembly hits the too-many-paths cap."""


def assemble_in_vivo(
    amplicons: list,
    limit: int,
    circular_only: bool = True,
    proceed_on_ambiguity: bool = True,
) -> tuple[object, list[dict]]:
    """Run pydna in_vivo_assembly, degrading gracefully on the path explosion.

    Repeated homology blocks (reused promoters/terminators) make the recombination graph
    dense; pydna caps simple-cycle enumeration at 10000 and raises
    ``ValueError: Too many possible paths``. Rather than crash, we fall back to the
    deterministic constructed assembly (concatenation of the amplicon bodies) and flag the
    result as truncated.

    Args:
        amplicons: pydna Amplicons in assembly order.
        limit: Assembly overlap limit (~20-25).
        circular_only: True for TAR (circular plasmid), False for linear integration.
        proceed_on_ambiguity: If False, re-raise the explosion instead of degrading.

    Returns:
        (product, notes). product is the chosen assembled Dseqrecord (or the constructed
        fallback). notes flags any degradation.
    """
    notes: list[dict] = []
    try:
        products = in_vivo_assembly(amplicons, limit=limit, circular_only=circular_only)
    except ValueError as exc:
        if "Too many possible paths" not in str(exc) or not proceed_on_ambiguity:
            raise
        notes.append({
            "kind": "truncated",
            "message": "assembly graph too dense (repeated homology); returning a "
                       "deterministically constructed product instead of enumerating",
        })
        return _fallback_product(amplicons, circular_only), notes

    if circular_only:
        candidates = [p for p in products if p.circular]
    else:
        candidates = [p for p in products if not p.circular] or list(products)

    if not candidates:
        notes.append({
            "kind": "no_product",
            "message": "in_vivo_assembly returned no matching product; returning a "
                       "deterministically constructed product",
        })
        return _fallback_product(amplicons, circular_only), notes

    product = max(candidates, key=len)
    return product, notes


def _fallback_product(amplicons: list, circular_only: bool) -> Dseqrecord:
    """Constructed product used when enumeration fails: concatenated amplicon bodies."""
    seq = "".join(str(a.seq) for a in amplicons)
    product = Dseqrecord(seq, circular=circular_only)
    return product


def _find_part_origin(product, first_part: SeqRecord, primers: dict, overhang: int) -> int:
    """Locate where the first part starts within the assembled product.

    in_vivo_assembly origins a circular product at the 5' end of the first primer, i.e. at
    the start of that primer's overhang (the last ``overhang`` bases of the previous part).
    The first part therefore starts one overhang in. We locate it robustly by searching for
    the forward primer's designed footprint (its 3' portion = the start of the first part),
    falling back to the fixed overhang offset if the search fails.
    """
    fp_seq = list(primers.values())[0]
    footprint = str(fp_seq)[overhang:]
    seq = str(product.seq)
    doubled = seq + seq[: len(footprint)]
    if footprint:
        idx = doubled.find(footprint)
        if idx != -1:
            return idx % len(seq)
    return overhang % len(seq) if len(seq) else 0


def annotate_circular_product(
    product,
    parts: list[SeqRecord],
    primers: dict,
    overhang: int,
    name: str = "assembly",
) -> SeqRecord:
    """Re-origin a circular pydna product to the first part and ensure component features.

    Handles the coordinate-origin gotcha: the product is rotated with ``.shifted(n)`` so the
    first part starts at position 0, matching PYEAST's convention. PYEAST component features
    normally arrive already on the product, propagated (and rotated by the shift) from the
    amplicons tagged in build_amplicons; primer_bind features pydna attached are likewise
    preserved. They are only synthesised here by cumulative part length as a fallback when a
    degraded (constructed/concatenated) product carries none.

    Returns the annotated product (a Dseqrecord, i.e. a SeqRecord subclass) ready for
    SeqIO.write(..., "genbank").
    """
    origin = _find_part_origin(product, parts[0], primers, overhang)
    shifted = product.shifted(origin) if origin else product

    # Fallback only: a degraded product lost the propagated features - rebuild them by
    # cumulative part length so the .gb still carries component annotations.
    if not any(is_pyeast_component(f) for f in shifted.features):
        current = 0
        for i, part in enumerate(parts):
            part_len = len(part.seq)
            shifted.features.append(SeqFeature(
                FeatureLocation(current, current + part_len),
                type="misc_feature",
                qualifiers={"label": [f"{part.id}_{i}"], "note": [PYEAST_COMPONENT_NOTE]},
            ))
            current += part_len

    _finalize(shifted, name, topology="circular", description=TAR_CONSTRUCT_DESCRIPTION)
    return shifted


def annotate_linear_product(
    product,
    parts: list[SeqRecord],
    name: str = "assembly",
) -> SeqRecord:
    """Ensure PYEAST component features on a linear integration product.

    The linear product's origin already sits at the first part (the upstream integration
    site's forward primer carries no overhang), so no re-origin is needed. Component features
    normally arrive already on the product, propagated from the amplicons tagged in
    build_amplicons (labelled f"{part.id}_{i}"; the up/downstream flanks are distinguished by
    their FASTA-derived part.id). They are only synthesised here by cumulative part length as
    a fallback when a degraded (constructed/concatenated) product carries none.
    """
    # Orient so the first part is at the 5' end (pick the matching strand).
    anchor = str(parts[0].seq)[:30]
    if anchor and anchor not in str(product.seq) and hasattr(product, "reverse_complement"):
        rc = product.reverse_complement()
        if anchor in str(rc.seq):
            product = rc

    # Fallback only: a degraded product lost the propagated features - rebuild them by
    # cumulative part length so the .gb still carries component annotations.
    if not any(is_pyeast_component(f) for f in product.features):
        current = 0
        for i, part in enumerate(parts):
            part_len = len(part.seq)
            product.features.append(SeqFeature(
                FeatureLocation(current, current + part_len),
                type="misc_feature",
                qualifiers={"label": [f"{part.id}_{i}"], "note": [PYEAST_COMPONENT_NOTE]},
            ))
            current += part_len

    _finalize(product, name, topology="linear", description=INTEGRATION_CONSTRUCT_DESCRIPTION)
    return product


def _finalize(record, name: str, topology: str, description: str) -> None:
    """Set the id/name/description and the GenBank annotations PYEAST outputs expect.

    ``description`` becomes the GenBank DEFINITION line; batch processing relies on it to
    recognise valid tar/integrate constructs, so it must not be left as pydna's default.
    """
    record.id = name
    record.name = name
    record.description = description
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = topology
    record.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()
