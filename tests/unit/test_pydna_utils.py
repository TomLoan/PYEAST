"""Unit tests for the pydna adapter and the TAR/integration engine swap."""

from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import json
import random

from pyeast.core.tar import TARDesigner
from pyeast.utils.primer_utils import design_circular_primers
from pyeast.utils.pydna_utils import (
    annotate_circular_product,
    assemble_in_vivo,
    build_amplicons,
    cloning_history_json,
    pcr_specificity,
)
from pyeast.utils.sequence_utils import is_pyeast_component, load_sequences, rationalize_templates

GOLDEN_TAR_INPUT = Path("tests/fixtures/golden_tar/input")
GOLDEN_TAR_EXPECTED = Path("tests/fixtures/golden_tar/expected/tar_test.gb")
GOLDEN_TAR_ORDER = ["pTEF1", "YeRFP", "tDIT1", "Ura3", "AmpR_ColE1", "2Micron"]


def _designer(order):
    seqs = load_sequences(GOLDEN_TAR_INPUT)
    d = TARDesigner(homology_length=25)
    d.available_sequences = seqs
    d.set_assembly_order(order)
    d.design_tar_primers()
    return d


def test_clean_assembly_diagnostic_clean_and_identity_preserved():
    """A clean design diagnoses clean and its assembled sequence matches the golden master."""
    d = _designer(GOLDEN_TAR_ORDER)

    diagnostics = d.diagnose()
    assert diagnostics == [], "clean assembly should have no ambiguity findings"

    asm = d.create_assembly()
    assert asm.annotations["topology"] == "circular"
    assert len(asm.seq) == 5722

    expected = SeqIO.read(GOLDEN_TAR_EXPECTED, "genbank")
    assert str(asm.seq).upper() == str(expected.seq).upper()

    components = [f for f in asm.features if is_pyeast_component(f)]
    primer_binds = [f for f in asm.features if f.type == "primer_bind"]
    assert len(components) == len(GOLDEN_TAR_ORDER)
    assert len(primer_binds) == 2 * len(GOLDEN_TAR_ORDER)


def test_reused_part_is_flagged_and_assembly_survives():
    """Reusing a part creates a repeated homology block that the diagnostic must flag."""
    order = ["pTEF1", "YeRFP", "tDIT1", "Ura3", "tDIT1", "AmpR_ColE1", "2Micron"]
    d = _designer(order)

    diagnostics = d.diagnose()
    assert diagnostics, "reused part should be flagged as ambiguous"
    assert "tDIT1" in {f["part"] for f in diagnostics}

    # Proceeding must not crash even on the ambiguous design.
    asm = d.create_assembly(proceed_on_ambiguity=True)
    assert len(asm.seq) > 0


def test_rationalize_templates_prefers_explicit_param():
    """An explicit preferred list steers selection over the shortest-match default."""
    template_dict = {"partA": ["some_long_genome_name", "pUC19"]}
    chosen = rationalize_templates(template_dict, preferred_templates=["pUC19"])
    assert chosen["partA"] == "pUC19"


def test_rationalize_templates_empty_config_uses_shortest_match():
    """With no preference (fixtures config is empty), the shorter template name wins on a tie."""
    chosen = rationalize_templates({"partA": ["AAAAAA", "BB"]})
    assert chosen["partA"] == "BB"


# Unique flanks so the part sits in genomic-style context with no accidental priming.
_FLANK_L = "GATTACAGGCATACGTAGCTAGGCATCGATCGTAGCATCGATCGGATCCATG" * 4
_FLANK_R = "CCTAGGCATCGGATCGATCGATGCATACGCATGCATCGATCGTAGCTAGCAT" * 4


def test_pcr_specificity_silent_on_single_copy_template():
    """A single-copy template (part embedded in unique flanks) yields no off-target warnings."""
    part = load_sequences(GOLDEN_TAR_INPUT)["pTEF1"]
    primers = design_circular_primers([part], target_tm=50, overhang_length=25)
    template = SeqRecord(Seq(_FLANK_L + str(part.seq) + _FLANK_R), name="T_single")
    warnings = pcr_specificity([part], primers, {part.name: "T_single"}, {"T_single": template}, 25)
    assert warnings == []


def test_pcr_specificity_flags_tandem_duplicate_offtarget():
    """A duplicated part in the template gives a real high-Tm secondary product -> warning fires."""
    part = load_sequences(GOLDEN_TAR_INPUT)["pTEF1"]
    primers = design_circular_primers([part], target_tm=50, overhang_length=25)
    dup = SeqRecord(
        Seq(_FLANK_L + str(part.seq) + _FLANK_R + str(part.seq) + _FLANK_L),
        name="T_dup",
    )
    warnings = pcr_specificity([part], primers, {part.name: "T_dup"}, {"T_dup": dup}, 25)
    assert warnings, "a duplicated part should produce an off-target warning"
    assert warnings[0]["part"] == "pTEF1"


def test_pcr_specificity_skips_parts_without_real_template():
    """No real template -> specificity is skipped (no spurious warnings on the dummy path)."""
    part = load_sequences(GOLDEN_TAR_INPUT)["pTEF1"]
    primers = design_circular_primers([part], target_tm=50, overhang_length=25)
    warnings = pcr_specificity([part], primers, {}, {}, 25)
    assert warnings == []


def _two_part_assembly_with_templates():
    """Build a clean 2-part circular assembly whose parts have real, annotated templates."""
    rng = random.Random(11)
    def rnd(n):
        return "".join(rng.choice("ACGT") for _ in range(n))
    flank = "GATTACAGGCATACGTAGCTAGGCATCGATCGTAGCATCGATCGGATCCATG"

    parts = []
    template_index = {}
    rationalized = {}
    for name in ("PartA", "PartB"):
        body = rnd(300)
        part = SeqRecord(Seq(body), id=name, name=name)
        parts.append(part)
        tname = f"T_{name}"
        template = SeqRecord(Seq(flank + body + flank), id=tname, name=tname)
        template.features.append(SeqFeature(
            FeatureLocation(len(flank) + 40, len(flank) + 120),
            type="CDS", qualifiers={"label": [f"{name}_marker"]},
        ))
        template_index[tname] = template
        rationalized[name] = tname

    primers = design_circular_primers(parts, target_tm=50, overhang_length=25)
    return parts, primers, rationalized, template_index


def test_cloning_history_json_captures_pcr_sources_and_primers():
    """The saved history is valid OpenCloning JSON with a PCR source per part + primers."""
    parts, primers, rationalized, template_index = _two_part_assembly_with_templates()
    amplicons, _ = build_amplicons(parts, primers, rationalized, template_index, 25)
    product, _ = assemble_in_vivo(amplicons, limit=25, circular_only=True)

    history = cloning_history_json(product)
    assert history, "expected a cloning history for a real-template assembly"
    doc = json.loads(history)

    source_types = [s.get("type") for s in doc["sources"]]
    assert "InVivoAssemblySource" in source_types
    assert source_types.count("PCRSource") == len(parts)
    assert len(doc["primers"]) == 2 * len(parts)


def test_cloning_history_json_includes_component_labels():
    """Component parts are tagged on the amplicons, so they appear in the history JSON."""
    parts, primers, rationalized, template_index = _two_part_assembly_with_templates()
    amplicons, _ = build_amplicons(parts, primers, rationalized, template_index, 25)
    product, _ = assemble_in_vivo(amplicons, limit=25, circular_only=True)

    history = cloning_history_json(product)
    assert history, "expected a cloning history for a real-template assembly"
    doc = json.loads(history)

    # The component marker + labels are serialised in the sequence nodes' GenBank content.
    all_content = "".join(s.get("file_content", "") for s in doc["sequences"])
    assert '/note="PYEAST_component"' in all_content
    for i, part in enumerate(parts):
        assert f"{part.id}_{i}" in all_content


def test_template_features_propagate_into_assembled_gb():
    """A template feature within the amplified region lands on the annotated assembly."""
    parts, primers, rationalized, template_index = _two_part_assembly_with_templates()
    amplicons, _ = build_amplicons(parts, primers, rationalized, template_index, 25)
    product, _ = assemble_in_vivo(amplicons, limit=25, circular_only=True)
    assembly = annotate_circular_product(product, parts, primers, 25)

    labels = {f.qualifiers.get("label", ["?"])[0] for f in assembly.features}
    assert "PartA_marker" in labels
    assert "PartB_marker" in labels
    # Exactly one PYEAST component feature per part (propagated from the amplicons; the
    # annotate fallback loop must not add a duplicate set).
    components = [f for f in assembly.features if is_pyeast_component(f)]
    assert len(components) == len(parts)
    # Components are standard misc_features carrying the /note marker (not a custom type).
    assert all(f.type == "misc_feature" for f in components)
    comp_labels = {f.qualifiers.get("label", ["?"])[0] for f in components}
    assert comp_labels == {f"{p.id}_{i}" for i, p in enumerate(parts)}
