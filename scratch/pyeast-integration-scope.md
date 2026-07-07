---
name: pyeast-integration-scope
description: How PYEAST is structured for integrating pydna in_vivo_assembly diagnostics
metadata: 
  node_type: memory
  type: project
  originSessionId: 64ba15d8-d3d4-44a6-9905-84c4fa367928
---

Scoped 2026-07-03 (read-only). Goal: add pydna in_vivo_assembly-based ambiguity diagnostics to PYEAST's TAR cloning. See [[pydna-ivassembly-failure]].

**Repo:** `../PYEAST Local` (i.e. C:/Users/u5390772/Desktop/Code/PYEAST Local), git repo, package `pyeast`. Test data: `../PYEAST_data_TomL/component_libraries/Saccharomyces_cerevisiae/` has the violacein parts (VioA-E + promoters/terminators). **pydna>=5.5.0 is already a dependency** — no new deps.

**TAR data model:**
- Parts = Biopython `SeqRecord` (.seq/.id/.name); `load_sequences` -> `available_sequences` dict; ordered into `assembly_sequences` list by user selection.
- Primers = `design_circular_primers()` -> dict{name: Seq}, 2 per part, keys `{id}_{i}F/R`. Each = annealing region (part first/last 50bp, Tm-trimmed ~50C) + 25bp neighbour overhang. 25bp/side = 50bp junction overlap.
- `assemble_parts_circular` (sequence_utils.py:370) just CONCATENATES parts in given order + annotates primers + does an ADJACENT-junction-only Bio.Align similarity>0.8 warning. Not a recombination sim; blind to non-adjacent/internal homology.

**Integration plan:**
1. Adapter: SeqRecord->`Dseqrecord(rec)`; primer Seq->`Primer(str(seq),id=name)`; `pcr(fp,rp,template)`->Amplicon with real designed tails. Pair primers to parts via existing `find_matching_primer`.
2. New `TARDesigner.diagnose_assembly(limit=20)`: build `Assembly(amps,limit,use_all_fragments=True).G`, flag fragments with >2 overlap partners, report partner + terminal/internal location from edge `locations`.
3. CLI hook: cli/main.py `run_tar_interactive_mode` at ~line 652-656, after `designer.design()` and at existing `click.confirm("Proceed with assembly?")`. Branch: re-pick order vs proceed-anyway. API: field on TARResult + `on_ambiguity` flag.
4. Phase 2 (optional): augment (not replace) assemble_parts_circular with bounded in_vivo_assembly for real product/misassembly enumeration + provenance.

**Key decision:** PYEAST current check = fuzzy alignment (mismatch-tolerant); pydna common_sub_strings = EXACT >=limit only. They flag different sets — consider keeping both. First slice = adapter + diagnose_assembly + CLI hook, no change to existing assembly output.

**Full handoff brief written to repo `PYEAST_INTEGRATION_HANDOFF.md`** (../PyDNA Assembly). User will do the integration in a fresh session inside ../PYEAST Local. Coordinate gotcha to remember: `in_vivo_assembly` indexes the circular plasmid with the 5' end of the FIRST PRIMER at position 1; PYEAST's `assemble_parts_circular` indexes from the first base of the first PART (offset ~overhang_length, inside the last<->first junction). So PYEAST part annotations (`PYEAST_component` features) still need re-adding to the pydna product — either adjust the annotation for-loop offsets, or `Dseqrecord.shifted(n)` to re-origin then reuse the existing annotation fn. Template annotations otherwise pass through template->PCR->in_vivo into the saved .gb. `Dseqrecord` IS a Bio SeqRecord subclass; product carries only primer_bind feats by default.

**SANDBOX VALIDATED (2026-07-03):** `sandbox_pyeast_pydna.py` runs PYEAST's real `design_circular_primers`/`get_templates` on the 20-part violacein set through pydna `pcr` -> `in_vivo_assembly` = 1 circular 18024bp plasmid with InVivoAssemblySource provenance. Confirmed:
- Adapter is trivial: SeqRecord->Dseqrecord, primer Seq->Primer(name=), `pcr(fp,rp,template)`; pair primers to parts by dict order (2/part). To import pyeast utils without the CLI/plotting deps, stub `sys.modules['pyeast.utils.visualisation']` before import and set PYEAST_DATA_DIR (needs pandas, primer3-py, pyyaml, matplotlib in the sandbox venv).
- Template-as-freezer model works: real templates link natively (`amp.template.name`) — promoters/terminators PCR from genome chromosomes, backbone from pUC19/pRG205MX; synthetic Vio genes -> dummy Dseqrecord named `<part>_missing_template`. Always something in history.
- Primers auto-attach as `primer_bind` features when amplicon comes from `pcr` (label = primer.name).
- **PCR anneal limit must be per-part = designed footprint length (primer minus 25bp overhang), NOT pydna's default 13.** The default 13bp seed is shorter than real footprints and gives FALSE "non-specific"/"no product" calls on large genome templates and repeat-containing parts. Corrected finding: at the true footprint length, AmpR (15/18-mer, pUC19), tFBA1 (41/29-mer, chr11), tGPM1, and even 2Micron-vs-its-own-part all give exactly ONE product — every earlier "specificity finding" was a seed artifact. `anneal_limit = max(12, min(len(fp),len(rp)) - 25)`.
- Only genuinely irreducible case: 2Micron amplified from its 6318bp *plasmid* template (extra repeat copies beyond the part) stays non-specific -> falls back to dummy and still assembles. This is a template-choice matter (pYES2/pESC-Trp have shorter origins that work), not a primer/design bug. User's principle: fix by making part EDGES unique (extend boundary into backbone), not by tweaking primers. Don't over-flag PCR non-specificity as errors — most amplify fine in the lab.
- Fallback ladder in sandbox: real template -> dummy(part) -> deterministically constructed amplicon (fwd_oh+part+rc(rev_oh)).

**PCR SPECIFICITY = practice-informed Tm analysis (not footprint length).** User's steer: keep seed=13 to FIND all candidate sites, then judge each by whether it would PRIME at the bench = melting temp. Impl (`pcr_specificity.py`): `Anneal([fp,rp],template,limit=13).products`; for each product min-Tm = min(footprint_tm(fwd), footprint_tm(rev)) via Bio `mt.Tm_NN` (matches PYEAST design scale). Intended product = closest to expected size (part+2*overhang); every other product with min-Tm >= cutoff is an off-target warning. **Cutoff = 40 C (user's choice, tunable param).** Validated on real templates: 2Micron warns (real 2-micron plasmid has repeated origin copies -> genuine 52.8 C secondary site -> ~4016bp mixed band, matches lab); tFBA1/tGPM1 genome-scale secondary products score Tm 19/25 C -> correctly silent. Use exact-3'-footprint Tm (robust/computable); Bio mismatch-aware Tm tables have gaps (errors on some neighbors) so don't rely on full mismatch Tm. Surface as a soft "you may see off-target amplification (~Nbp)" note, user can still proceed.
