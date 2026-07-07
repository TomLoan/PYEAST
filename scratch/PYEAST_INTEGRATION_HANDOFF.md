# PYEAST × pydna integration — handoff brief

Context for a fresh agent session working **inside `../PYEAST Local`**. This
work was scoped and de-risked in a sibling repo (`../PyDNA Assembly`) using
PYEAST's own functions against the real violacein data. Nothing in PYEAST was
changed yet. Read this before touching code.

## Goal

Add three pydna-powered capabilities to PYEAST's TAR flow, as **advisory /
validation layers** — none should hard-fail a design:

1. **Assembly-ambiguity diagnostic** — before assembling, tell the user which
   parts create ambiguous recombination (usually accidental part reuse or shared
   homology), *and where the offending overlap is*. Then offer: re-pick parts, or
   proceed anyway (intentional ambiguity / unavoidable repeats are legitimate).
2. **PCR-specificity note** — warn when a primer pair is likely to give a mixed
   product at the bench (off-target amplification), judged by **melting temp**,
   not match length.
3. **Provenance** — carry template → amplicon → plasmid history through the
   pydna objects (bonus; pydna gives it for free).

## Artifacts to copy from `../PyDNA Assembly` (reference implementations)

- `sandbox_pyeast_pydna.py` — the working adapter: PYEAST `SeqRecord`+primers →
  pydna `pcr` → `in_vivo_assembly`, with template resolution + dummy fallback +
  provenance. **This is the template for the adapter.**
- `pcr_specificity.py` — the working Tm-based specificity analysis.
- `reproduce_failure.py` — control vs reused-terminator; use as the basis for an
  ambiguity **test case**.
- `MEMORY.md` + `memory/pyeast-integration-scope.md` +
  `memory/pydna-ivassembly-failure.md` — the full findings log. The memory dir is
  per-project, so the new session won't see it automatically; copy the content
  or keep this repo open alongside.

## Object model & adapter (trivial)

- Parts: Bio `SeqRecord` from `load_sequences(dir)` → keyed by FASTA header id
  (NB: header id, not filename — e.g. `AmpR_ColE1`, `2Micron`, `tTef1`).
- Primers: `design_circular_primers(parts, target_tm, overhang_length)` →
  ordered dict `{name: Seq}`, **2 per part in list order** (part0 F, part0 R, …).
- Adapter: `Dseqrecord(str(seqrec.seq))`, `Primer(str(primer_seq), name=…)`,
  `pcr(fp, rp, template)`. `Dseqrecord` **is** a `SeqRecord` subclass, so pydna
  products flow straight into `SeqIO.write`.
- **Reuse PYEAST's own functions** (`design_circular_primers`, `get_templates`,
  `rationalize_templates`) — don't reimplement.
- **Do NOT copy the sandbox's `sys.path.insert`, `PYEAST_DATA_DIR`, or
  `sys.modules['pyeast.utils.visualisation']` stub** — those are sandbox-only
  hacks to import PYEAST from outside. Inside PYEAST none are needed, and pydna is
  already a declared dependency (`pydna>=5.5.0`).

## Feature 1 — ambiguity diagnostic

- Build the graph (cheap; only cycle enumeration is expensive):
  `from pydna.assembly2 import Assembly; G = Assembly(amplicons, limit=<20-25>, use_all_fragments=True).G`
- For each forward node, collect overlap partners; **>2 partners = ambiguous**.
  Report the offending part(s) **and the overlap location + whether it's terminal
  or internal** (from the edge `.locations` / source left/right locations) — the
  location is the high-value part; internal homology that truncates a part reads
  differently to a clean terminal reuse.
- Run this **before** `in_vivo_assembly`. If ambiguous, present the re-pick /
  proceed choice. If the user proceeds on a genuinely dense design, **wrap the
  `in_vivo_assembly` call** to catch `ValueError: Too many possible paths (more
  than 10000)` and degrade gracefully rather than crash (see
  `pydna-ivassembly-failure.md`).
- The violacein set is a clean control (every fragment overlaps exactly its 2
  neighbours); a reused-terminator variant is the positive case.

## Feature 2 — PCR specificity (Tm-based, not length-based)

- Keep pydna's permissive **seed `limit=13`** to *find* all candidate sites:
  `Anneal([fp, rp], template, limit=13).products`.
- Score each product by its **weaker** primer site:
  `min(Tm_NN(prod.forward_primer.footprint), Tm_NN(prod.reverse_primer.footprint))`
  using Bio `MeltingTemp.Tm_NN` (matches PYEAST's design-Tm scale).
- Intended product = the one closest to `len(part) + 2*overhang`; any *other*
  product whose min-Tm ≥ **off-target cutoff = 40 °C (user's call; make it a
  tunable param)** → surface a soft "you may see a ~N bp off-target band" note.
- Calibration that validated the cutoff: tFBA1/tGPM1 genomic secondary sites
  score Tm 19/25 °C → correctly silent; a repeat-rich template gives real 52.8 °C
  secondary sites → flagged. Use the **exact 3′ footprint** Tm (robust); Bio's
  mismatch-aware NN tables have gaps and error on some neighbours — don't rely on
  them.
- **Specificity runs against the rationalized template, and that's correct.**
  Score against `get_templates` → `rationalize_templates` output (the same
  template used for lab instructions). The `preferred_templates` list
  (`['pUC19','pYES2','pESC-TRP']`, sequence_utils.py ~L179) is the intended lever
  for the "sequence occurs in many places" case — prefer the specific, easy,
  artifact-free plasmid. Trust it; don't add a second disambiguation layer.
  - **Boundary condition to handle:** the preferred list can only act when a
    preferred plasmid is BOTH present in the templates folder AND an exact
    substring of the part (`get_templates` matches exactly). When neither holds,
    rationalize falls back to whatever did match — which can be the wild-type
    genome (`BY4741_Toronto_2012.fsa`). That happened to `2Micron`: pYES2/pESC-TRP
    aren't in the folder and no present plasmid contains this 2µ variant, so it
    fell to the genomic 2µ (inherently repeat-rich) and produced a spurious
    52.8 °C "off-target". Fixes: (a) keep real freezer plasmids in the templates
    data; (b) **label genome-fallback templates as low-confidence** so those
    warnings are flagged, not trusted. Do NOT treat a genome-sourced warning the
    same as a plasmid-sourced one.

### Scoped task: move `preferred_templates` to config

- Currently hardcoded in `sequence_utils.py` (~L179) as
  `['pUC19','pYES2','pESC-TRP']` — invisible to users and now stale.
- Move it into `~/PYEAST/config.yaml` (read via the existing `PyeastConfig`),
  **empty by default**. `rationalize_templates` reads the list from config
  instead of the constant. Empty ⇒ the shortest-matching-template sorter alone
  decides (the intended default; drops the stale steering with no behaviour change
  for typical parts).
- Add an explanatory comment in the scaffolded config. **Critical wording:** the
  list matches template **record names** (the sequence name inside the `.gb`,
  e.g. `pUC19`, `pRG205MX`) — NOT the filename (`puc19.gb`) and NOT the part name.
  A part name or filename here silently never matches. Suggested comment:
  "preferred_templates: when a part is found in several templates, prefer these.
  Use the template record name (inside the .gb, e.g. pUC19), not the filename or
  part name. Format: [\"pUC19\",\"pRG205MX\"]. Empty = use shortest match."
- This is the short-term answer to discoverability. Do NOT add template choice to
  the instructions table (already info-dense). Per-run `{part: template}` override
  and runtime transparency are deferred.
- Building amplicons (separate from the seed above): use a **per-primer anneal
  limit = the designed footprint length** (`primer minus the overhang`), not a
  fixed value. Fixed 13 over-flags on big templates; fixed 16 breaks legit short
  GC-rich primers (AmpR). Fallback ladder for non-specific/no-product: real
  template → dummy(part) → deterministically constructed amplicon
  (`fwd_oh + part + rc(rev_oh)`), flagging the finding.

## Feature 3 — provenance + annotations (the coordinate gotcha)

- pydna carries history for free: amplicon `.template`, `.forward_primer`,
  `.reverse_primer`; the assembled product's `.source` (`InVivoAssemblySource`,
  an OpenCloning LinkML model) serialises via `.model_dump_json()`. Template
  annotations pass through template → PCR → in-vivo and land in the saved `.gb`.
- **PYEAST's own part annotations (`PYEAST_component` features) still need to be
  added** to the assembled product — the pydna product only carries `primer_bind`
  features by default.
- **THE gotcha — coordinate origin differs:** `in_vivo_assembly` indexes the
  circular plasmid with the **5′ end of the first primer at position 1**, whereas
  PYEAST's existing `assemble_parts_circular` annotation logic indexes from the
  **first base of the first part** (offset by the first primer's overhang tail,
  ~`overhang_length` bp, sitting inside the last↔first junction). Two options:
  1. Adjust the position offsets in the annotation-adding for-loops (already done
     once before; "not that hard"), or
  2. Re-index the pydna product to the part origin first (pydna `Dseqrecord`
     circular sequences support `.shifted(n)`), then reuse PYEAST's existing
     annotation function unchanged.
  Either works; (2) may be cleaner if you want one annotation code path.

## Where to hook (cli/main.py)

- Interactive TAR: `run_tar_interactive_mode` (starts ~L623); slot the diagnostic
  + specificity report + re-pick/proceed choice at the **"Proceed with assembly?"**
  confirm (~L656). There's a parallel integration path at ~L729.
- Programmatic: `TARDesigner.design()` (core/tar.py) — add a diagnostic step and a
  `proceed_on_ambiguity: bool` flag so headless callers get the same behaviour.
- `assemble_parts_circular` already has an adjacent-junction fuzzy-alignment
  warning (`>0.8` similarity, mismatch-tolerant). Decide: keep it alongside the
  exact-match graph (they catch different things) or replace it. The graph is
  strictly more powerful for non-adjacent/internal homology.

## Two `limit`s — never conflate

- **PCR anneal limit** — per-footprint for building amplicons; seed 13 for
  specificity *finding*.
- **Assembly overlap limit** — ~20–25 for `in_vivo_assembly` and the diagnostic
  graph (20 = practical yeast HR minimum; keep it — don't raise it to make the
  computation easier, that hides real biology).

## Suggested tests

- Clean assembly (violacein set) → 1 circular product, diagnostic clean, no
  specificity warnings on real cloning templates.
- Reused-terminator variant → diagnostic flags the reused part with location;
  `in_vivo_assembly` wrapped so the >10000-paths case is handled gracefully.
- A part whose primer has a real off-target on its chosen template → specificity
  note fires at 40 °C; an A/T-rich low-Tm secondary site → stays silent.
