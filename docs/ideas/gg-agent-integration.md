# Idea: Agent integration for the Golden Gate (`gg`) workflow

Status: parked idea, not scheduled. Captured 2026-07-08.

## Summary

Add an optional LLM-agent layer that helps the user design Golden Gate (MoClo/YTK
level-1) assemblies by turning a natural-language goal into a validated, correctly
ordered part list, and by interpreting simulation failures. The deterministic
dnacauldron simulation remains the source of truth; the agent only proposes inputs
and explains outputs.

Guiding principle: **agent proposes, dnacauldron disposes.** The agent works the
front of the workflow (intent -> grammar -> ordered part list) and the failure-
interpretation step. It must never touch part-to-plasmid mapping, the simulation
itself, or instruction generation.

## Motivation

The current `gg` workflow (`ggDesigner` in `src/pyeast/core/gg.py`, driven by
`run_gg_interactive_mode` in `src/pyeast/cli/main.py`) is six steps:

1. Load sequences from a components directory.
2. User types the ordered part list (space-separated; `/` and `/allX` for multiplex).
3. Map each part to the plasmid file containing it (exact sequence match).
4. Run the dnacauldron `Type2sRestrictionAssembly` simulation.
5. Validate (errors, warnings, construct count).
6. Write report + liquid-handler instructions.

Steps 3-6 are deterministic and mechanical; an LLM adds nothing there except
hallucination risk. The genuine cognitive load is at step 2: parts must be in the
correct MoClo type order with compatible fusion sites, and the tool currently just
trusts the user to know that (the only safety net is the reminder line
"Remember that golden gate parts need to be assembled in order").

## Where an agent adds value

- **Intent -> valid ordered part list.** e.g. "strong constitutive GFP cassette" ->
  a correctly ordered, fusion-site-compatible selection, checked *before* a
  simulation is spent finding out it is wrong.
- **Combinatorial library design.** Translate "all promoters x all terminators"
  into the right `/all` / `/` selection and report the resulting construct count.
- **Error interpretation.** Turn dnacauldron's cryptic `"Mixtures found"` /
  construct-count-mismatch warnings into plain-language causes and fixes
  (e.g. "parts X and Y share a fusion site; swap the connector").
- **Multi-level orchestration (strongest long-term case).** lvl0 -> lvl1 -> lvl2 is
  bookkeeping-heavy; tracking which sub-assemblies feed which is exactly where a
  human loses the thread and an agent does not. (Current code is level-1 only.)

## Where an agent does NOT help (explicit non-goals)

- Single level-1 assembly by someone who already knows their parts: typing a few
  names is faster than conversing. Value only appears as combinatorial size,
  assembly levels, or user unfamiliarity increase.
- Part-to-plasmid mapping (exact sequence match) - deterministic, no agent.
- The assembly simulation - dnacauldron is the ground truth; never let the model
  guess assembly outcomes.
- Liquid-handler instruction generation - pure formatting, no agent.

## Feasibility notes

- `ggDesigner` is cleanly decoupled from the CLI prompts, so the agent layer should
  call the class methods directly (`load_and_get_sequences`, `get_plasmid_names`,
  `gg_assembly`, `gg_save_output`, `gg_instructions`) rather than try to drive the
  interactive TUI.
- The agent will NOT infer fusion-site compatibility from the Python source - that
  knowledge lives in the plasmid GenBank annotations and the MoClo/YTK standard.
  It must be given the YTK grammar explicitly (fixed part-type positions) as
  context or as a lookup tool. That grammar is small and well-defined, so this is
  tractable.

## Suggested first slice (de-risk before investing)

Build the narrowest thing that proves the value: a **MoClo-grammar
validator / ordering helper**.

- Input: available part types (+ a goal, optionally).
- Output: a validated, ordered part list, OR an explanation of why a given order is
  invalid (wrong type order / incompatible fusion sites).
- It feeds its output into the existing `gg` flow unchanged.

If that feels useful in practice, extend in order:

1. Error interpretation over `_validate_assembly_results` output.
2. Combinatorial library intent -> `/all` / `/` selection.
3. Multi-level orchestration.

## Routing: integrating with the existing agent

This should slot into the existing `run_agent` loop as a fifth design tool
(`design_gg`) alongside `design_tar` / `design_integration` / `design_replacement`
/ `design_deletion`. Routing between design methods is ALREADY how that agent works:
there is no dispatcher choosing methods, the model picks based on the "Design
Heuristics" block in the system prompt plus each tool's `description` field. So
adding gg is a new branch of a solved problem, not new machinery. Structurally it
mirrors deletion/replacement (its own instruction generation, not batched via
`run_batch`, so `batch_eligible: false`).

### Feasibility is not preference (the core distinction)

The choice of gg vs tar/integrate has two separable axes and they must not be
conflated:

- **Feasibility (technical, the agent may decide):** hard constraints that rule a
  method in or out. e.g. gg needs MoClo/YTK level-0 parts in plasmids with
  compatible fusion sites; a cassette that must land at a genomic locus is
  integration, not a plasmid method. These are gates, not opinions.
- **Preference (human/lab, the agent must NOT silently decide):** among the
  methods that are technically viable, which one to use depends on the lab's skill
  set, instruments, tooling, and individual researcher habit. These preferences
  are legitimate, they vary between labs and between researchers, and there is
  often no universally "correct" answer.

Inventory (what parts exist, in what form) is a **feasibility** signal, not a
preference signal. Having a YTK kit on the shelf bounds what is possible; it does
not mean gg is the right call - the lab might have the parts but no Type IIS
workflow, or prefer to PCR them into a TAR assembly. Earlier framing over-loaded
inventory into the routing decision; it belongs only to the feasibility half.

### The design: decide feasibility, elicit preference once, then apply silently

Avoid both failure modes: a tool that auto-picks as if one method were universally
best (dishonest / wrong for many labs), AND a tool that interrogates the user about
their general philosophy every session and rebuilds a heuristic from scratch (a
token furnace, and annoying). The middle path:

1. **Narrow to the feasible set** from hard constraints / inventory. The agent
   decides this.
2. **If exactly one method is viable,** use it (state which, briefly).
3. **If more than one is viable,** consult a durable, stored preference
   (see below). If a stored preference covers this case, apply it silently.
4. **Only if the set is still ambiguous and no stored preference applies,** ask -
   but ask concretely about THIS construct ("your parts support Golden Gate or
   TAR; which do you want?"), not abstractly about philosophy, and offer to
   remember the answer. One turn, not a survey.

### Storing preference so it is asked once, not every time

The anti-token-furnace mechanism is persistence: capture the lab/researcher
preference as data (e.g. an `assembly_preferences` block in `~/PYEAST/config.yaml`),
set once - either explicitly by the user or offered by the agent when it first
learns the preference - and reused silently thereafter. Externalising it as config
rather than baking it into the system prompt also respects that different agents /
models may weigh these trade-offs differently; the preference is the user's data,
not the model's opinion.

### Mis-routes should be cheap and self-correcting

Perfect routing is not required if a wrong pick bounces back. Two mechanisms the
existing agent already has:

- **Validating tool errors.** Mirror `_find_missing_parts`: `design_gg` returns a
  structured `{"success": false, ...}` when parts are not level-0 / lack fusion
  sites, and the agent loop feeds that back so the model retries with a different
  method.
- **The plan-approval checkpoint.** The multi-modification workflow already makes
  the model state its plan and wait for approval. That is the natural place for it
  to say "I'll use Golden Gate because ..." and for the user to redirect before
  anything runs - so ambiguity is cheap.

## Open questions

- Where does the authoritative YTK part-type / fusion-site grammar come from -
  hard-coded table, or parsed from the plasmid annotations?
- Concrete shape of the stored preference: what does an `assembly_preferences`
  schema look like, and how fine-grained is it (per construct type? per part
  class? a single default with overrides)?
- When the agent offers to remember a preference, at what scope does it save -
  per user, per lab/data-dir, per project?
- Guardrails: how to make it structurally impossible for the agent to bypass the
  simulation as the source of truth?
