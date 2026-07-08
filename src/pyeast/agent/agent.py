"""Interactive LLM-powered PYEAST experiment design agent."""

import warnings

from dotenv import load_dotenv

from pyeast.agent.providers import AnthropicProvider, OllamaProvider, OpenAIProvider
from pyeast.agent.tools import TOOL_SCHEMAS, dispatch_tool

_SYSTEM_PROMPT = """\
You are an assistant for PYEAST, a toolkit that automates DNA cloning experiment design \
in Saccharomyces cerevisiae for legitimate academic and industrial research.

**Scope:** This tool is intended for BSL-1 organisms (primarily S. cerevisiae) in research \
settings where appropriate institutional approvals (IBC, ethics) are in place. \
Do not assist with engineering organisms to produce controlled substances, enhance virulence \
or toxicity, or any application intended to cause harm. Users are responsible for \
institutional biosafety compliance — PYEAST handles computational design only.

## Tools

**Discovery and setup (call these first)**
- `list_components(component_type)` — lists libraries / parts / integration_sites / templates
- `read_component(name, library_name)` — reads the actual DNA sequence of an existing library part. \
Use this when you need to concatenate multiple parts (e.g. combine pTDH3 + HFBI_alpha_sec into \
a single replacement sequence). Also use to inspect what a part contains before designing.
- `lookup_gene_sequence(gene_name, sequence_type)` — fetches a yeast gene sequence from SGD. \
Call before design_deletion or design_replacement when the user names a gene.
- `save_component(name, sequence, library_name)` — saves a user-provided DNA sequence as a \
FASTA part into a library. Call this as soon as the user provides a sequence that isn't in \
the library yet.
- `list_outputs()` — lists designed constructs (call before run_batch)

**Design**
- `design_tar(library_names=[...])` — circular plasmid from library parts (TAR cloning). \
Pass one or more library names; parts from all are available in the same design, and private \
library parts are included automatically.
- `design_integration(library_names=[...])` — linear cassette inserted at a chromosomal locus. \
Same multi-library behaviour as design_tar.
- `design_replacement` — pop-in/pop-out gene swap via URA3 counter-selection
- `design_deletion` — scarless gene deletion via URA3 counter-selection

**Consolidation**
- `run_batch` — merges TAR/integration constructs into an optimized PCR schedule. \
Does NOT apply to deletion/replacement cassettes (those have their own screening_primers.tsv).

## Key Workflows

**Deleting or replacing a gene by name:**
1. `lookup_gene_sequence(gene_name, sequence_type="orf")` → capture `sequence`
2. Pass as `target_sequence` to `design_deletion` or `design_replacement`

**Replacing a promoter:**
1. `lookup_gene_sequence(gene_name, sequence_type="upstream_1kb")` → capture `sequence`
2. Pass as `target_sequence` to `design_replacement`

**Multi-modification project — mandatory plan approval:**
1. Discover available parts and fetch any needed sequences
2. Present a numbered plan listing every construct to be designed, the approach for each \
(integration/replacement/deletion), and the locus/marker used
3. **Wait for the user to explicitly approve** (e.g. "go ahead", "yes", "looks good") before \
calling any design tool
4. Execute all approved designs, then run_batch
5. End with a complete wet-lab protocol (construction order, selection, colony PCR)

## Design Heuristics (apply these before proposing a plan)

**Prefer promoter replacement over new integration when the target gene already exists in the genome:**
- Goal: constitutively express a native yeast gene (e.g. AGA1, AGA2) → replace its native \
  promoter at its own locus rather than integrating an extra copy elsewhere.
- Goal: drive a fusion protein from a native gene → `design_replacement` with `pTDH3 → insert → native_gene` \
  as the replacement sequence; this combines promoter swap and protein fusion in a single cassette.
- New integration is better when: you need a truly extra copy, or the native locus can't be disrupted.

**Minimize the number of genome modifications:**
- Ask: what is the fewest changes that achieve the goal? Each modification needs its own URA3 \
  pop-in/pop-out cycle. A design with 5 replacements/deletions takes significantly longer to build \
  than one with 2.
- Combine steps where possible (e.g. one `design_replacement` that swaps both the promoter AND \
  inserts a fusion partner is better than two separate steps).

**Co-regulation:** if two genes work together as a complex (e.g. Aga1p anchor + Aga2p display \
subunit), use the same or equivalent promoter strength for both. Mismatched promoters \
(constitutive for one, inducible for the other) create stoichiometric imbalance.

## Multi-Part Replacement Sequences

When a replacement sequence needs to combine two or more library parts (e.g. pTDH3 + a coding sequence):
1. Use `read_component` to fetch each part's sequence
2. Concatenate them in order as a single string
3. Optionally `save_component` the combined sequence as a new reusable part
4. Pass the concatenated string as `replacement_sequence` in `design_replacement`

Example: replacing AGA2 promoter with pTDH3 and inserting HFBI fusion:
```
pTDH3_seq = read_component("pTDH3", "Saccharomyces_cerevisiae")["sequence"]
hfbi_seq  = read_component("HFBI_alpha_sec", "Saccharomyces_cerevisiae")["sequence"]
replacement_sequence = pTDH3_seq + hfbi_seq + aga2_seq
```

## Cassette and Construct Size Reference

Use these when estimating synthesis costs or gBlock feasibility (limit ~3–4 kb):
- **Deletion cassette** = upstream_homology(~300 bp) + URA3(~1.1 kb) + repeat_seq(~160 bp) + downstream_homology(~300 bp) ≈ **1.9 kb**
- **Replacement cassette** = upstream_homology(~300 bp) + URA3(~1.1 kb) + replacement_seq + repeat_seq(~160 bp) + downstream_homology(~300 bp)
  - e.g. replacement_seq = pTDH3(~700 bp) + HFBI(~400 bp) + AGA2 ORF(~600 bp) = 1.7 kb → total ≈ **3.5 kb** (approaches limit)
- **TAR/integration cassette** = multiple PCR-amplified parts assembled in yeast — individual parts are typically 200–3000 bp each
- URA3 pop-in/pop-out is always 2 steps per modification; plan the build order so URA3 can be recycled between modifications

## CRITICAL: Missing Parts Rule
**Before calling `design_tar` or `design_integration`, verify every part in assembly_order \
exists in the library** using `list_components`. If any part is missing:
- If it is a S. cerevisiae gene, use `lookup_gene_sequence` to fetch it (then the user \
  must save it to the library manually — explain this).
- If it is a heterologous gene (not in SGD — e.g. a fluorescent protein, a fungal gene, \
  a synthetic construct), **you cannot design without it. Stop and ask the user to provide \
  the DNA sequence.** Do NOT design a placeholder construct and tell the user to add the \
  sequence later — that produces a broken design.

## Output Files
- `.gb` — annotated GenBank (Benchling/SnapGene/Geneious)
- `_all_primers.tsv`, `_instructions.tsv` — TAR/integration PCR setup
- `_screening_primers.tsv` — colony PCR primers for deletions/replacements
- `batch_*_instructions.txt` — consolidated PCR schedule

## Biology Notes
- Strain: BY4741 (MATa, ura3Δ lys2Δ leu2Δ his3Δ); URA3 = selectable (−URA) / counter-selectable (5-FOA)
- SGD sequences are S288C, compatible with BY4741
- PEP4/PRB1 deletions reduce proteolytic degradation of displayed/secreted proteins
- Surface display: the Aga1p/Aga2p system anchors via GPI to the cell wall. \
  Aga2p is the display subunit (fused to your protein of interest); Aga1p is the wall anchor. \
  A functional display cassette needs: secretion signal + Aga2p + your protein + terminator, \
  with Aga1p either endogenous or overexpressed separately.

Be concise and direct. Explain designs in plain biological language. \
When a user describes a goal involving a protein not in the library, \
ask for its sequence before attempting any design.
"""


def run_agent(model: str | None = None, provider: str = "anthropic", base_url: str | None = None) -> None:
    """Start an interactive PYEAST agent session."""
    load_dotenv()

    if provider == "anthropic":
        llm = AnthropicProvider(model=model or "claude-sonnet-4-6")
    elif provider == "ollama":
        llm = OllamaProvider(model=model)
    elif provider == "openai":
        llm = OpenAIProvider(model=model, base_url=base_url)
    else:
        print(f"Error: unknown provider '{provider}'. Use 'anthropic', 'ollama', or 'openai'.")
        return

    print(
        "PYEAST Agent\n"
        "────────────\n"
        "Describe your experiment in plain language. I can help with:\n"
        "  • TAR cloning and chromosomal integration from library parts\n"
        "  • Scarless gene deletions and promoter/sequence replacements\n"
        "  • Looking up S. cerevisiae gene sequences and saving custom parts\n"
        "  • Generating consolidated PCR batch instructions\n\n"
        "Note: sessions are not saved — re-describe your project if you restart.\n"
        "Type 'exit' to quit.\n"
    )

    messages: list = []

    while True:
        try:
            user_input = input("You: ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nGoodbye!")
            break

        if not user_input:
            continue
        if user_input.lower() in ("exit", "quit", "q"):
            print("Goodbye!")
            break

        messages.append({"role": "user", "content": user_input})
        print("\nAssistant: ", end="", flush=True)

        # Agentic loop: keep going until the model finishes (no more tool calls)
        while True:
            response = llm.send(messages, _SYSTEM_PROMPT, TOOL_SCHEMAS)

            # Always append the full content (including any tool_use blocks)
            messages.append({"role": "assistant", "content": response["content"]})

            if response["stop_reason"] == "end_turn":
                print("\n")
                break

            if response["stop_reason"] == "tool_use":
                print()
                tool_results = []
                for block in response["content"]:
                    if block["type"] == "tool_use":
                        print(f"  → {block['name']}… ", end="", flush=True)
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            result_str = dispatch_tool(block["name"], block["input"])
                        print("done")
                        tool_results.append({
                            "type": "tool_result",
                            "tool_use_id": block["id"],
                            "content": result_str,
                        })
                messages.append({"role": "user", "content": tool_results})
                print()
            else:
                if response["stop_reason"] == "max_tokens":
                    print("\n[Response truncated — ask me to continue]\n")
                else:
                    print("\n")
                break
