---
name: pydna-ivassembly-failure
description: "Why pydna assembly2.in_vivo_assembly fails on the user's yeast assemblies"
metadata: 
  node_type: memory
  type: project
  originSessionId: 64ba15d8-d3d4-44a6-9905-84c4fa367928
---

User is evaluating whether pydna's `assembly2.in_vivo_assembly` can replace bespoke yeast-assembly logic in another repo. They reported it "failing for large assemblies."

Root cause (reproduced 2026-07-03 in this repo, pydna 5.5.16): the blow-up is driven by **repeated homology blocks, not fragment count**. The full 20-part violacein plasmid assembles in ~1.3 s when every junction has unique homology. But reusing a part (e.g. the same terminator/promoter across cassettes — common in multi-gene yeast pathways) makes the duplicated ~250 bp body a valid recombination site in multiple places. The assembly graph then has many parallel simple cycles; pydna enumerates them capped at 10,000 (`assembly2.py:1789`, `limit_iterator(nx.cycles.simple_cycles(self.G), 10000)`) and raises `ValueError: Too many possible paths (more than 10000)`.

**Why:** the user believed part count was the problem; it's actually graph density from repeats. limit=20 is their deliberate choice — yeast really does recombine at ~20 bp even when longer homology is present, and the "alternative" products (e.g. ones that loop out a costly cassette like violacein) are biologically real and worth simulating. So DO NOT advise raising `limit` to make it compute — that hides real biology. The target tool is **PYEAST**; it needs to handle the explosion *gracefully*, not avoid it.

**How to apply (graceful-handling design):**
1. Diagnose cheaply BEFORE enumerating: build `Assembly(amps, limit).G` (fast) and count each fragment's forward overlap partners; any fragment with >2 partners = ambiguous junction. Pinpoints the offending repeat in ms even when full enumeration would blow up.
2. Bounded/ranked enumeration instead of the hard 10000 cap: `get_circular_assemblies(max_assemblies=N)` and `get_possible_assembly_number()` exist; enumerate lazily with a time+count budget, return partial results + a "truncated" flag, rank biologically-interesting products (smallest circles, cassette-dropouts) first.
3. Isolation: run enumeration in a killable subprocess (multiprocessing) with a wall-clock budget — threads can't be killed, which is why build_plasmid.py has to os._exit.

**Provenance/history:** pydna 5.x products carry `.source` (OpenCloning LinkML models: InVivoAssemblySource, PCRSource, primers, templates) serializable via `.model_dump_json()`. PYEAST's history store = serialize that source DAG; no custom format needed. See `reproduce_failure.py` and `build_plasmid.py` in the repo.
