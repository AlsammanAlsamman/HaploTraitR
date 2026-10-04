# Graph Report - HaploTraitR  (2026-10-04)

## Corpus Check
- 25 files · ~458,490 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 23 nodes · 21 edges · 4 communities
- Extraction: 100% EXTRACTED · 0% INFERRED · 0% AMBIGUOUS
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `9c1bf142`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- [[_COMMUNITY_Community 0|Community 0]]
- [[_COMMUNITY_Community 1|Community 1]]
- [[_COMMUNITY_Community 2|Community 2]]

## God Nodes (most connected - your core abstractions)
1. `HaploTraitR <img src="icon/logo.png" align="right" width="120"/>` - 10 edges
2. `Team Members` - 8 edges
3. `Example output` - 5 edges
4. `Installation` - 1 edges
5. `Quick start: one function call` - 1 edges
6. `What you get` - 1 edges
7. `How to read the results` - 1 edges
8. `Haplotype block: alleles, trait effect and LD` - 1 edges
9. `Trait by haplotype` - 1 edges
10. `Haplotype effects across blocks` - 1 edges

## Surprising Connections (you probably didn't know these)
- None detected - all connections are within the same source files.

## Import Cycles
- None detected.

## Communities (4 total, 0 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.22
Nodes (8): HaploTraitR <img src="icon/logo.png" align="right" width="120"/>, How to read the results, Installation, Method notes, Quick start: one function call, Settings, Step-by-step workflow (advanced), What you get

### Community 1 - "Community 1"
Cohesion: 0.25
Nodes (8): Dr. Alsamman Alsamman, Dr. Andrea Visioni, Dr. Outmane Bouhlal, Dr. Zakaria Kehel, Mr. Doaa Korkar, Mr. Khaled Helmy, Ms. Tamara Ortiz, Team Members

### Community 2 - "Community 2"
Cohesion: 0.40
Nodes (5): Example output, Genome overview, Haplotype block: alleles, trait effect and LD, Haplotype effects across blocks, Trait by haplotype

## Knowledge Gaps
- **18 isolated node(s):** `Installation`, `Quick start: one function call`, `What you get`, `How to read the results`, `Haplotype block: alleles, trait effect and LD` (+13 more)
  These have ≤1 connection - possible missing edges or undocumented components.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `HaploTraitR <img src="icon/logo.png" align="right" width="120"/>` connect `Community 0` to `Community 1`, `Community 2`?**
  _High betweenness centrality (0.745) - this node is a cross-community bridge._
- **Why does `Team Members` connect `Community 1` to `Community 0`?**
  _High betweenness centrality (0.515) - this node is a cross-community bridge._
- **Why does `Example output` connect `Community 2` to `Community 0`?**
  _High betweenness centrality (0.320) - this node is a cross-community bridge._
- **What connects `Installation`, `Quick start: one function call`, `What you get` to the rest of the system?**
  _18 weakly-connected nodes found - possible documentation gaps or missing edges._