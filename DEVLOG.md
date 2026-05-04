# pyPept Dev Log

Dense notes for session continuity. Not a tutorial — just enough to get up to speed fast.

---

## Current state (2026-05-04, updated)

**Branch:** master  
**Committed:** 100% library round-trip (1001/1001 monomers), C-attachment carboxyl convention, all slot renumbering.

### What the codebase is now

| File | What changed |
|------|--------------|
| `src/pyPept/sequence.py` | Bond validation (`_check_bond_chemistry`), CABILN inline attachment (`_expand_inline_caps`), old-BILN detection, attachment-point helpers |
| `src/pyPept/molecule.py` | `__restore_and_remove_rgroups` (Kekulize-first), `SanitizeMol` wrapped as `ValueError` |
| `src/pyPept/interfaces/monomer_pipeline.py` | Full SMARTS pipeline: `pre_activate`, `build_library_from_csv`, `import_helm_sdf`, `normalize_input`, `find_sidechain_slots` (CABILN slot order) |
| `src/pyPept/interfaces/map_monomers.py` | Deprecation shim only — redirects to `monomer_pipeline` |
| `src/pyPept/data/monomers.sdf` | **Replaced** — 1001-monomer CHUCKLES-format file, CABILN slot convention |
| `src/pyPept/data/monomers.csv` | Authoring source for the 52-monomer library |
| `tests/test_bond_validation_and_assembly.py` | 572 tests (571 pass, 1 xfail) — full coverage including library round-trip |

---

## CABILN (Chemistry Aware BILN)

**CABILN** = Chemistry Aware BILN. Extension of BILN (Boehringer Ingelheim Line Notation).

### What's new vs BILN

1. **R-group numbering** is pinned: R1=backbone N, R2=backbone C(=O), R3=backbone-N mod (α-N methylation etc.), R4+=sidechain slots. Old BILN had no convention; slot numbers varied by monomer.

2. **Inline attachment notation** `.Token(y,z)` — attaches cap/branch inline. `y`=R-group on host, `z`=R-group on the cap/branch partner. Auto-assigns a bond ID ≥ 100. The cap is appended as a pendant chain.

3. **Crosslink notation** `.!n(y,z)` — marks a crosslink endpoint. Two occurrences of `.!n` form a bond. Second occurrence must declare the *inverse* R-groups (swap y↔z). Conflict raises `ValueError`.

4. **Old BILN crosslinks** `Token(bid,rg)` (no `.` prefix) are **detected and rejected** with a helpful error showing the CABILN equivalent. No auto-converter.

### CABILN R-group convention

| Slot | Meaning |
|------|---------|
| R1 | Backbone N (N-terminal attachment) |
| R2 | Backbone C=O (C-terminal attachment) |
| R3 | Backbone-N modification (α-N methylation, formylation, etc.) |
| R4 | First sidechain reactive atom (thiol, ε-amine, guanidinium, sidechain COOH, etc.) |
| R5 | Second sidechain H or second reactive atom |
| R6 | Third sidechain slot (Arg only) |

**Slot summary:**

| AA  | R3 | R4 | R5 | R6 |
|-----|----|----|----|----|
| Gly/Ala | backbone-N mod | — | — | — |
| Pro | (none — ring N, only R1/R2) | — | — | — |
| Cys | backbone-N mod | thiol | — | — |
| Lys | backbone-N mod | ε-NH2 (1st H) | ε-NH2 (2nd H) | — |
| Arg | backbone-N mod | guanidinium terminal NH2 | guanidinium ε-N | terminal NH2 (2nd H) |
| Asp/Glu | backbone-N mod | sidechain COOH | — | — |
| His | backbone-N mod | imidazole NH | — | — |
| Ser/Thr/Tyr | backbone-N mod | sidechain OH | — | — |

### Inline attachment examples

```
ac-C.!1(4,4)-G-C.!1(4,4)-am       # disulfide crosslink via Cys thiol (R4)
fmoc-C.trt(4,2)-am                 # Cys thiol protected with Trt (trt R2)
fmoc-C.acm(4,2)-am                 # Cys thiol protected with Acm (acm R2)
fmoc-R.pbf(4,2)-am                 # Arg guanidinium protected with Pbf (pbf R2)
K.ac(2,1)                          # ac cap on Lys C-terminal R2
```

---

## Critical concepts

### CHUCKLES convention
Isotope-labelled dummy atoms encode attachment slots:
- `[1*]` = R1 = backbone N
- `[2*]` = R2 = backbone C=O
- `[3*]` = R3 = backbone-N mod
- `[4*]` = R4 = first sidechain slot
- etc.

The dummy atom's **neighbour** is the heavy atom where the inter-monomer bond forms.  
`_attachment_idx(mol, slot)` returns that neighbour's index.  
`_rgroup_atom_idx(mol, slot)` returns the dummy's own index.

### BILN separator gotcha
- `-` = monomer join within a chain → `A-G-A` is a tripeptide
- `.` = chain separator → `A.G.A` is three isolated single-residue chains (no backbone bonds)

`.Token(` = inline attachment (CABILN); `.Token` without `(` = chain separator (BILN).

### m_Rgroups format
Stored as a list of SMILES strings (or `None`), e.g. `['[H]', '[OH]', '[H]', None]`.  
**Not** a string of indices. The list index = slot index (0-based).  
`__restore_and_remove_rgroups` uses these to decide what to put back at unbound sites.

### Leaving group restoration (the Asp/His/Trp fix)
Old code used `DeleteSubstructs('[#0]')` — deleted ALL dummies including sidechain ones, turning Asp's COOH into an aldehyde.

New code (`__restore_and_remove_rgroups` in `molecule.py`):
1. Kekulize the combined mol **before** any atom removal (clears aromatic flags → explicit alternating bonds). Prevents `SanitizeMol` failing on Trp/His aromatic NH after dummy removal.
2. For each unbound slot: if leaving group is `[H]` → delete dummy; if `[OH]` or other → replace dummy with the leaving group's heavy atom.
3. Delete all remaining dummies (bonded slots).
4. `SanitizeMol` re-assigns aromaticity from the Kekulé form.

### Bond validation layers
`_check_bond_chemistry(mol1, at1, mol2, at2)` in `sequence.py`:
- **Silent**: amide N–C(=O), disulfide S–S, ester O–C(=O), diselenide Se–Se
- **UserWarning**: thioester S–C(=O), sulfenamide S–N, hydrazide N–N, non-carbonyl N–C, non-carbonyl O–C (ether), thioether S–C (aliphatic)
- **ValueError**: anything else (C–C, O–O, unknown element pairs) — raises before bond is added

`Molecule.__from_sequence` also wraps `SanitizeMol` — catches RDKit `AtomValenceException` and re-raises as `ValueError`.

Sulfenamide (S–N) and thioether (S–C aliphatic) are **intentionally warnings, not errors** — both arise in legitimate protecting-group and PTM chemistry.

### Monomer pipeline flow
```
CSV (token, input, chuckles, r1_leaving, r2_leaving, ...)
  ↓ build_library_from_csv()
  ├─ if chuckles column valid → use directly (caps: ac, am, fmoc, boc, tbu, otbu, trt, acm, pbf)
  └─ else → normalize_input(inp) → SMILES → pre_activate() → CHUCKLES + leaving groups
  ↓
SDF (one mol per token, properties: token, name, type, chuckles, r1_leaving, ...)
```

`normalize_input(raw)` auto-detects:
- CHUCKLES: regex `\[\d+\*\]`
- SMILES: RDKit parse attempt
- FASTA single-letter: `_AA_SMILES` dict lookup
- BILN token: same dict (keys: `'DAla'`, `'Hyp'`, `'Aib'`, etc. — **not** `'dA'`)

`build_library_from_csv` `rebuild` flag:
- `rebuild=False` (default): uses pre-filled CHUCKLES for all rows
- `rebuild=True`: re-derives `aa` entries via `pre_activate`; always preserves `cap` CHUCKLES; skips stale `r*_leaving` overrides for `aa` rows

---

## Protecting-group caps (2026-04-28)

Total library: **1001 monomers** (amino acids, caps, linkers, bioconjugation handles, etc.).

| Token | Description | Bond to residue | Bond type |
|-------|-------------|-----------------|-----------|
| `ac` | Acetyl N-terminal cap | R2=carbonyl C | amide (silent) |
| `am` | Amide C-terminal cap | R1=N | amide (silent) |
| `fmoc` | Fmoc N-terminal cap | R2=carbonyl C | amide (silent) |
| `boc` | Boc sidechain amine PG | R2=carbonyl C | N–C(=O) carbamate |
| `tbu` | tBu ester sidechain PG | R1=O | O–C ester (silent) |
| `otbu` | tBu ether sidechain PG | R1=C | O–C ether (UserWarning) |
| `trt` | Trityl Cys thiol PG | R2=C (LG=[Cl]) | S–C thioether (UserWarning) |
| `acm` | Acetamidomethyl Cys thiol PG | R2=C (LG=[Cl]), R4=amide N | S–C thioether (UserWarning) |
| `pbf` | Pbf Arg guanidinium PG | R2=sulfonyl S | S–N sulfenamide (UserWarning) |

Usage:
```
fmoc-C.trt(4,2)-am        # Cys thiol (R4) protected with Trt; trt R2→Cys R4
fmoc-C.acm(4,2)-am        # Cys thiol protected with Acm; acm R2→Cys R4
fmoc-R.pbf(4,2)-am        # Arg guanidinium (R4) protected with Pbf; pbf R2→Arg R4
```

---

## `find_sidechain_slots` three-pass order

**Pass 1 — backbone-N modification (always R3 for non-Pro):**
- If `backbone_n_idx` provided and that N has ≥ 2 explicit H → R3

**Pass 2 — sidechain reactive atoms (R4+), SMARTS priority order:**
- thiol > selenol > primary amine `[NX3;H2]` > guanidinium secondary N `[NX3;H1][CX3](=N)` > extra COOH > phenolic OH > aliphatic OH > aromatic NH > amide NH

**Pass 3 — second H on sidechain primary amines:**
- For every `[NX3;H2]` in `seen` that is NOT the backbone N → add a second slot

Pro is excluded naturally (ring N has only 1 H, never reaches the ≥ 2 threshold).

---

## Phase 2 CABILN parser (2026-04-28)

`_preprocess_cabiln`, `_handle_terminal_bond_markers`, and a rewritten `_expand_inline_caps` in `sequence.py`.

### New forms supported

| Form | Meaning |
|------|---------|
| `%` or `\n` | Segment separator — main chain first, branches follow |
| `.!1(4,2)` | First crosslink endpoint — explicit R-groups |
| `.!1` | Second endpoint — parens omitted, inverse inferred from `_seen` |
| `!1-A-B-C` | N-terminal branch marker — annotates first residue with `(!1,1)` |
| `A-B-C-!1` | C-terminal branch marker — annotates last residue with `(!1,2)` |

### Validation rules

- First `.!n(y,z)` sets `_seen[tok]`; second explicit `(y,z)` must be inverse or raises `ValueError`.
- No-parens `.!n` infers inverse from `_seen`; raises if no prior occurrence exists.
- Terminal markers cross-validate against `_seen[tok][1]` (partner rgroup): `!n-` implies R1, `-!n` implies R2.
- Exactly 2 endpoints per bond ID — 1 or 3+ raise `ValueError`.
- Preliminary scan over all segments before transformation, so branch-before-main-chain ordering works.

### Implementation note

`_INLINE_ATTACH_RE` is replaced by two separate regexes:
- `_INLINE_BOND_RE = re.compile(r'\.(!\w+)(?:\((\d+),(\d+)\))?')` — bond markers (optional parens)
- `_INLINE_CAP_RE  = re.compile(r'\.([A-Za-z]\w*)\((\d+),(\d+)\)')` — named caps (parens required)

---

## Pending / known issues

- **Mid-chain branch continuation `!x-` in main chain** — `!x-` as a leading marker on a sub-chain that is not a full pendant (e.g. `A-B-C-!1-D-E` meaning the branch grafts mid-chain). Not yet handled; would require representing the main chain as multiple BILN sub-chains. Current implementation handles start/end terminal markers on separate `%` segments only.
- **Legacy tests not audited** — `sequence_Test.py`, `molecule_Test.py`, `converter_Test.py` use old `.`-separated BILN and old SDF properties. They likely fail and need updating.
- **`describe_monomer` CLI utility** — show R-group slot map for any token. Not yet implemented.
- **Missing SMARTS for bioconjugation/click groups** — `_SC_PATTERNS` was written for the 20 standard AAs. The following functional groups are not auto-detected by `pre_activate` and require hand-authored CHUCKLES:
  - Azide (`[N-]=[N+]=N`)
  - Terminal alkyne (`C#CH`)
  - Maleimide (thiol-reactive alkene)
  - NHS ester (electrophilic, not a nucleophile — but marks the attachment C)
  - Tetrazine, DBCO, BCN, TCO (strained rings — not easily expressed as SMARTS)
  - Aldehyde/ketone (oxime ligation)
  - Boronic acid (covalent warheads, Pd coupling)
  - Halide leaving groups
  SMARTS for azide/alkyne/maleimide would be straightforward to add. Strained-ring bioorthogonal groups require hand-authoring regardless due to SMARTS complexity. Note: this gap only affects library *building* — joining itself uses pre-encoded `[n*]` atoms and is unaffected.

---

## Key file paths

```
src/pyPept/
  sequence.py          — Sequence class, _check_bond_chemistry, _expand_inline_caps, CABILN detection
  molecule.py          — Molecule class, __restore_and_remove_rgroups
  interfaces/
    monomer_pipeline.py — pre_activate, build_library_from_csv, find_sidechain_slots
    map_monomers.py     — deprecation shim only
  data/
    monomers.sdf        — 1001-monomer CABILN library (authoritative)
    monomers.csv        — CSV authoring source

tests/
  test_bond_validation_and_assembly.py  — 572 tests (571 pass, 1 xfail)
```
