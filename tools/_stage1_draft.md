Every bond in CABILN is validated before assembly. The system runs in two stages:

**Stage 1 — R-group labelling.** When a monomer is registered via `pre_activate()`, three detection passes assign a chemistry type to every attachment point:

1. **Backbone detection** (graph topology) — shortest path between amino N and carboxyl C → R1 (`backbone_n`), R2 (`backbone_c`), R3 (`backbone_n_mod` if backbone N has ≥ 2 H)
2. **Sidechain scan** (24 SMARTS patterns in priority order) — first match wins per atom; multi-atom patterns protect non-attachment atoms from re-matching
3. **Second-H pass** — sidechain primary amines (`amine_primary`) get an extra slot for the second N–H → `amine_secondary`

The complete type set, in priority order:

| # | Type | SMARTS | LG | label_only | Example |
|---|------|--------|----|------------|---------|
| — | `backbone_n` | *(graph topology)* | `[H]` | — | all amino acids R1 |
| — | `backbone_c` | *(graph topology)* | `[OH]` | — | all amino acids R2 |
| — | `backbone_n_mod` | *(≥ 2 H on backbone N)* | `[H]` | — | all except Pro R3 |
| 1 | `thiol` | `[SX2H1:1]` | `[H]` | no | Cys |
| 2 | `selenol` | `[SeX2H1:1]` | `[H]` | no | Sec |
| 3 | `alkyl_halide_c` | `[CX4;!H0:1][Cl,Br,I]` | `[H]` | no | ClAc, BrAc |
| 4 | `aminooxy` | `[NH2:1][OX2H0]` | `[H]` | no | Aoa |
| 5 | `hydrazide` | `[NX3H1:1][NX3H2]` | `[H]` | no | HzAla |
| 6 | `amine_primary` | `[NX3;H2:1]` | `[H]` | no | Lys, Orn |
| 7 | `guanidinium` | `[NX3;H1:1][CX3](=N)` | `[H]` | **yes** | Arg |
| 8 | `guanidinium_imine` | `[NX2H1:1]=[CX3]([NX3])[NX3]` | `[H]` | no | *(shadow)* |
| 9 | `carboxyl` | `[CX3:1](=O)[OX2H1]` | `[OH]` | no | Asp, Glu |
| 10 | `hydroxyl_phenolic` | `[OX2H1:1][c]` | `[H]` | **yes** | Tyr |
| 11 | `hydroxyl` | `[OX2H1:1][CX4]` | `[H]` | no | Ser, Thr |
| 12 | `aromatic_nh` | `[nH:1]` | `[H]` | **yes** | His, Trp |
| 13 | `amide_nh` | `[NX3;H1:1][CX3]=O` | `[H]` | **yes** | Lys_Boc, Cit |
| 14 | `phosphate_p` | `[P:1](=O)([OH])[OH]` | `[H]` | no | pSer, pTyr |
| 15 | `cyclooctyne_c` | `[CX4;!H0:1][C;r]#[C;r]` | `[H]` | no | CyoAla |
| 16 | `alkyne_c` | `[CX4;!H0:1]C#[CH]` | `[H]` | no | Pra |
| 17 | `azide_alpha_c` | `[CX4;!H0:1][N]=[N+]=[N-]` | `[H]` | no | AzAla |
| 18 | `terminal_alkene` | `[CH1:1]=[CH2]` | `[H]` | no | S5, R8 |
| 19 | `tetrazine_c` | `[CX4;!H0:1][c]1nnc(nn1)` | `[H]` | no | TzAla |
| 20 | `tco_c` | `[CX4;!H0;!r:1][C;r]=[C;r]` | `[H]` | no | TcoAla |
| 21 | `aldehyde` | `[CX3H1:1](=O)[!#7]` | `[H]` | no | Ald |
| 22 | `formamide_c` | `[CX3H1:1](=O)[NX3]` | `[H]` | no | *(shadow)* |
| 23 | `nhs_ester` | `[CX4;!H0:1]C(=O)ON1C(=O)CCC1=O` | `[H]` | no | NHSAla |
| 24 | `maleimide_c` | `[CH1:1]1=[C][C](=[O])[N][C]1=[O]` | `[H]` | no | MalAla |
| — | `amine_secondary` | *(Pass 3: second H)* | `[H]` | — | Lys R5 |

**`label_only`** types get a visible R-group slot but no entry in the bond reaction table. This prevents incorrect bonding of detectable-but-unreactive atoms (e.g. Arg's guanidinium NH, Tyr's phenolic OH, His's imidazole NH). The slot is still useful for cap attachment and for tools that inspect monomer connectivity.

**Shadowed types**: `formamide_c` (#22) and `guanidinium_imine` (#8) are shadowed at registration time. For formamide_c, amide_nh (#13) matches the formamide N first, and the `protected` set blocks the adjacent C. For guanidinium_imine, guanidinium (#7) matches the ε-NH first, and its `(=N)` branch puts the imine N in `protected`. Both types activate only at assembly time via `infer_chem_type`.

Each type is documented below with 5 walkthrough examples showing step-by-step detection logic.

---

#### Backbone detection — graph topology

R1/R2 are found by the shortest path between the amino N and carboxyl C. R3 is assigned if the backbone N has ≥ 2 hydrogens. No SMARTS pattern — pure graph search.

<details>
<summary>5 examples</summary>

```python
# 1. Alanine — standard α-amino acid
pre_activate('N[C@@H](C)C(=O)O')
# Backbone path: N → Cα → C(=O) → O. R1=backbone_n, R2=backbone_c.
# Backbone N has 2 H → R3=backbone_n_mod.
# No sidechain reactive group. Result: {1: 'backbone_n', 2: 'backbone_c', 3: 'backbone_n_mod'}

# 2. Proline — ring nitrogen, no R3
pre_activate('O=C(O)[C@@H]1CCCN1')
# Backbone path found through pyrrolidine ring.
# Ring N has only 1 H → R3 threshold (≥ 2 H) NOT met. No R3 assigned.
# Result: {1: 'backbone_n', 2: 'backbone_c'}

# 3. β-alanine (bAla) — β-amino acid
pre_activate('NCCC(=O)O')
# Backbone path: N → C → C → C(=O). Longer path than α-amino acid.
# Graph topology handles β-, γ-amino acids automatically.
# Result: {1: 'backbone_n', 2: 'backbone_c', 3: 'backbone_n_mod'}

# 4. Acetyl cap (ac) — no nitrogen
pre_activate('CC(=O)O')
# No amino N found → cap detection. COOH present → R2=backbone_c.
# No R1 or R3. Single-ended cap for N-terminal capping.
# Result: {2: 'backbone_c'}

# 5. Glycolic acid — depsipeptide building block
pre_activate('OCC(=O)O')
# No amino N → cap detection. Hydroxyl path to COOH.
# R2=backbone_c, R4=hydroxyl (the sidechain -OH).
# Result: {2: 'backbone_c', 4: 'hydroxyl'}
```

</details>

#### thiol — `[SX2H1:1]` (priority 1)

Sulfhydryl –SH. Highest sidechain priority. LG=`[H]` — the S–H hydrogen departs during disulfide or thioether bond formation.

<details>
<summary>5 examples</summary>

```python
# 1. Cysteine (C) — canonical thiol
pre_activate('N[C@@H](CS)C(=O)O')
# Sidechain scan: [SX2H1:1] matches the -SH sulfur (degree 2, 1 H). → R4=thiol
# CHUCKLES: [1*]N([3*])[C@@H](C[S][4*])C([2*])=O

# 2. Penicillamine (Pen) — β,β-dimethyl cysteine
pre_activate('NC(C)(CS)C(=O)O')
# Same [SX2H1:1] match on -SH. gem-dimethyl backbone doesn't affect sidechain.
# CHUCKLES: [1*]N([3*])C(C)(C[S][4*])C([2*])=O

# 3. Homocysteine (Hcy) — one extra CH₂ in sidechain
pre_activate('N[C@@H](CCS)C(=O)O')
# [SX2H1:1] matches regardless of chain length to sulfur.
# CHUCKLES: [1*]N([3*])[C@@H](CC[S][4*])C([2*])=O

# 4. D-Cysteine (DCys) — opposite chirality
pre_activate('N[C@H](CS)C(=O)O')
# SMARTS is chirality-agnostic. Same R4=thiol assignment.

# 5. Deamino-cysteine — cap with thiol sidechain
pre_activate('OC(=O)CS')
# No amino N → cap backbone. R2=backbone_c. Sidechain: R4=thiol.
# Result: {2: 'backbone_c', 4: 'thiol'}
```

</details>

#### selenol — `[SeX2H1:1]` (priority 2)

Selenohydryl –SeH. Diselenide chemistry analogue of thiol. LG=`[H]`.

<details>
<summary>5 examples</summary>

```python
# 1. Selenocysteine (Sec) — the 21st amino acid
pre_activate('N[C@@H](C[SeH])C(=O)O')
# [SeX2H1:1] matches -SeH (degree 2, 1 H). → R4=selenol

# 2. D-Selenocysteine (seC) — D-form
pre_activate('N[C@H](C[SeH])C(=O)O')
# Same match, opposite chirality.

# 3. Homo-selenocysteine (Se_Hcy) — longer sidechain
pre_activate('N[C@@H](CC[SeH])C(=O)O')
# Chain length doesn't affect [SeX2H1:1] match. → R4=selenol

# 4. D-Homo-selenocysteine (D_Se_Hcy)
pre_activate('N[C@H](CC[SeH])C(=O)O')
# D-form of Se_Hcy. Same detection.

# 5. Selenopenicillamine (Se_Pen) — gem-dimethyl
pre_activate('NC(C)(C[SeH])C(=O)O')
# [SeX2H1:1] on Se, unaffected by β,β-disubstitution. → R4=selenol
```

</details>

#### alkyl_halide_c — `[CX4;!H0:1][Cl,Br,I]` (priority 3)

sp3 carbon bearing a halide. The dummy replaces a C–H; the halide stays in the monomer and departs as part of the SMIRKS reaction (SN₂). LG=`[H]`.

<details>
<summary>5 examples</summary>

```python
# 1. Chloroacetyl (ClAc) — chloroacetic acid cap
pre_activate('OC(=O)CCl')
# [CX4;!H0:1] matches the CH₂ (sp3, has H); [Cl] matches chlorine. → R4=alkyl_halide_c
# No amino N → cap backbone. Result: {2: 'backbone_c', 4: 'alkyl_halide_c'}

# 2. Bromoacetyl (BrAc) — bromine variant
pre_activate('OC(=O)CBr')
# Same SMARTS, [Br] matches. → R4=alkyl_halide_c

# 3. β-Bromoalanine (Ala_3Br) — amino acid with halide sidechain
pre_activate('N[C@@H](CBr)C(=O)O')
# Full amino acid backbone + sidechain halide. R4=alkyl_halide_c

# 4. β-Chloroalanine (Ala_3Cl) — chlorine on Ala
pre_activate('N[C@@H](CCl)C(=O)O')
# Same pattern, Cl variant. R4=alkyl_halide_c

# 5. D-β-Bromoalanine (D_Ala_3Br) — D-form
pre_activate('N[C@H](CBr)C(=O)O')
# Chirality-agnostic. R4=alkyl_halide_c
```

</details>

#### aminooxy — `[NH2:1][OX2H0]` (priority 4)

Aminooxy –O–NH₂. Checked *before* amine_primary (#6) so `-ONH₂` is not misclassified as a generic amine. LG=`[H]`.

<details>
<summary>5 examples</summary>

```python
# 1. Aminooxy-serine (Aoa) — canonical aminooxy
pre_activate('N[C@@H](CON)C(=O)O')
# [NH2:1] matches the terminal -NH₂; [OX2H0] matches the O (2 bonds, no H).
# ⚠ PRIORITY: checked before amine_primary (#6), so -ONH₂ → aminooxy, not amine.
# CHUCKLES: [1*]N([3*])[C@@H](CON[4*])C([2*])=O → R4=aminooxy

# 2. Aminooxy-lysine (AoaLys) — Lys-length linker
pre_activate('N[C@@H](CCCCON)C(=O)O')
# Longer chain, same -ONH₂ terminus. → R4=aminooxy

# 3. Aminooxy-ornithine (AoaOrn) — Orn-length
pre_activate('N[C@@H](CCCON)C(=O)O')
# R4=aminooxy. Chain length irrelevant.

# 4. Aminooxy-Dab (AoaDab) — short linker
pre_activate('N[C@@H](CCON)C(=O)O')
# R4=aminooxy.

# 5. D-Aminooxy-serine (D_Aoa) — D-form
pre_activate('N[C@H](CON)C(=O)O')
# Same detection, opposite chirality. R4=aminooxy
```

</details>

#### hydrazide — `[NX3H1:1][NX3H2]` (priority 5)

Hydrazide –CO–NH–NH₂. The N–N bond distinguishes it from plain amine. The `:1` atom is the NH (not the terminal NH₂). Multi-atom match protects the terminal NH₂ from being re-classified as amine_primary. LG=`[H]`.

<details>
<summary>5 examples</summary>

```python
# 1. Hydrazide-alanine (HzAla) — canonical hydrazide
pre_activate('N[C@@H](CC(=O)NN)C(=O)O')
# [NX3H1:1] matches the -NH- (1 H); [NX3H2] matches the -NH₂.
# ⚠ PROTECTION: The terminal -NH₂ (match[1]) is added to `protected` set,
#   so amine_primary (#6) cannot re-match it. → R4=hydrazide only.
# CHUCKLES: [1*]N([3*])[C@@H](CC(=O)N([4*])N)C([2*])=O

# 2. Hydrazide-lysine (HzLys) — long linker
pre_activate('N[C@@H](CCCCC(=O)NN)C(=O)O')
# Same -CONHNH₂ terminus. → R4=hydrazide

# 3. Hydrazide-ornithine (HzOrn)
pre_activate('N[C@@H](CCCC(=O)NN)C(=O)O')
# R4=hydrazide.

# 4. Hydrazide-Dab (HzDab) — short linker
pre_activate('N[C@@H](CCC(=O)NN)C(=O)O')
# R4=hydrazide.

# 5. D-Hydrazide-alanine (D_HzAla) — D-form
pre_activate('N[C@H](CC(=O)NN)C(=O)O')
# R4=hydrazide. Chirality-agnostic.
```

</details>

#### amine_primary — `[NX3;H2:1]` (priority 6)

Generic sidechain –NH₂. Catches all primary amines not already claimed by aminooxy (#4) or hydrazide (#5). In Pass 3, a second slot is added for the second N–H → `amine_secondary`. LG=`[H]`.

<details>
<summary>5 examples</summary>

```python
# 1. Lysine (K) — ε-amine gets two slots
pre_activate('NCCCC[C@@H](N)C(=O)O')
# [NX3;H2:1] matches the ε-NH₂. → R4=amine_primary
# Pass 3: same atom gets R5=amine_secondary (second H).
# CHUCKLES: [1*]N([3*])[C@@H](CCCCN([4*])[5*])C([2*])=O

# 2. Ornithine (Orn) — one CH₂ shorter than Lys
pre_activate('N[C@@H](CCCN)C(=O)O')
# ε-NH₂ → R4=amine_primary, R5=amine_secondary.

# 3. 2,3-Diaminopropionic acid (Dap) — shortest amine sidechain
pre_activate('N[C@@H](CN)C(=O)O')
# β-NH₂ → R4=amine_primary, R5=amine_secondary.

# 4. D-Lysine (DLys) — D-form
pre_activate('NCCCC[C@H](N)C(=O)O')
# Same detection. R4=amine_primary, R5=amine_secondary.

# 5. Aminobutyric acid (Dab) — γ-amine
pre_activate('N[C@@H](CCN)C(=O)O')
# γ-NH₂ → R4=amine_primary, R5=amine_secondary.
```

</details>

#### guanidinium — `[NX3;H1:1][CX3](=N)` (priority 7, label_only)

Guanidinium NH on arginine. **label_only**: slot is visible but has no bond-table reaction — prevents the guanidinium N–H from participating in unintended bond formation. The multi-atom match protects the guanidinium C from later pattern re-matching.

<details>
<summary>5 examples</summary>

```python
# 1. Arginine (R) — multi-atom protection in action
pre_activate('N[C@@H](CCCNC(=N)N)C(=O)O')
# Sidechain scan (in order):
#   amine_primary (#6): terminal -NH₂ matches → R4=amine_primary
#   guanidinium (#7): the ε-NH matches [NX3;H1:1][CX3](=N) → R5=guanidinium (label_only)
#     ⚠ PROTECTION: match tuple includes [CX3] AND the (=N) branch nitrogen.
#     Both added to `protected`. The =NH nitrogen is now blocked.
#   guanidinium_imine (#8): =NH would match but is in `protected` → SKIPPED
# Pass 3: R4 amine_primary gets second H → R6=amine_secondary
# Result: {4: 'amine_primary', 5: 'guanidinium', 6: 'amine_secondary'}

# 2. D-Arginine (DArg)
pre_activate('N[C@H](CCCNC(=N)N)C(=O)O')
# Same protection mechanism. {4: amine_primary, 5: guanidinium, 6: amine_secondary}

# 3. Homoarginine (hArg) — one extra CH₂
pre_activate('N[C@@H](CCCCNC(=N)N)C(=O)O')
# Longer chain, same guanidinium detection and protection.

# 4. N-Methyl-arginine (meR) — slot numbering shift
pre_activate('CN[C@@H](CCCNC(=N)N)C(=O)O')
# Backbone N has 1 H → no R3. Sidechain starts at R3 (not R4).
# Result: {3: 'amine_primary', 4: 'guanidinium', 5: 'amine_secondary'}

# 5. Arginine aldehyde (Arg_al) — C-terminal aldehyde
pre_activate('N[C@@H](CCCNC(=N)N)C=O')
# Backbone C is aldehyde (not COOH). R2 still assigned.
# Sidechain guanidinium detection unaffected.
```

</details>

#### guanidinium_imine — `[NX2H1:1]=[CX3]([NX3])[NX3]` (priority 8, shadowed)

Guanidinium =NH imine nitrogen. In theory this matches the imine N in Arg-like guanidinium groups. In practice, **always shadowed**: guanidinium (#7) matches the ε-NH first, and its multi-atom match `[NX3;H1:1][CX3](=N)` includes the imine N via the `(=N)` branch — adding it to `protected`. Activates only at assembly time via `infer_chem_type`.

<details>
<summary>Why it's shadowed — step by step</summary>

```python
# Arginine (R): N[C@@H](CCCNC(=N)N)C(=O)O
# The guanidinium group -NHC(=NH)NH₂ has three matchable nitrogens:
#
# 1. Terminal -NH₂: amine_primary (#6) matches → R4
# 2. ε-NH: guanidinium (#7) matches [NX3;H1:1][CX3](=N)
#    → R5=guanidinium (label_only)
#    → match tuple = (ε-N, C, =NH-N) → C and =NH-N added to `protected`
# 3. =NH: guanidinium_imine (#8) would match [NX2H1:1]=[CX3]...
#    → BUT =NH-N is already in `protected` → SKIPPED
#
# Result: {4: 'amine_primary', 5: 'guanidinium', 6: 'amine_secondary'}
# guanidinium_imine is absent from pre_activate output.
# At assembly time, infer_chem_type() can detect it because it uses
# different SMARTS without the protection mechanism.

pre_activate('N[C@@H](CCCNC(=N)N)C(=O)O')
# → {4: 'amine_primary', 5: 'guanidinium', 6: 'amine_secondary'}
```

</details>

#### carboxyl — `[CX3:1](=O)[OX2H1]` (priority 9)

Sidechain –COOH. LG=`[OH]` (not `[H]`) — this distinguishes carboxyl from aldehyde at assembly time. Backbone COOH is excluded by graph topology claim before sidechain scan runs.

<details>
<summary>5 examples</summary>

```python
# 1. Aspartate (D) — two COOH groups, only sidechain detected
pre_activate('N[C@@H](CC(=O)O)C(=O)O')
# ⚠ EDGE CASE: Both backbone and sidechain have COOH.
# Backbone COOH claimed first by graph topology → excluded from scan.
# Sidechain COOH matches [CX3:1](=O)[OX2H1] → R4=carboxyl, LG=[OH]

# 2. Glutamate (E) — longer sidechain COOH
pre_activate('N[C@@H](CCC(=O)O)C(=O)O')
# Same logic. Backbone COOH excluded; sidechain → R4=carboxyl.

# 3. γ-Glutamate (gGlu) — γ-linked backbone
pre_activate('N[C@@H](CCC(=O)O)CC(=O)O')
# Graph topology picks the γ-carboxyl as backbone C.
# The α-carboxyl becomes sidechain → R4=carboxyl.

# 4. D-Aspartate (DAsp)
pre_activate('N[C@H](CC(=O)O)C(=O)O')
# Same as Asp, D-chirality. R4=carboxyl.

# 5. Glutaric acid derivative (Gla) — dicarboxylic sidechain
pre_activate('N[C@@H](CC(=O)CC(=O)O)C(=O)O')
# Backbone COOH excluded. Sidechain COOH matches → R4=carboxyl.
```

</details>

#### hydroxyl_phenolic — `[OX2H1:1][c]` (priority 10, label_only)

Phenolic –OH on aromatic ring. **label_only**: phenolic OH is unreactive under standard peptide conditions. Checked before aliphatic hydroxyl (#11) so Tyr doesn't get the wrong type.

<details>
<summary>5 examples</summary>

```python
# 1. Tyrosine (Y) — canonical phenol
pre_activate('N[C@@H](Cc1ccc(O)cc1)C(=O)O')
# [OX2H1:1][c] matches the phenolic -OH on the aromatic ring. → R4=hydroxyl_phenolic
# ⚠ PRIORITY: matched before hydroxyl (#11), so [OX2H1][c] wins over [OX2H1][CX4].

# 2. D-Tyrosine (DTyr)
pre_activate('N[C@H](Cc1ccc(O)cc1)C(=O)O')
# Same detection, D-form.

# 3. 3-Iodotyrosine (Tyr_3I) — halogenated ring
pre_activate('N[C@@H](Cc1cc(I)c(O)cc1)C(=O)O')
# Phenolic -OH still on aromatic ring. Iodo substituent doesn't block.

# 4. 3-Nitrotyrosine (Tyr_3NO2) — nitrated
pre_activate('N[C@@H](Cc1cc([N+](=O)[O-])c(O)cc1)C(=O)O')
# Still phenolic OH on aromatic carbon. → R4=hydroxyl_phenolic

# 5. 2,6-Dimethyltyrosine (Tyr_26diMe) — sterically hindered
pre_activate('N[C@@H](Cc1c(C)cc(O)cc1C)C(=O)O')
# Methyl groups don't change the O-H on aromatic ring match.
```

</details>

#### hydroxyl — `[OX2H1:1][CX4]` (priority 11)

Aliphatic –OH on sp3 carbon. LG=`[H]`. Matches Ser, Thr, hydroxyproline. Does not match phenolic OH (caught by hydroxyl_phenolic at #10).

<details>
<summary>5 examples</summary>

```python
# 1. Serine (S) — canonical aliphatic hydroxyl
pre_activate('N[C@@H](CO)C(=O)O')
# [OX2H1:1][CX4] matches the -OH on sp3 CH₂. → R4=hydroxyl

# 2. Threonine (T) — β-methyl serine
pre_activate('N[C@@H]([C@@H](C)O)C(=O)O')
# Same match on -OH. Methyl branch irrelevant.

# 3. Homoserine (Hser) — one extra CH₂
pre_activate('N[C@@H](CCO)C(=O)O')
# Chain length doesn't affect [OX2H1][CX4]. → R4=hydroxyl

# 4. Hydroxyproline (Hyp) — ring hydroxyl
pre_activate('O=C(O)[C@@H]1C[C@@H](O)CN1')
# Ring -OH on sp3 C. [CX4] is the ring carbon. → R4=hydroxyl

# 5. D-Serine (DSer) — D-form
pre_activate('N[C@H](CO)C(=O)O')
# Same hydroxyl detection, opposite chirality. R4=hydroxyl.
# Note: tert-butyl ether (COC(C)(C)C) does NOT match — ether O has no H.
```

</details>

#### aromatic_nh — `[nH:1]` (priority 12, label_only)

Aromatic ring NH (imidazole, indole). **label_only**: aromatic NH is unreactive under standard conditions but the slot is useful for monitoring and cap attachment.

<details>
<summary>5 examples</summary>

```python
# 1. Histidine (H) — imidazole NH
pre_activate('N[C@@H](Cc1c[nH]cn1)C(=O)O')
# [nH:1] matches the aromatic ring nitrogen bearing H. → R4=aromatic_nh

# 2. Tryptophan (W) — indole NH
pre_activate('N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O')
# Same [nH:1] match on the indole ring nitrogen. → R4=aromatic_nh

# 3. D-Histidine (DHis)
pre_activate('N[C@H](Cc1c[nH]cn1)C(=O)O')
# D-form. Same detection.

# 4. D-Tryptophan (DTrp)
pre_activate('N[C@H](Cc1c[nH]c2ccccc12)C(=O)O')
# D-form. R4=aromatic_nh.

# 5. 5-Hydroxytryptophan (Trp_5OH) — two sidechain slots
pre_activate('N[C@@H](Cc1c[nH]c2ccc(O)cc12)C(=O)O')
# R4=aromatic_nh (indole NH), R5=hydroxyl_phenolic (5-OH on aromatic ring).
# Both detected independently — different atoms, no priority conflict.
```

</details>

#### amide_nh — `[NX3;H1:1][CX3]=O` (priority 13, label_only)

Amide N–H (sidechains of Asn, Gln; Boc protecting group; formamide). **label_only**: prevents amide NH from participating in unintended bond formation. Multi-atom match: `[CX3]=O` is added to `protected`, blocking the adjacent carbonyl C from matching later patterns (notably formamide_c #22).

<details>
<summary>5 examples</summary>

```python
# 1. Boc-lysine (Lys_Boc) — carbamate NH
pre_activate('N[C@@H](CCCCNC(=O)OC(C)(C)C)C(=O)O')
# The carbamate -NHC(=O)O- has NH (1H). [NX3;H1:1][CX3]=O matches. → R4=amide_nh
# ⚠ WHY NOT amine_primary? amine_primary (#6) requires H2, but this N has only 1H.

# 2. Citrulline (Cit) — ureido sidechain, two N–H types
pre_activate('N[C@@H](CCCNC(=O)N)C(=O)O')
# -NHC(=O)NH₂ has BOTH N-H types:
#   amine_primary (#6) matches terminal -NH₂ (H2) → R4=amine_primary
#   amide_nh (#13) matches the ε-NH (H1) bonded to C=O → R5=amide_nh (label_only)
# Pass 3: R4 gets second H → R6=amine_secondary.
# Result: {4: 'amine_primary', 5: 'amide_nh', 6: 'amine_secondary'}

# 3. Formyl-lysine (Lys_For) — formamide sidechain
pre_activate('N[C@@H](CCCCNC=O)C(=O)O')
# ⚠ PRIORITY SHADOW: amide_nh matches the N-H → R4=amide_nh (label_only).
# The formamide C is added to `protected` by the multi-atom match [CX3]=O.
# formamide_c (#22) CANNOT match the C — it is blocked.

# 4. Acetyl-lysine (Lys_Ac) — acetamide sidechain
pre_activate('N[C@@H](CCCCNC(=O)C)C(=O)O')
# -NHC(=O)CH₃: the NH matches amide_nh. → R4=amide_nh (label_only)

# 5. Asparagine (N) — sidechain amide, NOT amide_nh
pre_activate('N[C@@H](CC(=O)N)C(=O)O')
# ⚠ SURPRISE: Asn's -C(=O)NH₂ has NH₂ (2H), not NH (1H).
# amine_primary (#6) matches first: [NX3;H2:1] → R4=amine_primary
# amide_nh (#13) requires H1 — does NOT match NH₂.
# Pass 3: R5=amine_secondary. Result: {4: 'amine_primary', 5: 'amine_secondary'}
```

</details>

#### phosphate_p — `[P:1](=O)([OH])[OH]` (priority 14)

Phosphate group. The dummy is placed on phosphorus. LG=`[H]`. Matched on raw SMILES; after CHUCKLES conversion the –OH groups may not be present, so `infer_chem_type` falls back to element-based heuristic (`element_15`) at assembly time.

<details>
<summary>5 examples</summary>

```python
# 1. Phosphoserine (pSer)
pre_activate('N[C@@H](COP(=O)(O)O)C(=O)O')
# [P:1](=O)([OH])[OH] matches the phosphate P. → R4=phosphate_p

# 2. Phosphothreonine (pThr)
pre_activate('N[C@@H]([C@@H](C)OP(=O)(O)O)C(=O)O')
# Same phosphate match. → R4=phosphate_p

# 3. Phosphotyrosine (pTyr)
pre_activate('N[C@@H](Cc1ccc(OP(=O)(O)O)cc1)C(=O)O')
# Phosphate on aromatic ring via ether. → R4=phosphate_p

# 4. Ser_PO3H2 — alternative naming
pre_activate('N[C@@H](COP(=O)(O)O)C(=O)O')
# Same as pSer. Naming convention varies.

# 5. Thr_PO3H2
pre_activate('N[C@@H]([C@@H](C)OP(=O)(O)O)C(=O)O')
# Same as pThr.
```

</details>

#### cyclooctyne_c — `[CX4;!H0:1][C;r]#[C;r]` (priority 15)

sp3 C–H adjacent to a ring-strained alkyne (cyclooctyne). Two ring carbons flank the triple bond → two slots per cyclooctyne. LG=`[H]`. Used in strain-promoted azide–alkyne cycloaddition (SPAAC).

<details>
<summary>5 examples</summary>

```python
# 1. Cyclooctyne-alanine (CyoAla) — two slots from one ring
pre_activate('N[C@@H](CC1CCCCC#CC1)C(=O)O')
# The 8-membered ring has C#C. Two ring carbons adjacent to the triple bond
# each match [CX4;!H0:1][C;r]#[C;r]. → R4=cyclooctyne_c, R5=cyclooctyne_c

# 2. Cyclooctyne-lysine (CyoLys) — long linker
pre_activate('N[C@@H](CCCCC1CCCCC#CC1)C(=O)O')
# Same ring, longer sidechain. Two cyclooctyne_c slots.

# 3. Cyclooctyne-ornithine (CyoOrn) — Orn-length
pre_activate('N[C@@H](CCCC1CCCCC#CC1)C(=O)O')
# R4=cyclooctyne_c, R5=cyclooctyne_c.

# 4. Cyclooctyne-Dab (CyoDab) — short linker
pre_activate('N[C@@H](CCC1CCCCC#CC1)C(=O)O')
# Two slots regardless of linker length.

# 5. Cyclooctyne cap (Cyo) — no amino group
pre_activate('OC(=O)CC1CCCCC#CC1')
# Cap backbone: R2=backbone_c. Ring → R4=cyclooctyne_c, R5=cyclooctyne_c.
```

</details>

#### alkyne_c — `[CX4;!H0:1]C#[CH]` (priority 16)

sp3 C–H adjacent to a terminal alkyne –C≡CH. LG=`[H]`. The dummy goes on the sp3 neighbor, not the alkyne C. Used in CuAAC click chemistry.

<details>
<summary>5 examples</summary>

```python
# 1. Propargylglycine (Pra) — canonical alkyne handle
pre_activate('N[C@@H](CC#C)C(=O)O')
# [CX4;!H0:1] matches the sp3 -CH₂-; C#[CH] matches the terminal alkyne.
# → R4=alkyne_c. Dummy on sp3 C, not on the alkyne carbon.

# 2. Homopropargylglycine (Hpg) — one extra CH₂
pre_activate('N[C@@H](CCC#C)C(=O)O')
# The sp3 C adjacent to C≡CH is the last CH₂ before the triple bond.
# → R4=alkyne_c

# 3. Phe-4-O-propargyl (Phe_4OPrg) — propargyl ether on Phe
pre_activate('N[C@@H](Cc1ccc(OCC#C)cc1)C(=O)O')
# The sp3 -OCH₂- matches [CX4;!H0:1]. → R4=alkyne_c

# 4. D-Homopropargylglycine (D_HPG) — D-form
pre_activate('N[C@H](CCC#C)C(=O)O')
# R4=alkyne_c.

# 5. Cysteine-propargyl (Cys_Prg) — thioether, not thiol
pre_activate('N[C@@H](CSCC#C)C(=O)O')
# ⚠ The S is a thioether (S bonded to two C's), NOT a thiol (no S-H).
# [SX2H1] does NOT match (H0). Only alkyne_c matches. → R4=alkyne_c
```

</details>

#### azide_alpha_c — `[CX4;!H0:1][N]=[N+]=[N-]` (priority 17)

sp3 C–H bearing an azide –N₃. LG=`[H]`. Dummy on the sp3 carbon, not on the azide nitrogen. Used in CuAAC and SPAAC.

<details>
<summary>5 examples</summary>

```python
# 1. Azidoalanine (AzAla) — canonical azide handle
pre_activate('N[C@@H](CN=[N+]=[N-])C(=O)O')
# [CX4;!H0:1] on the sp3 CH₂; azide [N]=[N+]=[N-] matched. → R4=azide_alpha_c

# 2. Azidolysine (AzK) — Lys scaffold
pre_activate('N[C@@H](CCCCN=[N+]=[N-])C(=O)O')
# Terminal azide on ε-carbon. → R4=azide_alpha_c

# 3. Azido-ornithine (AzOrn) — Orn scaffold
pre_activate('N[C@@H](CCCN=[N+]=[N-])C(=O)O')
# → R4=azide_alpha_c

# 4. D-Azidoalanine (D_AzAla) — D-form
pre_activate('N[C@H](CN=[N+]=[N-])C(=O)O')
# R4=azide_alpha_c.

# 5. Azidohomoalanine (AzHal) — methionine surrogate
pre_activate('N[C@@H](CCN=[N+]=[N-])C(=O)O')
# Extra CH₂. → R4=azide_alpha_c
```

</details>

#### terminal_alkene — `[CH1:1]=[CH2]` (priority 18)

Vinyl –CH=CH₂ for ring-closing metathesis (RCM). LG=`[H]`. The atom map `:1` is on the =CH–, not on =CH₂.

<details>
<summary>5 examples</summary>

```python
# 1. (S)-pentenylglycine (S5) — i,i+4 staple
pre_activate('N[C@@H](CCC=C)C(=O)O')
# [CH1:1]=[CH2] matches the terminal vinyl. → R4=terminal_alkene

# 2. (R)-octenylglycine (R8) — i,i+7 staple
pre_activate('N[C@@H](CCCCCC=C)C(=O)O')
# Same vinyl pattern, longer chain. → R4=terminal_alkene

# 3. Allylglycine (AllGly) — shortest vinyl sidechain
pre_activate('N[C@@H](CC=C)C(=O)O')
# CH₂-CH=CH₂. → R4=terminal_alkene

# 4. D-Allylglycine (D_AllGly) — D-form
pre_activate('N[C@H](CC=C)C(=O)O')
# R4=terminal_alkene.

# 5. Gly_allyl — glycine-backbone allyl cap
pre_activate('N[C@@H](C=C)C(=O)O')
# Vinyl directly on Cα. [CH1:1]=[CH2] still matches. → R4=terminal_alkene
```

</details>

#### tetrazine_c — `[CX4;!H0:1][c]1[n][n][c][n][n]1` (priority 19)

sp3 C–H adjacent to an s-tetrazine ring. LG=`[H]`. Used in inverse electron-demand Diels–Alder (IEDDA) with TCO.

<details>
<summary>5 examples</summary>

```python
# 1. Tetrazine-alanine (TzAla) — short linker
pre_activate('N[C@H](Cc1nncnn1)C(=O)O')
# [CX4;!H0:1] matches the CH₂ adjacent to the tetrazine ring.
# [c]1nnc(nn1) matches the aromatic tetrazine. → R4=tetrazine_c

# 2. Tetrazine-lysine (TzLys) — Lys-length linker
pre_activate('N[C@@H](CCCCCc1nncnn1)C(=O)O')
# Last sp3 CH₂ before the ring. → R4=tetrazine_c

# 3. Tetrazine cap (Tz) — no amino group
pre_activate('OC(=O)Cc1nncnn1')
# Cap backbone: R2=backbone_c. → R4=tetrazine_c

# 4. Tetrazine-ornithine (TzOrn) — Orn-length
pre_activate('N[C@@H](CCCCc1nncnn1)C(=O)O')
# → R4=tetrazine_c

# 5. Tetrazine-Dab (TzDab) — short linker
pre_activate('N[C@@H](CCCc1nncnn1)C(=O)O')
# → R4=tetrazine_c
```

</details>

#### tco_c — `[CX4;!H0;!r:1][C;r]=[C;r]` (priority 20)

Exocyclic sp3 C–H adjacent to a trans-cyclooctene ring C=C. LG=`[H]`. The `!r` constraint means the matched C must NOT be in the ring. Used in IEDDA with tetrazine.

<details>
<summary>5 examples</summary>

```python
# 1. TCO-alanine (TcoAla) — short linker
pre_activate('N[C@H](CC1=CCCCCCC1)C(=O)O')
# [CX4;!H0;!r:1] matches the exocyclic sp3 CH adjacent to ring C=C.
# [C;r]=[C;r] matches the ring double bond. → R4=tco_c

# 2. TCO-lysine (TcoLys) — Lys-length
pre_activate('N[C@@H](CCCCC1=CCCCCCC1)C(=O)O')
# Last exocyclic sp3 CH₂ before ring. → R4=tco_c

# 3. TCO cap (TCO) — no amino group
pre_activate('OC(=O)CC1=CCCCCCC1')
# Cap backbone. → R4=tco_c

# 4. TCO-ornithine (TcoOrn) — Orn-length
pre_activate('N[C@@H](CCCC1=CCCCCCC1)C(=O)O')
# → R4=tco_c

# 5. TCO-Dab (TcoDab) — short linker
pre_activate('N[C@@H](CCC1=CCCCCCC1)C(=O)O')
# → R4=tco_c
```

</details>

#### aldehyde — `[CX3H1:1](=O)[!#7]` (priority 21)

Aldehyde –CHO. LG=`[H]`. The `[!#7]` exclusion prevents matching formamide C–H (which is bonded to nitrogen). At assembly time, carboxyl and aldehyde are distinguished by their leaving groups: `[OH]` → carboxyl, `[H]` → aldehyde.

<details>
<summary>5 examples</summary>

```python
# 1. 4-Formyl-phenylalanine (Ald) — aromatic aldehyde handle
pre_activate('N[C@@H](Cc1ccc(C=O)cc1)C(=O)O')
# [CX3H1:1](=O) matches the aldehyde C-H; [!#7] confirms neighbor is not N.
# → R4=aldehyde

# 2. Alanine-aldehyde (Ala_al) — C-terminal aldehyde mimic
pre_activate('N[C@@H](C)C=O')
# Backbone C is now an aldehyde (not COOH). The graph-topology still finds it.
# The aldehyde C matches [CX3H1](=O)[!#7]. → R2=backbone_c (aldehyde-type)

# 3. Aspartate-aldehyde (Asp_al) — sidechain aldehyde
pre_activate('N[C@@H](CC=O)C(=O)O')
# Sidechain -CH=O matches. → R4=aldehyde

# 4. Glycine-aldehyde (Gly_al) — simplest aldehyde
pre_activate('NCC=O')
# Backbone C as aldehyde.

# 5. Arginine-aldehyde (Arg_al) — aldehyde + guanidinium
pre_activate('N[C@@H](CCCNC(=N)N)C=O')
# ⚠ TWO SYSTEMS: sidechain has guanidinium types; backbone C is aldehyde.
# Both detected independently — sidechain scan is separate from backbone.
```

</details>

#### formamide_c — `[CX3H1:1](=O)[NX3]` (priority 22, shadowed)

Formamide C–H. In theory matches formamide –N–CHO. In practice, **always shadowed** at registration time: amide_nh (#13) matches the adjacent N first, and the multi-atom protection adds the carbonyl C to `protected` — blocking formamide_c. Activates only at assembly time via `infer_chem_type`.

<details>
<summary>Why it's shadowed — step by step</summary>

```python
# Formyl-lysine (Lys_For): N[C@@H](CCCCNC=O)C(=O)O
# The formamide group -NH-CHO has two matchable atoms:
#
# Atom 1 (N): [NX3;H1:1][CX3]=O → amide_nh (#13) matches
#   → N gets R4=amide_nh (label_only)
#   → CX3 (the formamide C) added to `protected` set
#
# Atom 2 (C): [CX3H1:1](=O)[NX3] → formamide_c (#22) would match
#   → BUT the C is already in `protected` → SKIPPED
#
# Result: only R4=amide_nh. The formamide_c type is invisible at registration.
# At assembly time, infer_chem_type() can detect it because it uses different
# SMARTS without the protection mechanism.

pre_activate('N[C@@H](CCCCNC=O)C(=O)O')
# → {4: 'amide_nh'}  (formamide_c absent)
```

</details>

#### nhs_ester — `[CX4;!H0:1]C(=O)ON1C(=O)CCC1=O` (priority 23)

NHS-activated ester. The entire succinimide ring is encoded in the SMARTS. Dummy on the sp3 C adjacent to the ester carbonyl. LG=`[H]`.

<details>
<summary>5 examples</summary>

```python
# 1. NHS-alanine (NHSAla) — canonical NHS handle
pre_activate('N[C@H](CC(=O)ON1C(=O)CCC1=O)C(=O)O')
# The entire SMARTS matches as one unit. [CX4;!H0:1] on the sp3 CH₂.
# → R4=nhs_ester

# 2. NHS-lysine (NHSLys) — Lys-length linker
pre_activate('N[C@@H](CCCCC(=O)ON1C(=O)CCC1=O)C(=O)O')
# Longer chain, same NHS terminus. → R4=nhs_ester

# 3. NHS-ornithine (NHSOrn) — Orn-length
pre_activate('N[C@@H](CCCC(=O)ON1C(=O)CCC1=O)C(=O)O')
# → R4=nhs_ester

# 4. NHS-Dab (NHSDab) — short linker
pre_activate('N[C@@H](CCC(=O)ON1C(=O)CCC1=O)C(=O)O')
# → R4=nhs_ester

# 5. D-NHS-alanine (D_NHSAla) — D-form
pre_activate('N[C@H](CC(=O)ON1C(=O)CCC1=O)C(=O)O')
# → R4=nhs_ester
```

</details>

#### maleimide_c — `[CH1:1]1=[C][C](=[O])[N][C]1=[O]` (priority 24)

Maleimide ring C=C hydrogen. LG=`[H]`. Used in thiol–maleimide conjugation.

<details>
<summary>5 examples</summary>

```python
# 1. Maleimide-alanine (MalAla) — short linker
pre_activate('N[C@@H](CN1C(=O)C=CC1=O)C(=O)O')
# The maleimide ring pattern matches the =CH-. → R4=maleimide_c

# 2. Maleimide-lysine (MalLys) — Lys-length
pre_activate('N[C@@H](CCCCN1C(=O)C=CC1=O)C(=O)O')
# Same ring pattern. → R4=maleimide_c

# 3. Maleimide cap (Mal) — no amino group
pre_activate('OC(=O)CCN1C(=O)C=CC1=O')
# Cap backbone. → R4=maleimide_c

# 4. Maleimide-ornithine (MalOrn) — Orn-length
pre_activate('N[C@@H](CCCN1C(=O)C=CC1=O)C(=O)O')
# → R4=maleimide_c

# 5. Maleimide-Dab (MalDab) — short linker
pre_activate('N[C@@H](CCN1C(=O)C=CC1=O)C(=O)O')
# → R4=maleimide_c
```

</details>

#### amine_secondary — Pass 3 (after registry scan)

Not in the SMARTS registry. For every sidechain atom already assigned `amine_primary`, if the nitrogen had 2 hydrogens, a second slot is created for the second N–H. LG=`[H]`.

<details>
<summary>5 examples</summary>

```python
# 1. Lysine (K) — R4 + R5 on same nitrogen
pre_activate('NCCCC[C@@H](N)C(=O)O')
# Pass 2: ε-NH₂ → R4=amine_primary
# Pass 3: same ε-N had 2H, already amine_primary → R5=amine_secondary

# 2. Ornithine (Orn) — δ-amine
pre_activate('N[C@@H](CCCN)C(=O)O')
# R4=amine_primary, R5=amine_secondary

# 3. Diaminopropionic acid (Dap) — β-amine
pre_activate('N[C@@H](CN)C(=O)O')
# R4=amine_primary, R5=amine_secondary

# 4. Diaminobutyric acid (Dab) — γ-amine
pre_activate('N[C@@H](CCN)C(=O)O')
# R4=amine_primary, R5=amine_secondary

# 5. Arginine (R) — second H on the guanidinium -NH₂
pre_activate('N[C@@H](CCCNC(=N)N)C(=O)O')
# The terminal -NH₂ got R4=amine_primary in Pass 2.
# Pass 3: same atom gets R6=amine_secondary (R5 taken by guanidinium).
```

</details>
