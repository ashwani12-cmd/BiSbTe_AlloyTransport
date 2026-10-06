# AMM and DMM — what the two models are, and what they say about this interface

Written 2026-10-04. Numbers from `s19_models.json` (`s19_mismatch_models.py`), all on the
corrected R-3m cell at the quasi-harmonically expanded geometries.

---

## 1. The one equation both models share

A phonon hits the interface. What fraction of its energy gets across? Everything else is
bookkeeping:

```
G_K  =  ∫ ħω · Φ₁(ω) · α₁₂(ω) · (dn/dT) · dω
```

- **Φ₁(ω)** — the one-sided phonon flux from side 1: how many phonons per second per area
  arrive at the interface at frequency ω. Computed here from the full NEP dispersion
  (density of states × group velocity along z, summed over modes with v_z > 0).
- **dn/dT** — how much the Bose–Einstein occupation changes per kelvin; this is what makes
  G_K temperature dependent.
- **α₁₂(ω)** — the transmission coefficient. **This is the only thing AMM and DMM
  disagree about.** They are two opposite guesses for it.

Set α = 1 and you get the **phonon radiation limit**: every incident phonon crosses. That
is a hard upper bound on any harmonic model, and it is the right thing to measure both
models against.

---

## 2. AMM — the Acoustic Mismatch Model (Little, 1959)

**The picture:** the interface is a perfectly flat, perfectly bonded plane between two
elastic continua. Phonons are plane sound waves. They cross exactly the way light crosses
a glass boundary — they **refract**, and the transmitted fraction is set by the **acoustic
impedance** `Z = ρv`, exactly like impedance matching in a transmission line.

**The physics it assumes:** scattering is **specular**. The phonon keeps its polarisation
and its in-plane momentum; nothing is randomised.

At normal incidence

```
α = 4 Z₁Z₂ / (Z₁ + Z₂)²
```

and at angle θ₁, with Snell's law `sinθ₂ = (v₂/v₁) sinθ₁`,

```
α(θ₁) = 4 Z₁Z₂ cosθ₁ cosθ₂ / (Z₁cosθ₁ + Z₂cosθ₂)²
```

If `v₂ > v₁` there is a **critical angle** `θ_c = arcsin(v₁/v₂)` past which nothing is
transmitted at all — total internal reflection, just as in optics. The angle integral over
the transmission cone is what `amm_transmission()` does.

**When it is right:** low temperature, where only long-wavelength phonons are excited and
the interface genuinely looks flat to them.
**Its known weakness:** it uses one sound speed per polarisation — an elastic-continuum
(Debye) quantity — and applies it to the whole spectrum, including short-wavelength
phonons that see a rough, atomically structured interface. That is AMM's standard
criticism, and it applies here too.

---

## 3. DMM — the Diffuse Mismatch Model (Swartz & Pohl, 1989)

**The picture:** the exact opposite extreme. **Every** phonon that reaches the interface is
scattered so completely that it forgets which side it came from. It is then re-emitted
into side 1 or side 2 with a probability set purely by how many states each side has
available.

Detailed balance then fixes the transmission with no free parameters:

```
α₁₂(ω) = Φ₂(ω) / (Φ₁(ω) + Φ₂(ω))
```

No angle dependence. No polarisation memory. Only the two flux spectra matter.

**When it is right:** high temperature and/or a genuinely rough, disordered interface,
where diffuse scattering dominates.

**Its known pathology — and it bites hard here.** For two *identical* materials
Φ₁ = Φ₂ and DMM gives **α = 0.5**: a 50 % reflecting interface inside a perfect crystal,
where there should be no interface at all. So **DMM's transmission can never exceed ~0.5**,
no matter how well matched the materials are.

---

## 4. The two models are the two extremes, so they bracket reality

| | AMM | DMM |
|---|---|---|
| scattering | fully **specular** | fully **diffuse** |
| what sets α | acoustic impedance `Z = ρv` | phonon flux / density of states |
| angle dependence | yes, with a critical angle | none |
| polarisation | conserved | randomised |
| best at | low T, smooth interface | high T, rough interface |
| identical materials | α → 1 (correct) | α → 0.5 (**wrong**) |

A real interface is neither, so the true answer should sit between them. That is exactly
what we find.

---

## 5. TABLE A — DMM on the corrected cell

| T (K) | Φ-weighted α_DMM | G_DMM (MW m⁻² K⁻¹) | as a fraction of the radiation limit |
|---|---|---|---|
| 200 | 0.321 | **33.96** | 0.321 |
| 300 | 0.320 | **33.31** | 0.320 |
| 400 | 0.319 | **32.47** | 0.319 |
| 500 | 0.317 | **31.56** | 0.317 |

α_DMM ≈ 0.32 rather than 0.5 because the two flux spectra are not identical — Sb₂Te₃ runs
to 4.674 THz against Bi₂Te₃'s 4.081 THz, so the lighter, stiffer material has more
high-frequency states and detailed balance sends the flux preferentially the other way.

**G_DMM falls 7.1 % from 200 to 500 K.** Nothing in DMM makes it rise: the only
temperature dependence is dn/dT and the thermal expansion of the lattice, and the widening
van der Waals gaps soften the spectrum faster than the occupation grows.

---

## 6. TABLE B — AMM on the corrected cell

Per-polarisation ingredients (sound speeds along z, the transport direction, from the
slope of the NEP acoustic branches at Γ):

| T (K) | pol | v Bi₂Te₃ (m/s) | v Sb₂Te₃ (m/s) | Z Bi₂Te₃ (MRayl) | Z Sb₂Te₃ (MRayl) | Z₂/Z₁ | Γ (angle-averaged) | θ_c |
|---|---|---|---|---|---|---|---|---|
| 200 | TA (×2) | 2070.3 | 2136.1 | 16.489 | 13.941 | 0.845 | 0.8228 | 75.7° |
| 200 | LA | 2162.3 | 1884.5 | 17.221 | 12.299 | 0.714 | 0.9748 | none |
| 300 | TA (×2) | 2047.7 | 2113.1 | 16.284 | 13.776 | 0.846 | 0.8223 | 75.7° |
| 300 | LA | 2128.4 | 1874.4 | 16.926 | 12.220 | 0.722 | 0.9759 | none |
| 400 | TA (×2) | 2025.7 | 2090.2 | 16.086 | 13.613 | 0.846 | 0.8225 | 75.7° |
| 400 | LA | 2097.4 | 1864.4 | 16.655 | 12.142 | 0.729 | 0.9769 | none |
| 500 | TA (×2) | 2004.1 | 2067.5 | 15.891 | 13.450 | 0.846 | 0.8232 | 75.8° |
| 500 | LA | 2068.7 | 1854.6 | 16.404 | 12.065 | 0.736 | 0.9779 | none |

Densities: ρ(Bi₂Te₃) 7964.6 → 7929.5, ρ(Sb₂Te₃) 6526.3 → 6505.5 kg/m³ over 200→500 K.

Result:

| T (K) | α_AMM = (Γ_LA + 2Γ_TA)/3 | G_AMM (MW m⁻² K⁻¹) | fraction of the radiation limit |
|---|---|---|---|
| 200 | 0.873 | **92.45** | 0.873 |
| 300 | 0.874 | **90.95** | 0.874 |
| 400 | 0.874 | **89.03** | 0.874 |
| 500 | 0.875 | **87.09** | 0.875 |

**LA crosses almost freely (Γ = 0.975, no critical angle) because Sb₂Te₃'s LA is the
*slower* of the two — a wave entering a slower medium has no total internal reflection.**
TA is the lossy channel: Sb₂Te₃'s TA is faster, so there is a 75.7° cone outside which
transverse phonons are entirely reflected.

---

## 7. What the comparison actually says

| T (K) | DMM | AMM | radiation limit (Sb) | NEMD (seed 1) | paper |
|---|---|---|---|---|---|
| 200 | 33.96 | 92.45 | 95.75 | **58.09** | 34.75 |
| 300 | 33.31 | 90.95 | 95.13 | **57.75** | 56.94 |
| 400 | 32.47 | 89.03 | 93.50 | **69.63** | 60.76 |
| 500 | 31.56 | 87.09 | 91.56 | **75.57** | 98.60 |

| T (K) | NEMD / DMM | NEMD / AMM | NEMD / radiation limit |
|---|---|---|---|
| 200 | 1.71 | 0.63 | 0.549 |
| 300 | 1.73 | 0.63 | 0.555 |
| 400 | 2.14 | 0.78 | 0.684 |
| 500 | 2.39 | 0.87 | 0.759 |

1. **NEMD sits between DMM and AMM at every temperature** — 1.7–2.4× DMM and 0.63–0.87×
   AMM. A partly specular, partly diffuse interface, which is what a nearly lattice-matched
   commensurate junction should be.
2. **The interface becomes more specular as it heats**, in this ratio sense: NEMD climbs
   from 0.63 to 0.87 of AMM while all three models drift down. The rise is a real transport
   effect, not a change in the harmonic bracket.
3. **For THIS pair, DMM is structurally the wrong model.** Bi₂Te₃ and Sb₂Te₃ are nearly
   acoustically matched (Z₂/Z₁ = 0.71–0.85), so a smooth interface should transmit very
   well — AMM says 87 %. But DMM's α can never exceed ~0.5 by construction, so it is
   guaranteed to come out roughly 2× low. **Quote DMM as a lower bound, not as an
   estimate.** This matters for the manuscript, which currently leans on DMM alone.
4. **All three harmonic models FALL with T** (DMM −7.1 %, AMM −5.8 %, radiation limit
   −6.0 %) because the van der Waals gaps widen. The paper reports **+184 %**. Whatever
   drives that is not harmonic phonon transmission.
5. **The paper's 500 K point (98.60) exceeds the radiation limit (91.56) by 1.08×** — above
   a bound that no harmonic model can cross.

---

## 8. Where the numbers come from

- `s19_mismatch_models.py` → `s19_models.json` — all three models, same geometries.
- Sound speeds validated against the independent LAMMPS NEP elastic runs in
  `elastic_const/stress_strain/*/out.dat`: Sb₂Te₃ C33all 30.76 < C44all 35.37 GPa,
  Bi₂Te₃ C33all 46.85 > C44all 38.67 — confirming that **Sb₂Te₃ really does have LA slower
  than TA along z** in this potential.
- NEP caveat: against literature, Bi₂Te₃ C₃₃ is excellent (46.85 vs ~47.7 GPa) but C₄₄ is
  ~40 % too stiff (38.67 vs ~27.4), and Sb₂Te₃ C₃₃ is ~30 % too soft (30.76 vs ~44.6). The
  potential's weak spot is cross-plane stiffness in Sb₂Te₃ — the transport direction here.
- Figure: `s20_figure_GK.py` writes `Fig_GK_models.pdf/.png` in the published paper style.
  It **refuses to run** until ≥2 independent seeds exist at every temperature, so the NEMD
  points carry real seed-to-seed error bars; `--preview` renders the current state instead.
