# s19b — reviewer audit of the AMM/DMM temperature trend: what survived, what didn't

2026-10-05. An independent reviewer reinstalled calorine/phonopy, re-ran the pipeline and
challenged three supporting claims behind the "the decrease is quasi-harmonic" argument.
**All three challenges are correct.** Everything below is measured here, not accepted on
report: `s19b_mismatch_models.py` (order 2) and `s19b_cfit3.py` (order 3), logs `s19b.log`
and `s19b_cfit3.log`.

First, the control that makes the rest trustworthy: **the patched order-2 script reproduces
the backed-up `s19_models.json` to 1e-9 on all of G_DMM, G_AMM, G_rad(Bi/Sb), alpha_DMM,
alpha_AMM, a and cQL, at all four temperatures.** Backup kept as
`s19_models.json.bak-preReview`. So nothing below is a side effect of the patch.

---

## 1. The geometry/statistics split — my attribution was right, my arithmetic was loose

G_rad,Bi on a (geometry T) x (Bose T) grid, order-2 fit:

```
geom\Bose      200      300      400      500   classical
      200   105.85   106.86   107.22   107.38      107.68
      300   103.13   104.12   104.47   104.64      104.93
      400   100.55   101.52   101.87   102.03      102.31
      500    98.11    99.06    99.40    99.56       99.84
```

- statistics alone (top row):  **+1.45 %**
- geometry alone (left column): **-7.31 %**
- the diagonal, i.e. what the paper reports: **-5.95 %**

So the geometry costs **more** than the ~6 % I attributed to it; the diagonal is the two
effects partly cancelling. **The ratio is 5.0, not the "roughly six" I wrote.** Fix that.

Why the Bose factor is such a weak counterweight: the spectrum tops out at 5.2 THz, i.e.
hbar*w/k_B ~ 250 K, so most of the flux is already near-classical at 200 K. Hence the
classical (k_B) limit is only +1.4 % above the 500 K quantum value.

## 2. My mechanism sentence quotes the WRONG velocity — the Gamma slope understates it

Decomposing the flux integral (Bi2Te3, 200 -> 500 K, order 2):

```
sum(Phi)df acoustic : -8.59 %   (50.1 % of the total flux)
sum(Phi)df optical  : -5.98 %   (49.9 %)
mean positive v_z   : 173.6 -> 161.7 m/s   -6.87 %
band-mean freq, acoustic : -0.96 %      optical : -0.41 %
Gamma-slope LA      : 2162.3 -> 2068.7 m/s  -4.32 %
```

Three things follow, and the subsection currently gets the emphasis wrong:

- **The ACOUSTIC branches soften more than the optical ones** (-8.59 vs -5.98 %), and they
  are a 50/50 split of the flux. 
- **Frequencies barely move** (<1 %) while **v_z collapses by 6.9 %**. The vdW gap widening
  flattens the k_z dispersion *away from zone centre* without shifting frequencies much,
  and v_z is exactly what G_rad integrates.
- **The Gamma-slope sound speeds (-4.3 %) therefore UNDERSTATE the softening**, which is why
  quoting them alone cannot account for a -5.95 % drop in G_rad. The mesh-averaged v_z
  (-6.9 %) is the number that explains it. The subsection quotes only the Gamma slopes.

## 3. alpha_DMM_eff is partly circular — the real content is alpha(omega)

`alpha_DMM_eff` is *defined* in the script as `G_DMM/G_radB`, so "the decrease is inherited
from G_rad and not from alpha" is true by construction for that quantity. The defensible
version is that the **frequency-resolved** alpha(omega) is flat: re-weighting it classically
instead of with the Bose factor moves it from **0.3208 to 0.3211** at 200 K. Say that.

**And drop the "alpha_AMM rises" claim — it is fit-order dependent.** Order 2 gives
0.8734 -> 0.8747 (rising); order 3 gives 0.8589 -> 0.8462 (falling). What survives either
way is the magnitude: the transmissions move by <=1.5 % while the flux moves by 6-10 %.

## 4. The classical limit does NOT rescue the comparison — state it, don't leave it open

Removing the Bose factor entirely (the like-for-like comparison against classical NEMD):

```
order 2:  classical DMM_mod  69.16 -> 63.32  =  -8.45 %
order 3:  classical DMM_mod  76.30 -> 67.65  = -11.34 %
```

Still monotonically down, still nowhere near the NEMD rise. Worth one sentence in the paper
so a reader does not wonder whether the quantum factor was the problem. It buys +1.4 %.

## 5. THE BIG ONE: the QHA c-fit is not converged, and it is not a small effect

`Fmin_c` fits F(c) at fixed a with a **quadratic**. The sampled c grid spans **7.93 %** of
the mean, over which F(c) is visibly anharmonic. Fit residuals at 300 K, Bi2Te3:

```
order 2: max|residual| 0.4028 meV/atom, rms 0.2699
order 3: max|residual| 0.1007            rms 0.0626
order 4: max|residual| 0.0142            rms 0.0077
```

The quadratic is simply a bad fit. And the thermal expansion it returns is correspondingly
wrong:

```
Delta c_QL(Bi2Te3, 200->500 K):  order 2  +0.0269 A   order 3  +0.0534 A   order 4  +0.0593 A
```

**Orders 3 and 4 agree with each other (11 % apart) far better than either agrees with
order 2 (2x).** So order 2 is the outlier, and it UNDERESTIMATES the expansion by about a
factor of two.

### What that does to the whole table

```
 T    ---- order 2 (in the paper) ----   |   ---- order 3 (converged) ----
      cQL_Bi  DMMmod    AMM  rad_Sb      |   cQL_Bi  DMMmod    AMM  rad_Sb
 200 10.3665   67.92  92.45   95.75      |  10.3185   74.94  96.97  105.03
 300 10.3766   66.61  90.95   95.13      |  10.3367   72.46  93.96  103.76
 400 10.3855   64.94  89.03   93.50      |  10.3547   69.86  90.31  101.30
 500 10.3934   63.13  87.09   91.56      |  10.3719   67.45  86.50   98.59

 200->500 K     order 2     order 3      200 K absolute value shifts by
   DMM           -7.06 %    -10.00 %            +10.3 %
   AMM           -5.81 %    -10.80 %             +4.9 %
   rad (Bi)      -5.95 %     -9.46 %             +6.7 %
   rad (Sb)      -4.38 %     -6.13 %             +9.7 %
```

**This is not only a trend problem. Every absolute number in Table 1 moves by 5-10 %.**
The manuscript's "4--7 %" becomes "6--11 %", and 74.94 vs Chowdhury's 75 MW/m2K is a
coincidence worth noticing but not leaning on.

**The sign and the mechanism are fit-independent. The magnitude is uncertain by ~1.6-2x.**

### What to do

Either (a) refine the QHA grid — a denser, narrower c-grid around the minimum, which is
pure NEP compute and therefore cheap — and requote everything, or (b) quote the decrease
qualitatively and stop giving percentages to two significant figures. **Do not leave
order-2 percentages in the paper as if they were converged.** Option (a) is the right one:
the grid is the actual defect, and `s8_qha.py` can simply be rerun on a tighter mesh.

## 6. Mesh convergence is fine

The reviewer pushed 24x24x8 -> 32x32x16 and the trend was stable to +-0.1 %. The k_z
sampling is not the issue. (Not re-measured here; recorded as their result.)
