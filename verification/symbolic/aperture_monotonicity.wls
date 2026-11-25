#!/usr/bin/env wolframscript
(* Global monotonicity of the aperture of the uniform-density brachistochrones.

   Criterion. With J' > 0, J(0) = 0, J(1-) = oo, the optical radial distance l and j = J(x) = b Sec[theta]:
     Delta/2 = Int_0^{Pi/2} Phi(b Sec[theta]) d theta,   Phi = dl/dj = 1/L,   L = dJ/dl = y (1 + x n'/n),
   and dDelta/db = 2 Int_0^{Pi/2} Phi'(b Sec[theta]) Sec[theta] d theta. Since dy/dj < 0, L_y <= 0 with isolated zeros
   implies Phi strictly decreasing, hence Delta strictly decreasing in b and in q. Geometrically L_y < 0 means that the
   optical section dl^2 + J^2 dphi^2 has negative Gaussian curvature -J_ll/J.

   Explicit form (a = Sqrt[1 - mu], Q = (y - a)(5a - y)): L = y + (1 - y^2)[+-1/(3a - y) + (3a - y)/Q], + for tau, - for t;
   N := (Q (3a - y))^2 L_y is a polynomial with N(a, a) = 32 a^4 (a^2 - 1) and
   N(a, 1) = (a - 1)(3a - 1)(5a - 1) E(a), E_tau = 15a^3 - 15a^2 + 9a - 1, E_t = 15a^3 + 5a^2 - 15a + 3.

   Thresholds of the pointwise criterion (exact decisions; each threshold is the unique root of its polynomial in (1/3, 1)):
     tau: L_y < 0 on (a, 1) iff a > a_tau, a_tau root of 1638781a^8 - 2030112a^6 + 418014a^4 - 30132a^2 + 729;
          at a_tau a single interior double zero, so Delta_tau is strictly decreasing for 0 < mu <= mu_tau = 1 - a_tau^2;
     t:   L_y < 0 on (a, 1) iff a >= a_t, a_t root of 15a^3 + 5a^2 - 15a + 3 (mu_t = 1 - a_t^2).
   Second, independent proof on smaller ranges: with a = (a0 + u)/(1 + u), y = (a + v)/(1 + v), u, v >= 0, every non-zero
   coefficient of the numerator of -L_y is positive for a0 = 3/8 (tau, mu <= 55/64) and a0 = 7/10 (t, mu <= 51/100).

   Integral extension for t to mu <= 429/625 (a >= 14/25). With A = Delta/2 and H = Phi'(j)/j, A'(b) = Int_b^oo Phi' dj/
   Sqrt[j^2 - b^2] and A'(0) = C_t; if H has a single sign change (from + to -), is decreasing where positive, and C_t <= 0,
   then A'(b) < 0 for b > 0. Checked on 14/25 <= a <= 7/10 (which contains [14/25, a_t)): N_t < 0 for a < y <= 9/10,
   dN_t/dy > 0 for 9/10 <= y <= 1, and N_t(a, 1) > 0 for a < a_t, so N_t has exactly one zero in (a, 1), in y > 9/10;
   d/dy (H/mu) > 0 on [9/10, 1]; C_t(429/625) < 0 by a rational certificate (truncated binomial series of the closed-form
   integrand, error bound 1/64, exact arcsine moments), and also by the value of the closed form.
   Diagnostics: Delta monotone on 217 radii near the proved limits. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "common.wl"}]];
q = (y - a) (5 a - y);
nx["proper"][x_] := Module[{yy = Sqrt[1 - (1 - a^2) x^2]}, (3 a - yy)/Sqrt[(yy - a) (5 a - yy)]];
nx["coordinate"][x_] := Module[{yy = Sqrt[1 - (1 - a^2) x^2]}, 4/((3 a - yy) Sqrt[(yy - a) (5 a - yy)])];
lUnified["proper"] = y + (1 - y^2) (1/(3 a - y) + (3 a - y)/q);
lUnified["coordinate"] = y + (1 - y^2) (-1/(3 a - y) + (3 a - y)/q);
dom = 1/3 < a < 1 && a < y < 1;
uniqueRootIn[poly_] := Module[{rs = Select[Table[Root[poly, k], {k, Exponent[poly[x], x]}], Element[#, Reals] &]},
   rs = Select[rs, 1/3 < # < 1 &];
   If[Length[rs] != 1, tovAbort["threshold polynomial does not have a unique root in (1/3, 1)"]];
   First[rs]];

tovRecord["monotonicity_global.angular_substitution", "exact",
  "with j = b Sec[theta]: (b/j) dj/Sqrt[j^2 - b^2] = d theta",
  Simplify[((b/j) dj/Sqrt[j^2 - b^2]) /. {j -> b Sec[th], dj -> b Sec[th] Tan[th]}, b > 0 && 0 < th < Pi/2] === 1];

Do[Module[{lx, ly},
   lx = Sqrt[1 - (1 - a^2) x^2] (1 + x D[nx[m][x], x]/nx[m][x]);
   ly = lx /. x -> Sqrt[(1 - y^2)/(1 - a^2)];
   tovRecord["monotonicity_global.l_unified." <> m, "exact",
    "L = y (1 + x n'/n) equals the unified form y + (1 - y^2)[+-1/(3a - y) + (3a - y)/Q], " <> m <> " time",
    Simplify[ly - lUnified[m], dom] === 0]],
  {m, {"proper", "coordinate"}}];

nPoly[m_] := Factor[Together[(q (3 a - y))^2 D[lUnified[m], y]]];
eFactor["proper"] = 15 a^3 - 15 a^2 + 9 a - 1;
eFactor["coordinate"] = 15 a^3 + 5 a^2 - 15 a + 3;
Do[tovRecord["monotonicity_global.n_endpoints." <> m, "exact",
   "N is a polynomial with N(a, a) = 32 a^4 (a^2 - 1) and N(a, 1) = (a - 1)(3a - 1)(5a - 1) E(a), " <> m <> " time",
   PolynomialQ[nPoly[m], y] && Simplify[(nPoly[m] /. y -> a) - 32 a^4 (a^2 - 1)] === 0 &&
    Simplify[(nPoly[m] /. y -> 1) - (a - 1) (3 a - 1) (5 a - 1) eFactor[m]] === 0],
  {m, {"proper", "coordinate"}}];

(* exact decisions of the pointwise criterion *)
aTau = uniqueRootIn[729 - 30132 #^2 + 418014 #^4 - 2030112 #^6 + 1638781 #^8 &];
aT = uniqueRootIn[3 - 15 # + 5 #^2 + 15 #^3 &];
tovValue["mu_tau", 1 - aTau^2, 25];
tovValue["mu_t_pointwise", 1 - aT^2, 25];
(* The decision is made on the polynomial N, whose sign equals that of L_y because L_y = N/(Q (3a - y))^2 with a positive
   denominator on the domain. Note: Reduce with the quantifier applied to the rational function L_y (after Together)
   returned False for both indices in Mathematica 13.3, which is wrong; the polynomial form gives the thresholds that are
   also certified by Sturm sequences and by the positivity certificate below. *)
Do[tovRecord["monotonicity_global.sign_of_l_y." <> m, "exact",
   "L_y = N/(Q (3a - y))^2 with (Q (3a - y))^2 > 0 on a < y < 1, " <> m <> " time",
   Simplify[D[lUnified[m], y] - nPoly[m]/(q (3 a - y))^2] === 0 &&
    Reduce[dom && (q (3 a - y))^2 <= 0, {a, y}, Reals] === False],
  {m, {"proper", "coordinate"}}];
boundaryOf[r_] := Replace[r, {HoldPattern[Inequality[c_, op_, a, Less, 1]] :> {c, op}, _ :> {Indeterminate, None}}];
Module[{r = Reduce[1/3 < a < 1 && ForAll[y, a < y < 1, nPoly["proper"] < 0], a, Reals], bc},
  bc = boundaryOf[r];
  tovRecord["monotonicity_global.threshold.proper", "exact",
   "N_tau < 0 on (a, 1) iff a > a_tau (unique root in (1/3, 1) of the degree-8 polynomial)",
   bc[[2]] === Less && RootReduce[bc[[1]] - aTau] === 0]];
Module[{r = Reduce[1/3 < a < 1 && ForAll[y, a < y < 1, nPoly["coordinate"] < 0], a, Reals], bc},
  bc = boundaryOf[r];
  tovRecord["monotonicity_global.threshold.coordinate", "exact",
   "N_t < 0 on (a, 1) iff a >= a_t (unique root in (1/3, 1) of 15a^3 + 5a^2 - 15a + 3)",
   bc[[2]] === LessEqual && RootReduce[bc[[1]] - aT] === 0]];
Module[{g = PolynomialGCD[nPoly["proper"] /. a -> aTau, D[nPoly["proper"] /. a -> aTau, y], Extension -> Automatic], yc},
  If[Exponent[g, y] != 1, tovAbort["no double zero of N_tau at the threshold"]];
  yc = y /. First[Solve[g == 0, y]];
  tovValue["tau_threshold_double_zero_y", yc, 20];
  tovRecord["monotonicity_global.threshold_double_zero.proper", "exact",
   "at a = a_tau, N_tau has a single double zero, interior to (a, 1): Phi' vanishes only at an isolated point, so the threshold is included",
   Exponent[g, y] == 1 && aTau < yc < 1]];

(* second, independent proof: positivity of coefficients *)
Do[Module[{expr, num, coeffs},
   expr = Together[(-D[lUnified[c[[1]]], y]) /. y -> (a + v)/(1 + v) /. a -> (c[[2]] + u)/(1 + u)];
   num = Expand[Numerator[expr]];
   coeffs = Select[Flatten[CoefficientList[num, {u, v}]], # != 0 &];
   tovRecord["monotonicity_global.positivity_certificate." <> c[[1]], "exact",
    "after a = (a0 + u)/(1 + u), y = (a + v)/(1 + v), every non-zero coefficient of the numerator of -L_y is positive (a0 = " <>
     ToString[c[[2]], InputForm] <> ") and the denominator is a product of squares and positive linear forms",
    Min[coeffs] > 0 && Reduce[u >= 0 && v > 0 && Denominator[expr] <= 0, {u, v}, Reals] === False]],
  {c, {{"proper", 3/8}, {"coordinate", 7/10}}}];

(* negative control: at a0 = 1/3 (the Buchdahl limit, where Delta_tau is not monotone) the certificate must fail *)
Module[{expr, coeffs},
  expr = Together[(-D[lUnified["proper"], y]) /. y -> (a + v)/(1 + v) /. a -> (1/3 + u)/(1 + u)];
  coeffs = Select[Flatten[CoefficientList[Expand[Numerator[expr]], {u, v}]], # != 0 &];
  tovRecord["monotonicity_global.negative_control.buchdahl_certificate", "negative_control",
   "the positivity certificate for tau fails at a0 = 1/3 (some coefficient is negative)", Min[coeffs] < 0]];

(* integral extension for t to a >= 14/25, that is mu <= 429/625 *)
nT = nPoly["coordinate"];
extDom = 14/25 <= a <= 7/10;   (* contains [14/25, a_t), a_t = 0.6931... *)
tovRecord["monotonicity_global.extension.negative_core", "exact",
  "N_t < 0 for a < y <= 9/10 and a in [14/25, 7/10]: the positive part of Phi' lies in y > 9/10",
  Reduce[extDom && a < y <= 9/10 && nT >= 0, {a, y}, Reals] === False];
tovRecord["monotonicity_global.negative_control.extension_cut", "negative_control",
  "the cut at y = 9/10 fails if the range is extended to a = 1/2: N_t(1/2, 9/10) = 5317/31250 > 0 (rational point test)",
  With[{v = nT /. {a -> 1/2, y -> 9/10}}, Head[v] === Rational && v > 0]];
tovRecord["monotonicity_global.extension.increasing_tail", "exact",
  "dN_t/dy > 0 on [14/25, 7/10] x [9/10, 1]",
  Reduce[extDom && 9/10 <= y <= 1 && D[nT, y] <= 0, {a, y}, Reals] === False];
tovRecord["monotonicity_global.extension.single_sign_change", "exact",
  "N_t(a, 1) > 0 exactly for a < a_t on [14/25, 7/10]; with the two records above, N_t has exactly one zero in (a, 1) for a in [14/25, a_t), from negative to positive toward the center",
  Reduce[14/25 <= a < aT && (nT /. y -> 1) <= 0, a, Reals] === False && RootReduce[nT /. {a -> aT, y -> 1}] === 0];
Module[{fF = 16 (1 - y^2)/((3 a - y)^2 q), hOverMu},
  hOverMu = 2 D[1/lUnified["coordinate"], y]/D[fF, y];
  tovRecord["monotonicity_global.extension.h_monotone", "exact",
   "d/dy (H/mu) > 0 on [14/25, 7/10] x [9/10, 1], i.e. H = Phi'(j)/j decreasing in j where Phi' > 0",
   Reduce[extDom && 9/10 <= y <= 1 && Together[D[hOverMu, y]] <= 0, {a, y}, Reals] === False]];
(* C_t(429/625) < 0 by exact arithmetic (a = 14/25, d = 1/a - 1 = 11/14). The paper writes
   C_t = -(Sqrt[mu] a^2/2) Int_0^1 p_d W_d/Sqrt[s (1 - s)] ds, p_d = 2 - 4 d s + d^2 s^2, W_d = (1 + d s)/Sqrt[B],
   B = (4 - d s)(2 + d (1 + s)). With z = 1 - B/16, 4 p_d W_d = p_d (1 + d s)(1 - z)^(-1/2). On [0, 1]: 8 < B < 16, so
   0 < z < 1/2; |p_d| <= 2 and 0 < 1 + d s < 2. The binomial coefficients Binomial[2k, k]/4^k are at most 1 for every k
   (Binomial[2k, k] is one term of the row sum 4^k; the finite check below is only a sanity check), so the tail
   after degree 8 is at most z^9/(1 - z) <= 1/256 and the error of the product is at most 1/64. The arcsine mean of the
   truncation is computed exactly with the moments E[s^k] = Binomial[2k, k]/4^k; it exceeds 1/64, so the exact mean is
   positive and C_t < 0. *)
Module[{d = 11/14, s, pd, bb, z, trunc, mean, ok},
  pd = 2 - 4 d s + d^2 s^2; bb = (4 - d s) (2 + d (1 + s)); z = 1 - bb/16;
  trunc = Expand[pd (1 + d s) Sum[Binomial[2 k, k]/4^k z^k, {k, 0, 8}]];
  mean = Total[CoefficientList[trunc, s] (Binomial[2 #, #]/4^# & /@ Range[0, Exponent[trunc, s]])];
  ok = d === 1/(14/25) - 1 &&
    Reduce[0 <= s <= 1 && ! (8 < bb < 16), s, Reals] === False &&
    Reduce[0 <= s <= 1 && ! (-2 <= pd <= 2), s, Reals] === False &&
    Reduce[0 <= s <= 1 && ! (0 < 1 + d s < 2), s, Reals] === False &&
    AllTrue[Range[0, 40], Binomial[2 #, #]/4^# <= 1 &] && mean - 1/64 > 0;
  tovValue["C_t_429_625.truncated_arcsine_mean", N[mean, 20], 12];
  tovRecord["monotonicity_global.extension.centre_slope_bound", "exact",
   "C_t(429/625) < 0: 8 < B < 16, |p_d| <= 2, 1 + d s < 2 on [0, 1], and the exact arcsine mean of the degree-8 binomial truncation exceeds the error bound 1/64",
   TrueQ[ok]]];
Module[{av = 14/25, mm, nn, ct},
  mm = (1 - av) (1 + 5 av)/(8 av); nn = (1 - av)/2;
  ct = N[(Sqrt[2 (1 - av^2)/av]/8) ((1 + 3 av + 6 av^2) EllipticK[mm] - 12 av^2 EllipticE[mm] - (1 + av) EllipticPi[nn, mm]), 30];
  tovValue["C_t_at_429_625", ct, 20];
  tovRecord["monotonicity_global.extension.centre_slope", "numerical_check",
   "C_t(429/625) < 0 (closed form evaluated at 30 digits; independent of the rational certificate above)", ct < 0]];

(* diagnostics near the proved limits *)
nsqY["proper"][z_, av_] := (3 av - z)^2/((z - av) (5 av - z));
nsqY["coordinate"][z_, av_] := 16/((3 av - z)^2 (z - av) (5 av - z));
halfAp[m_, muv_, qv_] := Module[{av = Sqrt[1 - muv], yq, bv, v},
   yq = Sqrt[1 - muv qv^2]; bv = Sqrt[nsqY[m][yq, av] (1 - yq^2)/muv];
   v = NIntegrate[Module[{yy = av + (yq - av) Sin[w]^2},
      bv Sqrt[muv]/((1 - yy^2) Sqrt[nsqY[m][yy, av] (1 - yy^2) - muv bv^2]) (yq - av) Sin[2 w]],
     {w, 0, Pi/2}, WorkingPrecision -> 40, PrecisionGoal -> 20, MaxRecursion -> 60];
   If[! NumberQ[v] || Abs[Im[v]] > 10^-20 || Re[v] <= 0, tovAbort["an aperture is not a positive real number"]];
   Re[v]];
grid = Join[Table[k/2000, {k, 1, 19}], Table[k/200, {k, 2, 199}]];
Do[tovRecord["monotonicity_global.diagnostic." <> c[[1]] <> ".mu_" <> tovTag[c[[2]]], "diagnostic",
   "Delta strictly decreasing on 217 radii at a compactness inside the proved range",
   AllTrue[Differences[halfAp[c[[1]], c[[2]], #] & /@ grid], # < 0 &]],
  {c, {{"proper", 886/1000}, {"coordinate", 429/625}}}];

(* shape of Psi_t beyond the proved range (Discussion of the paper). With P = -Q_t, R = N_t (3a - y)^3 D^2 and
   num = (P R_y - 3 P_y R)/((5a - y)(y - a)(3a - y)^2), the sign of d/dy (Psi/mu) on the positive part is that of num. *)
Module[{aStar = Sqrt[Root[3645 #^3 + 3888 #^2 - 1080 # + 64 &, 3]], pT, rT, numT},
  pT = -(y^4 - 3 a y^3 - (2 + 9 a^2) y^2 + (12 a + 15 a^3) y - 14 a^2);
  rT = nT (3 a - y)^3 q^2;
  numT = Cancel[(pT D[rT, y] - 3 D[pT, y] rT)/((5 a - y) (y - a) (3 a - y)^2)];
  tovRecord["monotonicity_global.coordinate.single_sign_change_to_mu_star", "exact",
   "for mu_t < mu < mu_* (a_* < a < a_t): N_t(a, a) < 0, N_t(a, 1) > 0 and every zero of N_t in (a, 1) has dN_t/dy > 0, so Psi_t changes sign exactly once, from positive near the center to negative",
   Reduce[aStar < a < aT && ((nT /. y -> a) >= 0 || (nT /. y -> 1) <= 0), a, Reals] === False &&
    Reduce[aStar < a < aT && a < y < 1 && nT == 0 && D[nT, y] <= 0, {a, y}, Reals] === False];
  tovRecord["monotonicity_global.coordinate.psi_hump_witness", "exact",
   "at mu = 21/25 (a = 2/5), (a, y) = (2/5, 9/10): P > 0, N_t > 0 and the reduced numerator of d/dy (Psi/mu) is negative, so Psi_t is positive and increasing in j near that point (rational point test); the monotonicity hypothesis of the integral criterion fails there",
   With[{pt = {a -> 2/5, y -> 9/10}}, With[{vp = pT /. pt, vn = nT /. pt, vd = numT /. pt},
     AllTrue[{vp, vn, vd}, Head[#] === Rational &] && vp > 0 && vn > 0 && vd < 0]]]];
(* diagnostics beyond 429/625: decreasing below mu_Delta, a single maximum above it (sampled, not proved) *)
Module[{vals = halfAp["coordinate", 3/4, #] & /@ grid},
  tovRecord["monotonicity_global.diagnostic.coordinate.mu_3_4", "diagnostic",
   "mu = 3/4 (between 429/625 and mu_Delta): Delta_t strictly decreasing on 217 radii", AllTrue[Differences[vals], # < 0 &]]];
Do[Module[{vals = halfAp["coordinate", muv, #] & /@ grid, d},
   d = Sign[Differences[vals]];
   tovRecord["monotonicity_global.diagnostic.coordinate_unimodal.mu_" <> tovTag[muv], "diagnostic",
    "mu = " <> ToString[muv, InputForm] <> " (between mu_Delta and mu_*): on 217 radii Delta_t starts above pi, increases and then decreases, with a single change of direction",
    FreeQ[d, 0] && Length[Split[d]] == 2 && First[d] == 1 && 2 First[vals] > Pi]],
  {muv, {4/5, 21/25}}];

tovFinish["aperture_monotonicity"];
