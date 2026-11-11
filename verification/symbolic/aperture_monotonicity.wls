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

   Integral extension for t to mu <= 16/25 (a >= 3/5). With A = Delta/2 and H = Phi'(j)/j, A'(b) = Int_b^oo Phi' dj/
   Sqrt[j^2 - b^2] and A'(0) = C_t; if H has a single sign change (from + to -), is decreasing where positive, and C_t <= 0,
   then A'(b) < 0 for b > 0. Checked: one sign change of N_t on (a, 1) for a in [3/5, a_t); N_t(a, 9/10) < 0 on [3/5, 7/10];
   d/dy (H/mu) > 0 on [3/5, 7/10] x [9/10, 1]; C_t(16/25) < 0 by an analytic inequality between K, E and Pi (exact
   algebraic steps checked here), and also by the value of the closed form.
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

(* integral extension for t to a >= 3/5 *)
nT = nPoly["coordinate"];
(* N_t is monic in y up to a constant, so the number of its zeros in (a, 1) can change only where its discriminant vanishes
   or where N_t vanishes at y = a or y = 1. On [3/5, a_t): the discriminant has no zero, N_t(a, a) = 32 a^4 (a^2 - 1) < 0 and
   N_t(a, 1) vanishes only at a_t; one sample then fixes the count to one. *)
tovRecord["monotonicity_global.extension.single_sign_change", "exact",
  "on [3/5, a_t) the discriminant of N_t in y has no zero, N_t does not vanish at y = a or y = 1, and N_t has exactly one zero in (a, 1) at a = 3/5",
  Reduce[3/5 <= a < aT && Discriminant[nT, y] == 0, a, Reals] === False &&
   Reduce[3/5 <= a < aT && ((nT /. y -> a) == 0 || (nT /. y -> 1) == 0), a, Reals] === False &&
   CountRoots[nT /. a -> 3/5, {y, 3/5, 1}] == 1];
tovRecord["monotonicity_global.extension.negative_core", "exact",
  "N_t(a, 9/10) < 0 for a in [3/5, 7/10]: the positive part of Phi' lies in y > 9/10",
  Reduce[3/5 <= a <= 7/10 && (nT /. y -> 9/10) >= 0, a, Reals] === False];
Module[{fF = 16 (1 - y^2)/((3 a - y)^2 q), hOverMu},
  hOverMu = 2 D[1/lUnified["coordinate"], y]/D[fF, y];
  tovRecord["monotonicity_global.extension.h_monotone", "exact",
   "d/dy (H/mu) > 0 on [3/5, 7/10] x [9/10, 1], i.e. H = Phi'(j)/j decreasing in j where Phi' > 0",
   Reduce[3/5 <= a <= 7/10 && 9/10 <= y <= 1 && Together[D[hOverMu, y]] <= 0, {a, y}, Reals] === False]];
(* C_t(16/25) < 0 without numerical evaluation (a = 3/5). With s = Sin[theta]^2, Sqrt[1 - m s] = 1/Sqrt[1 - m s] -
   m s/Sqrt[1 - m s] gives E = K - m I, I = Int_0^{Pi/2} s d theta/Sqrt[1 - m s]. Analytic steps (not executed here):
   Pi(nu|m) >= K for nu > 0; K >= Pi/2 and I <= (Pi/4)/Sqrt[1 - m], so I/K <= 1/(2 Sqrt[1 - m]). Exact steps below: m = 1/3,
   nu = 1/5, the bracket is at most -(24/25) K + (36/25) I, and (36/25)/(2 Sqrt[1 - m]) < 24/25, so the bracket is negative. *)
Module[{av = 3/5, mm, nn, kk, ee, ii, pp, bracket, ok},
  mm = (1 - av) (1 + 5 av)/(8 av); nn = (1 - av)/2;
  bracket = (1 + 3 av + 6 av^2) kk - 12 av^2 ee - (1 + av) pp;
  ok = mm === 1/3 && nn === 1/5 &&
    Simplify[Sqrt[1 - mm s] - (1/Sqrt[1 - mm s] - mm s/Sqrt[1 - mm s]), 0 < s < 1] === 0 &&
    Expand[(bracket /. {ee -> kk - mm ii, pp -> kk}) - (-(24/25) kk + (36/25) ii)] === 0 &&
    Coefficient[bracket, pp] < 0 &&
    (36/25)/(2 Sqrt[1 - mm]) < 24/25;
  tovRecord["monotonicity_global.extension.centre_slope_bound", "exact",
   "C_t(16/25) < 0 analytically: m = 1/3, E = K - m I, Pi >= K, I/K <= 1/(2 Sqrt[1 - m]) < 2/3 make the bracket negative",
   TrueQ[ok]]];
Module[{av = 3/5, mm, nn, ct},
  mm = (1 - av) (1 + 5 av)/(8 av); nn = (1 - av)/2;
  ct = N[(Sqrt[2 (1 - av^2)/av]/8) ((1 + 3 av + 6 av^2) EllipticK[mm] - 12 av^2 EllipticE[mm] - (1 + av) EllipticPi[nn, mm]), 30];
  tovValue["C_t_at_16_25", ct, 20];
  tovRecord["monotonicity_global.extension.centre_slope", "numerical_check",
   "C_t(16/25) < 0 (closed form evaluated at 30 digits; independent of the analytic bound above)", ct < 0]];

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
  {c, {{"proper", 886/1000}, {"coordinate", 16/25}}}];

tovFinish["aperture_monotonicity"];
