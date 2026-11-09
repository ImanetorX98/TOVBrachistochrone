#!/usr/bin/env wolframscript
(* Centre expansion of the aperture and the threshold mu_Delta (uniform density).

   Claims checked:
     - C = lim_{eps->0}[Int_eps^1 Sqrt[h] dx/(x J) - 1/J(eps)] = -Sqrt[mu] Int_a^1 y (1/n)'(y) dy/Sqrt[1 - y^2]
       (integration by parts; numerical check of the two expressions) and Delta/2 = Pi/2 + C b + O(b^3)
       (numerical check of the slope at q = 1e-3, 1e-4, 1e-5);
     - (1/n_tau)' = 4a^2/((3a - y)^2 Sqrt[Q]) and (1/n_t)' = (y^2 - 6ay + 7a^2)/(2 Sqrt[Q]) (exact), hence
       C_tau < 0 for all mu < 8/9 and C_t < 0 for mu <= mu0 = 1 - 1/(11 - 6 Sqrt[2]) (exact sign decision);
     - existence and uniqueness of the zero of C_t (weighted-mean argument): with d = 1/a - 1 and s = (y - a)/(1 - a),
       C_t = -(Sqrt[mu] a^2/2) Int_0^1 p_d W_d ds/Sqrt[s(1 - s)], p_d = 2 - 4ds + d^2 s^2,
       W_d = (1 + ds)/Sqrt[(4 - ds)(2 + d(1 + s))]; the change of variables, the four sign conditions, the arcsine
       mean 2 - 2d + 3d^2/8 and its value -5/32 at d = 3/2 (mu = 21/25) are checked exactly. The covariance
       (Chebyshev) step is an analytic argument recorded in the paper, not executed here;
     - mu_Delta located by bisection on the closed form of C_t (the digits evaluate the exact equation);
     - L_tau = 1/F = y + 4a^2 (1 - y^2)/((3a - y) Q) > 1 (exact), hence 0 < Delta_tau < Pi for all mu < 8/9;
     - Delta_tau is not monotone near the Buchdahl limit (numerical check at mu = 8/9 - 1e-12). *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "common.wl"}]];
qq = (y - a) (5 a - y);
invN["proper"] = Sqrt[qq]/(3 a - y);
invN["coordinate"] = (3 a - y) Sqrt[qq]/4;
nsqY[m_][z_, av_] := 1/(invN[m] /. {y -> z, a -> av})^2;
wp = 40;

tovRecord["threshold.derivative_proper", "exact", "(1/n_tau)' = 4a^2/((3a - y)^2 Sqrt[Q])",
  Simplify[D[invN["proper"], y] - 4 a^2/((3 a - y)^2 Sqrt[qq]), 1/3 < a < 1 && a < y < 1] === 0];
tovRecord["threshold.derivative_coordinate", "exact", "(1/n_t)' = (y^2 - 6ay + 7a^2)/(2 Sqrt[Q])",
  Simplify[D[invN["coordinate"], y] - (y^2 - 6 a y + 7 a^2)/(2 Sqrt[qq]), 1/3 < a < 1 && a < y < 1] === 0];
tovRecord["threshold.sufficient_sign_bound", "exact",
  "y^2 - 6ay + 7a^2 > 0 on a < y < 1 whenever a >= 1/(3 - Sqrt[2]), i.e. mu <= 1 - 1/(11 - 6 Sqrt[2])",
  Reduce[1/(3 - Sqrt[2]) <= a < 1 && a < y < 1 && y^2 - 6 a y + 7 a^2 <= 0, {a, y}, Reals] === False &&
   Simplify[1 - 1/(3 - Sqrt[2])^2 - (1 - 1/(11 - 6 Sqrt[2]))] === 0];
tovValue["mu0", 1 - 1/(11 - 6 Sqrt[2])];

(* integration by parts and centre slope (numerical checks) *)
cIBP[m_, muv_] := Module[{av = Sqrt[1 - muv], dinv},
   dinv[yy_] = D[invN[m] /. a -> av, y] /. y -> yy;
   -Sqrt[muv] NIntegrate[Module[{yy = av + (1 - av) Sin[u]^2},
      yy dinv[yy]/Sqrt[1 - yy^2] 2 (1 - av) Sin[u] Cos[u]],
     {u, 0, Pi/2}, WorkingPrecision -> wp, PrecisionGoal -> 25, MaxRecursion -> 40]];
cLimit[m_, muv_, eps_] := Module[{av = Sqrt[1 - muv], jj, yx},
   yx[x_] := Sqrt[1 - muv x^2]; jj[x_] := x Sqrt[nsqY[m][yx[x], av]];
   NIntegrate[1/(yx[x] x jj[x]), {x, eps, 1}, WorkingPrecision -> wp, PrecisionGoal -> 20,
     MaxRecursion -> 40] - 1/jj[eps]];
halfAp[m_, muv_, qv_] := Module[{av = Sqrt[1 - muv], yq, bv},
   yq = Sqrt[1 - muv qv^2]; bv = Sqrt[nsqY[m][yq, av] (1 - yq^2)/muv];
   NIntegrate[Module[{yy = av + (yq - av) Sin[u]^2},
      bv Sqrt[muv]/((1 - yy^2) Sqrt[nsqY[m][yy, av] (1 - yy^2) - muv bv^2]) (yq - av) Sin[2 u]],
     {u, 0, Pi/2}, WorkingPrecision -> wp, PrecisionGoal -> 25, MaxRecursion -> 50]];
jOf[m_, muv_, qv_] := qv Sqrt[nsqY[m][Sqrt[1 - muv qv^2], Sqrt[1 - muv]]];
Do[Module[{c = cIBP[m, muv], tag = m <> ".mu_" <> tovTag[muv], ratios},
   tovRecord["threshold.by_parts." <> tag, "numerical_check",
    "the integrated-by-parts C equals the defining limit at eps = 1e-8 to 1e-6 relative",
    Abs[c - cLimit[m, muv, 10^-8]] < 10^-6 Abs[c]];
   ratios = Table[(halfAp[m, muv, qv] - Pi/2)/jOf[m, muv, qv], {qv, {10^-3, 10^-4, 10^-5}}];
   tovRecord["threshold.centre_slope." <> tag, "numerical_check",
    "(Delta/2 - Pi/2)/b tends to C as q = 1e-3, 1e-4, 1e-5, with a decreasing error",
    Abs[ratios[[3]] - c] < 10^-3 Abs[c] && Abs[ratios[[3]] - c] < Abs[ratios[[1]] - c]]],
  {m, {"proper", "coordinate"}}, {muv, {1/10, 1/2, 4/5}}];

(* uniqueness of the zero of C_t (exact parts of the weighted-mean argument) *)
pd = 2 - 4 d s + d^2 s^2;
wd = (1 + d s)/Sqrt[(4 - d s) (2 + d (1 + s))];
Module[{orig, sub, rhs},
  orig = (y (y^2 - 6 a y + 7 a^2))/Sqrt[(y - a) (5 a - y) (1 - y^2)];
  sub = orig /. y -> a + (1 - a) s;
  rhs = a^2 pd wd/Sqrt[s (1 - s)] /. d -> 1/a - 1;
  (* with y = a + (1 - a) s = a (1 + d s), d = 1/a - 1: the numerator is a^3 (1 + d s) p_d exactly, and the radicand
     Q4(y) equals a^2 (1 - a)^2 s (1 - s)(4 - d s)(2 + d (1 + s)) exactly; all radicands are positive on the domain,
     so the integrand (1 - a) N/Sqrt[Q4] equals a^2 p_d W_d/Sqrt[s (1 - s)] with no sign ambiguity *)
  tovRecord["threshold.weighted_form", "exact",
   "numerator y(y^2 - 6ay + 7a^2) = a^3 (1 + ds) p_d and radicand Q4 = a^2 (1 - a)^2 s(1 - s)(4 - ds)(2 + d(1 + s)) at y = a(1 + ds)",
   Expand[((y (y^2 - 6 a y + 7 a^2)) /. y -> a (1 + d s)) - a^3 (1 + d s) pd] === 0 &&
    Simplify[(((y - a) (5 a - y) (1 - y^2)) /. y -> a + (1 - a) s) -
       (a^2 (1 - a)^2 s (1 - s) (4 - d s) (2 + d (1 + s)) /. d -> 1/a - 1)] === 0 &&
    Simplify[(a + (1 - a) s) - a (1 + d s) /. d -> 1/a - 1] === 0]];
dom = 0 < d < 2 && 0 < s < 1;
lw = Log[1 + d s] - Log[4 - d s]/2 - Log[2 + d (1 + s)]/2;
tovRecord["threshold.p_decreasing_in_d", "exact", "d_d p_d = -2s(2 - ds) < 0",
  Simplify[D[pd, d] + 2 s (2 - d s)] === 0 && Reduce[dom && 2 s (2 - d s) <= 0, {d, s}, Reals] === False];
tovRecord["threshold.p_decreasing_in_s", "exact", "d_s p_d < 0", Reduce[dom && D[pd, s] >= 0, {d, s}, Reals] === False];
tovRecord["threshold.weight_mixed_derivative", "exact", "d_s d_d log W_d > 0",
  Reduce[dom && Together[D[lw, s, d]] <= 0, {d, s}, Reals] === False];
tovRecord["threshold.weight_increasing_in_s", "exact", "d_s log W_d > 0",
  Reduce[dom && Together[D[lw, s]] <= 0, {d, s}, Reals] === False];
tovRecord["threshold.arcsine_mean", "exact", "arcsine mean of p_d = 2 - 2d + 3d^2/8, equal to -5/32 at d = 3/2 (mu = 21/25)",
  Simplify[Integrate[pd/(Pi Sqrt[s (1 - s)]), {s, 0, 1}] - (2 - 2 d + 3 d^2/8)] === 0 &&
   (2 - 2 d + 3 d^2/8 /. d -> 3/2) == -5/32 && (1 - (1/(1 + d))^2 /. d -> 3/2) == 21/25];

(* mu_Delta: bisection on the closed form of C_t *)
ctClosed[muv_] := Module[{av = Sqrt[1 - muv], mm, nn},
   mm = (1 - av) (1 + 5 av)/(8 av); nn = (1 - av)/2;
   N[(Sqrt[2 muv/av]/8) ((1 + 3 av + 6 av^2) EllipticK[mm] - 12 av^2 EllipticE[mm] - (1 + av) EllipticPi[nn, mm]), 70]];
Module[{lo = 7/10, hi = 4/5, mid},
  If[! (ctClosed[lo] < 0 < ctClosed[hi]), tovAbort["C_t does not change sign on [0.7, 0.8]"]];
  Do[mid = (lo + hi)/2; If[ctClosed[mid] < 0, lo = mid, hi = mid], {80}];
  tovValue["mu_Delta_bracket_lo", lo, 40];
  tovValue["mu_Delta_bracket_hi", hi, 40];
  tovValue["C_t_at_bracket_lo", ctClosed[lo], 10];
  tovValue["C_t_at_bracket_hi", ctClosed[hi], 10];
  tovValue["mu_Delta_evaluation_precision_digits", "70"];
  tovValue["mu_Delta_midpoint_not_certified", (lo + hi)/2, 26];
  tovRecord["threshold.mu_delta_bracket", "numerical_check",
   "the closed form of C_t is negative at the lower and positive at the upper end of the final bisection bracket (width 2^-80/10, closed form evaluated at 70 digits)",
   ctClosed[lo] < 0 < ctClosed[hi]]];

(* proper time: L_tau > 1, hence 0 < Delta_tau < Pi *)
Module[{nt, lt, claim},
  nt[x_] := With[{yy = Sqrt[1 - mu x^2], aa = Sqrt[1 - mu]}, (3 aa - yy)/Sqrt[(yy - aa) (5 aa - yy)]];
  lt = Sqrt[1 - mu x^2] (1 + x D[nt[x], x]/nt[x]);
  claim = With[{yy = Sqrt[1 - mu x^2], aa = Sqrt[1 - mu]}, yy + 4 aa^2 (1 - yy^2)/((3 aa - yy) (yy - aa) (5 aa - yy))];
  tovRecord["threshold.proper_l_identity", "exact", "L_tau = y (1 + x n'/n) = y + 4a^2 (1 - y^2)/((3a - y) Q)",
   Simplify[lt - claim, 0 < mu < 8/9 && 0 < x < 1] === 0]];
tovRecord["threshold.proper_l_bound", "exact", "L_tau > 1 for a < y < 1 and 1/3 < a < 1",
  Reduce[1/3 < a < 1 && a < y < 1 && y + 4 a^2 (1 - y^2)/((3 a - y) (y - a) (5 a - y)) <= 1, {a, y}, Reals] === False];
Module[{muB = 8/9 - 10^-12, h1, h2},
  h1 = halfAp["proper", muB, 1/100]; h2 = halfAp["proper", muB, 2/100];
  tovValue["half_aperture_proper.near_buchdahl.q_0.01", h1, 15];
  tovValue["half_aperture_proper.near_buchdahl.q_0.02", h2, 15];
  tovRecord["threshold.proper_not_monotone_near_buchdahl", "numerical_check",
   "at mu = 8/9 - 1e-12 the proper-time half-aperture increases from q = 0.01 to q = 0.02", h2 > h1]];

tovFinish["centre_slope_threshold"];
