#!/usr/bin/env wolframscript
(* Hamilton-Jacobi identity and the antipodal comparison for the uniform-density brachistochrones.

   Along a boundary-to-boundary extremal with Clairaut constant b = J(q):
     T(q) = 2 Int_q^1 n^2 x Sqrt[h] dx/Sqrt[J^2 - b^2]   (optical length),
     Delta(q) = 2 Int_q^1 b Sqrt[h] dx/(x Sqrt[J^2 - b^2]),
     W = T - b Delta = 2 Int_q^1 Sqrt[h] Sqrt[J^2 - b^2]/x dx,  dW/db = -Delta,  hence dT/dq = b dDelta/dq.
   The physical times are tau = T for the proper-time index and t = a T for the coordinate-time index.

   Claims checked:
     - the integrand identity behind W (exact);
     - dT/dq = b dDelta/dq by central differences at sample points, both indices (numerical check);
     - antipodal points, coordinate time at mu = 4/5 (mu_Delta < mu < mu_star): a radius q with Delta(q) = Pi, found by
       a scan on a grid of step 1/100 and bisection, gives a coordinate time below the diameter time (numerical check).
       The scan samples the aperture; it does not certify that q is the smallest root. The analytic statement (Theorem B)
       does not rest on this check. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "common.wl"}]];

tovRecord["hj.w_integrand_identity", "exact",
  "n^2 x Sqrt[h]/Sqrt[J^2 - b^2] - b^2 Sqrt[h]/(x Sqrt[J^2 - b^2]) = Sqrt[h] Sqrt[J^2 - b^2]/x with J = n x",
  Simplify[(nn^2 x sh/Sqrt[nn^2 x^2 - bb^2] - bb^2 sh/(x Sqrt[nn^2 x^2 - bb^2])) - sh Sqrt[nn^2 x^2 - bb^2]/x,
    nn > 0 && x > 0 && nn x > bb > 0 && sh > 0] === 0];

nsqY["proper"][z_, a_] := (3 a - z)^2/((z - a) (5 a - z));
nsqY["coordinate"][z_, a_] := 16/((3 a - z)^2 (z - a) (5 a - z));
(* half-aperture and half optical length in y = Sqrt[1 - mu x^2]; every value must be a positive real number *)
quadY[m_, muv_, qv_, g_] := Module[{a = Sqrt[1 - muv], yq, bv, v},
   yq = Sqrt[1 - muv qv^2]; bv = Sqrt[nsqY[m][yq, a] (1 - yq^2)/muv];
   v = NIntegrate[Module[{yy = a + (yq - a) Sin[u]^2, n2},
       n2 = nsqY[m][yy, a];
       g[n2, bv, yy] (yq - a) Sin[2 u]/Sqrt[n2 (1 - yy^2) - muv bv^2]], {u, 0, Pi/2},
      WorkingPrecision -> 50, PrecisionGoal -> 30, MaxRecursion -> 60];
   If[! NumberQ[v] || Abs[Im[v]] > 10^-30 Abs[v] || Re[v] <= 0, tovAbort["a quadrature is not a positive real number"]];
   Re[v]];
delta[m_, muv_, qv_?NumericQ] := 2 quadY[m, muv, qv, (#2 Sqrt[muv]/(1 - #3^2)) &];
tOpt[m_, muv_, qv_?NumericQ] := 2 quadY[m, muv, qv, (#1/Sqrt[muv]) &];
bOf[m_, muv_, qv_] := qv Sqrt[nsqY[m][Sqrt[1 - muv qv^2], Sqrt[1 - muv]]];

Do[Module[{h = 10^-6, dT, dD, rel},
   dT = (tOpt[m, s[[1]], s[[2]] + h] - tOpt[m, s[[1]], s[[2]] - h])/(2 h);
   dD = (delta[m, s[[1]], s[[2]] + h] - delta[m, s[[1]], s[[2]] - h])/(2 h);
   rel = Abs[dT - bOf[m, s[[1]], s[[2]]] dD]/Abs[dT];
   tovValue["hj_relative_residual." <> m <> ".mu_" <> tovTag[s[[1]]] <> ".q_" <> tovTag[s[[2]]], rel, 5];
   tovRecord["hj.derivative_identity." <> m <> ".mu_" <> tovTag[s[[1]]] <> ".q_" <> tovTag[s[[2]]], "numerical_check",
    "dT/dq = b dDelta/dq by central differences (step 1e-6), relative residual below 1e-9", rel < 10^-9]],
  {m, {"proper", "coordinate"}}, {s, {{1/10, 1/2}, {1/2, 3/10}, {4/5, 7/10}, {4/5, 1/10}}}];

Module[{muv = 4/5, a = Sqrt[1 - 4/5], grid, vals, k, lo, hi, mid, q1, res, tq1, tDiam},
  grid = Range[1/100, 99/100, 1/100];
  vals = (delta["coordinate", muv, #] - Pi) & /@ grid;
  k = FirstPosition[vals, _?Negative, {None}][[1]];
  If[k === None || k < 2 || ! AllTrue[vals[[;; k - 1]], Positive],
   tovAbort["Delta - Pi does not change sign once from positive to negative on the left of the grid"]];
  {lo, hi} = {grid[[k - 1]], grid[[k]]};
  Do[mid = (lo + hi)/2; If[delta["coordinate", muv, mid] > Pi, lo = mid, hi = mid], {50}];
  q1 = (lo + hi)/2; res = delta["coordinate", muv, q1] - Pi;
  tq1 = a tOpt["coordinate", muv, q1];
  tDiam = a 2 NIntegrate[Sqrt[nsqY["coordinate"][Sqrt[1 - muv x^2], a]]/Sqrt[1 - muv x^2], {x, 0, 1},
     WorkingPrecision -> 40, PrecisionGoal -> 25, MaxRecursion -> 60];
  tovValue["antipodal.mu_4_5.q_bracket_lower", N[lo, 25], 20];
  tovValue["antipodal.mu_4_5.q_bracket_upper", N[hi, 25], 20];
  tovValue["antipodal.mu_4_5.aperture_residual", res, 5];
  tovValue["antipodal.mu_4_5.coordinate_time_extremal_R_over_c", tq1, 15];
  tovValue["antipodal.mu_4_5.coordinate_time_diameter_R_over_c", tDiam, 15];
  tovValue["antipodal.mu_4_5.saving_percent_of_diameter_time", 100 (tDiam - tq1)/tDiam, 8];
  tovRecord["hj.antipodal.aperture_pi", "numerical_check",
   "mu = 4/5, coordinate time: the bisected radius has |Delta(q) - Pi| < 1e-12 (grid step 1/100, 50 bisections)",
   Abs[res] < 10^-12];
  tovRecord["hj.antipodal.faster_than_diameter", "numerical_check",
   "mu = 4/5: the coordinate time along that antipodal regular extremal is below the diameter time (units R/c)",
   tq1 < tDiam]];

tovFinish["hamilton_jacobi"];
