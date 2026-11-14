#!/usr/bin/env wolframscript
(* A numerical candidate for a coexistence (Maxwell) separation for the proper time of the uniform-density interior near the Buchdahl bound.

   At mu = 8/9 - 10^-4 the proper-time aperture Delta_tau(q) falls from Pi to a local minimum, rises to a local maximum
   and then falls to 0. For an angle theta between these two values the equation Delta_tau(q) = theta has three roots: an
   inner one, one on the rising branch (never minimal, Jacobi exclusion) and an outer one. The separation theta_M at which
   the inner and the outer extremal take the same proper time is located by regula falsi (Illinois variant). By the
   Hamilton-Jacobi identity this is the equal-area condition Int (theta - Delta) db = 0 between the two roots.

   Claims checked (numerical checks and one diagnostic; the existence of a coexistence separation for compactness close
   to 8/9 is proved analytically in the paper, but this script does not certify that theta_M is that separation):
     - the aperture is not monotone: the located local minimum is below the located local maximum, at a larger radius;
     - at theta_M the inner and outer extremals (aperture residual below 1e-9) have equal proper time within 1e-10
       (units R/c). At 36 working digits the roots are reached to about 1e-11 (at 50 digits to about 1e-16, see the
       paper's remark); the time difference changes by about 5e-8 per 1e-6 in theta, and it has opposite signs at
       theta_M -/+ 1e-6 (a bracket of the coexistence angle);
     - the extremal on the rising branch (aperture residual below 1e-9, positive sampled slope) is slower than the
       inner and outer ones by more than 1e-6;
     - diagnostic: on a grid of 205 radii, Delta - theta_M changes sign in exactly three cells, those of the three
       located roots. A grid can miss tangential roots or two crossings in one cell, so this is not a root count.
   The extrema of Delta are located as zeros of a centred difference in brackets read off a sampled profile; they are
   not certified to be the only critical points. The digits are a numerical localization at one compactness. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "common.wl"}]];

nsq[z_, a_] := (3 a - z)^2/((z - a) (5 a - z));      (* proper-time index squared, as a function of y *)
muv = 8/9 - 1/10^4; a0 = Sqrt[1 - muv];
(* half-quantities in y = Sqrt[1 - mu x^2], substitution y = a + (y_q - a) sin^2 u; every value must be real *)
quadY[qv_?NumericQ, g_] := Module[{qe = Rationalize[qv, 10^-28], yq, bv, v},
   yq = Sqrt[1 - muv qe^2]; bv = Sqrt[nsq[yq, a0] (1 - yq^2)/muv];
   v = NIntegrate[Module[{yy = a0 + (yq - a0) Sin[u]^2, n2}, n2 = nsq[yy, a0];
       g[n2, bv, yy] (yq - a0) Sin[2 u]/Sqrt[n2 (1 - yy^2) - muv bv^2]], {u, 0, Pi/2},
      WorkingPrecision -> 36, PrecisionGoal -> 20, MaxRecursion -> 50];
   If[! NumberQ[v] || Abs[Im[v]] > 10^-25 Abs[v], tovAbort["a quadrature is not real"]];
   Re[v]];
delta[q_?NumericQ] := 2 quadY[q, (#2 Sqrt[muv]/(1 - #3^2)) &];
tau[q_?NumericQ] := 2 quadY[q, (#1/Sqrt[muv]) &];       (* proper time in units R/c *)
bOf[q_] := q Sqrt[nsq[Sqrt[1 - muv q^2], a0]];

illinois[f_, x1_, x2_, n_] := Module[{a = N[x1, 30], b = N[x2, 30], fa, fb, c, fc, side = 0},
   fa = f[a]; fb = f[b];
   If[fa fb > 0, tovAbort["regula falsi: the initial bracket has no sign change"]];
   Do[c = (a fb - b fa)/(fb - fa); fc = f[c];
     If[fc fb > 0, b = c; fb = fc; If[side == -1, fa /= 2]; side = -1,
                   a = c; fa = fc; If[side == 1, fb /= 2]; side = 1];
     If[Abs[b - a] < 10^-22, Break[]], {n}];
   c];

(* local extrema of Delta as zeros of a centred difference, in brackets read off a sampled profile *)
dd[q_?NumericQ] := (delta[q + 10^-8] - delta[q - 10^-8]) 10^8/2;
qMin = illinois[dd, 15/100, 26/100, 60]; qMax = illinois[dd, 30/100, 45/100, 60];
dMin = delta[qMin]; dMax = delta[qMax];
tovValue["mu", "8/9 - 1/10^4"];
tovValue["aperture_local_min.q", qMin, 10]; tovValue["aperture_local_min.Delta", dMin, 12];
tovValue["aperture_local_max.q", qMax, 10]; tovValue["aperture_local_max.Delta", dMax, 12];
tovRecord["maxwell.aperture_not_monotone", "numerical_check",
  "mu = 8/9 - 1e-4, proper time: the located local minimum of Delta lies at a smaller radius and below the located local maximum",
  qMin < qMax && dMin < dMax];

rootIn[th_] := illinois[delta[#] - th &, 5/100, qMin, 80];
rootOut[th_] := illinois[delta[#] - th &, qMax, 7/10, 80];
gap[th_?NumericQ] := tau[rootOut[th]] - tau[rootIn[th]];
thM = illinois[gap, dMin + (dMax - dMin)/4, dMax - (dMax - dMin)/4, 40];
qa = rootIn[thM]; qc = rootOut[thM]; qm = illinois[delta[#] - thM &, qMin, qMax, 80];
{ta, tc, tm} = tau /@ {qa, qc, qm};
tovValue["theta_M_rad", thM, 11]; tovValue["theta_M_deg", thM 180/Pi, 10];
tovValue["inner.q", qa, 10]; tovValue["inner.b", bOf[qa], 10]; tovValue["inner.tau_R_over_c", ta, 13];
tovValue["outer.q", qc, 10]; tovValue["outer.b", bOf[qc], 10]; tovValue["outer.tau_R_over_c", tc, 13];
tovValue["rising.q", qm, 10]; tovValue["rising.tau_R_over_c", tm, 13];
tovValue["aperture_residual.inner_outer_max", Max[Abs[delta[#] - thM] & /@ {qa, qc}], 3];
tovValue["aperture_residual.rising", Abs[delta[qm] - thM], 3];
tovValue["time_difference.outer_minus_inner", tc - ta, 3];
tovRecord["maxwell.equal_times", "numerical_check",
  "at theta_M the inner and outer extremals (Delta = theta_M within 1e-9) have equal proper time within 1e-10",
  Abs[tc - ta] < 10^-10 && Max[Abs[delta[#] - thM] & /@ {qa, qc}] < 10^-9];
Module[{res = {}, gaps = {}},
  Do[Module[{ri = rootIn[th], ro = rootOut[th]},
     AppendTo[res, Max[Abs[delta[ri] - th], Abs[delta[ro] - th]]]; AppendTo[gaps, tau[ro] - tau[ri]]],
   {th, {thM - 10^-6, thM + 10^-6}}];
  tovValue["time_difference.theta_M_minus_1e-6", gaps[[1]], 5]; tovValue["time_difference.theta_M_plus_1e-6", gaps[[2]], 5];
  tovValue["aperture_residual.bracket_max", Max[res], 3];
  tovRecord["maxwell.bracket", "numerical_check",
   "the outer-minus-inner proper-time difference is negative at theta_M - 1e-6 and positive at theta_M + 1e-6, with aperture residuals below 1e-9 at both angles",
   gaps[[1]] < 0 < gaps[[2]] && Max[res] < 10^-9]];
tovValue["rising.slope_dDelta_dq", dd[qm], 5];
tovRecord["maxwell.rising_branch_slower", "numerical_check",
  "the extremal on the rising branch (Delta = theta_M within 1e-9, positive sampled slope) is slower than the inner and outer ones by more than 1e-6",
  Abs[delta[qm] - thM] < 10^-9 && dd[qm] > 0 && tm - Max[ta, tc] > 10^-6];

Module[{grid, vals, cells},
  grid = Union[{1/1000, 2/1000, 3/1000}, Range[5/1000, 995/1000, 5/1000], {996/1000, 997/1000, 998/1000, 999/1000}];
  vals = (delta[#] - thM) & /@ grid;
  cells = Select[Range[Length[grid] - 1], vals[[#]] vals[[# + 1]] < 0 &];
  tovValue["scan.grid_points", Length[grid]];
  tovValue["scan.sign_change_cells", ToString[N[{grid[[#]], grid[[# + 1]]} & /@ cells, 4], InputForm]];
  tovRecord["maxwell.scan_three_crossings", "diagnostic",
   "on a grid of 205 radii, Delta - theta_M changes sign in exactly three cells, which contain the inner, rising and outer roots",
   Length[cells] == 3 && And @@ MapThread[grid[[#1]] < #2 < grid[[#1 + 1]] &, {cells, {qa, qm, qc}}]]];

tovFinish["maxwell_angle"];
