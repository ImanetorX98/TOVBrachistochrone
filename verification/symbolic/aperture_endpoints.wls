#!/usr/bin/env wolframscript
(* Aperture of the uniform-density brachistochrones at the ends of the regular set (Theorem A, part 2).

   Surface (q -> 1), eps = 1 - q, u = 1 - x = eps v, J(1 - u)^2 = (c/u) M(u), M(0) = 1:
     Delta/2 = eps Int_0^1 Sqrt[v/(1 - v)] Sqrt[M(eps)/B(eps, v)]/((1 - eps v) y(1 - eps v)) dv,  B = 1 + O(eps^2),
   so Delta = (Pi/a) eps [1 + kappa eps + O(eps^2)], kappa = M'(0)/2 + (3/4)(1 - mu/a^2) (3/4 is the mean of v for the
   weight Sqrt[v/(1 - v)]). Checked exactly: M'(0) from the series of J^2 and kappa_tau = -7 mu/(8 a^2),
   kappa_t = mu/(8 a^2). Checked numerically: the ratio and the slope at eps = 1e-4, 1e-5, 1e-6, and (ratio - 1)/eps -> kappa.

   Above mu_star (coordinate time, mu = 22/25), with F = J_t^2, x_m the local minimum of J_t (smaller root y1 in (a, 1) of
   the quartic Qt), y_m = y(x_m), q0 < x_M with J(q0) = J(x_m):
     Delta(q) = [2 J(x_m)/(x_m y_m Sqrt[F''(x_m)/2])] Log[1/d] + O(1)   on both sides (q = x_m + d, q = q0 - d).
   Checked numerically: the increase of Delta per decade of d tends to Log[10] times this coefficient.

   At mu_star, with q_star the double zero of J_t', B_star = F'''(q_star)/6, W_star = J(q_star)/(q_star y(q_star)):
     Delta(q) ~ (2 W_star/Sqrt[B_star]) |q - q_star|^(-1/2) I_side,
     I_right = Beta(1/6, 1/2)/3,  I_left = [Beta(1/3, 1/2) + Beta(1/3, 1/6)]/3.
   Checked numerically: Delta Sqrt[|q - q_star|] approaches the coefficient on both sides. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "common.wl"}]];
nsqY["proper"][z_, a_] := (3 a - z)^2/((z - a) (5 a - z));
nsqY["coordinate"][z_, a_] := 16/((3 a - z)^2 (z - a) (5 a - z));
(* every aperture must be a positive real number: an unevaluated or complex quadrature aborts the job *)
halfAp[m_, muv_, qv_] := Module[{a = Sqrt[1 - muv], yq, bv, v},
   yq = Sqrt[1 - muv qv^2]; bv = Sqrt[nsqY[m][yq, a] (1 - yq^2)/muv];
   v = NIntegrate[Module[{yy = a + (yq - a) Sin[u]^2},
      bv Sqrt[muv]/((1 - yy^2) Sqrt[nsqY[m][yy, a] (1 - yy^2) - muv bv^2]) (yq - a) Sin[2 u]],
     {u, 0, Pi/2}, WorkingPrecision -> 50, PrecisionGoal -> 30, MaxRecursion -> 60];
   If[! NumberQ[v], tovAbort["an aperture quadrature did not evaluate to a number"]];
   If[Abs[Im[v]] > 10^-20 Abs[v] || Re[v] <= 0, tovAbort["an aperture is not a positive real number"]];
   Re[v]];
strictlyDecreasing[l_List] := And @@ Thread[Rest[l] < Most[l]];

(* --- surface: exact second-order coefficient --- *)
Module[{u, mm, yU, fU, jsq, mSeries, kappa},
  yU = Sqrt[1 - mm (1 - u)^2];
  fU = ((3 Sqrt[1 - mm] - yU)/2)^2;
  jsq["proper"] = (1 - u)^2 fU/((1 - mm) - fU);
  jsq["coordinate"] = (1 - u)^2/(fU ((1 - mm) - fU));
  Do[
   mSeries = Normal[Series[u jsq[m], {u, 0, 1}]];
   kappa = Simplify[(Coefficient[mSeries, u, 1]/Coefficient[mSeries, u, 0])/2 + (3/4) (1 - mm/(1 - mm)),
     0 < mm < 8/9];
   tovRecord["endpoints.surface_coefficient." <> m, "exact",
    "kappa = M'(0)/2 + (3/4)(1 - mu/a^2) equals " <> If[m == "proper", "-7 mu/(8 a^2)", "mu/(8 a^2)"],
    Simplify[kappa - If[m == "proper", -7 mm/(8 (1 - mm)), mm/(8 (1 - mm))], 0 < mm < 8/9] === 0],
   {m, {"proper", "coordinate"}}]];

(* --- surface: numerical checks of the expansion --- *)
epsList = {10^-4, 10^-5, 10^-6};
Do[Module[{a = Sqrt[1 - muv], kappa, ratios, slopes, second},
   kappa = If[m == "proper", -7 muv/(8 a^2), muv/(8 a^2)];
   ratios = Table[2 halfAp[m, muv, 1 - eps]/(eps Pi/a), {eps, epsList}];
   slopes = Table[(2 halfAp[m, muv, 1 - eps] - 2 halfAp[m, muv, 1 - 2 eps])/eps/(Pi/a), {eps, epsList}];
   second = MapThread[(#1 - 1)/#2 &, {ratios, epsList}];
   tovValue["surface_ratio." <> m <> ".mu_" <> tovTag[muv] <> ".eps_1e-6", ratios[[3]], 15];
   tovRecord["endpoints.surface." <> m <> ".mu_" <> tovTag[muv], "numerical_check",
    "Delta/(Pi (1 - q)/a) -> 1 and -Delta'(q)/(Pi/a) -> 1 at 1 - q = 1e-4, 1e-5, 1e-6, errors decreasing at every step, below 1e-4 at 1e-6",
    Abs[ratios[[3]] - 1] < 10^-4 && Abs[slopes[[3]] + 1] < 10^-4 &&
     strictlyDecreasing[Abs[ratios - 1]] && strictlyDecreasing[Abs[slopes + 1]]];
   tovRecord["endpoints.surface_second_order." <> m <> ".mu_" <> tovTag[muv], "numerical_check",
    "(Delta/(Pi (1 - q)/a) - 1)/(1 - q) tends to kappa: relative error below 1e-4 at 1 - q = 1e-6, decreasing at every step",
    Abs[second[[3]]/kappa - 1] < 10^-4 && strictlyDecreasing[Abs[second/kappa - 1]]]],
  {m, {"proper", "coordinate"}}, {muv, {1/10, 1/2, 4/5, 22/25}}];

(* --- above mu_star: logarithmic divergences with the analytic coefficient --- *)
Module[{muv = 22/25, a, qt, roots, y1, y2, xm, xM, jt, fsq, q0, coeff, perDecade, dels, nearXm, nearQ0, dXm, dQ0},
  a = Sqrt[1 - muv];
  qt = y^4 - 3 a y^3 - (2 + 9 a^2) y^2 + (12 a + 15 a^3) y - 14 a^2;
  roots = Sort[Select[y /. Solve[qt == 0, y, Reals], a < # < 1 &]];
  If[Length[roots] != 2, tovAbort["expected two critical radii above mu_star"]];
  tovRecord["endpoints.two_critical_points", "exact",
   "at mu = 22/25 the quartic Qt has exactly two roots in (a, 1)", Length[roots] == 2];
  {y1, y2} = roots;
  xm = Sqrt[(1 - y1^2)/muv]; xM = Sqrt[(1 - y2^2)/muv];
  jt[x_] := x Sqrt[nsqY["coordinate"][Sqrt[1 - muv x^2], a]];
  fsq[x_] := jt[x]^2;
  (* q0 is rationalized to 1e-40 so that the quadratures below receive exact input *)
  q0 = Rationalize[x /. FindRoot[jt[x] == jt[xm], {x, xM/2, 10^-3, xM}, WorkingPrecision -> 60], 10^-40];
  tovValue["x_m", xm, 20]; tovValue["x_M", xM, 20]; tovValue["q0", q0, 20];
  tovRecord["endpoints.order_of_critical_points", "numerical_check",
   "0 < q0 < x_M < x_m < 1 and J_t(q0) = J_t(x_m) to 1e-30",
   0 < q0 < xM < xm < 1 && Abs[N[jt[q0] - jt[xm], 60]] < 10^-30];
  coeff = N[2 jt[xm]/(xm Sqrt[1 - muv xm^2] Sqrt[(D[fsq[x], {x, 2}] /. x -> xm)/2]), 40];
  perDecade = Log[10] coeff;
  tovValue["log_coefficient", coeff, 15];
  tovValue["log_coefficient_per_decade", perDecade, 15];
  dels = {10^-3, 10^-4, 10^-5, 10^-6};
  nearXm = Table[2 halfAp["coordinate", muv, xm (1 + del)], {del, dels}];
  nearQ0 = Table[2 halfAp["coordinate", muv, q0 (1 - del)], {del, dels}];
  dXm = Differences[nearXm]; dQ0 = Differences[nearQ0];
  MapThread[tovValue["decade_increment.x_m." <> #1, #2, 12] &, {{"3_4", "4_5", "5_6"}, dXm}];
  MapThread[tovValue["decade_increment.q0." <> #1, #2, 12] &, {{"3_4", "4_5", "5_6"}, dQ0}];
  tovRecord["endpoints.log_coefficient_at_x_m", "numerical_check",
   "the increase of Delta_t per decade of d (q = x_m (1 + d)) tends to Log[10] times the analytic coefficient: error below 1e-3 at the last decade and decreasing",
   Abs[Last[dXm] - perDecade] < 10^-3 && strictlyDecreasing[Abs[dXm - perDecade]]];
  tovRecord["endpoints.log_coefficient_at_q0", "numerical_check",
   "the increase of Delta_t per decade of d (q = q0 (1 - d)) tends to Log[10] times the same coefficient: error below 1e-3 at the last decade and decreasing",
   Abs[Last[dQ0] - perDecade] < 10^-3 && strictlyDecreasing[Abs[dQ0 - perDecade]]]];

(* --- at mu_star: power-law divergence with the Beta-function coefficients --- *)
Module[{aStarSq, muStar, a, qt, g, yc, qStar, jt, fsq, bStar, wStar, base, iRight, iLeft, sides},
  aStarSq = SelectFirst[Table[Root[3645 #^3 + 3888 #^2 - 1080 # + 64 &, k], {k, 3}], Element[#, Reals] && 1/9 < # < 1 &];
  muStar = 1 - aStarSq; a = Sqrt[aStarSq];
  qt = y^4 - 3 a y^3 - (2 + 9 a^2) y^2 + (12 a + 15 a^3) y - 14 a^2;
  g = PolynomialGCD[qt, D[qt, y], Extension -> Automatic];
  If[Exponent[g, y] != 1, tovAbort["no double critical radius at mu_star"]];
  yc = y /. First[Solve[g == 0, y]];
  qStar = Sqrt[(1 - yc^2)/muStar];
  jt[x_] := x Sqrt[nsqY["coordinate"][Sqrt[1 - muStar x^2], a]];
  fsq[x_] := jt[x]^2;
  bStar = N[(D[fsq[x], {x, 3}] /. x -> qStar)/6, 40];
  wStar = N[jt[qStar]/(qStar Sqrt[1 - muStar qStar^2]), 40];
  base = 2 wStar/Sqrt[bStar];
  iRight = N[Beta[1/6, 1/2]/3, 40]; iLeft = N[(Beta[1/3, 1/2] + Beta[1/3, 1/6])/3, 40];
  tovValue["q_star", qStar, 20];
  tovValue["power_coefficient.right", base iRight, 15];
  tovValue["power_coefficient.left", base iLeft, 15];
  (* exact parameters for the quadratures: mu_star and q_star rationalized to 1e-40 *)
  sides = Table[Module[{qs = Rationalize[N[qStar, 60], 10^-40], ms = Rationalize[N[muStar, 60], 10^-40], vals},
      vals = Table[2 halfAp["coordinate", ms, qs (1 + sgn del)] Sqrt[qs del], {del, {10^-5, 10^-6, 10^-7}}];
      {sgn, vals, If[sgn > 0, base iRight, base iLeft]}], {sgn, {1, -1}}];
  Do[Module[{sgn = side[[1]], vals = side[[2]], target = side[[3]]},
     tovValue["power_scaled.decade_7." <> If[sgn > 0, "right", "left"], Last[vals], 12];
     tovRecord["endpoints.power_law." <> If[sgn > 0, "right", "left"], "numerical_check",
      "at mu_star, Delta_t Sqrt[|q - q_star|] tends to the Beta-function coefficient on the " <>
       If[sgn > 0, "outer", "inner"] <> " side: relative error below 1e-2 at 1e-7 and decreasing",
      Abs[Last[vals]/target - 1] < 10^-2 && strictlyDecreasing[Abs[vals/target - 1]]]],
   {side, sides}]];

tovFinish["aperture_endpoints"];
