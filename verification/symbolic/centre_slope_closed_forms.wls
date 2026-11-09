#!/usr/bin/env wolframscript
(* Closed forms of the centre slopes C_t and C_tau (uniform density), Delta/2 = Pi/2 + C b + O(b^3) as q -> 0.

   With a = Sqrt[1 - mu], Q4 = (y - a)(5a - y)(1 - y^2) (roots -1 < a < 1 < 5a):
     C_t   = -(Sqrt[mu]/2) Int_a^1 y (y^2 - 6 a y + 7 a^2) dy/Sqrt[Q4],
     C_tau = -4 a^2 Sqrt[mu] Int_a^1 y dy/((3a - y)^2 Sqrt[Q4]).
   Derivation (all steps exact):
     - Hermite: y (y^2 - 6ay + 7a^2) = Q4'/4 - 3a/2 + (1 + 9a^2) y/2 - 3a y^2/2, and [Sqrt[Q4]/2]_a^1 = 0;
     - Legendre: sn^2 u = 2 (y - a)/((1 - a)(y + 1)), y = -1 + (1 + a)/(1 - n sn^2), dy/Sqrt[Q4] = du/Sqrt[2a],
       m = k^2 = (1 - a)(1 + 5a)/(8a), n = (1 - a)/2 (0 < n < m < 1: hyperbolic case in DLMF 19.2(ii));
     - 3a - y = 2a (1 - N sn^2)/(1 - n sn^2), N = (1 + 3a)(1 - a)/(4a) (m < N < 1: circular case in DLMF);
     - W_2(nn) = Int_0^K du/(1 - nn sn^2)^2 from the exact identity for d/du[sn cn dn/(1 - nn sn^2)].
   Results (parameter convention m = k^2, Pi(nn|m) with denominator 1 - nn Sin[theta]^2):
     C_t   = (Sqrt[2 mu/a]/8) [ (1 + 3a + 6a^2) K(m) - 12 a^2 E(m) - (1 + a) Pi(n|m) ],
     C_tau = Sqrt[2 a mu]/((3a - 1)(1 + 3a)^2) [ a (1 + 9a) K(m) - 6a (1 + 3a) E(m) - (1 + a) Pi(N|m) ].
   The script derives both and checks that they equal these displayed formulas exactly. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "common.wl"}]];
assume = 1/3 < a < 1;
q4 = (y - a) (5 a - y) (1 - y^2);
m = (1 - a) (1 + 5 a)/(8 a);
n = (1 - a)/2;
bigN = (1 + 3 a) (1 - a)/(4 a);
yOfS = -1 + (1 + a)/(1 - n s);
ellK = EllipticK[m]; ellE = EllipticE[m];

(* Legendre substitution *)
tovRecord["centre_slope.legendre_map", "exact", "s = 2(y - a)/((1 - a)(y + 1)) inverts to y = -1 + (1 + a)/(1 - n s)",
  Simplify[(2 (yOfS - a)/((1 - a) (yOfS + 1))) - s] === 0];
tovRecord["centre_slope.legendre_jacobian", "exact", "(dy/du)^2 = Q4(y)/(2a) with s = sn^2(u|m)",
  Simplify[(D[yOfS, s] 2)^2 s (1 - s) (1 - m s) - (q4 /. y -> yOfS)/(2 a), assume && 0 < s < 1] === 0];
tovRecord["centre_slope.hyperbolic_case", "exact", "0 < n < m < 1 on 1/3 < a < 1",
  Reduce[assume && ! (0 < n < m < 1), a, Reals] === False];
tovRecord["centre_slope.circular_case", "exact", "m < N < 1 on 1/3 < a < 1, N -> 1 at a = 1/3",
  Reduce[assume && ! (m < bigN < 1), a, Reals] === False && (bigN /. a -> 1/3) == 1];
tovRecord["centre_slope.pole_factor", "exact", "3a - y = 2a (1 - N s)/(1 - n s)",
  Simplify[(3 a - yOfS) - 2 a (1 - bigN s)/(1 - n s)] === 0];

(* W_2 identity for a generic characteristic nn *)
dNum = (1 - s) (1 - m s) - s (1 - m s) - m s (1 - s);
dQuot = (dNum (1 - nn s) + 2 nn s (1 - s) (1 - m s))/(1 - nn s)^2;
r2 = Simplify[Limit[dQuot (1 - nn s)^2, s -> 1/nn]];
r1 = Simplify[Limit[D[dQuot (1 - nn s)^2, s], s -> 1/nn]/(-nn)];
rest = Simplify[dQuot - r1/(1 - nn s) - r2/(1 - nn s)^2];
{p0, p1} = {Coefficient[Expand[rest], s, 0], Coefficient[Expand[rest], s, 1]};
tovRecord["centre_slope.w2_identity", "exact",
  "d/du[sn cn dn/(1 - nn sn^2)] = p0 + p1 sn^2 + r1/(1 - nn sn^2) + r2/(1 - nn sn^2)^2",
  Simplify[dQuot - (p0 + p1 s + r1/(1 - nn s) + r2/(1 - nn s)^2)] === 0 && FreeQ[{p0, p1}, s]];
w2[nv_] := Simplify[-(p0 ellK + p1 (ellK - ellE)/m + r1 EllipticPi[nn, m])/r2 /. nn -> nv];

(* C_t: Hermite reduction and moments *)
numT = y (y^2 - 6 a y + 7 a^2);
reducedT = -3 a/2 + (1 + 9 a^2) y/2 - 3 a y^2/2;
tovRecord["centre_slope.hermite_coordinate", "exact",
  "y (y^2 - 6ay + 7a^2) = Q4'/4 + reduced quadratic, and Sqrt[Q4] vanishes at y = a, 1",
  Expand[numT - (D[q4, y]/4 + reducedT)] === 0 && (q4 /. y -> a) === 0 && Expand[q4 /. y -> 1] === 0];
c0 = -1; c1 = 1 + a;
mom0 = ellK; mom1 = c0 ellK + c1 EllipticPi[n, m];
mom2 = c0^2 ellK + 2 c0 c1 EllipticPi[n, m] + c1^2 w2[n];
ctDerived = -(Sqrt[1 - a^2]/2) (1/Sqrt[2 a]) ({Coefficient[reducedT, y, 0], Coefficient[reducedT, y, 1],
      Coefficient[reducedT, y, 2]}.{mom0, mom1, mom2});
ctDisplayed = (Sqrt[2 (1 - a^2)/a]/8) ((1 + 3 a + 6 a^2) ellK - 12 a^2 ellE - (1 + a) EllipticPi[n, m]);
tovRecord["centre_slope.coordinate_formula", "exact", "the derived C_t equals the displayed K, E, Pi formula",
  Simplify[ctDerived - ctDisplayed, assume] === 0];

(* C_tau: partial fractions in s and moments *)
rat = Together[yOfS/(3 a - yOfS)^2];
{d0, d1, d2} = {-1/(1 + 3 a)^2, (3 a - 1) (1 + a)/(2 a (1 + 3 a)^2), 3 (1 + a)^2/(4 a (1 + 3 a)^2)};
tovRecord["centre_slope.partial_fractions_proper", "exact",
  "y/(3a - y)^2 = d0 + d1/(1 - N s) + d2/(1 - N s)^2 with the displayed coefficients",
  Simplify[rat - (d0 + d1/(1 - bigN s) + d2/(1 - bigN s)^2)] === 0];
ctauDerived = -4 a^2 Sqrt[1 - a^2] (1/Sqrt[2 a]) (d0 ellK + d1 EllipticPi[bigN, m] + d2 w2[bigN]);
ctauDisplayed = Sqrt[2 a (1 - a^2)]/((3 a - 1) (1 + 3 a)^2) (a (1 + 9 a) ellK - 6 a (1 + 3 a) ellE -
     (1 + a) EllipticPi[bigN, m]);
tovRecord["centre_slope.proper_formula", "exact", "the derived C_tau equals the displayed K, E, Pi formula",
  Simplify[ctauDerived - ctauDisplayed, assume] === 0];
tovRecord["centre_slope.proper_bracket_bound", "exact",
  "a(1 + 9a) - 6a(1 + 3a)(1 - m) = -(3a - 1)(15a^2 + 10a + 3)/4 < 0 (sign of the C_tau bracket)",
  Simplify[a (1 + 9 a) - 6 a (1 + 3 a) (1 - m) + (3 a - 1) (15 a^2 + 10 a + 3)/4] === 0];

(* numerical checks of the displayed formulas against direct quadrature in y *)
quadT[av_, muv_] := -(Sqrt[muv]/2) NIntegrate[Module[{yy = av + (1 - av) Sin[t]^2},
     (numT /. {a -> av, y -> yy})/Sqrt[q4 /. {a -> av, y -> yy}] 2 (1 - av) Sin[t] Cos[t]],
    {t, 0, Pi/2}, WorkingPrecision -> 40, PrecisionGoal -> 25, MaxRecursion -> 40];
quadTau[av_, muv_] := -4 av^2 Sqrt[muv] NIntegrate[Module[{yy = av + (1 - av) Sin[t]^2},
     yy/((3 av - yy)^2 Sqrt[q4 /. {a -> av, y -> yy}]) 2 (1 - av) Sin[t] Cos[t]],
    {t, 0, Pi/2}, WorkingPrecision -> 40, PrecisionGoal -> 25, MaxRecursion -> 40];
Do[Module[{av = Sqrt[1 - muv], ct, ctau, tag = tovTag[muv]},
   ct = N[ctDisplayed /. a -> av, 30]; ctau = N[ctauDisplayed /. a -> av, 30];
   tovValue["C_t.mu_" <> tag, ct]; tovValue["C_tau.mu_" <> tag, ctau];
   tovRecord["centre_slope.coordinate_vs_quadrature.mu_" <> tag, "numerical_check",
    "displayed C_t equals the direct quadrature to 1e-20 relative", Abs[ct - quadT[av, muv]] < 10^-20 Abs[ct]];
   tovRecord["centre_slope.proper_vs_quadrature.mu_" <> tag, "numerical_check",
    "displayed C_tau equals the direct quadrature to 1e-20 relative",
    Abs[ctau - quadTau[av, muv]] < 10^-20 Abs[ctau]]],
  {muv, {1/10, 1/2, 4/5, 87/100}}];

(* Buchdahl asymptotics of C_tau: C_tau = -(2 Pi/9) delta^(-3/2) (1 + delta/8 + O(delta^(3/2))), delta = 3a - 1 *)
Module[{ratios = Table[Module[{av = (1 + del)/3, ctau},
       ctau = N[ctauDisplayed /. a -> av, 40];
       (ctau/(-(2 Pi/9) del^(-3/2) (1 + del/8)) - 1)/del^(3/2)], {del, {10^-4, 10^-6, 10^-8}}]},
  tovValue["buchdahl_remainder_over_delta_3_2", Last[ratios], 8];
  tovRecord["centre_slope.proper_buchdahl_asymptotics", "numerical_check",
   "(C_tau/(-(2Pi/9) delta^(-3/2)(1 + delta/8)) - 1)/delta^(3/2) stays bounded for delta = 1e-4, 1e-6, 1e-8",
   Max[Abs[ratios]] < 10 && Abs[ratios[[3]] - ratios[[2]]] < Abs[ratios[[2]]]]];

(* negative control: drop the double-pole term of C_tau *)
Module[{bad = -4 a^2 Sqrt[1 - a^2] (1/Sqrt[2 a]) (d0 ellK + d1 EllipticPi[bigN, m]), av = Sqrt[9/10]},
  tovRecord["centre_slope.negative_control.drop_double_pole", "negative_control",
   "dropping the double-pole term of C_tau changes the value by more than 1e-6 relative",
   Abs[N[(bad - ctauDisplayed) /. a -> av, 30]] > 10^-6 Abs[N[ctauDisplayed /. a -> av, 30]]]];

tovFinish["centre_slope_closed_forms"];
