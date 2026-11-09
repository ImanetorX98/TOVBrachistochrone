#!/usr/bin/env wolframscript
(* Uniform-density (Schwarzschild interior) brachistochrones: reduction to a genus-two curve.

   Units R = c = G = 1, mu = 2M/R < 8/9, x = r/R, y = Sqrt[1 - mu x^2], yR = Sqrt[1 - mu].
   Lapse squared f = ((3 yR - y)/2)^2, radial metric h = 1/y^2, energy E^2 = f(1) = yR^2
   (release from rest at the surface). Indices (as in tov_brachistochrone.py):
     proper time      n^2 = f/(E^2 - f),
     coordinate time  n^2 = 1/(f (E^2 - f)).
   Clairaut constant b = n(xs) xs at the turning radius xs = q; dphi/dx = b Sqrt[h]/(x Sqrt[n^2 x^2 - b^2]).

   Claims checked:
     - E^2 - f = (y - yR)(5 yR - y)/4;
     - dphi/dy = -b Sqrt[mu] N(y)/((1 - y^2) w), w^2 = S(y) a sextic, square-free for generic parameters;
     - P(yq) = 0 at the Clairaut constant (the turning point is a known root);
     - N/(1 - y^2) = c.(1, y) + aP/(y - 1) + aM/(y + 1) with exact coefficients (first plus third kind);
     - the half-aperture by the reduced y-integral equals the x-quadrature (numerical check);
     - no Moebius involution permutes the six roots on the samples (diagnostic). *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "common.wl"}]];

f[y_] := ((3 yR - y)/2)^2;
e2 = yR^2;
nsq["proper"][y_] := f[y]/(e2 - f[y]);
nsq["coordinate"][y_] := 1/(f[y] (e2 - f[y]));
xOfY[y_] := Sqrt[(1 - y^2)/mu];
dPoly = (y - yR) (5 yR - y);
pPoly["proper"] = (3 yR - y)^2 (1 - y^2) - mu b^2 dPoly;
pPoly["coordinate"] = 16 (1 - y^2) - mu b^2 (3 yR - y)^2 dPoly;
numer["proper"] = dPoly;
numer["coordinate"] = (3 yR - y) dPoly;
sextic[m_] := Expand[dPoly pPoly[m]];
decomposition["proper"] = {{1, 0}, (yR - 1) (5 yR - 1)/2, -(1 + yR) (1 + 5 yR)/2};
decomposition["coordinate"] = {{9 yR, -1}, (yR - 1) (3 yR - 1) (5 yR - 1)/2,
   -(1 + yR) (1 + 3 yR) (1 + 5 yR)/2};
modes = {"proper", "coordinate"};

tovRecord["reduction.energy_gap", "exact", "E^2 - f = (y - yR)(5 yR - y)/4",
  Simplify[e2 - f[y] - dPoly/4] === 0];

Do[
  Module[{lhs, rhs},
   lhs = b^2 (1/y^2)/(xOfY[y]^2 (nsq[m][y] xOfY[y]^2 - b^2)) (D[xOfY[z], z] /. z -> y)^2;
   rhs = b^2 mu numer[m]^2/((1 - y^2)^2 sextic[m]);
   tovRecord["reduction.dphi_dy." <> m, "exact",
    "(dphi/dy)^2 = b^2 mu N^2/((1 - y^2)^2 S), " <> m <> " time",
    Simplify[Together[lhs - rhs], Assumptions -> 0 < mu < 8/9 && yR > 0 && y > 0] === 0]],
  {m, modes}];

Do[tovRecord["reduction.sextic_squarefree." <> m, "exact",
   "S has degree 6 and gcd(S, S') = 1 over the parameter field, " <> m <> " time",
   Exponent[sextic[m], y] == 6 && Exponent[PolynomialGCD[sextic[m], D[sextic[m], y]], y] == 0],
  {m, modes}];

Do[tovRecord["reduction.turning_root." <> m, "exact",
   "P(yq) = 0 when b = n(yq) x(yq), " <> m <> " time",
   Simplify[(pPoly[m] /. {b^2 -> nsq[m][ys] xOfY[ys]^2, y -> ys})] === 0],
  {m, modes}];

Do[tovRecord["reduction.decomposition." <> m, "exact",
   "N/(1 - y^2) = c.(1, y) + aP/(y - 1) + aM/(y + 1), " <> m <> " time",
   Simplify[decomposition[m][[1]].{1, y} + decomposition[m][[2]]/(y - 1) +
      decomposition[m][[3]]/(y + 1) - numer[m]/(1 - y^2)] === 0],
  {m, modes}];
tovRecord["reduction.buchdahl_factor", "exact",
  "the coefficient aP of the coordinate-time curve carries the factor (3 yR - 1)",
  Simplify[decomposition["coordinate"][[2]]/(3 yR - 1)] === (yR - 1) (5 yR - 1)/2];

tovRecord["reduction.buchdahl_branch_point", "exact",
  "at yR = 1/3: P_t(1) = 0 and P_t'(1) = -32 for every b, so the two points over y = 1 merge into a simple branch point",
  Simplify[(pPoly["coordinate"] /. {yR -> 1/3, y -> 1})] === 0 &&
   Simplify[(D[pPoly["coordinate"], y] /. {yR -> 1/3, y -> 1, mu -> 8/9})] === -32];

(* numerical check: half-aperture in x (metric written directly) against the reduced y-integral *)
halfApertureX[m_, muv_, qv_] := Module[{a = Sqrt[1 - muv], fx, nx, bv},
   fx[x_] := ((3 a - Sqrt[1 - muv x^2])/2)^2;
   nx[x_] := If[m == "proper", Sqrt[fx[x]/(a^2 - fx[x])], 1/Sqrt[fx[x] (a^2 - fx[x])]];
   bv = nx[qv] qv;
   NIntegrate[Module[{x = qv + (1 - qv) Sin[u]^2},
     bv/Sqrt[1 - muv x^2]/(x Sqrt[nx[x]^2 x^2 - bv^2]) (1 - qv) Sin[2 u]],
    {u, 0, Pi/2}, WorkingPrecision -> 40, PrecisionGoal -> 25, MaxRecursion -> 30]];
halfApertureY[m_, muv_, qv_] := Module[{a = Sqrt[1 - muv], ys, bv, pp, nn},
   ys = Sqrt[1 - muv qv^2];
   bv = Sqrt[nsq[m][ys] xOfY[ys]^2 /. {yR -> a, mu -> muv}];
   pp = sextic[m] /. {yR -> a, mu -> muv, b -> bv};
   nn = numer[m] /. yR -> a;
   NIntegrate[Module[{yy = a + (ys - a) Sin[u]^2},
     bv Sqrt[muv] (nn /. y -> yy)/((1 - yy^2) Sqrt[pp /. y -> yy]) (ys - a) Sin[2 u]],
    {u, 0, Pi/2}, WorkingPrecision -> 40, PrecisionGoal -> 25, MaxRecursion -> 30]];
samples = {{1/10, 1/2}, {3/10, 3/10}, {1/2, 7/10}, {4/5, 1/5}};
Do[Module[{vx = halfApertureX[m, smp[[1]], smp[[2]]], vy = halfApertureY[m, smp[[1]], smp[[2]]],
    tag = m <> ".mu_" <> tovTag[smp[[1]]] <> ".q_" <> tovTag[smp[[2]]]},
   tovValue["half_aperture." <> tag, vy];
   tovRecord["reduction.half_aperture_routes." <> tag, "numerical_check",
    "half-aperture: reduced y-integral equals the x-quadrature to 1e-20 relative",
    Abs[vx - vy] < 10^-20 Abs[vx]]],
  {m, modes}, {smp, samples}];

(* diagnostic: no Moebius involution permutes the six roots (three pairs with Det[{1, a+b, ab}] = 0) *)
pairings[l_] := If[Length[l] == 0, {{}},
   Flatten[Table[Prepend[#, {l[[1]], l[[k]]}] & /@ pairings[Delete[l, {{1}, {k}}]], {k, 2, Length[l]}], 1]];
minInvolutionDet[m_, muv_, qv_] := Module[{a = Sqrt[1 - muv], ys, bv, rts},
   ys = Sqrt[1 - muv qv^2];
   bv = Sqrt[nsq[m][ys] xOfY[ys]^2 /. {yR -> a, mu -> muv}];
   rts = y /. NSolve[(sextic[m] /. {yR -> a, mu -> muv, b -> bv}) == 0, y, WorkingPrecision -> 40];
   Min[Abs[Det[{{1, #[[1, 1]] + #[[1, 2]], #[[1, 1]] #[[1, 2]]},
          {1, #[[2, 1]] + #[[2, 2]], #[[2, 1]] #[[2, 2]]},
          {1, #[[3, 1]] + #[[3, 2]], #[[3, 1]] #[[3, 2]]}}]] & /@ pairings[rts]]];
Module[{dets = Flatten[Table[minInvolutionDet[m, smp[[1]], smp[[2]]], {m, modes}, {smp, samples}]]},
  tovValue["min_involution_determinant", Min[dets], 6];
  tovRecord["reduction.no_obvious_elliptic_splitting", "diagnostic",
   "no Moebius involution permutes the six roots on the 8 samples (min determinant > 1e-3)",
   Min[dets] > 10^-3]];

tovFinish["uniform_density_reduction"];
