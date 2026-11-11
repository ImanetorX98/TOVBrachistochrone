#!/usr/bin/env wolframscript
(* An anisotropic interior beyond the Buchdahl bound: the constant-density interior with vanishing radial pressure
   (the Florides solution, in the form of Boehmer and Lobo, arXiv:gr-qc/0703024, eqs. (10)-(11)).

   Units R = c = G = 1, 0 < mu < 1 (no Buchdahl bound), x = r/R, y = Sqrt[1 - mu x^2], a = Sqrt[1 - mu].
   Lapse squared f = a^3/y, radial metric h = 1/y^2, energy E^2 = f(1) = a^2 (release from rest at the surface).
   Indices: proper time n^2 = f/(E^2 - f) = a/(y - a); coordinate time n^2 = 1/(f (E^2 - f)) = y^2/(a^5 (y - a)).

   Claims checked:
     - Einstein equations: rho = 3 mu/(8 Pi), p_r = 0, p_T/rho = mu x^2/(4 (1 - mu x^2)); matching data at x = 1;
       the dominant energy condition p_T <= rho holds on the whole sphere iff mu <= 4/5;
     - indices; (dphi/dy)^2 = b^2 mu N^2/((1 - y^2)^2 w^2) with a cubic w^2 for tau (square-free) and a quintic for t;
     - J_tau increasing for every mu < 1; J_t increasing iff mu <= 125/128 (discriminant 4(a^2 - 1)(128 a^2 - 3));
     - curvature: L_tau,y < 0 on (a, 1) for every mu < 1 (numerator (y - a)^2 - (1 - a^2) < -2a(1 - a));
       L_t,y <= 0 iff a >= 2/3 (mu <= 5/9), by convexity in a and the end values N(y, y), N(2/3, y);
     - central slopes: C_tau < 0 for every mu; C_t = -(Sqrt[mu] a^(5/2)/(2 Sqrt[1 + a])) M(a) with M strictly increasing
       (exact derivative identity, positive term by term), M(0+) < 0 and M(1/2) > 0, so C_t has exactly one zero; its
       location mu_Delta^F is bracketed by bisection (numerical localisation);
     - negative control: the isotropic uniform-density lapse ((3a - y)/2)^2 fails the tau curvature inequality near the
       Buchdahl bound, so the inequality is not automatic;
     - diagnostic: Delta_tau decreases on a grid of radii at mu = 19/20. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "common.wl"}]];

dom = 0 < a < 1 && a < y < 1;
fY = a^3/y; e2 = a^2;
ntau2 = a/(y - a); nt2 = y^2/(a^5 (y - a));

(* --- Einstein equations for -f dt^2 + h dr^2 + r^2 dOmega^2 with f = a^3/y(r), h = 1/y(r)^2 --- *)
Module[{yr, fr, hr, crd, g, gi, gam, ric, rs, gm, rho, pr, pt, th, ph, t, r, mu},
  mu = 1 - a^2;
  yr = Sqrt[1 - mu r^2]; fr = a^3/yr; hr = 1/yr^2;
  crd = {t, r, th, ph}; g = DiagonalMatrix[{-fr, hr, r^2, r^2 Sin[th]^2}]; gi = Inverse[g];
  gam = Table[Sum[gi[[i, l]] (D[g[[l, j]], crd[[k]]] + D[g[[l, k]], crd[[j]]] - D[g[[j, k]], crd[[l]]])/2, {l, 4}],
    {i, 4}, {j, 4}, {k, 4}];
  ric = Table[Sum[D[gam[[i, j, k]], crd[[i]]], {i, 4}] - Sum[D[gam[[i, j, i]], crd[[k]]], {i, 4}] +
     Sum[gam[[i, i, l]] gam[[l, j, k]] - gam[[i, k, l]] gam[[l, j, i]], {i, 4}, {l, 4}], {j, 4}, {k, 4}];
  rs = Simplify[Sum[gi[[i, j]] ric[[i, j]], {i, 4}, {j, 4}]];
  gm = Simplify[gi.ric - rs/2 IdentityMatrix[4]];
  rho = Simplify[-gm[[1, 1]]/(8 Pi)]; pr = Simplify[gm[[2, 2]]/(8 Pi)]; pt = Simplify[gm[[3, 3]]/(8 Pi)];
  tovRecord["florides.einstein", "exact",
   "rho = 3 mu/(8 Pi), p_r = 0, p_T/rho = mu x^2/(4 (1 - mu x^2)) for f = a^3/y, h = 1/y^2",
   Simplify[rho - 3 mu/(8 Pi)] === 0 && Simplify[pr] === 0 &&
    Simplify[pt/rho - mu r^2/(4 (1 - mu r^2)), 0 < a < 1 && 0 < r < 1] === 0];
  tovRecord["florides.matching", "exact",
   "at x = 1: f = a^2, h = 1/a^2 and f' = mu, the values of the exterior Schwarzschild metric",
   Simplify[(fr /. r -> 1) - a^2, 0 < a < 1] === 0 && Simplify[(hr /. r -> 1) - 1/a^2, 0 < a < 1] === 0 &&
    Simplify[(D[fr, r] /. r -> 1) - mu, 0 < a < 1] === 0];
  tovRecord["florides.dominant_energy", "exact",
   "p_T <= rho on 0 <= x <= 1 iff mu <= 4/5",
   Reduce[0 < m < 1 && ForAll[xx, 0 <= xx <= 1, m xx^2/(4 (1 - m xx^2)) <= 1], m, Reals] ===
    Inequality[0, Less, m, LessEqual, 4/5]]];

tovRecord["florides.indices", "exact", "n_tau^2 = a/(y - a) and n_t^2 = y^2/(a^5 (y - a)) with E^2 = a^2",
  Simplify[fY/(e2 - fY) - ntau2, dom] === 0 && Simplify[1/(fY (e2 - fY)) - nt2, dom] === 0];

(* --- reduction: (dphi/dy)^2 = b^2 h (dx/dy)^2/(x^2 (J^2 - b^2)), x^2 = (1 - y^2)/mu --- *)
Module[{mu = 1 - a^2, x2, lhs, pTau, pT, sTau, sT},
  x2 = (1 - y^2)/mu;
  lhs[n2_] := b^2 (1/y^2) (y^2/(mu^2 x2))/(x2 (x2 n2 - b^2));
  pTau = a (1 - y^2) - mu b^2 (y - a); pT = y^2 (1 - y^2) - mu a^5 b^2 (y - a);
  sTau = Expand[(y - a) pTau]; sT = Expand[(y - a) pT];
  tovRecord["florides.reduction.proper", "exact",
   "tau: N = y - a, w^2 = (y - a)(a (1 - y^2) - mu b^2 (y - a)), a cubic, square-free (distinct zeros of the quadratic factor, which is a(1 - a^2) > 0 at y = a)",
   Simplify[lhs[ntau2] - b^2 mu (y - a)^2/((1 - y^2)^2 sTau)] === 0 && Exponent[sTau, y] == 3 &&
    Reduce[0 < a < 1 && b > 0 && Discriminant[pTau, y] <= 0, {a, b}, Reals] === False &&
    Simplify[(pTau /. y -> a) - a (1 - a^2)] === 0];
  tovRecord["florides.reduction.coordinate", "exact",
   "t: N = a^(5/2)(y - a), w^2 = (y - a)(y^2 (1 - y^2) - mu a^5 b^2 (y - a)), a quintic, square-free at the witness (a, b) = (1/2, 1)",
   Simplify[lhs[nt2] - b^2 mu a^5 (y - a)^2/((1 - y^2)^2 sT)] === 0 && Exponent[sT, y] == 5 &&
    Exponent[PolynomialGCD[sT /. {a -> 1/2, b -> 1}, D[sT /. {a -> 1/2, b -> 1}, y]], y] == 0]];

(* --- Clairaut functions; J increasing in x iff d(J^2)/dy < 0 since dy/dx < 0 --- *)
Module[{mu = 1 - a^2, jTau2, jT2, pp},
  jTau2 = (1 - y^2)/mu ntau2; jT2 = (1 - y^2)/mu nt2;
  pp = 3 y^3 - 4 a y^2 - y + 2 a;
  tovRecord["florides.clairaut.proper", "exact", "d(J_tau^2)/dy < 0 on (a, 1) for every 0 < a < 1",
   Reduce[dom && Numerator[Together[D[jTau2, y]]] (Denominator[Together[D[jTau2, y]]]) >= 0, {a, y}, Reals] === False];
  tovRecord["florides.clairaut.coordinate", "exact",
   "d(J_t^2)/dy = -y p/(mu a^5 (y - a)^2), p = 3y^3 - 4ay^2 - y + 2a, Disc_y p = 4(a^2 - 1)(128 a^2 - 3), p(a), p(1) > 0; at a^2 = 3/128, p = 3 (y - Sqrt[6]/6)^2 (y + Sqrt[6]/4); p > 0 on (a, 1) for a^2 > 3/128 (decided on the polynomial)",
   Simplify[D[jT2, y] + y pp/(mu a^5 (y - a)^2)] === 0 &&
    Factor[Discriminant[pp, y] - 4 (a^2 - 1) (128 a^2 - 3)] === 0 &&
    Simplify[(pp /. a -> Sqrt[6]/16) - 3 (y - Sqrt[6]/6)^2 (y + Sqrt[6]/4)] === 0 &&
    Reduce[Sqrt[3/128] < a < 1 && a < y < 1 && pp <= 0, {a, y}, Reals] === False &&
    Reduce[0 < a < Sqrt[3/128] && ForAll[y, a < y < 1, pp > 0], a, Reals] === False]];

(* --- curvature criterion: L = y (1 + x n'/n), x n'/n = -(1 - y^2)/y d(log n)/dy; decided on numerators --- *)
Module[{lOf, lTau, lT, nTau, nT},
  lOf[n2_] := Simplify[y (1 - (1 - y^2)/y D[Log[n2]/2, y]), dom];
  lTau = lOf[ntau2]; lT = lOf[nt2];
  nTau = (y - a)^2 - (1 - a^2); nT = 3 y^4 - 6 a y^3 + (1 + 4 a^2) y^2 - 4 a y + 2 a^2;
  tovRecord["florides.curvature.proper", "exact",
   "L_tau,y = ((y - a)^2 - (1 - a^2))/(2 (y - a)^2) < 0 on (a, 1) for every 0 < a < 1",
   Simplify[D[lTau, y] - nTau/(2 (y - a)^2), dom] === 0 &&
    Reduce[dom && nTau >= 0, {a, y}, Reals] === False];
  tovRecord["florides.curvature.coordinate", "exact",
   "L_t,y = N/(2 y^2 (y - a)^2); N(a, 1) = 2(3a - 2)(a - 1) > 0 for a < 2/3; N convex in a, N(y, y) = y^2 (y^2 - 1) < 0 and N(2/3, y) = (y - 1)(27y^3 - 9y^2 + 16y - 8)/9 < 0 on (2/3, 1)",
   Simplify[D[lT, y] - nT/(2 y^2 (y - a)^2), dom] === 0 &&
    Factor[(nT /. y -> 1) - 2 (3 a - 2) (a - 1)] === 0 && Coefficient[nT, a, 2] === 2 + 4 y^2 &&
    Expand[(nT /. a -> y) - y^2 (y^2 - 1)] === 0 &&
    Expand[(nT /. a -> 2/3) - (y - 1) (27 y^3 - 9 y^2 + 16 y - 8)/9] === 0 &&
    Reduce[2/3 < y < 1 && 27 y^3 - 9 y^2 + 16 y - 8 <= 0, y, Reals] === False]];

(* --- central slopes: C = -Sqrt[mu] Int_a^1 y (1/n)'(y) dy/Sqrt[1 - y^2] --- *)
Module[{mu = 1 - a^2, dTau, dT, yS, pP, wW, integrand, mOf, ctOf, lo, hi, mid},
  dTau = D[Sqrt[(y - a)/a], y]; dT = D[a^(5/2) Sqrt[y - a]/y, y];
  tovRecord["florides.centre.proper", "exact", "(1/n_tau)' = 1/(2 Sqrt[a (y - a)]) > 0, so C_tau < 0 for every mu",
   Simplify[dTau - 1/(2 Sqrt[a (y - a)]), dom] === 0];
  yS = a + (1 - a) s; pP = 2 a/yS - 1; wW = Sqrt[(1 + a)/(1 + yS)];
  integrand = (y dT/Sqrt[1 - y^2]) /. y -> yS;
  tovRecord["florides.centre.coordinate_form", "exact",
   "with y = a + (1 - a) s: the C_t integrand is (a^(5/2)/(2 Sqrt[1 + a])) (2a/Y - 1) Sqrt[(1 + a)/(1 + Y)]/Sqrt[s (1 - s)]",
   Simplify[integrand (1 - a) - a^(5/2)/(2 Sqrt[1 + a]) pP wW/Sqrt[s (1 - s)], 0 < a < 1 && 0 < s < 1] === 0];
  tovRecord["florides.centre.coordinate_monotone", "exact",
   "d_a(P W)/W = 2s/Y^2 + P s/((1 + a)(1 + Y)) > 0 on 0 < a < 1, 0 < s < 1, so M is strictly increasing; P >= 0 at a = 1/2",
   Simplify[D[pP wW, a]/wW - (2 s/yS^2 + pP s/((1 + a) (1 + yS))), 0 < a < 1 && 0 < s < 1] === 0 &&
    Reduce[0 < a < 1 && 0 < s < 1 && 2 s/yS^2 + pP s/((1 + a) (1 + yS)) <= 0, {a, s}, Reals] === False &&
    Reduce[0 < s < 1 && (pP /. a -> 1/2) < 0, s, Reals] === False];
  mOf[av_?NumericQ] := Module[{g = (pP wW) /. a -> av},
     2 NIntegrate[g /. s -> Sin[u]^2, {u, 0, Pi/2}, WorkingPrecision -> 40, PrecisionGoal -> 30]];
  {lo, hi} = {1/10, 1/2};
  If[! (mOf[lo] < 0 < mOf[hi]), tovAbort["M does not change sign between a = 1/10 and 1/2"]];
  Do[mid = (lo + hi)/2; If[mOf[mid] < 0, lo = mid, hi = mid], {60}];
  tovValue["mu_delta_F_lower", N[1 - hi^2, 25], 20];
  tovValue["mu_delta_F_upper", N[1 - lo^2, 25], 20];
  tovRecord["florides.centre.mu_delta_bracket", "numerical_check",
   "M changes sign on the final bisection bracket in a (width 2^-60 (4/10)); mu_Delta^F = 1 - a^2 lies in the recorded bracket",
   mOf[lo] < 0 < mOf[hi] && 947/1000 < 1 - hi^2 < 1 - lo^2 < 948/1000]];

(* --- negative control: the isotropic uniform-density lapse does not satisfy the tau curvature inequality --- *)
Module[{nIso, lIso, aB = 335/1000},
  nIso = ((3 a - y)/2)^2/(a^2 - ((3 a - y)/2)^2);
  lIso = Simplify[y (1 - (1 - y^2)/y D[Log[nIso]/2, y]), dom];
  tovRecord["florides.negative_control.isotropic", "negative_control",
   "for the isotropic interior at a = 335/1000 (mu above mu_tau) the proper-time L_y is positive somewhere in (a, 1)",
   Reduce[aB < y < 1 && (D[lIso, y] /. a -> aB) > 0, y, Reals] =!= False]];

(* --- diagnostic: half-aperture of tau decreasing on a grid at mu = 19/20 --- *)
Module[{muv = 19/20, av, half, qs, vals},
  av = Sqrt[1 - muv];
  half[qv_] := Module[{yq = Sqrt[1 - muv qv^2], bv, v},
    bv = Sqrt[(1 - yq^2)/muv av/(yq - av)];
    v = NIntegrate[Module[{yy = av + (yq - av) Sin[u]^2},
       bv Sqrt[muv] (yy - av)/((1 - yy^2) Sqrt[(yy - av) (av (1 - yy^2) - muv bv^2 (yy - av))]) (yq - av) Sin[2 u]],
      {u, 0, Pi/2}, WorkingPrecision -> 40, PrecisionGoal -> 25, MaxRecursion -> 40];
    If[! NumberQ[v] || Abs[Im[v]] > 10^-20, tovAbort["an aperture quadrature failed"]]; Re[v]];
  qs = Table[k/20, {k, 1, 19}]; vals = half /@ qs;
  tovValue["half_aperture_proper_mu_19_20_q_1_20", First[vals], 15];
  tovRecord["florides.diagnostic.proper_decreasing", "diagnostic",
   "Delta_tau/2 strictly decreasing on q = 1/20, ..., 19/20 at mu = 19/20, between Pi/2 and 0",
   And @@ Thread[Rest[vals] < Most[vals]] && First[vals] < N[Pi/2] && Last[vals] > 0]];

tovFinish["florides_interior"];
