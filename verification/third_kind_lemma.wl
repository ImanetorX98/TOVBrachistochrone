(* Complete third-kind integrals on a genus-two sextic in closed form (shared library).

   Lemma (paper, label lem:third-kind). Let w^2 = S(y) be a square-free sextic, lo < hi consecutive real
   roots with w = +Sqrt[S] > 0 on (lo, hi), U = (1, y) and Baker's second-kind numerators
     R1(y) = lam3 y/4 + lam4 y^2/2 + 3 lam5 y^3/4 + lam6 y^4,   R2(y) = lam5 y^2/4 + lam6 y^3/2.
   With omegahat_i = Int_lo^hi U_i dy/w and etahat_i = Int_lo^hi R_i dy/w, for a finite branch point e and a point
   p off [lo, hi] with S(p) != 0, and u_i, r_i the integrals of U_i dy/w, R_i dy/w from e to p along one lifted
   path whose projection avoids [lo, hi]:
     Ihat(p) := Int_lo^hi dy/((y - p) w) = (2/w(p)) Sum_i [etahat_i u_i(e -> p) - omegahat_i r_i(e -> p)].

   Hypotheses of the lemma: S square-free of degree six; lo < hi consecutive real roots with S > 0 between them;
   e a finite simple root outside [lo, hi]; p finite, outside [lo, hi], S(p) != 0; one lifted affine path.
   This implementation covers REAL e and p only. Path: e -> e + i H -> p + i H -> p, with 0 < H at most half the
   smallest positive imaginary part of the other roots (H = 1/2 if all are real). Then the horizontal segment lies
   strictly between the real axis and every non-real root, the vertical segments meet the real axis only at e and p,
   and the projection misses [lo, hi]: the path avoids every root and the cycle. On each straight segment
   (z - z_j)^(1/2) = c_j Sqrt[(z - z_j)/(z_s - z_j)] with the principal root is the exact continuation, because the
   ratio runs on a segment from 1 that misses 0.
   Enforced here (exit status 2 otherwise): e is a root, p is not, e and p are real and outside [lo, hi].
   NOT certified here: the degree, square-freeness, that lo and hi are consecutive roots and the sign of S between
   them; the callers check these as diagnostics. Root locations come from NSolve at the working precision.

   The numerical integrals evaluate the closed form; they are not a fit. *)

tovLams = Array[tovLam, 7, 0];
tovBakerF = 2 tovLam[0] + tovLam[1] (tovX + tovP) + 2 tovLam[2] tovX tovP + tovLam[3] tovX tovP (tovX + tovP) +
   2 tovLam[4] tovX^2 tovP^2 + tovLam[5] tovX^2 tovP^2 (tovX + tovP) + 2 tovLam[6] tovX^3 tovP^3;
tovGenericS[z_] := Sum[tovLam[k] z^k, {k, 0, 6}];
tovBakerQuotient = PolynomialQuotient[
   Expand[tovBakerF - D[tovGenericS[tovP], tovP] (tovX - tovP) - 2 tovGenericS[tovP]], (tovX - tovP)^2, tovX];
tovBakerRGeneric = Expand[{Coefficient[tovBakerQuotient, tovX, 0], Coefficient[tovBakerQuotient, tovX, 1]}/4];

(* R1, R2 as polynomials in the variable var for a concrete sextic sx in var *)
tovBakerR[sx_, var_] := tovBakerRGeneric /. Thread[tovLams -> PadRight[CoefficientList[sx, var], 7]] /. tovP -> var;

tovIntervalIntegral[g_, sx_, var_, lo_, hi_, wp_] := NIntegrate[
   Module[{yy = lo + (hi - lo) Sin[tovU]^2}, (g /. var -> yy)/Sqrt[sx /. var -> yy] (hi - lo) Sin[2 tovU]],
   {tovU, 0, Pi/2}, WorkingPrecision -> wp, PrecisionGoal -> wp - 20, MaxRecursion -> 40];

(* mutation: None, "sign" (eta.u + omega.r) or "R1" (drop lam3 y/4 from R1); used only by negative controls *)
tovThirdKindComplete[sx_, var_, lo_, hi_, e_, pv_, wp_, mutation_: None] := Module[
   {rr, roots, others, hgt, lead, acc, cs, mult, om, et, wEnd},
   rr = tovBakerR[sx, var] + If[mutation === "R1", {-Coefficient[sx, var, 3] var/4, 0}, {0, 0}];
   roots = var /. NSolve[sx == 0, var, WorkingPrecision -> wp];
   If[Min[Abs[roots - e]] > 10^-30, tovAbort["the base point is not a root of S"]];
   If[Min[Abs[roots - pv]] < 10^-40, tovAbort["the pole p is a root of S"]];
   If[! (Element[e, Reals] && Element[pv, Reals]), tovAbort["this implementation needs a real base and pole"]];
   If[lo <= e <= hi || lo <= pv <= hi, tovAbort["the base or the pole lies on the projection of the cycle"]];
   others = Delete[roots, First[Ordering[Abs[roots - e], 1]]];
   hgt = With[{im = Select[Abs[Im[others]], # > 10^-30 &]}, If[im == {}, 1/2, Min[im]/2]];
   lead = Coefficient[sx, var, 6];
   acc = {{0, 0}, {0, 0}};
   (* first segment z = e + I hgt s^2: the base factor is Sqrt[I hgt] s, regular after dz = 2 I hgt s ds *)
   Module[{c0 = Sqrt[e - #] & /@ others, prod, g},
    prod[sv_] := Times @@ (c0 Sqrt[(e + I hgt sv^2 - others)/(e - others)]);
    g[h_] := NIntegrate[(h /. var -> e + I hgt sv^2) 2 I hgt/(Sqrt[lead] Sqrt[I hgt] prod[sv]),
      {sv, 0, 1}, WorkingPrecision -> wp, PrecisionGoal -> wp - 20, MaxRecursion -> 40];
    acc += {g /@ {1, var}, g /@ rr};
    cs = c0 Sqrt[(e + I hgt - others)/(e - others)];
    mult = Sqrt[I hgt]];
   Do[Module[{zs = seg[[1]], ze = seg[[2]], prod, g},
     (* numerical guard only; the construction of the path already excludes the roots *)
     If[Min[Abs[(zs + (ze - zs) Range[0, 1, 1/200]) - #] & /@ others] < 10^-20,
      tovAbort["a path segment passes through a root"]];
     prod[t_] := Times @@ (cs Sqrt[(zs + (ze - zs) t - others)/(zs - others)]);
     g[h_] := NIntegrate[(h /. var -> zs + (ze - zs) t) (ze - zs)/
         (Sqrt[lead] mult Sqrt[(zs + (ze - zs) t - e)/(zs - e)] prod[t]),
       {t, 0, 1}, WorkingPrecision -> wp, PrecisionGoal -> wp - 20, MaxRecursion -> 40];
     acc += {g /@ {1, var}, g /@ rr};
     cs = cs Sqrt[(ze - others)/(zs - others)];
     mult = mult Sqrt[(ze - e)/(zs - e)]],
    {seg, {{e + I hgt, pv + I hgt}, {pv + I hgt, pv}}}];
   wEnd = Sqrt[lead] mult Times @@ cs;
   om = tovIntervalIntegral[#, sx, var, lo, hi, wp] & /@ {1, var};
   et = tovIntervalIntegral[#, sx, var, lo, hi, wp] & /@ rr;
   If[mutation === "sign", (2/wEnd) (et.acc[[1]] + om.acc[[2]]), (2/wEnd) (et.acc[[1]] - om.acc[[2]])]];
