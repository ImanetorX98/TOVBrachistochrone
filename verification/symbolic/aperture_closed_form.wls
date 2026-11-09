#!/usr/bin/env wolframscript
(* The aperture of the uniform-density brachistochrones in closed form.

   Claims checked:
     - Baker's symmetric polynomial F(x,p) for a general sextic, the division
       F - S'(p)(x - p) - 2 S(p) = 4 (x - p)^2 Sum_i U_i(x) R_i(p), and the bidifferential identity
       Omega = d_p[(w(x) + w(p)) dx/(2 (x - p) w(x))] + Sum_i U_i(x) R_i(p) dx dp/(w(x) w(p))  (exact);
     - the complete third-kind integrals Ihat(+1), Ihat(-1) by the lemma (third_kind_lemma.wl), with the uniform
       base point e = 5 yR, against direct quadrature (numerical check);
     - Delta/2 = b Sqrt[mu] [ c.omegahat + aP Ihat(1) + aM Ihat(-1) ] against the x-quadrature (numerical check);
     - no root of S inside (yR, yq) on the samples (diagnostic; the hypothesis of the lemma);
     - two mutations of the lemma are detected (negative controls). *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "common.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "third_kind_lemma.wl"}]];
wp = 50;

(* exact identities for a general sextic *)
tovRecord["aperture.baker_symmetric", "exact", "Baker's F(x,p) is symmetric in x <-> p",
  Expand[tovBakerF - (tovBakerF /. {tovX -> tovP, tovP -> tovX})] === 0];
tovRecord["aperture.baker_division", "exact",
  "F - S'(p)(x - p) - 2 S(p) is divisible by (x - p)^2 with a quotient linear in x",
  PolynomialRemainder[Expand[tovBakerF - D[tovGenericS[tovP], tovP] (tovX - tovP) - 2 tovGenericS[tovP]],
     (tovX - tovP)^2, tovX] === 0 && Exponent[tovBakerQuotient, tovX] == 1];
Block[{wx, wq, lhs, rhs, numr},
  lhs = (2 wx wq + tovBakerF)/(4 (tovX - tovP)^2 wx wq);
  rhs = D[(wx + wpf[tovP])/(2 (tovX - tovP) wx), tovP] /.
     {Derivative[1][wpf][tovP] -> D[tovGenericS[tovP], tovP]/(2 wq), wpf[tovP] -> wq};
  rhs = rhs + {1, tovX}.tovBakerRGeneric/(wx wq);
  numr = Numerator[Together[lhs - rhs]];
  numr = PolynomialRemainder[PolynomialRemainder[Expand[numr], wq^2 - tovGenericS[tovP], wq],
    wx^2 - tovGenericS[tovX], wx];
  tovRecord["aperture.bidifferential_identity", "exact",
   "Omega = d_p[(w(x)+w(p)) dx/(2(x-p)w(x))] + Sum U_i(x) R_i(p) dx dp/(w(x) w(p))", Expand[numr] === 0]];
tovValue["baker_R1", ToString[tovBakerRGeneric[[1]], InputForm]];
tovValue["baker_R2", ToString[tovBakerRGeneric[[2]], InputForm]];

(* the two brachistochrone sextics *)
dPoly = (y - yR) (5 yR - y);
pPoly["proper"] = (3 yR - y)^2 (1 - y^2) - mu b^2 dPoly;
pPoly["coordinate"] = 16 (1 - y^2) - mu b^2 (3 yR - y)^2 dPoly;
nsq["proper"][z_] := ((3 yR - z)/2)^2/((z - yR) (5 yR - z)/4);
nsq["coordinate"][z_] := 1/(((3 yR - z)/2)^2 (z - yR) (5 yR - z)/4);
decomposition["proper"] = {{1, 0}, (yR - 1) (5 yR - 1)/2, -(1 + yR) (1 + 5 yR)/2};
decomposition["coordinate"] = {{9 yR, -1}, (yR - 1) (3 yR - 1) (5 yR - 1)/2,
   -(1 + yR) (1 + 3 yR) (1 + 5 yR)/2};
data[m_, muv_, qv_] := Module[{a = Sqrt[1 - muv], yq, bv},
   yq = Sqrt[1 - muv qv^2];
   bv = Sqrt[nsq[m][yq] (1 - yq^2)/muv /. yR -> a];
   <|"a" -> a, "yq" -> yq, "b" -> bv, "S" -> Expand[dPoly pPoly[m] /. {yR -> a, mu -> muv, b -> bv}]|>];
ihat[d_, pv_, mutation_: None] := tovThirdKindComplete[d["S"], y, d["a"], d["yq"], 5 d["a"], pv, wp, mutation];
ihatDirect[d_, pv_] := tovIntervalIntegral[1/(y - pv), d["S"], y, d["a"], d["yq"], wp];
inDomain[d_] := Module[{rts = y /. NSolve[d["S"] == 0, y, WorkingPrecision -> 60]},
   ! Or @@ (Abs[Im[#]] < 10^-25 && d["a"] + 10^-25 < Re[#] < d["yq"] - 10^-25 & /@ rts)];
halfApertureClosed[m_, muv_, qv_] := Module[{d = data[m, muv, qv], pc, aP, aM, om},
   {pc, aP, aM} = decomposition[m] /. yR -> d["a"];
   om = tovIntervalIntegral[#, d["S"], y, d["a"], d["yq"], wp] & /@ {1, y};
   d["b"] Sqrt[muv] (pc.om + aP ihat[d, 1] + aM ihat[d, -1])];
halfApertureX[m_, muv_, qv_] := Module[{a = Sqrt[1 - muv], fx, nx, bv},
   fx[xx_] := ((3 a - Sqrt[1 - muv xx^2])/2)^2;
   nx[xx_] := If[m == "proper", Sqrt[fx[xx]/(a^2 - fx[xx])], 1/Sqrt[fx[xx] (a^2 - fx[xx])]];
   bv = nx[qv] qv;
   NIntegrate[Module[{xx = qv + (1 - qv) Sin[u]^2},
     bv/Sqrt[1 - muv xx^2]/(xx Sqrt[nx[xx]^2 xx^2 - bv^2]) (1 - qv) Sin[2 u]],
    {u, 0, Pi/2}, WorkingPrecision -> wp, PrecisionGoal -> 30, MaxRecursion -> 40]];

samples = {{1/10, 1/2}, {3/10, 3/10}, {1/2, 7/10}, {4/5, 1/5},
   {8888/10000, 9/10}, {8888/10000, 9999/10000}, {2/10, 1/10000}, {8888/10000, 1/10000}};
tag[m_, smp_] := m <> ".mu_" <> tovTag[smp[[1]]] <> ".q_" <> tovTag[smp[[2]]];
Do[Module[{d = data[m, smp[[1]], smp[[2]]]},
   tovRecord["aperture.domain." <> tag[m, smp], "diagnostic",
    "no root of S inside (yR, yq), roots located at 60 digits", inDomain[d]];
   Do[Module[{cl = ihat[d, pv], di = ihatDirect[d, pv]},
     tovRecord["aperture.third_kind.p" <> If[pv > 0, "_plus1.", "_minus1."] <> tag[m, smp], "numerical_check",
      "Ihat(" <> ToString[pv] <> ") by the lemma equals the direct quadrature to 1e-25 relative",
      Abs[cl - di] < 10^-25 Abs[di]]], {pv, {1, -1}}]],
  {m, {"proper", "coordinate"}}, {smp, samples}];
Do[Module[{cl = halfApertureClosed[m, smp[[1]], smp[[2]]], dx = halfApertureX[m, smp[[1]], smp[[2]]]},
   tovValue["half_aperture." <> tag[m, smp], Re[cl]];
   tovRecord["aperture.half_aperture." <> tag[m, smp], "numerical_check",
    "Delta/2 in closed form equals the x-quadrature to 1e-25 relative, imaginary part below 1e-25",
    Abs[Im[cl]] < 10^-25 && Abs[cl - dx] < 10^-25 Abs[dx]]],
  {m, {"proper", "coordinate"}}, {smp, samples}];

(* negative controls on the coordinate-time case (0.8888, 0.9), where all six roots are real *)
Do[Module[{d = data["coordinate", 8888/10000, 9/10], cl, di},
   cl = ihat[d, -1, mut]; di = ihatDirect[d, -1];
   tovRecord["aperture.negative_control." <> mut, "negative_control",
    "mutation '" <> mut <> "' of the lemma changes Ihat(-1) by more than 1e-10 relative",
    Abs[cl - di] > 10^-10 Abs[di]]], {mut, {"sign", "R1"}}];

tovFinish["aperture_closed_form"];
