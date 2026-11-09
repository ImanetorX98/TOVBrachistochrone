#!/usr/bin/env wolframscript
(* Diagnostics (sampled, not proofs): the profile q -> Delta(q) of the uniform-density brachistochrones.

   Statements sampled here, labelled as diagnostics in the paper:
     - proper time: Delta_tau(q) decreases on a 39-point grid for mu = 0.01 ... 0.88 (the decrease is FALSE near the
       Buchdahl limit, see centre_slope_threshold.wls; only the bound Delta_tau < Pi is proved);
     - coordinate time: Delta_t(q) decreases on 108 values of q at mu = 0.70, 0.75, 0.76 and mu_Delta - 0.001,
       and is not monotone at mu_Delta + 0.001 and 0.78 (supports the conjecture that mu_Delta is the global onset). *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "common.wl"}]];
nsqY["proper"][z_, a_] := (3 a - z)^2/((z - a) (5 a - z));
nsqY["coordinate"][z_, a_] := 16/((3 a - z)^2 (z - a) (5 a - z));
halfAp[m_, muv_, qv_] := Module[{a = Sqrt[1 - muv], yq, bv},
   yq = Sqrt[1 - muv qv^2]; bv = Sqrt[nsqY[m][yq, a] (1 - yq^2)/muv];
   NIntegrate[Module[{yy = a + (yq - a) Sin[u]^2},
      bv Sqrt[muv]/((1 - yy^2) Sqrt[nsqY[m][yy, a] (1 - yy^2) - muv bv^2]) (yq - a) Sin[2 u]],
     {u, 0, Pi/2}, WorkingPrecision -> 30, PrecisionGoal -> 15, MaxRecursion -> 40]];

coarse = Table[k/40, {k, 1, 39}];
Module[{ok = True},
  Do[Module[{v = halfAp["proper", muv, #] & /@ coarse},
     ok = ok && AllTrue[Differences[v], # < 0 &] && Max[v] < Pi/2],
   {muv, {1/100, 1/10, 3/10, 1/2, 7/10, 8/10, 85/100, 88/100}}];
  tovRecord["profile.proper_decreasing_on_grid", "diagnostic",
   "Delta_tau(q) decreases and stays below Pi on a 39-point grid for mu in {0.01, ..., 0.88}", ok]];

muDelta = 7685752721110638606517046/10^25;
fine = Join[Table[k/1000, {k, 1, 9}], Table[k/100, {k, 1, 99}]];
onsetLabel[muv_] := Which[muv === muDelta - 1/1000, "mu_delta_minus_1_1000", muv === muDelta + 1/1000,
   "mu_delta_plus_1_1000", True, "mu_" <> tovTag[muv]];
Do[Module[{v = halfAp["coordinate", muv, #] & /@ fine, mono, expected},
   mono = AllTrue[Differences[v], # < 0 &];
   expected = muv < muDelta;
   tovRecord["profile.coordinate_onset." <> onsetLabel[muv], "diagnostic",
    "Delta_t(q) on 108 values of q is " <> If[expected, "monotone (below mu_Delta)", "not monotone (above mu_Delta)"],
    mono === expected]],
  {muv, {70/100, 75/100, 76/100, muDelta - 1/1000, muDelta + 1/1000, 78/100}}];

tovFinish["aperture_profile"];
