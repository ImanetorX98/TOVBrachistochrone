#!/usr/bin/env wolframscript
(* The whole uniform-density brachistochrone in closed form.

   For a <= y <= yq (descent from the surface to the turning point), w = +Sqrt[S] on the physical sheet:
     phi(y) = b Sqrt[mu] [ c0 u0(y) + c1 u1(y) + aP PiP(y) + aM PiM(y) ],
     u_k(y) = Int_a^y z^k dz/w (first kind),  Pi+-(y) = Int_a^y dz/((z -+ 1) w) (third kind),
   with the exact coefficients of uniform_density_reduction.wls; the ascent branch is Delta - phi(y).
   In physical units y = Sqrt[1 - mu r^2/R^2] = Sqrt[1 - 2GM r^2/(R^3 c^2)], and the Clairaut constant with the
   dimension of a length is R b.

   Independent check (different method): geodesics of the optical metric n^2 (h dx^2 + x^2 dphi^2) in Hamiltonian
   form, H = [P^2/(n^2 h) + b^2/(n^2 x^2)]/2, integrated from the turning point (P = 0, regular for a simple
   turning point because dP/ds = J'(q)/b > 0) to the surface; phi compared at five radii. The tolerance 1e-10
   reflects the integrator (its drift of 2H = 1 is reported), not the closed form. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "common.wl"}]];
wp = 40;
dPoly = (y - yR) (5 yR - y);
pPoly["proper"] = (3 yR - y)^2 (1 - y^2) - mu b^2 dPoly;
pPoly["coordinate"] = 16 (1 - y^2) - mu b^2 (3 yR - y)^2 dPoly;
nsqY["proper"][z_, a_] := (3 a - z)^2/((z - a) (5 a - z));
nsqY["coordinate"][z_, a_] := 16/((3 a - z)^2 (z - a) (5 a - z));
decomposition["proper"] = {{1, 0}, (yR - 1) (5 yR - 1)/2, -(1 + yR) (1 + 5 yR)/2};
decomposition["coordinate"] = {{9 yR, -1}, (yR - 1) (3 yR - 1) (5 yR - 1)/2,
   -(1 + yR) (1 + 3 yR) (1 + 5 yR)/2};

setup[m_, muv_, qv_] := Module[{a = Sqrt[1 - muv], yq, bv},
   yq = Sqrt[1 - muv qv^2]; bv = Sqrt[nsqY[m][yq, a] (1 - yq^2)/muv];
   <|"a" -> a, "yq" -> yq, "b" -> bv, "S" -> Expand[dPoly pPoly[m] /. {yR -> a, mu -> muv, b -> bv}]|>];
incomplete[g_, st_, y1_] := NIntegrate[
   Module[{zz = st["a"] + (y1 - st["a"]) Sin[t]^2}, (g /. y -> zz)/Sqrt[st["S"] /. y -> zz] (y1 - st["a"]) Sin[2 t]],
   {t, 0, Pi/2}, WorkingPrecision -> wp, PrecisionGoal -> 25, MaxRecursion -> 40];
phiClosed[m_, st_, y1_, muv_, swap_: False] := Module[{pc, aP, aM},
   {pc, aP, aM} = decomposition[m] /. yR -> st["a"];
   If[swap, {aP, aM} = {aM, aP}];
   st["b"] Sqrt[muv] (pc.(incomplete[#, st, y1] & /@ {1, y}) + aP incomplete[1/(y - 1), st, y1] +
      aM incomplete[1/(y + 1), st, y1])];

odeSolution[m_, muv_, qv_] := Module[{a = Sqrt[1 - muv], nsqx, bv, ham, dHdX, sol, xs, px, ps, s, XX, PP},
   nsqx[X_] := nsqY[m][Sqrt[1 - muv X^2], a];
   bv = Sqrt[nsqx[qv]] qv;
   ham = (PP^2 (1 - muv XX^2)/nsqx[XX] + bv^2/(nsqx[XX] XX^2))/2;
   dHdX = D[ham, XX];
   sol = NDSolve[{xs'[s] == PP (1 - muv XX^2)/nsqx[XX] /. {XX -> xs[s], PP -> px[s]},
      px'[s] == -dHdX /. {XX -> xs[s], PP -> px[s]}, ps'[s] == bv/(nsqx[xs[s]] xs[s]^2),
      xs[0] == qv, px[0] == 0, ps[0] == 0, WhenEvent[xs[s] > 1 - 10^-6, "StopIntegration"]},
     {xs, px, ps}, {s, 0, 100}, WorkingPrecision -> 45, PrecisionGoal -> 28, AccuracyGoal -> 28,
     Method -> "ExplicitRungeKutta", MaxSteps -> 10^7];
   <|"x" -> (xs /. First[sol]), "phi" -> (ps /. First[sol]),
     "twoH" -> Function[{sv}, Evaluate[2 ham /. {XX -> xs[sv], PP -> px[sv]} /. First[sol]]]|>];

samples = {{"proper", 3/10, 3/10}, {"proper", 4/5, 1/5}, {"coordinate", 1/2, 7/10}, {"coordinate", 7/10, 2/5}};
Do[Module[{m = smp[[1]], muv = smp[[2]], qv = smp[[3]], st, half, ode, smax, errs, drift, tag},
   tag = m <> ".mu_" <> tovTag[muv] <> ".q_" <> tovTag[qv];
   st = setup[m, muv, qv];
   half = phiClosed[m, st, st["yq"], muv];
   ode = odeSolution[m, muv, qv];
   smax = ode["x"]["Domain"][[1, 2]];
   If[Abs[ode["x"][smax] - (1 - 10^-6)] > 10^-15, tovAbort["the ODE stopped before reaching x = 1 - 1e-6"]];
   errs = Table[Abs[(half - phiClosed[m, st, Sqrt[1 - muv ode["x"][sv]^2], muv]) - ode["phi"][sv]],
     {sv, smax {1/10, 3/10, 1/2, 7/10, 9/10}}];
   drift = Max[Table[Abs[ode["twoH"][sv] - 1], {sv, smax Range[1/10, 9/10, 1/10]}]];
   tovValue["half_aperture." <> tag, half];
   tovValue["max_curve_deviation." <> tag, Max[errs], 3];
   tovValue["ode_energy_drift." <> tag, drift, 3];
   tovRecord["whole_curve.vs_ode." <> tag, "numerical_check",
    "closed-form phi equals the Hamiltonian ODE at five radii within 1e-10, with 2H = 1 conserved within 1e-10",
    Max[errs] < 10^-10 && drift < 10^-10]],
  {smp, samples}];

Module[{st = setup["proper", 3/10, 3/10], yv = Sqrt[1 - 3/10 (1/2)^2]},
  tovRecord["whole_curve.negative_control.swapped_residues", "negative_control",
   "swapping the residues aP and aM changes phi at y(x = 1/2) by more than 1e-6",
   Abs[phiClosed["proper", st, yv, 3/10, True] - phiClosed["proper", st, yv, 3/10]] > 10^-6]];

tovFinish["whole_curve"];
