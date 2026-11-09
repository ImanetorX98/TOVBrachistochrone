#!/usr/bin/env wolframscript
(* Travel times of the uniform-density brachistochrones in closed form (units R = c = 1; multiply by R/c).

   Along an extremal, with the sextic S = D P and w = +Sqrt[S] on (a, yq):
     n_tau dl = (3a - y)^2 dy/(Sqrt[mu] w),   n_t dl = 16 dy/(Sqrt[mu] (3a - y) w)   (descent, y increasing),
   and d tau = n_tau dl, dt = a n_t dl (E = f dt/d tau = a). Hence
     t_half   = (16 a/Sqrt[mu]) Int_a^yq dy/((3a - y) w) = -(16 a/Sqrt[mu]) Ihat_S(3a),
     tau_half = (1/Sqrt[mu]) [9 a^2 omegahat_0 - 6 a omegahat_1 + m2],
     m2 = Int_a^yq y^2 dy/w = Ihat_Stilde(0) on Stilde(v) = v^6 S(1/v), interval [1/yq, 1/a]   (v = 1/y),
   with Ihat from the lemma in third_kind_lemma.wl (bases 5a and 1/(5a)).

   Checks: the time elements (exact); the closed forms against direct quadrature; the optical length of the
   Hamiltonian ODE (the ODE stops at x = 1 - 1e-12, whose missing tails 2 Sqrt[eps]/Sqrt[mu] for tau and
   2 Sqrt[eps]/(a Sqrt[mu]) for t are added back); tau < t on the samples; the Newtonian limit
   t_half ~ tau_half ~ (Pi/2) Sqrt[2/mu] Sqrt[1 - q^2] at mu = 1e-8; a misplaced pole (negative control). *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "common.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "third_kind_lemma.wl"}]];
wp = 50;
dPoly = (y - yR) (5 yR - y);
pPoly["proper"] = (3 yR - y)^2 (1 - y^2) - mu b^2 dPoly;
pPoly["coordinate"] = 16 (1 - y^2) - mu b^2 (3 yR - y)^2 dPoly;
nsqY["proper"][z_, a_] := (3 a - z)^2/((z - a) (5 a - z));
nsqY["coordinate"][z_, a_] := 16/((3 a - z)^2 (z - a) (5 a - z));

Do[Module[{nn = nsqY[m][y, yR], lhs, rhs},
   lhs = nn^2/(mu (nn (1 - y^2) - mu b^2));
   rhs = If[m == "proper", (3 yR - y)^4/(mu dPoly pPoly[m]), 256/(mu (3 yR - y)^2 dPoly pPoly[m])];
   tovRecord["travel_times.time_element." <> m, "exact",
    "(n dl/dy)^2 = n^4/(mu (n^2 (1 - y^2) - mu b^2)) reduces to the displayed rational form, " <> m <> " time",
    Simplify[Together[lhs - rhs]] === 0]],
  {m, {"proper", "coordinate"}}];

setup[m_, muv_, qv_] := Module[{a = Sqrt[1 - muv], yq, bv},
   yq = Sqrt[1 - muv qv^2]; bv = Sqrt[nsqY[m][yq, a] (1 - yq^2)/muv];
   <|"a" -> a, "yq" -> yq, "b" -> bv, "S" -> Expand[dPoly pPoly[m] /. {yR -> a, mu -> muv, b -> bv}]|>];
tHalfClosed[muv_, qv_, pole_: 3] := Module[{st = setup["coordinate", muv, qv]},
   (16 st["a"]/Sqrt[muv]) (-tovThirdKindComplete[st["S"], y, st["a"], st["yq"], 5 st["a"], pole st["a"], wp])];
tauHalfClosed[muv_, qv_] := Module[{st = setup["proper", muv, qv], a, stil, om, m2},
   a = st["a"];
   stil = Expand[y^6 (st["S"] /. y -> 1/y)];
   om = tovIntervalIntegral[#, st["S"], y, a, st["yq"], wp] & /@ {1, y};
   m2 = tovThirdKindComplete[stil, y, 1/st["yq"], 1/a, 1/(5 a), 0, wp];
   (1/Sqrt[muv]) (9 a^2 om[[1]] - 6 a om[[2]] + m2)];
tHalfDirect[muv_, qv_] := Module[{st = setup["coordinate", muv, qv]},
   (16 st["a"]/Sqrt[muv]) tovIntervalIntegral[1/(3 st["a"] - y), st["S"], y, st["a"], st["yq"], wp]];
tauHalfDirect[muv_, qv_] := Module[{st = setup["proper", muv, qv]},
   (1/Sqrt[muv]) tovIntervalIntegral[(3 st["a"] - y)^2, st["S"], y, st["a"], st["yq"], wp]];

samples = {{1/10, 1/2}, {3/10, 3/10}, {1/2, 7/10}, {4/5, 1/5}};
tag[smp_] := "mu_" <> tovTag[smp[[1]]] <> ".q_" <> tovTag[smp[[2]]];
Do[Module[{tc = tHalfClosed[smp[[1]], smp[[2]]], td = tHalfDirect[smp[[1]], smp[[2]]],
    tauc = tauHalfClosed[smp[[1]], smp[[2]]], taud = tauHalfDirect[smp[[1]], smp[[2]]]},
   tovValue["t_half." <> tag[smp], Re[tc]]; tovValue["tau_half." <> tag[smp], Re[tauc]];
   tovRecord["travel_times.coordinate." <> tag[smp], "numerical_check",
    "t_half by the lemma at p = 3a equals the direct quadrature to 1e-25 relative (imaginary part < 1e-25)",
    Abs[Im[tc]] < 10^-25 && Abs[tc - td] < 10^-25 Abs[td]];
   tovRecord["travel_times.proper." <> tag[smp], "numerical_check",
    "tau_half by the lemma on the inverted sextic equals the direct quadrature to 1e-25 relative",
    Abs[Im[tauc]] < 10^-25 && Abs[tauc - taud] < 10^-25 Abs[taud]];
   tovRecord["travel_times.proper_less_than_coordinate." <> tag[smp], "numerical_check",
    "tau_half < t_half", Re[tauc] < Re[tc]]],
  {smp, samples}];

(* optical length of the Hamiltonian ODE, with the analytic tail added back *)
odeLength[m_, muv_, qv_] := Module[{a = Sqrt[1 - muv], nsqx, bv, ham, dHdX, sol, xs, px, s, XX, PP, xf, smax},
   nsqx[X_] := nsqY[m][Sqrt[1 - muv X^2], a];
   bv = Sqrt[nsqx[qv]] qv;
   ham = (PP^2 (1 - muv XX^2)/nsqx[XX] + bv^2/(nsqx[XX] XX^2))/2;
   dHdX = D[ham, XX];
   sol = NDSolve[{xs'[s] == PP (1 - muv XX^2)/nsqx[XX] /. {XX -> xs[s], PP -> px[s]},
      px'[s] == -dHdX /. {XX -> xs[s], PP -> px[s]}, xs[0] == qv, px[0] == 0,
      WhenEvent[xs[s] == 1 - 10^-12, "StopIntegration"]},
     {xs, px}, {s, 0, 100}, WorkingPrecision -> 45, PrecisionGoal -> 28, AccuracyGoal -> 28,
     Method -> "ExplicitRungeKutta", MaxSteps -> 10^7];
   xf = xs /. First[sol]; smax = xf["Domain"][[1, 2]];
   If[Abs[xf[smax] - (1 - 10^-12)] > 10^-20, tovAbort["the ODE stopped before reaching x = 1 - 1e-12"]];
   smax];
Do[Module[{muv = smp[[1]], qv = smp[[2]], a = Sqrt[1 - smp[[1]]], eps = 10^-12, lt, ltau},
   lt = a odeLength["coordinate", muv, qv] + 2 Sqrt[eps]/(a Sqrt[muv]);
   ltau = odeLength["proper", muv, qv] + 2 Sqrt[eps]/Sqrt[muv];
   tovRecord["travel_times.coordinate_vs_ode." <> tag[smp], "numerical_check",
    "a times the ODE optical length plus the analytic tail equals t_half to 1e-10 relative",
    Abs[lt - tHalfDirect[muv, qv]] < 10^-10 Abs[lt]];
   tovRecord["travel_times.proper_vs_ode." <> tag[smp], "numerical_check",
    "the ODE optical length plus the analytic tail equals tau_half to 1e-10 relative",
    Abs[ltau - tauHalfDirect[muv, qv]] < 10^-10 Abs[ltau]]],
  {smp, {{3/10, 3/10}, {1/2, 7/10}}}];

(* Newtonian limit *)
Do[Module[{muv = 10^-8, newton, rt, rtau},
   newton = (Pi/2) Sqrt[2/muv] Sqrt[1 - qv^2];
   rt = tHalfDirect[muv, qv]/newton; rtau = tauHalfDirect[muv, qv]/newton;
   tovRecord["travel_times.newtonian_limit.q_" <> tovTag[qv], "numerical_check",
    "at mu = 1e-8, t_half and tau_half equal (Pi/2) Sqrt[2/mu] Sqrt[1 - q^2] within 1e-6 relative",
    Abs[rt - 1] < 10^-6 && Abs[rtau - 1] < 10^-6]],
  {qv, {1/5, 1/2, 9/10}}];

tovRecord["travel_times.negative_control.misplaced_pole", "negative_control",
  "the closed form with the pole at -3a differs from the direct quadrature by more than 1e-6 relative",
  Abs[tHalfClosed[1/2, 7/10, -3] - tHalfDirect[1/2, 7/10]] > 10^-6 Abs[tHalfDirect[1/2, 7/10]]];

tovFinish["travel_times"];
