#!/usr/bin/env wolframscript
(* Monotonicity of the Clairaut function J(x) = x n(x) for the uniform-density star (all checks exact).

   With a = yR = Sqrt[1 - mu], Q = (y - a)(5a - y):
     coordinate time  J_t^2 = 16 H_t/mu,  H_t = (1 - y^2)/((3a - y)^2 Q),
     proper time      J_tau^2 = H_tau/mu, H_tau = (1 - y^2)(3a - y)^2/Q.
   Since dy/dx < 0, J is increasing in x where H' < 0.

   Claims checked:
     - H_t' = 2 Qt(y)/((3a - y)^3 Q^2) with the monic quartic Qt = y^4 - 3a y^3 - (2 + 9a^2) y^2 + (12a + 15a^3) y - 14a^2;
     - Disc_y(Qt) = 16 a^2 (a^2 - 1)(25 a^2 - 1) C(a^2), C(A) = 3645 A^3 + 3888 A^2 - 1080 A + 64;
     - C has exactly one root A_star in (1/9, 1), simple; mu_star = 1 - A_star;
     - Qt(a) = 4 a^2 (a^2 - 1) < 0 and Qt(1) = (a - 1)(3a - 1)(5a - 1) < 0 on 1/3 < a < 1;
     - Qt has two roots in (a, 1) at a = 17/50 and none at a = 2/5 (Sturm counts); hence two simple interior zeros for
       mu > mu_star, none for mu < mu_star;
     - at a = Sqrt[A_star], gcd(Qt, Qt') has degree one: a double root (critical turning radius q_star);
     - proper time: H_tau' < 0 on a < y < 1 for every 1/3 < a < 1 (Reduce), so J_tau is increasing for all mu < 8/9. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "common.wl"}]];
q = (y - a) (5 a - y);
ht = (1 - y^2)/((3 a - y)^2 q);
htau = (1 - y^2) (3 a - y)^2/q;
qt = y^4 - 3 a y^3 - (2 + 9 a^2) y^2 + (12 a + 15 a^3) y - 14 a^2;
cubic[aa_] := 3645 aa^3 + 3888 aa^2 - 1080 aa + 64;

tovRecord["monotonicity.coordinate_derivative", "exact", "H_t' = 2 Qt/((3a - y)^3 Q^2)",
  Simplify[D[ht, y] - 2 qt/((3 a - y)^3 q^2)] === 0];
tovRecord["monotonicity.discriminant", "exact", "Disc_y(Qt) = 16 a^2 (a^2 - 1)(25 a^2 - 1) C(a^2)",
  Expand[Discriminant[qt, y] - 16 a^2 (a^2 - 1) (25 a^2 - 1) cubic[a^2]] === 0];
tovRecord["monotonicity.threshold_unique", "exact", "C has exactly one root in (1/9, 1), and it is simple",
  CountRoots[cubic[A], {A, 1/9, 1}] == 1 && Exponent[PolynomialGCD[cubic[A], D[cubic[A], A]], A] == 0];
tovRecord["monotonicity.endpoints", "exact",
  "Qt(a) = 4a^2(a^2 - 1) and Qt(1) = (a - 1)(3a - 1)(5a - 1), both negative on 1/3 < a < 1",
  Expand[(qt /. y -> a) - 4 a^2 (a^2 - 1)] === 0 && Expand[(qt /. y -> 1) - (a - 1) (3 a - 1) (5 a - 1)] === 0 &&
   Reduce[1/3 < a < 1 && ((qt /. y -> a) >= 0 || (qt /. y -> 1) >= 0), a, Reals] === False];
tovRecord["monotonicity.root_counts", "exact", "Qt has 2 roots in (a, 1) at a = 17/50 and 0 at a = 2/5",
  CountRoots[qt /. a -> 17/50, {y, 17/50, 1}] == 2 && CountRoots[qt /. a -> 2/5, {y, 2/5, 1}] == 0];
Module[{rootA, g},
  rootA = SelectFirst[Table[Root[cubic[#] &, k], {k, 3}], Element[#, Reals] && 1/9 < # < 1 &];
  tovRecord["monotonicity.threshold_bracket", "exact", "the root A* of C in (1/9, 1) lies in ((17/50)^2, (2/5)^2)",
   (17/50)^2 < rootA < (2/5)^2 && 1/9 < rootA < 1];
  g = PolynomialGCD[qt /. a -> Sqrt[rootA], D[qt, y] /. a -> Sqrt[rootA], Extension -> Automatic];
  tovRecord["monotonicity.double_root_at_threshold", "exact", "gcd(Qt, Qt') has degree one at a = Sqrt[A*]",
   Exponent[g, y] == 1];
  Module[{yc = y /. First[Solve[g == 0, y]]},
   tovValue["mu_star", 1 - rootA];
   tovValue["q_star", Sqrt[(1 - yc^2)/(1 - rootA)]]]];
tovRecord["monotonicity.proper_time", "exact", "H_tau' < 0 for a < y < 1 and every 1/3 < a < 1",
  Reduce[1/3 < a < 1 && a < y < 1 && Together[D[htau, y]] >= 0, {a, y}, Reals] === False];
tovRecord["monotonicity.non_monotone_witness", "exact", "a = 167/500: H_t' > 0 at y = 99/100 and < 0 at y = 9999/10000",
  (D[ht, y] /. {a -> 167/500, y -> 99/100}) > 0 && (D[ht, y] /. {a -> 167/500, y -> 9999/10000}) < 0];

tovFinish["clairaut_monotonicity"];
