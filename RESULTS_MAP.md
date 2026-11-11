# Result map

The results implemented so far, with the script that checks them, the kind of evidence and the profile that runs
them. Claim labels are stable identifiers; the LaTeX labels of the manuscript will be added when the paper is
written, and only then can coverage of every claim of the paper be stated. Every check identifier in
`verification/manifests/jobs.json` names the claim it supports; the manifest is a reviewed contract and is never
regenerated from run outputs. Sample parameters in identifiers are exact rationals (`mu_1_10` is mu = 1/10).

Notation: units G = c = 1 and R = 1, `mu = 2M/R < 8/9`, `x = r/R`, `a = yR = sqrt(1 - mu)`,
`y = sqrt(1 - mu x^2)` (in physical units `y = sqrt(1 - mu r^2/R^2)`), `Q = (y - a)(5a - y)`,
`Q4 = Q (1 - y^2)`, `b = J(q)` the Clairaut constant at the turning radius `q`.

Evidence kinds: **exact** (symbolic identity or decision), **numerical check** (a quantitative check of a stated
result at given parameters and tolerance: a closed form against an independent quadrature or ODE, or the evaluation
of an analytic statement such as a sign, a bracket or an asymptotic ratio), **diagnostic** (sampled statement, not a proof),
**negative control** (a deliberate mutation that must be detected), **analytic** (a proof step reviewed outside
the executable suite).

## Uniform density: reduction to a genus-two curve

| Claim | Statement | Script | Checks | Evidence | Profile |
|---|---|---|---|---|---|
| `uniform.energy_gap` | E^2 - f = Q/4 | `symbolic/uniform_density_reduction.wls` | `reduction.energy_gap` | exact | quick |
| `uniform.genus_two` | dphi/dy = -b sqrt(mu) N/((1 - y^2) w) on the outgoing branch, w^2 = S sextic, square-free for generic (a, b) after substituting mu = 1 - a^2 (exact witness (a, b) = (1/2, 1)); Disc(S) != 0 is a hypothesis for specific parameters | same | `reduction.dphi_dy.*`, `reduction.sextic_squarefree.*` | exact (squared form, specialised gcd) + analytic (sign) | quick |
| `uniform.turning_root` | the turning point y_q is a root of S | same | `reduction.turning_root.*` | exact | quick |
| `uniform.decomposition` | N/(1 - y^2) = c.(1, y) + aP/(y - 1) + aM/(y + 1); aP of the t-curve carries (3a - 1) | same | `reduction.decomposition.*`, `reduction.buchdahl_factor` | exact | quick |
| `uniform.buchdahl_branch_point` | at mu = 8/9, P_t(1) = 0 and P_t'(1) = -32: the two points over the centre merge into a simple branch point, where dy/w is regular in the local parameter sqrt(y - 1) and N_t/(1 - y^2) = D/(1 + y), so the angular differential is regular there (at fixed q); the vanishing of aP alone would not suffice | same | `reduction.buchdahl_branch_point`, `reduction.buchdahl_angular_regular` | exact + analytic | quick |
| `uniform.reduced_aperture` | the reduced y-integral equals the x-quadrature | same | `reduction.half_aperture_routes.*` | numerical check | quick |
| `uniform.no_obvious_splitting` | no Moebius involution permutes the six roots (samples) | same | `reduction.no_obvious_elliptic_splitting` | diagnostic | quick |

## Clairaut function and accessibility

| Claim | Statement | Script | Checks | Evidence | Profile |
|---|---|---|---|---|---|
| `monotonicity.threshold` | J_t' > 0 inside the star for mu < mu_star; at mu_star J_t is still strictly increasing, with a single double zero of J_t' at q_star; two simple zeros (a local maximum x_M, then a local minimum x_m) for mu > mu_star. mu_star = 1 - A, A the root in (1/9, 1) of 3645A^3 + 3888A^2 - 1080A + 64 | `symbolic/clairaut_monotonicity.wls` | `monotonicity.*` | exact (counts, discriminant, gcd) + analytic (continuity of the root count between thresholds) | quick |
| `monotonicity.accessibility` | a turning radius q is regular iff K_q(x) = J(x)^2 - b^2 > 0 for q < x < 1 and K_q'(q) > 0 (pointwise criterion, distinct from Disc(S) != 0); hence every q is regular for tau and for t below mu_star; at mu_star every q except q_star (where the aperture diverges, K_q has a cubic zero); above mu_star exactly q in (0, q0) U (x_m, 1), with J(q0) = J(x_m) | — | — | analytic consequence of `monotonicity.threshold` | — |
| `monotonicity.proper` | J_tau increasing for all mu < 8/9 | same | `monotonicity.proper_time` | exact | quick |

## Aperture

| Claim | Statement | Script | Checks | Evidence | Profile |
|---|---|---|---|---|---|
| `lemma.third_kind` | Ihat(p) = (2/w(p)) Sum_i [etahat_i u_i(e->p) - omegahat_i r_i(e->p)] with Baker's R_1, R_2 | `symbolic/aperture_closed_form.wls`, library `third_kind_lemma.wl` | `aperture.baker_*`, `aperture.bidifferential_identity` (exact); `aperture.third_kind.*` (numerical check); `aperture.negative_control.*` | exact + numerical check + analytic (cycle integration step) | full |
| `aperture.closed_form` | Delta/2 = b sqrt(mu)[c.omegahat + aP Ihat(1) + aM Ihat(-1)] | same | `aperture.half_aperture.*` | numerical check | full |
| `aperture.domain` | no root of S inside (yR, yq) on the samples | same | `aperture.domain.*` | diagnostic | full |
| `aperture.surface_asymptotics` | Delta(q) = (pi/a)(1 - q)[1 + kappa (1 - q) + O((1 - q)^2)], kappa_tau = -7mu/(8a^2), kappa_t = mu/(8a^2); Delta extends analytically to q = 1 and Delta'(q) -> -pi/a; Newtonian limit Delta = pi(1 - q) | `symbolic/aperture_endpoints.wls` | `endpoints.surface_coefficient.*` (exact); `endpoints.surface.*`, `endpoints.surface_second_order.*` (numerical check) | exact + numerical check + analytic (factorised integrand analytic up to both ends) | full |
| `aperture.interior_divergence` | above mu_star the coordinate-time aperture diverges logarithmically at both inner ends x_m and q0 of the regular set, with the same coefficient 2 J(x_m)/(x_m y_m sqrt(F''(x_m)/2)), F = J_t^2 | same | `endpoints.two_critical_points` (exact); `endpoints.order_of_critical_points`, `endpoints.log_coefficient_at_*` (numerical check) | exact + numerical check + analytic | full |
| `aperture.threshold_divergence` | at mu_star, Delta ~ (2 W/sqrt(B)) |q - q_star|^(-1/2) I, W = J(q_star)/(q_star y(q_star)), B = F'''(q_star)/6, I = Beta(1/6, 1/2)/3 on the outer side and [Beta(1/3, 1/2) + Beta(1/3, 1/6)]/3 on the inner side | same | `endpoints.power_law.*` | numerical check + analytic | full |

## Centre expansion and thresholds

| Claim | Statement | Script | Checks | Evidence | Profile |
|---|---|---|---|---|---|
| `centre.expansion` | Delta/2 = pi/2 + C b + O(b^3), C = -sqrt(mu) Int_a^1 y (1/n)' dy/sqrt(1 - y^2) | `symbolic/centre_slope_threshold.wls` | `threshold.by_parts.*`, `threshold.centre_slope.*` | numerical check + analytic | quick |
| `centre.signs` | (1/n_tau)' = 4a^2/((3a - y)^2 sqrt(Q)), so C_tau < 0; (1/n_t)' = (y^2 - 6ay + 7a^2)/(2 sqrt(Q)), so C_t < 0 for mu <= mu0 = 1 - 1/(11 - 6 sqrt 2) | same | `threshold.derivative_*`, `threshold.sufficient_sign_bound` | exact | quick |
| `centre.coordinate_closed_form` | C_t = (sqrt(2mu/a)/8)[(1 + 3a + 6a^2)K(m) - 12a^2 E(m) - (1 + a)Pi(n|m)], m = (1 - a)(1 + 5a)/(8a), n = (1 - a)/2 | `symbolic/centre_slope_closed_forms.wls` | `centre_slope.legendre_*`, `centre_slope.w2_identity`, `centre_slope.hermite_coordinate`, `centre_slope.coordinate_formula` (exact); `centre_slope.coordinate_vs_quadrature.*` | exact + numerical check | quick |
| `centre.proper_closed_form` | C_tau = sqrt(2a mu)/((3a - 1)(1 + 3a)^2)[a(1 + 9a)K(m) - 6a(1 + 3a)E(m) - (1 + a)Pi(N|m)], N = (1 + 3a)(1 - a)/(4a); bracket negative | same | `centre_slope.partial_fractions_proper`, `centre_slope.proper_formula`, `centre_slope.proper_bracket_bound` (exact); `centre_slope.proper_vs_quadrature.*`; `centre_slope.negative_control.*` | exact + numerical check | quick |
| `centre.proper_buchdahl` | C_tau = -(2pi/9) delta^(-3/2)(1 + delta/8 + O(delta^(3/2))), delta = 3a - 1 | same | `centre_slope.proper_buchdahl_asymptotics` | numerical check + analytic | quick |
| `centre.mu_delta` | C_t has a unique simple zero mu_Delta in (mu0, 21/25); a bisection on the closed form brackets it in an interval of width 2^-80/10 (the JSON records the bracket and the values of C_t at its ends; the digits of the midpoint are not certified) | `symbolic/centre_slope_threshold.wls` | `threshold.weighted_form`, `threshold.p_*`, `threshold.weight_*`, `threshold.arcsine_mean` (exact); `threshold.mu_delta_bracket` (numerical check) | exact + numerical check + **analytic: covariance (Chebyshev) step** | quick |
| `centre.proper_bound` | 0 < Delta_tau < pi for all q and mu < 8/9 (L_tau > 1) | same | `threshold.proper_l_identity`, `threshold.proper_l_bound` | exact | quick |
| `centre.proper_not_monotone` | Delta_tau is not monotone near the Buchdahl limit | same | `threshold.proper_not_monotone_near_buchdahl` | numerical check + analytic | quick |
| `profile.diagnostics` | Delta_tau decreasing on a grid for mu <= 0.88; Delta_t monotone below and not monotone above mu_Delta on 108 radii | `diagnostics/aperture_profile.wls` | `profile.*` | diagnostic | full |
| `aperture.monotonicity_criterion` | with J' > 0, J(0) = 0, J(1-) = oo: Delta/2 = Int_0^{pi/2} Phi(b sec theta) d theta, Phi = 1/L, L = dJ/dl = y(1 + x n'/n); L_y <= 0 with isolated zeros (negative Gaussian curvature of the optical section) implies Delta strictly decreasing; L = y + (1 - y^2)[+-1/(3a - y) + (3a - y)/Q]; sign(L_y) = sign(N), N polynomial | `symbolic/aperture_monotonicity.wls` | `monotonicity_global.angular_substitution`, `.l_unified.*`, `.n_endpoints.*`, `.sign_of_l_y.*` | exact + analytic (differentiation under the integral) | full |
| `aperture.global_monotonicity_proper` | Delta_tau strictly decreasing from pi to 0 for 0 < mu <= mu_tau = 0.887085615449..., a_tau the unique root in (1/3, 1) of 1638781a^8 - 2030112a^6 + 418014a^4 - 30132a^2 + 729 (interior double zero at the threshold); independent positivity certificate for mu <= 55/64 | same | `.threshold.proper`, `.threshold_double_zero.proper`, `.positivity_certificate.proper` (exact); `.negative_control.*`; `.diagnostic.proper.*` | exact + diagnostic | full |
| `aperture.global_monotonicity_coordinate` | Delta_t strictly decreasing from pi to 0 for 0 < mu <= 16/25: pointwise criterion up to mu_t = 0.519590683381... (a_t the unique root in (1/3, 1) of 15a^3 + 5a^2 - 15a + 3), integral criterion (single sign change of H = Phi'/j, H decreasing where positive, C_t <= 0) up to 16/25 | same | `.threshold.coordinate`, `.positivity_certificate.coordinate`, `.extension.*` (exact); `.extension.centre_slope` (numerical check); `.diagnostic.coordinate.*` | exact + numerical check + analytic | full |
| `profile.coordinate_conjecture` | Delta_t monotone on (0, 1) for 16/25 < mu < mu_Delta | — | — | **conjecture** | — |

## Whole curve and travel times

| Claim | Statement | Script | Checks | Evidence | Profile |
|---|---|---|---|---|---|
| `curve.closed_form` | phi(y) = b sqrt(mu)[c0 u0(y) + c1 u1(y) + aP Pi+(y) + aM Pi-(y)]; ascent branch Delta - phi(y) | `symbolic/whole_curve.wls` | `whole_curve.vs_ode.*`, `whole_curve.negative_control.*` | numerical check (independent Hamiltonian ODE) | full |
| `times.elements` | on the descent (y increasing, integration from a to y_q, w > 0): n_tau dl = (3a - y)^2 dy/(sqrt(mu) w), n_t dl = 16 dy/(sqrt(mu)(3a - y) w); d tau = n_tau dl, dt = a n_t dl | `symbolic/travel_times.wls` | `travel_times.time_element.*` | exact (squared form) + analytic (orientation) | full |
| `times.closed_form` | t_half = -(16a/sqrt(mu)) Ihat(3a); tau_half = (1/sqrt(mu))[9a^2 omegahat_0 - 6a omegahat_1 + Ihat_tilde(0)] on v = 1/y | same | `travel_times.coordinate.*`, `travel_times.proper.*`, `*_vs_ode.*`, `travel_times.negative_control.*` | numerical check | full |
| `times.ordering` | tau_half < t_half at equal q | same | `travel_times.proper_less_than_coordinate.*` | numerical check + analytic | full |
| `times.newtonian_limit` | t_half ~ tau_half ~ (pi/2) sqrt(2/mu) sqrt(1 - q^2); total time ~ pi sqrt(R^3/(GM)) sqrt(1 - q^2) | same | `travel_times.newtonian_limit.*` | numerical check + analytic | full |

## Analytic steps not executed by the suite

- The cycle-integration step of `lemma.third_kind` (symmetry of the bidifferential, vanishing of exact terms on the
  cycle, G(E) = 0 at a branch point).
- The covariance (Chebyshev) inequality in the uniqueness proof of `centre.mu_delta`.
- The asymptotic arguments behind `centre.expansion` (remainder O(b^3)), `centre.proper_buchdahl`,
  `centre.proper_not_monotone`, `times.ordering` and `times.newtonian_limit`.

These will be written out in the paper; the suite checks the identities and the numbers they rest on.
