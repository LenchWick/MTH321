# Financial SDE Presentation Script

The script follows 17 main slides. The detailed confidence-interval and weak-bias derivation slides have been removed. Rehearse the remaining talk against the course time limit. The slide text is English because the project brief requires an English presentation. The explanations below are written as a continuous talk for one presenter. The group can divide the slides after reviewing the full story.

## Slide 1 — Financial SDEs

Good morning. Our project studies two financial stochastic differential equations. GBM gives us an exact pathwise benchmark, so it tells us whether the numerical methods are implemented correctly. The non-affine stochastic-volatility model removes that exact solution and forces us to build and test a numerical reference. That change is the main story of the presentation.

## Slide 2 — Anatomy of an SDE

An SDE has two parts. The drift is the conditional mean rate of change. The diffusion determines the conditional variance over a short interval. Brownian increments are proportional to the square root of the time step, rather than the time step itself. This scaling is why stochastic methods can have strong orders such as one half. It also explains why we need Itô calculus instead of ordinary calculus.

## Slide 3 — Itô’s formula and the exact GBM solution

The boxed equation at the top is Itô’s formula. It includes a second-derivative term because Brownian motion has nonzero quadratic variation. For GBM, a is mu S and b is sigma S.

On the left, choose f equal to log S. Its first derivative is one over S and its second derivative is minus one over S squared. Substituting them into Itô’s formula gives drift mu minus one half sigma squared and diffusion sigma.

On the right, integrate this equation and exponentiate to obtain the exact GBM solution. This also gives the lognormal terminal law and positive prices. Here mu is the expected growth parameter under our chosen probability measure.

## Slide 4 — Itô–Taylor expansion: which terms are retained?

The top equation shows the scalar Itô–Taylor terms used by these methods. EM retains the drift term and the first Brownian term. Milstein retains the same terms plus the next stochastic correction, which contains b times its derivative. The remainder collects omitted terms. The two boxed formulas make this difference explicit.

## Slide 5 — Log-price GBM: why EM and Milstein coincide

For log-price GBM, b equals sigma, a constant. Therefore, its derivative is zero. Substituting this into Milstein makes its additional correction zero, so EM and Milstein give the same exact grid update.

For the original price equation, b equals sigma times S and its derivative is sigma, not zero. The correction remains. We therefore use the original price equation to compare convergence orders.

## Slide 6 — Strong and weak convergence

Strong and weak convergence answer different questions. Strong error compares two terminal values driven by the same Brownian path. It matters when the path itself matters. Euler–Maruyama has strong order one half, while scalar Milstein has strong order one under the usual assumptions.

Weak error compares expectations of an observable. Both methods have weak order one for sufficiently smooth observables. For the mean of GBM, we can prove this directly. Both schemes have the discrete expectation S zero times one plus mu h to the power T over h. Expanding its logarithm gives a first-order time-step bias.

The bias can become smaller than Monte Carlo noise. We therefore use paired differences, confidence intervals, and enough independent paths before making a slope claim.

## Slide 7 — Strong convergence: from paths to the order

Strong convergence asks whether the numerical solution approaches the exact solution on the same Brownian path. Fix the final time T and choose a step h equal to T divided by N. Generate M independent paths, with each Brownian increment equal to square root h times a standard normal draw.

Starting from S zero, iterate EM or Milstein to the final time. Use the sum of those very same increments in the exact GBM formula. For each path, take the absolute difference between numerical and exact terminal prices. Average these absolute differences to estimate the terminal L1 strong error.

Repeat for h over two, h over four, and further grids. Coarse increments must be sums of fine increments so that grids share the same trajectory. If strong error behaves like C times h to the power p, dividing errors at h and h over two gives two to the power p. The logarithm of that ratio divided by log two estimates the strong convergence order. We fit several grids in the experiment. An empirical slope supports a theoretical rate; it is not a convergence proof.

## Slide 8 — Weak convergence: from expectations to the order

Weak convergence asks whether expectations of test functions approach their exact values. We choose the mean-price observable, phi of s equal to s. Using the endpoints already simulated, average signed numerical-minus-exact differences first. The absolute value of that average estimates weak error. Unlike strong error, positive and negative path differences can cancel.

For GBM we can also derive the mean bias analytically. The Brownian increment has mean zero, and delta W squared minus h has mean zero. Thus both EM and Milstein have expected one-step multiplier one plus mu h. Iterating N steps gives S zero times one plus mu h raised to T over h.

Subtract the exact mean S zero times exp of mu T. Expanding the logarithm gives mu T minus one half mu squared T h plus terms of order h squared. Expanding the exponential then gives a leading bias of minus one half S zero exp of mu T times mu squared T h. For our nonzero drift, its absolute value is proportional to h. Therefore both schemes have weak order one for the mean observable. If mu is zero this mean bias vanishes exactly.

We check the simulated biases against this analytic curve and fit their log–log slope across step sizes. This establishes a mean-price test; other observables, especially nonsmooth option payoffs, need their own justification.

## Slide 9 — GBM convergence and uncertainty

The left panel shows terminal pathwise error against the time step on log–log axes. The fitted slope is 0.4975 for Euler–Maruyama and 0.9972 for Milstein, which matches the theoretical orders.

The right panel shows the signed mean bias. The analytic curve confirms first-order weak bias for both schemes. Our empirical weak slopes are 0.9610 and 0.9948. On the finest Euler–Maruyama grid, the confidence interval crosses zero. We therefore exclude that unresolved point from the slope fit. The uncertainty bars determine which claims the experiment supports.

## Slide 10 — Distribution of S(T) versus the exact lognormal law

We simulate 65,536 terminal prices with EM and Milstein at 32 steps, alongside exact GBM samples. The left plot overlays normalised histograms with the analytic lognormal density. The samples track its centre, spread and right-skewed shape. The right plot places exact quantiles on the horizontal axis and sample quantiles on the vertical axis. Agreement puts points on the diagonal. Small tail differences are easier to see in the next table. Distribution agreement does not establish accuracy on the same individual path.

## Slide 11 — Quantile comparison: centre and tails

The log price has an exact Gaussian law, so its density and quantiles can be calculated analytically. Histogram heights are counts divided by total sample size and bin width. The table compares one, five, fifty, ninety-five and ninety-nine percent quantiles. Even exact simulation differs from the population quantiles because it is a finite random sample. These checks are descriptive rather than a formal simultaneous goodness-of-fit test.

## Slide 12 — Positivity: when can an update become negative?

Starting from a positive price, the update becomes negative precisely when its multiplier is negative. For EM, solve the linear inequality for the normal draw. Its unbounded lower tail means negative updates remain possible at any positive step size with nonzero volatility.

For Milstein, complete the square. Its minimum is one half plus mu minus one half sigma squared, multiplied by h. If that minimum is nonnegative, negative updates are impossible. Otherwise, normal draws between the two roots produce negative prices. Our baseline Milstein setting is positive almost surely, but this property does not hold for all parameters and steps.

## Slide 13 — When are negative prices possible?

For EM, an update becomes negative when the standard normal draw is below minus one plus mu h, divided by sigma times square root h. Avoiding a negative update requires the draw to stay at or above this threshold. Since normal draws have no lower bound, no finite positive step guarantees nonnegative price-space EM values when volatility is nonzero.

For Milstein, negative updates are possible when its quadratic minimum, one half plus mu minus one half sigma squared times h, is negative. They occur for draws between the two roots on the preceding slide. To guarantee no negative updates, require one plus two mu minus sigma squared times h to be nonnegative. If sigma squared exceeds twice mu, choose h no greater than one divided by sigma squared minus twice mu. If sigma squared is no greater than twice mu, every positive h satisfies this condition.

With a strict inequality, every Milstein multiplier is positive. At equality, zero is possible only at one exact normal draw, which has probability zero. Simulating log price and exponentiating is a direct way to guarantee positive prices in exact arithmetic.

## Slide 14 — Non-affine stochastic volatility

The extension does not make mu random. Mu remains a constant expected return parameter. The new state Y makes the volatility random.

Kappa pulls Y toward theta, while xi controls the size of factor shocks. The logistic function g maps every real value of Y to a smooth positive volatility between 0.10 and 0.50. This avoids negative volatility and keeps the price diffusion bounded.

The Brownian correlation is minus 0.7. A negative return shock therefore tends to coincide with a positive factor shock, which raises g of Y. This represents the leverage-style relationship observed in equity markets. The model is non-affine because the logistic volatility and the factor diffusion are nonlinear functions of the state.

## Slide 15 — Why the scalar Milstein formula does not transfer

*指向顶部公式，再指向两个向量。*

The first equation constructs correlated Brownian motions from two independent sources. The diffusion vectors show how each source affects the system: B one affects both X and Y, while B two directly affects Y. Since Y determines g, a change in Y also changes the price volatility.

*指向底部 commutator。*

The nonzero commutator shows that the diffusion fields do not commute. A full Milstein method therefore needs cross stochastic integrals that cannot be obtained from the endpoint increments alone. We cannot simply apply scalar Milstein to each equation and claim order one. We use EM for the baseline.

## Slide 16 — Nested Brownian grids

*指向右侧公式，再看左图。*

Without an exact solution, we use a finer EM solution as a reference. The formula sums fine Brownian increments into coarse increments, so both solutions follow the same random path.

The left plot shows terminal strong error against step size. The three reference grids give similar curves, and the finest fitted slope is 0.5003. This supports empirical strong order one half relative to the numerical reference.

*指向右图。*

The right plot shows signed call-payoff differences. Most error bars cross zero, so we do not claim a resolved weak convergence order.

## Slide 17 — Beyond the baseline

*指向右侧公式。*

The control variate subtracts a multiple of the centred stock price from the call payoff. These quantities are positively correlated, so the adjustment reduces sampling variance. We estimate beta on an independent pilot sample and obtain a variance reduction of 6.23 times. It does not remove discretisation bias.

*指向左边两幅图。*

The balanced factor update gives similar price errors to EM in the left plot. Its benefit appears in the right plot: the long-time factor variance stays close to the exact target as the step increases.

*指向右下角步数。*

The GBM cost comparison uses the same error target: EM needs 256 steps and Milstein needs 16. This efficiency result concerns GBM.

## Technical corrections to the original outline

- The prescribed non-affine model keeps \(\mu\) constant. It makes volatility stochastic through \(Y_t\) and \(g(Y_t)\).
- \(W^{(1)}\) and \(W^{(2)}\) are correlated with \(\rho=-0.7\). Two independent base Brownian motions generate them.
- Cross iterated stochastic integrals are mathematically meaningful. The obstacle is that ordinary endpoint increments do not determine Lévy area when diffusion fields do not commute.
- The scalar Milstein formula cannot simply be applied to each component. This would omit cross terms. A full multidimensional Milstein method remains possible with extra simulation work.
- Randomness is part of the model. Numerical error means the difference between an approximation and the exact or reference solution. It should not be described as all simulated values outside a normal range.
- The control variate does not reconstruct GBM. It uses the known value \(\mathbb E[S_T]=S_0e^{\mu T}\) to reduce the sampling variance of the expected call payoff.
- \(\mathbb E[(S_T-K)^+]\) is an expected payoff. It becomes an option price only after specifying the risk-neutral measure and discounting.

## Optional four-person split

After the group reviews the complete talk, a natural split is:

1. Presenter 1: Slides 1–3, SDE definition and exact GBM solution.
2. Presenter 2: Slides 4–8, Itô–Taylor methods and convergence definitions.
3. Presenter 3: Slides 9–14, GBM experiments, distribution, positivity, and non-affine model.
4. Presenter 4: Slides 15–17, multidimensional noise, nested validation, and extensions.

Each person should understand the full presentation because the rubric uses cross-questioning during the viva.
