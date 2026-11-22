# Problem 6 – Kalman filters (self-contained)

This folder is self-contained. It does not change the group repository (`MPC-me-dev-ch4-5`).

## How to run (MATLAB R2024b)
1. Open this folder in MATLAB.
2. Run `testKalmanFilter.m`. It runs 12 checks of the generic Kalman filter, and all should report PASS.
3. Run `main_problem6.m`. It creates 18 figures and `problem6_results.csv` in `figures/`.
   - Optional: to use the Problem 4 fit directly instead of the Table 5.3 values, run the repository's `problem4_main.m` first. `main_problem6.m` then takes `estimated_params` from the workspace.

## The generic Kalman filter
`kalmanFilter(model, y, u, 'dynamic' | 'static', x0, P0)` works for **any** linear model given as a struct:

```matlab
model.A, model.C, model.Q, model.R      % required
model.B, model.G, model.S, model.P0     % optional (default: no input, I, 0, I)
```
It implements the general form from Lecture 6, slides 19–23, including correlated noise S. `kalmanStatic(model)` solves the Riccati equation (slide 29). It uses `idare`/`dare` from the Control System Toolbox and always also runs the fixed-point iteration as a check.

## Files
| File | Content |
|---|---|
| `main_problem6.m` | Problem 6.1–6.4: models, simulations, figures, results table |
| `testKalmanFilter.m` | Verification of the generic filter (Riccati residual, idare vs iteration, dynamic → static, S = 0 vs textbook filter, Monte Carlo consistency, four-tank models) |
| `functions/kalmanFilter.m`, `kalmanStatic.m` | The generic static and dynamic Kalman filter |
| `functions/buildKalmanModelP5.m` | Problem 5 model: 6 states, or 8 states for the offset-free version (6.12–6.13) |
| `functions/buildKalmanModelP4.m` | Problem 4 model: 4 states with shared T1 and T2, or 6 states for the offset-free version |
| `functions/sharedTimeConstants.m` | Converts `estimated_params` from `problem4_main.m` to shared T1 and T2 |
| `functions/simulateKalmanPlant.m`, `createKalmanNoise.m` | True linear or nonlinear plant with SDE disturbances. Both plants use the same noise. |
| `functions/truthInModelStates.m` | True trajectory expressed in each model's states, for plotting |
| `functions/innovationStats.m`, `isDetectable.m`, `findResult.m` | Innovation checks (mean, NIS, whiteness), PBH detectability test, result lookup |
| `functions/plotKalmanStates.m`, `plotKalmanGains.m`, `plotInnovations.m`, `saveFigure.m`, `graphicsArray.m` | Figures |
| `functions/createParameters.m`, `createSimSettings.m`, `computeSteadyState.m`, `computeJacobians.m`, `computeAugmentedModel.m`, `c2dsde.m`, `c2dzoh.m`, `quadrupleTankProcess.m`, `systemSense.m`, `systemOutput.m` | **Unchanged copies** from the group repository (`MPC-me-dev-ch4-5/functions`, 21 Nov 2026) |

## Tuning (design choices – state them in the report)
- Problem 4 noise on x3 and x4: σw = σF/√(2aF)/(1−γ), the stationary standard deviation of the SDE disturbance in pump-flow units.
- Integrators: σN = 0.3 for the Problem 5 model and σN = 3 for the Problem 4 model. A larger σN follows steps faster, but the estimates are noisier.
- P0: 5 cm on the levels, 20–50 cm³/s on the flow-type states. The plant starts 3 cm off the steady state, while the filter starts at x̂ = 0.
