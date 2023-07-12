# Muscle-tendon-Wrapping
An effective muscle wrapping method based on tangency constraints that guarantee continuity along a geodesic path over parametric surfaces.

## Simulations

A simple but representative muscle wrapping was simulated under two cases: (1) The muscle-tendon unit remains in contact with the cylindrical surface and (2) The muscle-tendon unit transitions from a wrapping state to a non-wrapping state.

The numerical efficiency is investigated considering the classical Newton-Raphson and two Quasi-Newton (Broyden–Fletcher–Goldfarb–Shanno and Broyden) methods

- sim_NewtonRaphson_Case1.m: Newton-Raphson method implementation to solve the muscle wrapping problem under case 1
- sim_NewtonRaphson_Case2.m: Newton-Raphson method implementation to solve the muscle wrapping problem under case 2
- 
- sim_BFGS_Case1.m: BFGS method implementation to solve the muscle wrapping problem under case 1
- sim_BFGS_Case2.m: BFGS method implementation to solve the muscle wrapping problem under case 2
- 
- sim_BR_Case1.m: Broyden's method implementation to solve the muscle wrapping problem under case 1
- sim_BR_Case2.m: Broyden's method implementation to solve the muscle wrapping problem under case 2

> **Note:** Media files include videos only for Broyden's method solutions, for case 1, case 2, and an additional video of a simulation considering the surface translational and rotational motion.
