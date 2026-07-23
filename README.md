# CopperLP

This repository contains the SageMath research scripts accompanying our work on computing asymptotic bounds for the automated Coppersmith method via linear programming.

The repository includes an implementation of Algorithm 1 for computing asymptotic coefficients, together with research scripts for inspecting and reproducing the experiments reported in the paper. Each experiment is executed independently, and its parameters are configured manually near the beginning of the corresponding SageMath script.

## Algorithm 1

- `Algorithm_1.sage`: implementation of Algorithm 1 for computing the asymptotic coefficients and the resulting asymptotic small-root bound from a polynomial system and a rational polytope.

Run Algorithm 1 with:

```bash
sage Algorithm_1.sage
```

## Included Experiments

- `CIHNP_CSURF.sage`: CIHNP-CSURF experiments.
- `MIHNP.sage`: MIHNP experiments with different numbers of samples.
- `ECHNP.sage`: ECHNP experiments.
- `LCG.sage`: truncated linear congruential generator experiments.
- `LIPH_POKE.sage`: LIPH-POKE experiments.

## Requirements

The scripts require:

- SageMath 10.8 or later;
- sufficient memory for polyhedral computation, lattice reduction, and Gröbner-basis computation.

The experimental scripts use SageMath's built-in lattice-reduction and Gröbner-basis implementations by default.

The optional `msolve` backend can be used to accelerate Assumption 1 verification when it is available in the SageMath environment. It is disabled by default and is not required to run the scripts.

## Running the Scripts

Clone or download the repository and enter its directory.

Run Algorithm 1 with:

```bash
sage Algorithm_1.sage
```

Run an individual experiment with SageMath:

```bash
sage CIHNP_CSURF.sage
```

The other experimental scripts are executed in the same way:

```bash
sage MIHNP.sage
sage ECHNP.sage
sage LCG.sage
sage LIPH_POKE.sage
```

## Adjusting Parameters

Before running an experimental script, modify the configuration block near the beginning of that script.

For example:

```python
N_RUNS = 1
PBITS = 256
UBITS = 29
M = 2
SEED = None
VERIFY_ASSUMPTION1 = True
USE_MSOLVE_ASSUMPTION1 = False
```

The available parameters depend on the experiment. Common parameters include:

- `N_RUNS`: number of independent runs;
- `PBITS`: modulus size in bits;
- `UBITS`, `X1_BITS`, or `X2_BITS`: bounds for the unknown values;
- `NSAMPLES`: number of samples;
- `M`: lattice-construction or polytope scale;
- `TH`: polytope parameters;
- `DELTA`: LLL reduction parameter;
- `SEED`: optional random seed;
- `VERIFY_ASSUMPTION1`: whether to verify Assumption 1;
- `USE_MSOLVE_ASSUMPTION1`: whether to use the optional `msolve` backend.

To reproduce a particular row of Table 2, manually enter the corresponding parameters from the paper into the configuration block of the relevant script.

## Notes

The experimental scripts generate synthetic experimental instances. Running times may vary depending on the processor, available memory, SageMath version, lattice-reduction backend, and Gröbner-basis implementation.

Version information is recorded in the `RELEASE` file.
