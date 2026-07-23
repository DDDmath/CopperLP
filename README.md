# CopperLP

This repository contains the SageMath research scripts accompanying our work on computing asymptotic bounds for the automated Coppersmith method via linear programming.

The repository is intended to provide the experimental implementations used in the paper. Each script can be inspected and executed independently. Experimental parameters are configured manually near the beginning of the corresponding script.

## Included Experiments

- `CIHNP_CSURF.sage`: CIHNP-CSURF experiments.
- `MIHNP.sage`: MIHNP experiments with different numbers of samples.
- `ECHNP.sage`: ECHNP experiments.
- `LCG.sage`: experiments for truncated linear congruential generators.
- `LIPH_POKE.sage`: LIPH-POKE experiments.

## Requirements

The experiments require:

- SageMath 10.8 or later;
- sufficient memory for lattice reduction and Gröbner-basis computation.

The scripts use SageMath's built-in lattice-reduction and Gröbner-basis implementations by default.

The optional `msolve` backend can be used to accelerate Assumption 1 verification when it is available in the SageMath environment. It is disabled by default and is not required to run the scripts.

## Running an Experiment

Download or clone the repository and enter the repository directory.

Run an individual experiment with SageMath. For example:

```bash
sage CIHNP_CSURF.sage
```

The other scripts can be executed in the same way:

```bash
sage MIHNP.sage
sage ECHNP.sage
sage LCG.sage
sage LIPH_POKE.sage
```

## Adjusting Parameters

Before running an experiment, open the corresponding `.sage` file and modify the configuration block near the beginning of the script.

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
- `M`: lattice-construction or polytope scale;
- `TH`: polytope parameters;
- `DELTA`: LLL reduction parameter;
- `SEED`: optional random seed;
- `VERIFY_ASSUMPTION1`: whether to verify Assumption 1;
- `USE_MSOLVE_ASSUMPTION1`: whether to use the optional `msolve` backend.

After modifying the parameters, save the file and execute it with SageMath.

## Experimental Parameters

The parameter combinations corresponding to the experiments reported in Table 2 are recorded separately in:

```text
profiles/table2_profiles.json
```

The SageMath scripts do not automatically read this file. To reproduce a particular parameter setting, manually copy the desired values into the configuration block of the corresponding script.

## Notes

The scripts generate synthetic experimental instances. Running times may vary depending on the processor, available memory, SageMath version, lattice-reduction backend, and Gröbner-basis implementation.

The current version information is recorded in the `RELEASE` file.
