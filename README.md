# CopperLP

This repository contains the SageMath implementations and reproducibility
artifacts for the Coppersmith experiments reported in our work.

## Repository Contents

- `CIHNP_CSURF.sage`: experiments for CIHNP-CSURF.
- `MIHNP.sage`: experiments for MIHNP.
- `ECHNP.sage`: experiments for ECHNP.
- `LCG.sage`: experiments for truncated linear congruential generators.
- `LIPH_POKE.sage`: experiments for LIPH-POKE.
- `run_all_table2.sh`: runs all five experiment scripts sequentially.
- `profiles/table2_profiles.json`: parameter profiles for the Table 2 experiments.
- `instances/`: fixed test instances used for reproducibility.
- `logs/table2_results.csv`: results corresponding to the Table 2 experiments.
- `logs/stress_test_results.csv`: additional stress-test results.
- `RELEASE`: release and version information.

## Requirements

The experiments require:

- SageMath;
- Bash;
- sufficient memory for lattice reduction and Gröbner-basis computation.

Some optional verification procedures may require additional external
software when they are enabled in the corresponding SageMath script.

## Running All Experiments

Download or clone the repository, enter the repository directory, and run:

```bash
bash run_all_table2.sh
```

The scripts are executed in the following order:

1. CIHNP-CSURF;
2. MIHNP;
3. ECHNP;
4. LCG;
5. LIPH-POKE.

The runner prints the experimental results directly to the terminal.

## Running an Individual Experiment

Each experiment can also be executed separately. For example:

```bash
sage CIHNP_CSURF.sage
```

The other experiments can be run with:

```bash
sage MIHNP.sage
sage ECHNP.sage
sage LCG.sage
sage LIPH_POKE.sage
```

## Parameter Configuration

The configurable parameters are defined at the beginning of each SageMath
script and can be adjusted manually before execution.

The file:

```text
profiles/table2_profiles.json
```

is intended to record the exact parameter combinations used for the
experiments reported in Table 2. The current SageMath scripts do not
automatically load this JSON file. Users should manually copy the desired
parameter values into the configuration section at the beginning of the
corresponding script before running it.

## Fixed Instances

The `instances/` directory contains fixed synthetic test instances for the
five application families:

```text
instances/
├── fixed_instance_cihnp_csurf.json
├── fixed_instance_mihnp.json
├── fixed_instance_echnp.json
├── fixed_instance_lcg.json
└── fixed_instance_liph_poke.json
```

The fixed instances are included to make the experiments deterministic and
independently reproducible.

## Experimental Logs

Reference results are stored in:

```text
logs/table2_results.csv
logs/stress_test_results.csv
```

Running times depend on the processor, available memory, SageMath version,
lattice-reduction backend, and Gröbner-basis implementation. Therefore,
runtime values may differ across machines even when the recovered roots and
success indicators remain unchanged.

## Repository Structure

```text
CopperLP/
├── README.md
├── RELEASE
├── run_all_table2.sh
├── CIHNP_CSURF.sage
├── MIHNP.sage
├── ECHNP.sage
├── LCG.sage
├── LIPH_POKE.sage
├── profiles/
│   └── table2_profiles.json
├── instances/
│   ├── fixed_instance_cihnp_csurf.json
│   ├── fixed_instance_mihnp.json
│   ├── fixed_instance_echnp.json
│   ├── fixed_instance_lcg.json
│   └── fixed_instance_liph_poke.json
└── logs/
    ├── table2_results.csv
    └── stress_test_results.csv
```

## Release Information

The current release information is recorded in the `RELEASE` file.
