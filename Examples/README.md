# RNA–RNA interaction examples

This directory contains a set of biologically motivated example runs for
RIKin, alongside the scripts used to reproduce them and the plots shown in
the top-level [`README.md`](../README.md).

## Contents

* **FASTA files** — input sequences for each example (`exampleA.fasta`,
  `exampleB.fasta`, `HP1.fasta`, `HP2.fasta`, `HP3.fasta`, `MicA.fasta`,
  `OmpA.fasta`).
* **Config files** (`*.cfg`) — per-example JSON overrides, deep-merged on
  top of `rikin_pipeline.cfg` (see the main README's
  [Configuration file](../README.md#configuration-file) section for the
  schema).
* **`run-all.sh`** — runs `rikin_pipeline.py` for every example below, in
  one go.
* **`plot_kinetics.py`** — a [Jupytext](https://jupytext.readthedocs.io/)
  notebook (light-format `.py`, pairs with a `.ipynb`) that loads each
  example's results with `rikinplotlib.RikinRun` and regenerates the
  figures used throughout this repository.
* **`Plots/`** — the figures (SVG + PDF) produced by `plot_kinetics.py`,
  including the ones embedded in the top-level README.

## Examples

### Toy example

A short synthetic RNA pair, used as the minimal illustrative example in the
main README:
```bash
rikin_pipeline.py -o example --seqA AAAGGGGGGAAAAAAAGGGUGGGAAAAAAAGGGCGGGAAA --seqB CCCGCCC
```
or equivalently, from the FASTA files in this directory:
```bash
rikin_pipeline.py -o example --fastaA exampleA.fasta --fastaB exampleB.fasta
```

### Kissing-hairpin (KHP) system

Two designed hairpins that interact via a "kissing" loop–loop interaction,
run against two different partner hairpins:

```bash
rikin_pipeline.py --fastaA HP1.fasta --fastaB HP2.fasta -c khp.cfg -o HP1-HP2
rikin_pipeline.py --fastaA HP3.fasta --fastaB HP2.fasta -c khp.cfg -o HP3-HP2
```

| HP1:HP2 | HP3:HP2 |
|---|---|
| ![HP1-HP2 kinetics](Plots/HP1-HP2_kinetics.svg) | ![HP3-HP2 kinetics](Plots/HP3-HP2_kinetics.svg) |

### MicA–MicA homodimer

*E. coli* MicA is a small regulatory sRNA; here its self-interaction
(homodimerization) is analyzed:
```bash
rikin_pipeline.py --fastaA MicA.fasta --fastaB MicA.fasta -c MicA-MicA.cfg -o MicA-MicA
```
![MicA-MicA kinetics](Plots/micA_micA_kinetics.svg)

### MicA–OmpA heterodimer

MicA's biological target: the 5′UTR of the *ompA* mRNA (translation of
which MicA represses upon binding):
```bash
rikin_pipeline.py --fastaA MicA.fasta --fastaB OmpA.fasta -c MicA-ompA.cfg -o MicA-OmpA
```
![MicA-OmpA kinetics](Plots/micA_ompA_kinetics.svg)

## Reproducing everything

`run-all.sh` runs all of the above in sequence:
```bash
bash run-all.sh
```

> **Note:** at the time of writing, `run-all.sh` in this directory contains
> the example commands as plain text rather than a runnable script (it's
> effectively a copy of the command list above). Until it's fixed, either
> run the commands shown above directly, or use the corrected version:
> ```bash
> #!/bin/bash
> set -euo pipefail
>
> rikin_pipeline.py -o example --fastaA exampleA.fasta --fastaB exampleB.fasta
>
> rikin_pipeline.py --fastaA HP1.fasta --fastaB HP2.fasta -c khp.cfg -o HP1-HP2
> rikin_pipeline.py --fastaA HP3.fasta --fastaB HP2.fasta -c khp.cfg -o HP3-HP2
>
> rikin_pipeline.py --fastaA MicA.fasta --fastaB MicA.fasta -c MicA-MicA.cfg -o MicA-MicA
> rikin_pipeline.py --fastaA MicA.fasta --fastaB OmpA.fasta -c MicA-ompA.cfg -o MicA-OmpA
> ```

## Regenerating the plots

Once the runs above have produced their output directories, `plot_kinetics.py`
loads each one and renders the full set of figures (state-probability
kinetics, top-state structures, interaction dotplots, and per-nucleotide
paired-probability kinetics) into `Plots/`.

Since it's kept in [Jupytext](https://jupytext.readthedocs.io/) light
format rather than as a plain script, open it as a notebook:
```bash
pip install jupytext   # if not already available
jupytext --to notebook plot_kinetics.py   # produces plot_kinetics.ipynb
jupyter notebook plot_kinetics.ipynb
```
or, if you have the Jupytext extension enabled in Jupyter/JupyterLab, you
can open `plot_kinetics.py` directly as a notebook without a separate
conversion step.

> We're considering committing the paired `.ipynb` alongside `plot_kinetics.py`
> directly, to make this a one-step "open and run" without requiring Jupytext
> as a prerequisite — check whether `plot_kinetics.ipynb` is present in this
> directory, in which case you can skip straight to `jupyter notebook plot_kinetics.ipynb`.

The notebook currently has some machine-specific absolute paths for
locating each example's `input_directory` (e.g.
`/home/will/Research/Projects/...`) that will need adjusting to point at
wherever you ran the corresponding `rikin_pipeline.py` command above.
