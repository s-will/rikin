#!/usr/bin/env python3
"""
rikin_plot.py
==================

Stand-alone command-line tool that generates kinetics analysis plots for RNAInterKin.

It reads pre-computed RNA-RNA interaction kinetics results from a result
directory and writes out the standard set of plots (state-probability
kinetics, per-state structure panels, interaction dotplot, full dotplot,
paired-probability kinetics) as SVG and PDF.

Every plot method's real keyword arguments (figsize, seqfs, wspace,
hspace, row_gap, rna_names, state_names, show_rnas, show_pairing,
swap_axes, vmax, ...) are simply forwarded through from an optional JSON
`--config` file, keyed by plot name -- so this stays in sync with the
library's actual signatures instead of guessing at a fixed CLI flag set.

Usage
-----

Single run, sequences given directly:

    python rikin_plot_cli.py \
        --input-dir /path/to/results/HP1-HP2 \
        --name HP1-HP2 \
        --seqA GGACGAGGCAUUUCCCCUUGU \
        --seqB GGACAAGGGGAAAUGCCUUGU \
        --shown-state-threshold 0.02 \
        --output-dir Figs

Single run, sequences from a FASTA file with exactly two records:

    python rikin_plot_cli.py --input-dir results/example --name example \
        --fasta example.fasta

Only generate some plots:

    python rikin_plot_cli.py ... --plots kinetics dotplot

Batch mode (reproduces the notebook's HP1-HP2 / HP2-HP3 loop) via a JSON
file describing several runs at once:

    python rikin_plot_cli.py --runs-config runs.json --output-dir Figs

  where runs.json looks like:

    [
      {
        "name": "HP1-HP2",
        "input_dir": "/path/HP1-HP2",
        "seqA": "GGACGAGGCAUUUCCCCUUGU",
        "seqB": "GGACAAGGGGAAAUGCCUUGU",
        "shown_state_threshold": 0.02
      },
      {
        "name": "HP2-HP3",
        "input_dir": "/path/HP2-HP3",
        "seqA": "GGACAAGGGGAAAUGCCUUGU",
        "seqB": "GGACGAUCAGCAUUUCCCUGAUGU",
        "shown_state_threshold": 0.02
      }
    ]

Per-plot styling overrides go in an optional `--config config.json`,
keyed by plot name, and are forwarded verbatim as kwargs to the
matching RikinRun method (so anything the method accepts can go here):

    {
      "states":  {"seqfs": 6, "wspace": 0.15, "hspace": 0.6, "row_gap": 1.2},
      "dotplot": {"figsize": [13.5, 12], "barsize": 0.5},
      "paired_kinetics": {"show_rnas": [1], "suffix": "paired_kinetics_only_srna"}
    }
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt

try:
    import rikinplotlib as rkplt
    from rikinplotlib import RikinRun
except ImportError as exc:  # pragma: no cover
    sys.exit(
        "Error: could not import 'rikinplotlib' (and its dependency 'RNA' / "
        "ViennaRNA python bindings). Make sure it is installed / on your "
        f"PYTHONPATH.\n({exc})"
    )


# --------------------------------------------------------------------------
# Plot registry: name -> function(run, cfg)
# --------------------------------------------------------------------------

PLOT_REGISTRY = {}


def register(name):
    def deco(fn):
        PLOT_REGISTRY[name] = fn
        return fn
    return deco


def _kwargs_from_cfg(cfg: dict, defaults: dict) -> dict:
    """Merge user config on top of built-in defaults; tuple-ify figsize."""
    kwargs = dict(defaults)
    kwargs.update(cfg)
    if "figsize" in kwargs:
        kwargs["figsize"] = tuple(kwargs["figsize"])
    return kwargs


@register("kinetics")
def plot_kinetics(run: RikinRun, cfg: dict):
    figsize = tuple(cfg.get("figsize", (6, 6)))
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    run.plot_state_probability_kinetics(ax=ax)
    run.optsave("kinetics")  # library saves the current figure itself
    plt.close(fig)


@register("states")
def plot_states(run: RikinRun, cfg: dict):
    kwargs = _kwargs_from_cfg(cfg, {"figsize": (10, 4), "seqfs": 6})
    fig, ax = run.plot_states(**kwargs)  # saves internally (tag "states")
    plt.close(fig)


@register("interaction_dotplot")
def plot_interaction_dotplot(run: RikinRun, cfg: dict):
    figsize = tuple(cfg.get("figsize", (8, 1.75)))
    width_ratios = cfg.get("width_ratios", [30, 1])
    fig, ax = plt.subplots(1, 2, figsize=figsize, width_ratios=width_ratios)

    call_kwargs = _kwargs_from_cfg(
        {k: v for k, v in cfg.items() if k not in ("figsize", "width_ratios")},
        {"swap_axes": True, "vmax": 1},
    )
    run.interaction_dotplot(ax=ax[0], cbar_ax=ax[1], **call_kwargs)
    fig.tight_layout()
    run.optsave("interaction")
    plt.close(fig)


@register("dotplot")
def plot_dotplot(run: RikinRun, cfg: dict):
    kwargs = _kwargs_from_cfg(cfg, {"figsize": (11, 10), "seqfs": 8})
    fig, ax = run.dotplot(**kwargs)  # saves internally (tag "dotplot")
    plt.close(fig)


@register("paired_kinetics")
def plot_paired_kinetics(run: RikinRun, cfg: dict):
    kwargs = _kwargs_from_cfg(cfg, {"figsize": (12, 6), "seqfs": 7})
    fig, ax = run.plot_paired_probability_kinetics(**kwargs)  # saves internally
    plt.close(fig)


# --------------------------------------------------------------------------
# Sequence loading
# --------------------------------------------------------------------------

def read_fasta(filename: Path):
    seq = ""
    name = None
    with open(filename) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is None:
                    name = line[1:].strip().split()[0]
                else: break  # stop reading after the first sequence
            else:
                seq += line
    if name is None or name == "":
        raise ValueError(f"FASTA file {filename} empty or does not contain a valid sequence name")
    return name, seq

# --------------------------------------------------------------------------
# Output configuration (mirrors notebook's `rkplt.output_directory = "Figs"`)
# --------------------------------------------------------------------------

def configure_output(out_dir: Path, formats):
    out_dir.mkdir(parents=True, exist_ok=True)
    rkplt.output_directory = str(out_dir)

    dotted = [f.lstrip(".") for f in formats]
    dotted = [f".{f}" for f in dotted]
    # `RikinRun.save` binds `rikinplotlib.suffix` by reference as a default
    # argument, so we must mutate the existing list, not rebind the name.
    rkplt.suffix[:] = dotted


# --------------------------------------------------------------------------
# Core: run one RikinRun + its plots
# --------------------------------------------------------------------------

def run_one(run_spec: dict, plots_to_run, cfg_all: dict):
    name = run_spec["name"]
    print(f"=== {name} ===")

    run = RikinRun(
        run_spec["seqA"],
        run_spec["seqB"],
        input_directory=run_spec["input_dir"],
        shown_state_threshold=run_spec.get("shown_state_threshold", 0.02),
        autosave=name,
    )

    for plot_name in plots_to_run:
        fn = PLOT_REGISTRY[plot_name]
        cfg = cfg_all.get(plot_name, {})
        print(f"  generating '{plot_name}' ...")
        try:
            fn(run, cfg)
        except Exception as exc:
            print(f"    WARNING: '{plot_name}' failed: {exc}", file=sys.stderr)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Generate RNA-RNA interaction kinetics plots from RikinRun result directories.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    single = parser.add_argument_group("single-run options (ignored if --runs-config is given)")
    single.add_argument("--input-dir", help="Directory with RikinRun result files")
    single.add_argument("--name", help="Prefix used for output plot filenames (RikinRun's autosave name)")
    parser.add_argument("--fastaA", help="Fasta containing first sequence")
    parser.add_argument("--fastaB", help="Fasta containing second sequence")
    parser.add_argument("--seqA", help="First sequence")
    parser.add_argument("--seqB", help="Second sequence")
    single.add_argument("--shown-state-threshold", type=float, default=0.02)

    parser.add_argument(
        "--runs-config", type=Path,
        help="JSON file describing multiple runs (batch mode). See module docstring for the schema.",
    )
    parser.add_argument("--output-dir", default=".", help="Directory to write plots to (default: .)")
    parser.add_argument(
        "--formats", nargs="+", choices=["svg", "pdf", "png"], default=["svg", "pdf"],
        help="Output file formats (default: svg pdf, matching rikinplotlib's own default)",
    )
    parser.add_argument(
        "--plots", nargs="+",
        choices=list(PLOT_REGISTRY.keys()) + ["all"],
        default=["all"],
        help="Which plots to generate (default: all)",
    )
    parser.add_argument(
        "--config", type=Path,
        help="Optional JSON file with per-plot kwargs, keyed by plot name, forwarded straight "
             "to the matching RikinRun method.",
    )
    return parser


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    plots_to_run = list(PLOT_REGISTRY.keys()) if "all" in args.plots else args.plots
    cfg_all = json.loads(args.config.read_text()) if args.config else {}

    configure_output(Path(args.output_dir), args.formats)

    if args.runs_config:
        run_specs = json.loads(args.runs_config.read_text())
        if isinstance(run_specs, dict):
            run_specs = [run_specs]
    else:
        if not (args.input_dir and args.name):
            parser.error("either --runs-config, or both --input-dir and --name, are required")

    seqA = None
    seqB = None
    nameA = "seqA"
    nameB = "seqB"

    if args.seqA and args.fastaA or args.seqB and args.fastaB:
        parser.error("provide sequences via --seqA/--seqB or --fastaA/--fastaB")

    if args.fastaA:
        nameA,seqA = read_fasta(args.fastaA)
    else:
        seqA = args.seqA

    if args.fastaB:
        nameB,seqB = read_fasta(args.fastaB)
    else:
        seqB = args.seqB

    if not seqA or not seqB:
        parser.error("must provide both sequences via --seqA/--seqB or --fastaA/--fastaB")
    if set(seqA) - set("ACGUacgu"):
        parser.error(f"seqA contains invalid characters: {set(seqA) - set('ACGUacgu')}")
    if set(seqB) - set("ACGUacgu"):
        parser.error(f"seqB contains invalid characters: {set(seqB) - set('ACGUacgu')}")
        
    run_specs = [{
        "name": args.name,
        "input_dir": args.input_dir,
        "seqA": seqA,
        "seqB": seqB,
        "nameA": nameA,
        "nameB": nameB,
        "shown_state_threshold": args.shown_state_threshold,
    }]

    for run_spec in run_specs:
        run_one(run_spec, plots_to_run, cfg_all)

    print(f"\nDone. Plots written to {Path(args.output_dir).resolve()}/")


if __name__ == "__main__":
    main()
