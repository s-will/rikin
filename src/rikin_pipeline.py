#!/usr/bin/env python3
"""
rikin_pipeline.py
==================

Runs the RNAInterKin pipeline enumerate -> barriers -> prune -> solve master equation and
finally calls rikin_plot.py to produce the full set of kinetics plots for the run.

Usage
-----
    rikin_pipeline.py [-h] [-o OUTDIR] [-c CONFIG] [--dryrun] [--reuse] SEQA SEQB

Example:
    rikin_pipeline.py -o example AAAGGGGGGAAAAAAAGGGUGGGAAAAAAAGGGCGGGAAA CCCGCCC

Config file
-----------
On startup the script looks for `rikin_pipeline.cfg` next to itself
(override with --global-config) and, if present, loads it; `--config`
optionally supplies a second JSON file that is deep-merged on top
(dicts are merged key-by-key, everything else -- including lists -- is
fully replaced, matching how a second `source`d bash file would
overwrite a variable). See the accompanying rikin_pipeline.cfg for the
schema.
"""

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

VERSION = "0.9.5"

DEFAULT_SEQA_NAME = "seqA"
DEFAULT_SEQB_NAME = "seqB"
DEFAULT_OUTDIR = "Rikin-results"
DEFAULT_GLOBAL_CONFIG_NAME = "rikin_pipeline.cfg"


# --------------------------------------------------------------------------
# Config handling
# --------------------------------------------------------------------------

def deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge `override` into `base`. Dicts recurse; everything
    else (including lists) is fully replaced -- mirroring what a second
    `source`d bash config file would do to a plain variable."""
    result = dict(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def load_json_config(path: Path) -> dict:
    with open(path) as fh:
        return json.load(fh)


# --------------------------------------------------------------------------
# Helpers mirroring the bash script's call / tcall / call_redirect, plus
# dryrun-aware equivalents of the gzip/gunzip/sort housekeeping.
# --------------------------------------------------------------------------

def _fmt_cmd(cmd):
    return " ".join(str(c) for c in cmd)


def call(cmd, dryrun: bool) -> int:
    print()
    print(_fmt_cmd(cmd))
    if dryrun:
        return 0
    return subprocess.run(cmd).returncode


def tcall(cmd, dryrun: bool) -> int:
    print()
    print(_fmt_cmd(cmd))
    if dryrun:
        return 0
    start = time.time()
    rc = subprocess.run(cmd).returncode
    print(f"\nreal\t{time.time() - start:.3f}s")
    return rc


def call_redirect(target: Path, cmd, dryrun: bool) -> int:
    print()
    print(f"{_fmt_cmd(cmd)} > {target}")
    if dryrun:
        return 0
    with open(target, "w") as out:
        return subprocess.run(cmd, stdout=out).returncode


def gzip_file(path: Path, dryrun: bool):
    """Equivalent of `gzip -f path`."""
    if not path.exists():
        return
    print(f"gzip -f {path}")
    if dryrun:
        return
    gz_path = Path(str(path) + ".gz")
    with open(path, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    path.unlink()


def gunzip_file(gz_path: Path, dryrun: bool):
    """Equivalent of `gunzip -f path.gz`."""
    if not gz_path.exists():
        return
    print(f"gunzip -f {gz_path}")
    if dryrun:
        return
    assert str(gz_path).endswith(".gz")
    target = Path(str(gz_path)[:-3])
    with gzip.open(gz_path, "rb") as f_in, open(target, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    gz_path.unlink()


def sort_states(states_path: Path, sorted_path: Path, dryrun: bool):
    """Equivalent of `LC_ALL=C sort -k1,1n states > sorted_states; rm states`."""
    cmd = ["sort", "-k1,1n", str(states_path)]
    print()
    print(f"LC_ALL=C {_fmt_cmd(cmd)} > {sorted_path}")
    if dryrun:
        return
    env = dict(os.environ)
    env["LC_ALL"] = "C"
    with open(sorted_path, "w") as out:
        subprocess.run(cmd, stdout=out, env=env, check=True)
    states_path.unlink()


def remove_if_exists(path: Path):
    if path.exists():
        path.unlink()


def get_tool_version(bindir: Path) -> str:
    try:
        result = subprocess.run(
            [str(bindir / "rikin_enum"), "--version"],
            capture_output=True, text=True, timeout=10,
        )
        return result.stdout.strip() or result.stderr.strip() or "unknown"
    except Exception:
        return "unknown"


# --------------------------------------------------------------------------
# Argument parsing
# --------------------------------------------------------------------------

def build_arg_parser():
    parser = argparse.ArgumentParser(
        prog="rikin_pipeline.py",
        description="Run the RNAInterKin pipeline.",
        epilog=(
            "EXAMPLE CALL:\n"
            "  rikin_pipeline.py -o example AAAGGGGGGAAAAAAAGGGUGGGAAAAAAAGGGCGGGAAA CCCGCCC\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-o", "--outdir", default=DEFAULT_OUTDIR, help="output directory")
    parser.add_argument(
        "-c", "--config", type=Path,
        help="JSON file with pipeline configuration, deep-merged on top of the global config "
             "(rikin_pipeline.cfg next to this script). Comprises association_prefactor, "
             "common_opts, enum_opts, barriers_opts, prune_opts, xrates_opts, plot_opts.",
    )
    parser.add_argument(
        "--global-config", type=Path, default=None,
        help=f"Override path to the global config (default: {DEFAULT_GLOBAL_CONFIG_NAME} next to this script)",
    )
    parser.add_argument("--bindir", type=Path, default=None,
                         help="Directory containing the rikin_* tools (default: this script's directory)")
    parser.add_argument("--dryrun", action="store_true", help="don't run commands and/or write files")
    parser.add_argument("--reuse", action="store_true", help="reuse existing partial results in the output directory")
    parser.add_argument("--seqA-name", dest="seqA_name", default=DEFAULT_SEQA_NAME)
    parser.add_argument("--seqB-name", dest="seqB_name", default=DEFAULT_SEQB_NAME)
    parser.add_argument("seqA", help="RNA sequence A, as a word over A,C,G,U")
    parser.add_argument("seqB", help="RNA sequence B, as a word over A,C,G,U")
    return parser


# --------------------------------------------------------------------------
# Pipeline stages
# --------------------------------------------------------------------------

def stage_enumerate_and_sort(outdir: Path, bindir: Path, seqA, seqB, cfg, reuse, dryrun) -> bool:
    sorted_states = outdir / "sorted_states"
    sorted_states_gz = outdir / "sorted_states.gz"

    gzip_file(sorted_states, dryrun)

    if not reuse or not sorted_states_gz.exists():
        states = outdir / "states"
        cmd = [str(bindir / "rikin_enum"), seqA, seqB, *cfg.get("common_opts", []), *cfg.get("enum_opts", [])]
        rc = call_redirect(states, cmd, dryrun)
        if rc != 0:
            print("ERROR: rikin_enum died", file=sys.stderr)
            remove_if_exists(sorted_states_gz)
            return False

        sort_states(states, sorted_states, dryrun)

        return True  # caller sets reuse=False

    return None  # stage skipped, reuse unaffected


def stage_barriers(outdir: Path, bindir: Path, seqA, seqB, cfg, reuse, dryrun) -> bool:
    bg = outdir / "bg"
    barriers_track = outdir / "barriers_track.gz"
    track_ipps_barriers = outdir / "track-ipps-barriers.gz"
    sorted_states = outdir / "sorted_states"
    sorted_states_gz = outdir / "sorted_states.gz"

    ran = None
    if not reuse or not bg.exists() or not barriers_track.exists() or not track_ipps_barriers.exists():
        gunzip_file(sorted_states_gz, dryrun)

        cmd = [
            str(bindir / "rikin_barriers"), seqA, seqB,
            "-i", str(sorted_states),
            "-o", str(bg),
            "--compress-track", "--track", str(barriers_track),
            f"--track-ipps={track_ipps_barriers}",
            *cfg.get("common_opts", []), *cfg.get("barriers_opts", []),
        ]
        rc = call(cmd, dryrun)
        if rc != 0:
            print("ERROR: rikin_barriers died", file=sys.stderr)
            remove_if_exists(bg)
            return False
        ran = True

    # Unconditional cleanup, exactly as in the original script: whether or
    # not this stage ran, re-compress sorted_states if it was left uncompressed.
    gzip_file(sorted_states, dryrun)

    return ran


def stage_prune(outdir: Path, bindir: Path, cfg, reuse, dryrun) -> bool:
    out = outdir / "out"
    pf = outdir / "pf"
    bar = outdir / "bar"
    rates = outdir / "rates.gz"
    track_ipps_prune = outdir / "track-ipps-prune.gz"
    bg = outdir / "bg"
    prune_track = outdir / "prune_track.gz"
    track_ipps_barriers = outdir / "track-ipps-barriers.gz"

    if not reuse or not (out.exists() and pf.exists() and bar.exists()
                          and rates.exists() and track_ipps_prune.exists()):
        prefactor = cfg.get("association_prefactor", 1.0)
        prune_opts = ["--preexpf-first", str(prefactor), *cfg.get("prune_opts", [])]

        cmd = [
            str(bindir / "rikin_prune"), str(bg),
            *prune_opts,
            "--pffile", str(pf), "--verbose", "--barfile", str(bar),
            "--ratesfile", str(rates),
            "--compress-track", "--track", str(prune_track),
            f"--track-pps-out={track_ipps_prune}",
            f"--track-pps-in={track_ipps_barriers}",
        ]
        rc = call_redirect(out, cmd, dryrun)
        if rc != 0:
            print("ERROR: rikin_prune died", file=sys.stderr)
            remove_if_exists(rates)
            return False
        return True

    return None


def stage_solve_master_equation(outdir: Path, bindir: Path, cfg, reuse, dryrun) -> bool:
    kin = outdir / "kin"
    pf = outdir / "pf"
    error_file = outdir / "error"
    xrates_done = outdir / "xrates.done"

    if not reuse or not kin.exists():
        cmd = [str(bindir / "rikin_xrates.m"), str(pf), "--out", str(kin), *cfg.get("xrates_opts", [])]
        rc = tcall(cmd, dryrun)
        if rc != 0:
            if not dryrun:
                error_file.write_text(_fmt_cmd(cmd))
            print("ERROR: xrates died", file=sys.stderr)
            remove_if_exists(kin)
            return False
        else:
            if not dryrun:
                xrates_done.touch()
            remove_if_exists(error_file)
        return True

    return None


def stage_generate_all_plots(outdir: Path, bindir: Path, seqA, seqB, cfg, reuse, dryrun) -> bool:
    """New stage: call rikin_plot.py to produce the full set of kinetics plots."""
    plot_opts = cfg.get("plot_opts", {})
    plots_done = outdir / "plots.done"

    if reuse and plots_done.exists():
        print("\nSkipping plot generation (--reuse and plots.done already present).")
        return None

    plots_outdir = Path(plot_opts["output_dir"]) if plot_opts.get("output_dir") else outdir

    cmd = [
        str(bindir / "rikin_plot.py"),
        "--input-dir", str(outdir),
        "--name", outdir.name,
        "--seqA", seqA,
        "--seqB", seqB,
        "--output-dir", str(plots_outdir),
        "--shown-state-threshold", str(plot_opts.get("shown_state_threshold", 0.02)),
        "--formats", *plot_opts.get("formats", ["svg", "pdf"]),
        "--plots", *plot_opts.get("plots", ["all"]),
    ]
    if plot_opts.get("config"):
        cmd += ["--config", str(plot_opts["config"])]

    rc = call(cmd, dryrun)
    if rc != 0:
        print("ERROR: rikin_plot.py died", file=sys.stderr)
        remove_if_exists(plots_done)
        return False

    if not dryrun:
        plots_done.touch()
    return True


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    bindir = (args.bindir or Path(__file__).resolve().parent).resolve()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Copy the (optional) run-specific config into the output dir, as the
    # original script did, before it's parsed.
    if args.config is not None:
        if not args.config.exists():
            print(f"ERROR: configuration file {args.config} does not exist.", file=sys.stderr)
            sys.exit(1)
        shutil.copy(args.config, outdir)

    # a) global config
    global_config_path = args.global_config or (bindir / DEFAULT_GLOBAL_CONFIG_NAME)
    if global_config_path.exists():
        cfg = load_json_config(global_config_path)
        shutil.copy(global_config_path, outdir)
    else:
        print(f'WARNING: Cannot find global config at "{global_config_path}".')
        cfg = {}

    # b) optional run-specific config, deep-merged on top
    if args.config is not None:
        cfg = deep_merge(cfg, load_json_config(args.config))

    seqA, seqB = args.seqA, args.seqB

    version = get_tool_version(bindir)
    print("=" * 60)
    print(f"RNAInterKin Pipeline ver {VERSION}")
    print()
    print("Input:")
    print(f"  SeqA:   {seqA}")
    print(f"  SeqB:   {seqB}")
    print(f"  Outdir: {outdir}")
    print()
    print("Environment:")
    print(f"  Working dir: {Path.cwd()}")
    print(f"  BINDIR={bindir}")
    print(f"  RNAInterKin version {version}")
    print()
    print("Options:")
    print(f'  common_opts = {cfg.get("common_opts", [])}  #(controls rikin_enum and rikin_barriers)')
    print(f'  enum_opts = {cfg.get("enum_opts", [])}  #(controls rikin_enum)')
    print(f'  barriers_opts = {cfg.get("barriers_opts", [])}  #(controls rikin_barriers)')
    print(f'  prune_opts = {cfg.get("prune_opts", [])}  #(controls rikin_prune; '
          f'association_prefactor={cfg.get("association_prefactor", 1.0)})')
    print(f'  xrates_opts = {cfg.get("xrates_opts", [])}  #(controls solving of master equation by rikin_xrates.m)')
    print(f'  plot_opts = {cfg.get("plot_opts", {})}  #(controls plotting by rikin_plot.py)')
    print()
    print(f"Start date:  {datetime.now()}")
    print("=" * 60)

    reuse = args.reuse
    dryrun = args.dryrun

    ran = stage_enumerate_and_sort(outdir, bindir, seqA, seqB, cfg, reuse, dryrun)
    if ran is False:
        sys.exit(1)
    if ran:
        reuse = False

    ran = stage_barriers(outdir, bindir, seqA, seqB, cfg, reuse, dryrun)
    if ran is False:
        sys.exit(1)
    if ran:
        reuse = False

    ran = stage_prune(outdir, bindir, cfg, reuse, dryrun)
    if ran is False:
        sys.exit(1)
    if ran:
        reuse = False

    ran = stage_solve_master_equation(outdir, bindir, cfg, reuse, dryrun)
    if ran is False:
        sys.exit(1)
    if ran:
        reuse = False

    ran = stage_generate_all_plots(outdir, bindir, seqA, seqB, cfg, reuse, dryrun)
    if ran is False:
        sys.exit(1)

    print()
    print()
    print("=" * 60)
    print(f"RNAInterKin pipeline finished at {datetime.now()}.")
    print()
    print(f"Output files written to directory {outdir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
