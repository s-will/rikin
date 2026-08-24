# CLI reference

Full `--help` output for each RIKin command-line tool. See the [main README](https://github.com/s-will/rikin#readme) for a narrative overview and usage examples.

<a id="rikin_enum"></a>
## rikin_enum

Enumerate candidate interaction/structure states for the two input RNAs.

```
Usage: rikin_enum [options] seq1 seq2
Calculate kinetics of RNA-RNA-interaction

Enumerates states of RNA-RNA-interaction in the RRI hybrid ensemble model.
States are modelled as equilibrated partitions of the whole interaction
ensemble. Supports states with single and double hybridization sites.

  -h, --help                    Print help and exit
  -V, --version                 Print version and exit
      --max-hyb-energy=DOUBLE   Threshold for single hybridization site energy
                                  (default=`1e12')
      --max-total-energy=DOUBLE Threshold for single hybridization site energy
                                  (default=`1e12')
      --max-hyb-length=INT      Maximum length of each subsequence in
                                  hybridization site  (default=`-1')
      --max-hyb-length-diff=INT Maximum difference of lengths in hybridization
                                  sites  (default=`-1')
      --double-sites            Include double hybridization site
                                  configurations
      --binary                  Binary encoding of states in output
      --homodimer               Enumerate conformations of homodimer (only one
                                  sequence accepted)
      --antisense               Enumerate conformations of sense/antisense pair
                                  (only one sequence accepted)
      --regionA=STRING          Restrict region for interaction <n> bases at 5'
                                  head (>0) or 3' tail (<0) of sequence A; use
                                  i-j for region (0=no restriction).
                                  (default=`0')
      --regionB=STRING          Restrict region for interaction <n> bases at 5'
                                  head (>0) or 3' tail (<0) of sequence B; use
                                  i-j for region (0=no restriction).
                                  (default=`0')
      --span=INT                Base pair span for local folding. Set this to
                                  turn on local folding.  (default=`-1')
      --temperature=DOUBLE      Folding temperature  (default=`37.0')
      --verbose                 Turn on verbose output

Input sequences:
  seq1 first sequence in 5'->3' orientation
  seq2 second sequence in 5'->3' orientation
```

<a id="rikin_barriers"></a>
## rikin_barriers

Construct the discretely coarse-grained basin/state system.

```
Usage: rikin_barriers [options] seq1 seq2
Construct basin state system for kinetics of RNA-RNA interaction

Reads energy-sorted list of states of RNA-RNA-interaction in the RRI hybrid
ensemble model and combines them into basin macrostates. Furthermore, computes
rates between macro-states. States are modelled as equilibrated partitions of
the whole interaction ensemble.

  -h, --help                    Print help and exit
  -V, --version                 Print version and exit
      --max-hyb-length=INT      Maximum length hybridization site length to be
                                  precomputed  (default=`-1')
      --max-hyb-length-diff=INT Maximum difference of lengths in hybridization
                                  sites  (default=`-1')
      --max-recover-energy=DOUBLE
                                Recover states with smaller or equal energy;
                                  def=effectively turn off recovery
                                  (default=`-1e6')
      --min-rate=DOUBLE         Minimum transition rate; otherwise ignore
                                  (sparsifies edges)  (default=`1e-12')
      --double-sites            Include double hybridization site
                                  configurations
      --no-special-open-state   Treat open state like other states
      --binary                  Assume binary encoding of states in input
      --no-gradient             Don't combine into basins by gradient walk
      --homodimer               Compute rates for homodimer (only one sequence
                                  accepted)
      --antisense               Computes rates for sense/antisense pair (only
                                  one sequence accepted)
      --compress-track          Gzip-compress all written track or ipp track
                                  information
      --track=STRING            Turn on basin tracking, providing output
                                  filename
      --track-ipps=STRING       Track interaction pair probabilities, providing
                                  output filename
      --ipp-min-prob=DOUBLE     Minimum interaction pair probabilities for
                                  output  (default=`1e-2')
      --regionA=STRING          Restrict region for interaction <n> bases at 5'
                                  head (>0) or 3' tail (<0) of sequence A; use
                                  i-j for region (0=no restriction).
                                  (default=`0')
      --regionB=STRING          Restrict region for interaction <n> bases at 5'
                                  head (>0) or 3' tail (<0) of sequence B; use
                                  i-j for region (0=no restriction).
                                  (default=`0')
      --span=INT                Base pair span for local folding. Set this to
                                  turn on local folding.  (default=`-1')
      --temperature=DOUBLE      Folding temperature  (default=`37.0')
  -o, --output=STRING           Output filename
  -i, --input=STRING            Input filename
      --verbose                 Turn on verbose output
      --debug                   Turn on debugging output

Input sequences:
  seq1 first sequence in 5'->3' orientation
  seq2 second sequence in 5'->3' orientation
```

<a id="rikin_prune"></a>
## rikin_prune

Prune/continuously coarse-grain the state system and compute rates.

```
Usage: rikin_prune [options] bgfile
Prune state system for kinetics of RNA-RNA interaction

Prune barrier graph by various criteria (outflow, rate, and equilibrium
probability).

  -h, --help                    Print help and exit
  -V, --version                 Print version and exit
      --max-outflow=DOUBLE      Maximal outflow of a state; otherwise dissolve
                                  (default=`1e6')
      --min-pequ=DOUBLE         Minimum equilibrium probability of a state;
                                  otherwise dissolve  (default=`0')
      --qu-outflow=DOUBLE       State quantile (in percent) to keep due to
                                  max-outflow criterion (approximate)
      --qu-pequ=DOUBLE          State quantile (in percent) to keep due to
                                  min-pequ criterion
      --num-outflow=DOUBLE      Number of states to keep due to max-outflow
                                  criterion (approximate)
      --num-pequ=DOUBLE         Number of states to keep due to min-pequ
                                  criterion
      --min-rate=DOUBLE         Minimum transition rate; otherwise ignore
                                  (sparsifies edges)  (default=`0')
      --dont-simplify-graph     Turn of all barrier graph simplifications
      --no-special-first-state  Treat first state like any other state;
                                  def=never dissolve first state
      --preexpf-first=DOUBLE    Pre-exponential factor for rates from and to
                                  first (usually "open") state (default/other
                                  rates: factor 1; only with special first
                                  state)  (default=`1.0')
      --barfile=STRING          Write treekin compatible pseudo-barrier file
      --ratesfile=STRING        Write rates matrix to a treekin compatible file
      --pffile=STRING           Write partition functions of states and
                                  transition states
      --rxns=STRING             Write reactions to file in 'rxns' format
                                  (Stochastirator-compatible)
      --spcs=STRING             Write species to file in 'spcs' format
                                  (Stochastirator-compatible)
      --nameA=STRING            Specify name of molecule A in rxns/spcs (def=A)
      --nameB=STRING            Specify name of molecule B in rxns/spcs (def=B,
                                  in case of homodimer, defaults to nameA)
      --nameAB=STRING           Specify name of molecule AB in rxns/spcs
                                  (def=nameA+nameB)
      --to-keep=INT             Keep only the specified states (separate
                                  indices by comma)
      --only-first-component    If barrier graph is not connected, keep only
                                  first component [default=heuristically
                                  connect the components]
      --verbose                 Turn on verbose output
      --debug                   Turn on debugging output
      --compress-track          Gzip-compress all written track or pp track
                                  information
      --track=STRING            Turn on tracking of pruning, providing output
                                  filename
      --track-pps-out=STRING    Track interaction pair probabilities (specifies
                                  output file)
      --track-pps-in=STRING     Track interaction pair probabilities (specifies
                                  input file)
      --pp-min-prob=DOUBLE      Minimum interaction pair probabilities for
                                  output  (default=`1e-2')

Input file: bgfile = file as generated by rikin_barriers -o <filename>
If more than one specification of outflow or pequ is provided, the strictest
filter is applied.
```

<a id="rikin_pipeline_py"></a>
## rikin_pipeline.py

Run the complete pipeline (enumeration through plotting) for a pair of sequences.

```
usage: rikin_pipeline.py [-h] [-o OUTDIR] [-c CONFIG]
                         [--global-config GLOBAL_CONFIG] [--bindir BINDIR]
                         [--dryrun] [--reuse] [--fastaA FASTAA]
                         [--fastaB FASTAB] [--seqA SEQA] [--seqB SEQB]

Run the RNAInterKin pipeline.

options:
  -h, --help            show this help message and exit
  -o, --outdir OUTDIR   output directory
  -c, --config CONFIG   JSON file with pipeline configuration, deep-merged on
                        top of the global config (rikin_pipeline.cfg next to
                        this script). Comprises association_prefactor,
                        common_opts, enum_opts, barriers_opts, prune_opts,
                        xrates_opts, plot_opts.
  --global-config GLOBAL_CONFIG
                        Override path to the global config (default:
                        rikin_pipeline.cfg next to this script)
  --bindir BINDIR       Directory containing the rikin_* tools (default: this
                        script's directory)
  --dryrun              don't run commands and/or write files
  --reuse               reuse existing partial results in the output directory
  --fastaA FASTAA       Fasta containing first sequence
  --fastaB FASTAB       Fasta containing second sequence
  --seqA SEQA           First sequence
  --seqB SEQB           Second sequence

EXAMPLE CALLS:
   rikin_pipeline.py -o example --seqA AAAGGGGGGAAAAAAAGGGUGGGAAAAAAAGGGCGGGAAA --seqB CCCGCCC
   rikin_pipeline.py -o example --fastaA exampleA.fasta --fastaB exampleB.fasta -c example_config.cfg -o example
```

<a id="rikin_plot_py"></a>
## rikin_plot.py

Render kinetics plots from a completed run's result files.

```
usage: rikin_plot.py [-h] [--input-dir INPUT_DIR] [--name NAME]
                     [--fastaA FASTAA] [--fastaB FASTAB] [--seqA SEQA]
                     [--seqB SEQB]
                     [--shown-state-threshold SHOWN_STATE_THRESHOLD]
                     [--runs-config RUNS_CONFIG] [--output-dir OUTPUT_DIR]
                     [--formats {svg,pdf,png} [{svg,pdf,png} ...]]
                     [--plots {kinetics,states,dotplot,paired_kinetics,all} [{kinetics,states,dotplot,paired_kinetics,all} ...]]
                     [--config CONFIG]

Generate RNA-RNA interaction kinetics plots from RikinRun result directories.

options:
  -h, --help            show this help message and exit
  --fastaA FASTAA       Fasta containing first sequence
  --fastaB FASTAB       Fasta containing second sequence
  --seqA SEQA           First sequence
  --seqB SEQB           Second sequence
  --runs-config RUNS_CONFIG
                        JSON file describing multiple runs (batch mode). See
                        module docstring for the schema.
  --output-dir OUTPUT_DIR
                        Directory to write plots to (default: .)
  --formats {svg,pdf,png} [{svg,pdf,png} ...]
                        Output file formats (default: svg pdf, matching
                        rikinplotlib's own default)
  --plots {kinetics,states,dotplot,paired_kinetics,all} [{kinetics,states,dotplot,paired_kinetics,all} ...]
                        Which plots to generate (default: all)
  --config CONFIG       Optional JSON file with per-plot kwargs, keyed by plot
                        name, forwarded straight to the matching RikinRun
                        method.

single-run options (ignored if --runs-config is given):
  --input-dir INPUT_DIR
                        Directory with RikinRun result files
  --name NAME           Prefix used for output plot filenames (RikinRun's
                        autosave name)
  --shown-state-threshold SHOWN_STATE_THRESHOLD
```

<a id="rikin_xrates_m"></a>
## rikin_xrates.m

Solve the Master Equation via Padé matrix exponentiation. Has its own shebang; run directly, not via `octave rikin_xrates.m`.

```
XRates 0.2 -- Calculate kinetics from weights/pfs of states and transitions by solving
  the master equation for Arrhenius rates.
  --help,-h                     Print this help.
  --verbose,-v                  Be verbose.
  --binary,-b                   Binary input
  --nonbinary                   Non-binary input
  --out <string>                Name of output file. The distributions at
                                each time point until convergence are written
                                as a table to this file. (mandatory)
  --t0 <double>                 Start time
  --t8 <double>                 End time
  --tinc <double>               Time increment
  --mode <string>               Mode (expm or diag)
  --p0 <int>                    State with initial probability 1
  --absorb <int>                Absorbing state
  --mfpts                       Compute mean first passage times
  in <string>                   Input file. Matrix of transition
                                weights/partition functions; the diagonal
                                contains state weights. (mandatory)
```

