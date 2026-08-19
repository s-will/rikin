ChangeLog
=========

Post release
------------------
* Update/fix doxygen Documentation

0.9.5 (Aug, 2026)
-----------------
* Add rikin_pipeline.py that runs entire pipeline
  and generates plots
* Add Python code for new kinetic plots (matplotlib/seaborn)
* Add rikin_plot.py plotting CLI and rikinplotlib library
* Produce kinetics plots from Jupyter notebook
* Add/update Bioconda recipe (meta.yaml/build.sh) for the rikin package
* License the project under AGPL-3.0-or-later; add COPYING
* Rewrite README.md: updated installation instructions (Conda and
  from-source), usage/CLI documentation, configuration file schema, and
  example figures
* Add Examples directory with biologically relevant examples
  (kissing-hairpin system, MicA-MicA homodimer, MicA-OmpA heterodimer),
  an accompanying README, a run-all.sh script, and a plotting notebook
* Update convergence criterion in xrates

0.9.4 (Oct, 2023)
-----------------
* Update continuous coarse graining
* Refactorize code for basin transitions in pruning
* Add rikin_pipeline script
* Rename tool to RIkin
* Doxygen awesome style for generated documentation
* Extend documentation

0.9.3 (June, 2022)
------------------
* Update / fix reading of binary from octave

0.9.2 (Feb, 2017)
-------------------

* make rikin_prune argument only-first-component optional

0.9.1 (Feb, 2017)
-------------------

* code refactorization
* bug fixes

0.9 (Apr, 2016)
----------------

* connect disconnected components heuristically
* add pre-exponential factor for transitions from and to first state
* [=open state] in rrikin_prune
