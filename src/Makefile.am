SUBDIRS=Tests

bin_PROGRAMS=rikin_enum rikin_barriers rikin_prune test_moves

dist_bin_SCRIPTS=rikin_xrates.m rikin_plot.py rikinplotlib.py\
	rikin_pipeline.py rikin_pipeline.cfg

# base sources, common to all programs
BASE_SRCS=hybrid_pf.cc hybrid_pf.hh unpaired_pf.cc unpaired_pf.hh	\
hybrid_ensemble_model.hh hybrid_ensemble_model.cc moves.hh moves.cc	\
util.hh util.cc const_sparse_matrix.hh ordered_sparse_matrix.hh		\
rri_enumeration.hh rri_enumeration.cc

# extended base sources, common to programs operating on barrier graph
EXT_BASE_SRCS=$(BASE_SRCS) barrier_graph.hh barrier_graph.cc basin.hh	\
basin.cc basin_info.hh basin_info.cc basin_pruning_info.hh		\
basin_pruning_info.cc basin_transition.hh basin_transition.cc		\
pair_pfs.hh pair_pfs.cc


test_moves_SOURCES=test_moves.cc $(BASE_SRCS)

rikin_enum_SOURCES = rikin_enum.cc rikin_enum_cmdline.c	\
rikin_enum_cmdline.h $(BASE_SRCS)

rikin_barriers_SOURCES = rikin_barriers.cc			\
rikin_barriers_cmdline.c rikin_barriers_cmdline.h		\
rri_barrier_graph.hh rri_barrier_graph.cc $(EXT_BASE_SRCS)

rikin_prune_SOURCES = rikin_prune.cc rikin_prune_cmdline.c	\
rikin_prune_cmdline.h $(EXT_BASE_SRCS)


############################################################
## gengetopt command line parsers
##

## register parsers and config.h as built sources
BUILT_SOURCES = rikin_enum_cmdline.c rikin_enum_cmdline.h	\
rikin_barriers_cmdline.c rikin_barriers_cmdline.h		\
rikin_prune_cmdline.c rikin_prune_cmdline.h config.h

EXTRA_DIST=rikin_barriers.ggo rikin_enum.ggo rikin_prune.ggo

## register built sources as cleanfiles (otherwise they would not be
## removed by 'make clean')
CLEANFILES = $(BUILT_SOURCES)

## gen get opt rules (rules are expanded to avoid GNU-only pattern rule)
rikin_barriers_cmdline.c rikin_barriers_cmdline.h: $(srcdir)/rikin_barriers.ggo 
	gengetopt < $<

rikin_enum_cmdline.c rikin_enum_cmdline.h: $(srcdir)/rikin_enum.ggo 
	gengetopt < $<

rikin_prune_cmdline.c rikin_prune_cmdline.h: $(srcdir)/rikin_prune.ggo 
	gengetopt < $<

## get a copy of config.h in directory src 
config.h: ../config.h
	if [ -e "$@" ] ; then $(RM) $@ ; fi
	ln -s $< $@
