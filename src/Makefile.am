SUBDIRS=Tests

bin_PROGRAMS=rrikin_enum rrikin_barriers rrikin_prune test_moves

bin_SCRIPTS=xrates.m kinetics.R kinetics_lib.R

# base sources, common to all programs
BASE_SRCS=hybrid_pf.cc hybrid_pf.hh unpaired_pf.cc unpaired_pf.hh	\
LPfold.h LPfold.c hybrid_ensemble_model.hh hybrid_ensemble_model.cc	\
moves.hh moves.cc util.hh util.cc const_sparse_matrix.hh		\
ordered_sparse_matrix.hh rri_enumeration.hh rri_enumeration.cc

# extended base sources, common to programs operating on barrier graph
EXT_BASE_SRCS=$(BASE_SRCS) barrier_graph.hh barrier_graph.cc basin.hh	\
basin.cc basin_info.hh basin_info.cc basin_pruning_info.hh		\
basin_pruning_info.cc basin_transition.hh basin_transition.cc		\
pair_pfs.hh pair_pfs.cc


test_moves_SOURCES=test_moves.cc $(BASE_SRCS)

rrikin_enum_SOURCES = rrikin_enum.cc rrikin_enum_cmdline.c	\
rrikin_enum_cmdline.h $(BASE_SRCS)

rrikin_barriers_SOURCES = rrikin_barriers.cc			\
rrikin_barriers_cmdline.c rrikin_barriers_cmdline.h		\
rri_barrier_graph.hh rri_barrier_graph.cc $(EXT_BASE_SRCS)

rrikin_prune_SOURCES = rrikin_prune.cc rrikin_prune_cmdline.c	\
rrikin_prune_cmdline.h $(EXT_BASE_SRCS)


############################################################
## gengetopt command line parsers
##

## register parsers and config.h as built sources
BUILT_SOURCES = rrikin_enum_cmdline.c rrikin_enum_cmdline.h	\
rrikin_barriers_cmdline.c rrikin_barriers_cmdline.h		\
rrikin_prune_cmdline.c rrikin_prune_cmdline.h config.h

EXTRA_DIST=rrikin_barriers.ggo rrikin_enum.ggo rrikin_prune.ggo

## register built sources as cleanfiles (otherwise they would not be
## removed by 'make clean')
CLEANFILES = $(BUILT_SOURCES)

## gen get opt rules (rules are expanded to avoid GNU-only pattern rule)
rrikin_barriers_cmdline.c rrikin_barriers_cmdline.h: $(srcdir)/rrikin_barriers.ggo 
	gengetopt < $<

rrikin_enum_cmdline.c rrikin_enum_cmdline.h: $(srcdir)/rrikin_enum.ggo 
	gengetopt < $<

rrikin_prune_cmdline.c rrikin_prune_cmdline.h: $(srcdir)/rrikin_prune.ggo 
	gengetopt < $<

## get a copy of config.h in directory src 
config.h: ../config.h
	if [ -e "$@" ] ; then $(RM) $@ ; fi
	ln -s $< $@
