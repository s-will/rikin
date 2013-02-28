bin_PROGRAMS=rrikin_enum rrikin_barriers test_moves

test_moves_SOURCES=test_moves.cc hybrid_pf.cc hybrid_pf.hh		\
unpaired_pf.cc unpaired_pf.hh LocARNA/matrices.hh LPfold.h LPfold.c	\
hybrid_ensemble_model.hh hybrid_ensemble_model.cc moves.hh moves.cc

rrikin_enum_SOURCES = rrikin_enum.cc hybrid_pf.cc hybrid_pf.hh		\
unpaired_pf.cc unpaired_pf.hh LocARNA/matrices.hh LPfold.h LPfold.c	\
hybrid_ensemble_model.hh hybrid_ensemble_model.cc moves.hh moves.cc	\
rrikin_enum_cmdline.c rrikin_enum_cmdline.h

rrikin_barriers_SOURCES = rrikin_barriers.cc barrier_graph.hh		\
barrier_graph.cc basin_transition.hh hybrid_pf.cc hybrid_pf.hh		\
unpaired_pf.cc unpaired_pf.hh LocARNA/matrices.hh LPfold.h LPfold.c	\
hybrid_ensemble_model.hh hybrid_ensemble_model.cc moves.cc		\
rrikin_barriers_cmdline.c rrikin_barriers_cmdline.h


############################################################
## gengetopt command line parsers
##

## register parsers as built sources
BUILT_SOURCES = rrikin_enum_cmdline.c rrikin_enum_cmdline.h	\
rrikin_barriers_cmdline.c rrikin_barriers_cmdline.h config.h

## register built sources as cleanfiles (otherwise they would not be
## removed by 'make clean')
CLEANFILES = $(BUILT_SOURCES)

rrikin_enum_cmdline.c: $(srcdir)/rrikin_enum.ggo
	gengetopt < $<

rrikin_barriers_cmdline.c: $(srcdir)/rrikin_barriers.ggo
	gengetopt < $<


## get a copy of config.h in directory src 
config.h: ../config.h
	if [ -e "$@" ] ; then $(RM) $@ ; fi
	ln -s $< $@
