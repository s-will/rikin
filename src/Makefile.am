bin_PROGRAMS=rrikin_enum rrikin_barriers hpf

hpf_SOURCES=hpf.cc hybrid_pf.cc hybrid_pf.hh unpaired_pf.cc	\
unpaired_pf.hh LocARNA/matrices.hh LPfold.h LPfold.c		\
hybrid_ensemble_model.hh hybrid_ensemble_model.cc moves.cc
# hybrid_ensemble_state.hh hybrid_ensemble_state.cc

BUILT_SOURCES = rrikin_enum_cmdline.c rrikin_enum_cmdline.h	\
rrikin_barriers_cmdline.c rrikin_barriers_cmdline.h config.h

CLEANFILES = $(BUILT_SOURCES)

rrikin_enum_SOURCES = rrikin_enum.cc hybrid_pf.cc hybrid_pf.hh		\
unpaired_pf.cc unpaired_pf.hh LocARNA/matrices.hh LPfold.h LPfold.c	\
hybrid_ensemble_model.hh hybrid_ensemble_model.cc moves.cc		\
rrikin_enum_cmdline.c rrikin_enum_cmdline.h

rrikin_enum_cmdline.c: rrikin_enum.ggo
	gengetopt <rrikin_enum.ggo

rrikin_barriers_SOURCES = rrikin_barriers.cc hybrid_pf.cc hybrid_pf.hh		\
unpaired_pf.cc unpaired_pf.hh LocARNA/matrices.hh LPfold.h LPfold.c	\
hybrid_ensemble_model.hh hybrid_ensemble_model.cc moves.cc		\
rrikin_barriers_cmdline.c rrikin_barriers_cmdline.h


rrikin_barriers_cmdline.c: rrikin_barriers.ggo
	gengetopt <rrikin_barriers.ggo

config.h: ../config.h
	if [ -e "$@" ] ; then $(RM) $@ ; fi
	ln -s $< $@
