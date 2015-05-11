# $Id: makefile,v 1.1.6.1 2002/04/24 03:25:55 erik Exp $ 
ALL: lib 

CFLAGS    =
FFLAGS    =
SOURCEC	  = ESMC_Timer.c ESMC_TimerF.c
SOURCEF	  = ESMF_TimerMod.F
SOURCEH	  = ESMC_Timer.h
OBJSC     = $(addsuffix .o, $(basename $(SOURCEC)))
OBJSF     = $(addsuffix .o, $(basename $(SOURCEF)))
TEXFILES  = $(addsuffix .tex, $(basename $(SOURCEF)))
MANSEC	  = Timeobj
LIBBASE	  = libesmf
DIRS	  = examples doc
LOCDIR	  = src/Infrastructure/Timer/
CLEANFILES = $(LIBNAME) $(addprefix doc/, $(addsuffix .tex, $(basename $(SOURCEF))))

include ${ESMF_DIR}/build/${ESMF_ARCH}/base






