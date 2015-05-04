# $Id: makefile,v 1.1.6.1 2002/04/24 03:25:33 erik Exp $ 
ALL: lib 

CFLAGS    =
FFLAGS    =
SOURCEC	  = ESMC_Error.c ESMC_ErrorF.c
SOURCEF	  = ESMF_ErrorMod.F
SOURCEH	  = ESMC_Error.h
OBJSC     = $(addsuffix .o, $(basename $(SOURCEC)))
OBJSF     = $(addsuffix .o, $(basename $(SOURCEF)))
TEXFILES  = $(addsuffix .tex, $(basename $(SOURCEF)))
MANSEC	  = Error
LIBBASE	  = libesmf
DIRS	  = 
LOCDIR	  = src/Infrastructure/Error/
CLEANFILES = $(LIBNAME) $(addprefix doc/, $(addsuffix .tex, $(basename $(SOURCEF))))

include ${ESMF_DIR}/build/${ESMF_ARCH}/base






