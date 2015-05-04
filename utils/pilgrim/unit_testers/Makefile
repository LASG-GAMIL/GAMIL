
include ../Makefile.conf

PILGRIMLIB      =  ../libpilgrim.a
TESTSRCS        =  decomptest.F90 parutilitiestest.F90 ghosttest.F90 \
                   redistributetest.F90 parpatterntest.F90 unstructured.F90 
TESTS           =  decomptest parutilitiestest ghosttest \
                   redistributetest parpatterntest
TESTOBJS        = ${TESTSRCS:.F90=.o}

#OPTIONS         = -DSTAND_ALONE -DTIMING -I. -I..
OPTIONS         = -I..

.SUFFIXES:	.F90 .c .o

# -------------------------------------------------------------------

All:	tests

all:	tests

tests:	clean $(TESTS)

metis:
	(cd metis-4.0; make; cp libmetis.a ..)

decomptest:     decomptest.o
	$(LD) $(LDFLAGS) -o decomptest decomptest.o -L.. -lpilgrim

parutilitiestest:     parutilitiestest.o
	$(LD) $(LDFLAGS) -o parutilitiestest parutilitiestest.o -L.. -lpilgrim $(_LMPI) $(_lMPI) 

ghosttest:     ghosttest.o
	$(LD) $(LDFLAGS) -o ghosttest ghosttest.o -L.. -lpilgrim $(_LMPI) $(_lMPI)

redistributetest:     redistributetest.o
	$(LD) $(LDFLAGS) -o redistributetest redistributetest.o -L.. -lpilgrim $(_LMPI) $(_lMPI)

parpatterntest:       unstructured.o parpatterntest.o
	$(LD) $(LDFLAGS) -o parpatterntest parpatterntest.o unstructured.o -L.. -lpilgrim $(_LMPI) $(_lMPI)

.F90.o:
	$(FC) $(FFLAGS) $(DFLAGS) $(OPTIONS) -c $<

clean:
	-@${RM} -f $(TESTOBJS) $(TESTS) *~ *.mod *.MOD

libclean:
	-@${RM} -f $(PILGRIMLIB)

