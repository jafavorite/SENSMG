EXE             = ../bin/sensmg
FCOMPILATOR     = gfortran
FOPTIONS        = $(DEBUG) -Wunused -fcheck=bounds -fbacktrace -openmp -DGNU
FCOMPILE        = $(FCOMPILATOR) -c $(FOPTIONS)

FOBJS           = $(patsubst %.F, %.o, $(wildcard *.F))
F90OBJS         = $(patsubst %.F90, %.o, $(wildcard *.F90))

$(EXE): $(F90OBJS) $(FOBJS) $(OBJS)
	$(FCOMPILATOR) $(FOPTIONS) -o $(EXE) $(FOBJS) $(F90OBJS) $(OBJS) $(CPP_LIBS) -shared-libgcc
	@rm -f *__genmod.f90 *__genmod.mod

%.o: %.F90
	$(FCOMPILE) $<

%.o: %.F
	$(FCOMPILE) $<

$(FOBJS): coms.F90

clean:
	@echo "   cleaning obj files"
	@rm -f *.o
	@echo "   cleaning mod files"
	@rm -f *.mod
	@echo "   cleaning generated files"
	@rm -f *__genmod.f90 *__genmod.mod
	@echo "   cleaning trash files"
	@rm -f core *~

dist_clean:
	@echo "   cleaning exe files"
	@rm -f $(EXE)
	@echo "   cleaning obj files"
	@rm -f *.o
	@echo "   cleaning mod files"
	@rm -f *.mod
	@echo "   cleaning generated files"
	@rm -f *__genmod.f90 *__genmod.mod
	@echo "   cleaning trash files"
	@rm -f core *~	

.PHONY: clean dist_clean

