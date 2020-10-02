OBJS = mitoutil.o keys.o readkey.o genutil.o mt19937.o inputparams.o dynamics.o

FC = gfortran
FLAGS = -O3 -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -fcheck=all
#FLAGS = -O3 -fopenmp
DEBUGFLAGS =
LIBS =  -lblas -llapack
OSGFLAGS = -fopenmp -static

# executable name
EXE = ../ssphago.exe
OSGEXE=../ssphago_osg.exe
TESTEXE = ../testrun.exe

.SUFFIXES: .o .f90 

# general suffix rules
%.o: %.f90 
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

%.o: tools/%.f90
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

default: $(OBJS) main.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) main.o -o $(EXE) $(LIBS)

osg: $(OBJS) main.o
	$(FC) $(OSGFLAGS) $(DEBUGFLAGS) $(OBJS) main.o -o $(OSGEXE) $(LIBS)

test: $(OBJS) testmain.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) testmain.o -o $(TESTEXE) $(LIBS)

static: $(OBJS) main.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) -static $(OBJS) main.o -o $(EXE) $(LIBS) 



clean: 
	rm *.o *.mod

veryclean: clean
	rm *~ $(EXE)


# ------------
# individual file dependencies
# ------------
mitoutil.o: keys.o mt19937.o genutil.o
readkey.o: inputparams.o keys.o genutil.o mt19937.o
genutil.o: mt19937.o 
dynamics.o: mt19937.o mitoutil.o
