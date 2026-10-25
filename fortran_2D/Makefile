FC      = gfortran
FFLAGS  = -O3 -fdefault-real-8 -funroll-loops -ffast-math -march=native -flto
#FFLAGS = -O0 -g -fcheck=all -fbacktrace -Wall -Wextra
TARGET  = solver2d

# Object files in strict dependency order
OBJS = GaussQuadrature1D.o \
       LinAlg.o            \
       WENO.o              \
       PDE2D.o             \
       Prob2D.o            \
       Mesh2D.o            \
       Solver2D.o          \
       main_Solver2D.o

# Default target
all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

# Module compilation rules (order matters for .mod generation)
GaussQuadrature1D.o: GaussQuadrature1D.f90
	$(FC) $(FFLAGS) -c $<

LinAlg.o: LinAlg.f90
	$(FC) $(FFLAGS) -c $<

WENO.o: WENO.f90
	$(FC) $(FFLAGS) -c $<

PDE2D.o: PDE2D.f90
	$(FC) $(FFLAGS) -c $<

Prob2D.o: Prob2D.f90 PDE2D.o
	$(FC) $(FFLAGS) -c $<

Mesh2D.o: Mesh2D.f90 Prob2D.o GaussQuadrature1D.o LinAlg.o
	$(FC) $(FFLAGS) -c $<s

Solver2D.o: Solver2D.f90 PDE2D.o Prob2D.o Mesh2D.o WENO.o GaussQuadrature1D.o
	$(FC) $(FFLAGS) -c $<

main_Solver2D.o: main_Solver2D.f90 Prob2D.o Mesh2D.o Solver2D.o
	$(FC) $(FFLAGS) -c $<

# Utility targets
clean:
	del /f /q $(OBJS) *.mod $(TARGET).exe 2>nul || true

.PHONY: all clean
