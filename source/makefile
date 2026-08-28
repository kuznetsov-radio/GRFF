# -D LINUX in sources means "not Windows" (OpenMP vs MSVC PPL). It is set on
# Darwin and Linux only; Windows builds use xThermal.vcxproj, not this makefile.
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
BINDIR = ../binaries

CXX ?= g++
COMMON_COMPILE = -c -O3 -fPIC

ifeq ($(UNAME_S),Darwin)
LIBOMP := $(shell brew --prefix libomp 2>/dev/null)
ifeq ($(wildcard $(LIBOMP)/include/omp.h),)
LIBOMP := /opt/homebrew/opt/libomp
endif
ifeq ($(wildcard $(LIBOMP)/include/omp.h),)
LIBOMP := /usr/local/opt/libomp
endif
ifeq ($(wildcard $(LIBOMP)/include/omp.h),)
$(error Homebrew libomp not found (need omp.h). Install with: brew install libomp)
endif
ifeq ($(UNAME_M),arm64)
ONAME = $(BINDIR)/GRFF_DEM_Transfer_arm64.so
else
ONAME = $(BINDIR)/GRFF_DEM_Transfer_x86_64.so
endif
# Apple clang rejects -fopenmp; Homebrew libomp needs the preprocessor flag.
CXXFLAGS = $(COMMON_COMPILE) -D LINUX -Xpreprocessor -fopenmp -I$(LIBOMP)/include
LDFLAGS = -shared -L$(LIBOMP)/lib -lomp -Wl,-rpath,$(LIBOMP)/lib
else ifeq ($(UNAME_S),Linux)
ONAME = $(BINDIR)/GRFF_DEM_Transfer.so
CXXFLAGS = $(COMMON_COMPILE) -D LINUX -fopenmp
LDFLAGS = -shared -fopenmp
else ifeq ($(OS),Windows_NT)
ONAME = GRFF_DEM_Transfer_64.dll
CXXFLAGS = $(COMMON_COMPILE) -fopenmp
LDFLAGS = -shared -fopenmp
else
ONAME = GRFF_DEM_Transfer.so
CXXFLAGS = $(COMMON_COMPILE) -D LINUX -fopenmp
LDFLAGS = -shared -fopenmp
endif

OBJS = Coulomb.o DEM.o ExtMath.o FF.o getparms.o GR.o ExtInterface.o MWtransfer.o \
       Neutrals.o Plasma.o Zeta.o GXcode.o RFactorKappa.o RFactorN.o

GRFF_DEM_Transfer: $(OBJS)
	$(CXX) $(LDFLAGS) -o $(ONAME) $(OBJS)

Coulomb.o: Coulomb.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CXX) $(CXXFLAGS) Coulomb.cpp

DEM.o: DEM.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CXX) $(CXXFLAGS) DEM.cpp

ExtMath.o: ExtMath.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CXX) $(CXXFLAGS) ExtMath.cpp

FF.o: FF.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CXX) $(CXXFLAGS) FF.cpp

getparms.o: getparms.cpp
	$(CXX) $(CXXFLAGS) getparms.cpp

GR.o: GR.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CXX) $(CXXFLAGS) GR.cpp

ExtInterface.o: ExtInterface.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CXX) $(CXXFLAGS) ExtInterface.cpp

MWtransfer.o: MWtransfer.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CXX) $(CXXFLAGS) MWtransfer.cpp

Neutrals.o: Neutrals.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CXX) $(CXXFLAGS) Neutrals.cpp

Plasma.o: Plasma.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CXX) $(CXXFLAGS) Plasma.cpp

Zeta.o: Zeta.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CXX) $(CXXFLAGS) Zeta.cpp

GXcode.o: GXcode.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CXX) $(CXXFLAGS) GXcode.cpp

RFactorKappa.o: RFactorKappa.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CXX) $(CXXFLAGS) RFactorKappa.cpp

RFactorN.o: RFactorN.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CXX) $(CXXFLAGS) RFactorN.cpp

.PHONY: install clean

install: GRFF_DEM_Transfer

clean:
	rm -f *.o
