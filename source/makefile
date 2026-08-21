CC=g++
FLAGS=-fopenmp

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    ONAME = GRFF_DEM_Transfer.so
endif
ifeq ($(UNAME_S),Darwin)
    FLAGS=-arch arm64 -lomp 
    ONAME = GRFF_DEM_Transfer_arm64.so
    CC=gcc
    LDFLAGS=-L$(brew --prefix libomp)/lib -lomp
    CPPFLAGS=-I$(brew --prefix libomp)/include -Xpreprocessor -fopenmp
endif
ifeq ($(OS),Windows_NT)
    ONAME = GRFF_DEM_Transfer_64.dll
else
    ONAME = GRFF_DEM_Transfer.so
endif

GRFF_DEM_Transfer: Coulomb.o DEM.o ExtMath.o FF.o getparms.o GR.o ExtInterface.o MWtransfer.o Neutrals.o Plasma.o Zeta.o GXcode.o RFactorKappa.o RFactorN.o
	$(CC) $(LDFLAGS) -shared -o $(ONAME) Coulomb.o DEM.o ExtMath.o FF.o getparms.o GR.o ExtInterface.o MWtransfer.o Neutrals.o Plasma.o Zeta.o GXcode.o RFactorKappa.o RFactorN.o

Coulomb.o: Coulomb.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX Coulomb.cpp

DEM.o: DEM.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX DEM.cpp

ExtMath.o: ExtMath.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX ExtMath.cpp

FF.o: FF.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX FF.cpp

getparms.o: getparms.cpp
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX getparms.cpp

GR.o: GR.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX GR.cpp

ExtInterface.o: ExtInterface.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX ExtInterface.cpp

MWtransfer.o: MWtransfer.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX MWtransfer.cpp

Neutrals.o: Neutrals.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX Neutrals.cpp
 
Plasma.o: Plasma.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX Plasma.cpp

Zeta.o: Zeta.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX Zeta.cpp

GXcode.o: GXcode.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX GXcode.cpp

RFactorKappa.o: GXcode.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX RFactorKappa.cpp

RFactorN.o: GXcode.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h RFactorKappa.h RFactorN.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX RFactorN.cpp

install:
	cp $(ONAME) ../binaries

clean:
	rm -f *.o $(ONAME)