CC=gcc
FLAGS=-lomp 
LDFLAGS=-L/opt/homebrew/lib -L/opt/homebrew/opt/llvm/lib
CPPFLAGS=-I/opt/homebrew/include -I/opt/homebrew/opt/llvm/include

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    EXTENSION = .so
endif
ifeq ($(UNAME_S),Darwin)
	FLAGS=-arch arm64 -lomp 
    EXTENSION = .mac.so
endif
ifeq ($(OS),Windows_NT)
    EXTENSION = .dll
else
    EXTENSION = .so
endif

GRFF_DEM_Transfer: Coulomb.o DEM.o ExtMath.o FF.o getparms.o GR.o ExtInterface.o MWtransfer.o Neutrals.o Plasma.o Zeta.o GXcode.o
	$(CC) $(LDFLAGS) -shared -o GRFF_DEM_Transfer.so Coulomb.o DEM.o ExtMath.o FF.o getparms.o GR.o ExtInterface.o MWtransfer.o Neutrals.o Plasma.o Zeta.o GXcode.o

Coulomb.o: Coulomb.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX Coulomb.cpp

DEM.o: DEM.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX DEM.cpp

ExtMath.o: ExtMath.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX ExtMath.cpp

FF.o: FF.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX FF.cpp

getparms.o: getparms.cpp
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX getparms.cpp

GR.o: GR.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX GR.cpp

ExtInterface.o: ExtInterface.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX ExtInterface.cpp

MWtransfer.o: MWtransfer.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX MWtransfer.cpp

Neutrals.o: Neutrals.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX Neutrals.cpp
 
Plasma.o: Plasma.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX Plasma.cpp

Zeta.o: Zeta.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX Zeta.cpp

GXcode.o: GXcode.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h GXcode.h ExtInterface.h
	$(CC) $(FLAGS) $(CPPFLAGS) -c -O3 -fPIC -D  LINUX GXcode.cpp

install:
	cp *$(EXTENSION) ../binaries

clean:
	rm -f *.o *$(EXTENSION)