GRFF_DEM_Transfer	:	Coulomb.o DEM.o ExtMath.o FF.o getparms.o GR.o ExtInterface.o MWtransfer.o Neutrals.o Plasma.o Zeta.o
				g++ -shared -o GRFF_DEM_Transfer.so Coulomb.o DEM.o ExtMath.o FF.o getparms.o GR.o ExtInterface.o MWtransfer.o Neutrals.o Plasma.o Zeta.o
Coulomb.o		:	Coulomb.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h
				g++ -c -O3 -fPIC -D LINUX Coulomb.cpp
DEM.o			:	DEM.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h
				g++ -c -O3 -fPIC -D LINUX DEM.cpp
ExtMath.o		:	ExtMath.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h
				g++ -c -O3 -fPIC -D LINUX ExtMath.cpp
FF.o			:	FF.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h
				g++ -c -O3 -fPIC -D LINUX FF.cpp
getparms.o		:	getparms.cpp
				g++ -c -O3 -fPIC -D LINUX getparms.cpp
GR.o			:	GR.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h
				g++ -c -O3 -fPIC -D LINUX GR.cpp
ExtInterface.o		:	ExtInterface.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h
				g++ -c -O3 -fPIC -D LINUX ExtInterface.cpp
MWtransfer.o		:	MWtransfer.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h
				g++ -c -O3 -fPIC -D LINUX MWtransfer.cpp
Neutrals.o		:	Neutrals.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h
				g++ -c -O3 -fPIC -D LINUX Neutrals.cpp
Plasma.o		:	Plasma.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h
				g++ -c -O3 -fPIC -D LINUX Plasma.cpp
Zeta.o			:	Zeta.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h MWtransfer.h Neutrals.h Plasma.h Zeta.h
				g++ -c -O3 -fPIC -D LINUX Zeta.cpp
