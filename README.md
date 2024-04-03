Codes for computing the solar gyroresonance and free-free radio emissions; both the isothermal plasma and the sources described by the differential emission measure (DEM) and differential density metric (DDM) are supported.

Alexey Kuznetsov, February 2021.

The codes are implemented as Windows/Linux libraries callable from IDL (via call_external function) or Python (via the Ctypes interface). See the files CallingConventions.pdf and Diagram.pdf in the "doc" folder for more details and calling conventions, and the folder "examples" for the usage examples. 

The folder "binaries" contains the compiled Windows DLL and Linux SO libraries.

[![DOI](https://zenodo.org/badge/292502263.svg)](https://zenodo.org/badge/latestdoi/292502263)


#### for MacOS users 

need to install Homebrew first, then
```bash
brew install gcc llwv libomp

cd GRFF/source

make
make install
make clean
```