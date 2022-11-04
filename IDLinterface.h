#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW(int argc, void **argv);
extern "C" __declspec(dllexport) int GET_MW_SLICE(int argc, void **argv);
#else
extern "C" double GET_MW(int argc, void **argv);
extern "C" double GET_MW_SLICE(int argc, void **argv);
#endif
