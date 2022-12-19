#! /bin/bash
# remove current mod and object file
cd ..
rm -f *mod *a
cd source

gfortran -ffixed-line-length-none -ffree-line-length-none  -fPIC  -O  -c \
SingPrec.f90 \
SysGnuLinux.f90 \
NWTC_IO.f90 \
NWTC_Num.f90 \
NWTC_Aero.f90 \
ModMesh.f90 \
NWTC_Library.f90 \
SharedInflowDefs.f90 \
HHWind.f90 \
FFWind.f90 \
HAWCWind.f90 \
FDWind.f90 \
CTWind.f90 \
UserWind.f90 \
InflowWindMod.f90 \
SharedTypes.f90 \
AeroMods.f90 \
GenSubs.f90 \
AeroSubs.f90 \
AeroDyn.f90 \
DISCON.f90 \

ar ru libAeroDyn.a *.o

mv -f *mod libAeroDyn.a ../ 
rm *o 