cd ..
rm -f *a
cd moorsys
rm -f *a
#gfortran -ffixed-line-length-none -ffree-line-length-none -fPIC -O -c TU11-00075_moorsys.f
gfortran -fPIC -O -c TU11-00075_moorsys.f
ar ru libMoorsys.a TU11-00075_moorsys.o

cp libMoorsys.a ../