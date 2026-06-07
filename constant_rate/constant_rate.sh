gfortran -O3 pp2c.f randgen.f
./a.out
cp cat.dat c1.dat
cp cat2.dat c3.dat
gfortran legge_01.f90
./a.out
gfortran -O3 pp2d.f randgen.f
./a.out
cp cat.dat c2.dat
gfortran -O3 trasforma_file_ingresso.f90
./a.out
cp c4.dat cat2.dat
rm catq2.dat
gfortran -O3 legge_qe_fit_exp3.f90
./a.out
gfortran legge_02.f90
./a.out
gfortran -O3 legge_refilling6_runningthreshold.f90
./a.out
cp catq2.dat c1.dat
cp cat2.dat c2.dat
cp cat_filled.dat c3.dat
gfortran -O3 trasforma_file_ingresso.f90
./a.out
gfortran legge_2.f90
./a.out
gfortran legge_1.f90
./a.out
gfortran legge_1un.f90
./a.out
