make:
	gfortran Matrix.f90 -llapack
	./a.out
	gnuplot  -persist -e "plot 'eigen.dat' u 1:2  title '1' w lines, 'eigen.dat' u 1:3 w l, 'eigen.dat' u 1:4 w l, 'eigen.dat' u 1:5 w l, 'eigen.dat' u 1:6 w l, 'eigen.dat' u 1:7 w l ; q"
