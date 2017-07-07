make:
	gfortran Matrix.f90 -llapack
	./a.out
	gnuplot  -persist -e "set xlabel 'g'; set ylabel 'E'; f(x)=0;plot 'eigen.dat' u 1:2  title '1st' w lines, 'eigen.dat' u 1:3 title '2nd' w l, 'eigen.dat' u 1:4 title '3rd' w l, 'eigen.dat'  u 1:5 title '4th' w l, 'eigen.dat' u 1:6 title '5th' w l, 'eigen.dat' u 1:7 title '6th'  w l,f(x) title 'E=0'; q"
