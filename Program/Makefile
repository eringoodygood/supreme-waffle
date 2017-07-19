run:
	@reset
	@echo "Please choose the program."
	@echo "Type make first for the first part program and make second for the second part program."

first:
	@gfortran Matrix3.f90 -llapack
	@./a.out
	@gnuplot  -persist -e "set xlabel 'g'; set ylabel 'E'; f(x)=0;plot 'eigen.dat' u 1:2  title '1st' w lines, 'eigen.dat' u 1:3 title '2nd' w l, 'eigen.dat' u 1:4 title '3rd' w l, 'eigen.dat'  u 1:5 title '4th' w l, 'eigen.dat' u 1:6 title '5th' w l, 'eigen.dat' u 1:7 title '6th'  w l,f(x) title 'E=0'; q"

second:
	@gcc -c get_filled_ar.c
	@gfortran get_filled_ar.o Matrix2.f90 -llapack
	@echo ""
	@echo ""
	@echo "Welcome to the second program."
	@echo "The 2-body interaction (in M-scheme) is in the file data.dat, the 1-body is in data2.dat."
	@echo ""
	@echo ""
	@echo ""
	@./a.out
	@gnuplot  -persist -e "set xlabel 'A'; set ylabel 'E [MeV]'; set xrange [16:29]; plot 'energies.dat' u 1:2  title 'Ground state energy' w linespoints ; q"
