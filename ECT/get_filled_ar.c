#include <stdio.h>
FILE * pFile;
	
void get_filled_ar (float ar []) {


	float f;
	int ret = fscanf(pFile, "%f %f %f %f %f", &ar[0], &ar[1], &ar[2], &ar[3], &ar[4]);
	return ;

}

void open() {

	char status;
	float f;
	pFile = fopen ("data.dat","r+");
	return ;

}
