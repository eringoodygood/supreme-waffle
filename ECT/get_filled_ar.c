#include <stdio.h>
void get_filled_ar (int ar []) {

	char status;
	float f;
	FILE * pFile;
	pFile = fopen ("data.dat","r+");
	do {
	        int ret = fscanf(pFile, "%i %i %i %i %c", &ar[0], &ar[1], &ar[2], &ar[3], &status);
		if(ret == EOF)
			break;
	}while(1);
	return ;

}
