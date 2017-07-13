#include <stdio.h>
FILE * pFile;
	
void get_filled_ar (int i, float ar []) {
	//printf("hi %i",i);
	float gabage;
	if(i==1){
		int ret = fscanf(pFile, "%f %f %f %f %f", &ar[0], &ar[1], &ar[2], &ar[3], &ar[4]);
	}	
	else if(i==2){
		int ret = fscanf(pFile, "%f %f %f %f %f %f", &ar[0], &gabage, &gabage, &gabage, &gabage,&ar[1]);
	}
	return ;

}



void open(char* name) {
	//printf("hi %s",name);
	pFile = fopen (name,"r+");
	return ;

}
void close() {
	fclose(pFile);
	return ;

}
