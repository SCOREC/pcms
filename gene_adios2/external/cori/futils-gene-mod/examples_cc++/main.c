#include <stdio.h>
#include <string.h>

void readf_(char*, char*, double*, int*, int*, unsigned int, unsigned int);

int main(int argc, char *argv[])
{
    char file[]="para.h5";
    char name[]="/array_col";
    double a[10][10];
    int nx=10, ny=10;

    readf_(file, name, &a[0][0], &nx, &ny, strlen(file), strlen(name));
    for( int i=0; i<nx; i++ ) {
	for( int j=0; j<ny; j++ )
	    printf("%5.1f ", a[i][j]);
	printf("\n");
    }

    return 0;
}
