#include <iostream>
#include <string>

using namespace std;

extern "C" {
  void readf_(char*, char*, double*, int*, int*, unsigned int, unsigned int);
}

int main(int argc, char *argv[])
{
    char file[]="para.h5";
    char name[]="/array_col";
    double a[10][10];
    int nx=10, ny=10;

    readf_(file, name, &a[0][0], &nx, &ny, strlen(file), strlen(name));
    for( int i=0; i<nx; i++ ) {
	for( int j=0; j<ny; j++ )
	  cout << a[i][j] << " ";
	cout << "\n";
    }

    return 0;
}
