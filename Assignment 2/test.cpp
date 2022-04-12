#include "utility.cpp"

int main(){
	int n = 4;
	double arr[4][4] = {{21, -2.22045e-16, -5.6, -1.06581e-14}, {-2.22045e-16, 7.7, 0, -2.8336}, {-5.6, 0, 10.4664, 2.22045e-15}, {-1.06581e-14, -2.8336, 2.22045e-15, 11.0106}};
	//double b[4] = {24.118, 13.2345, 9.46837, 7.55944};
	//gj((double*) arr, (double*) b, 4);
	//print((double*) b, 4, 1);
	//int n = 3;
    double inv[n][n], arr1[n][n], b[n];
    assign((double*) arr1, (double*) arr, n, n);
    // print((double*) arr, n,n);
    for(int j = 0; j < n; j++){
        assign((double*) arr, (double*) arr1, n, n);
        for(int k = 0; k < n; k++){
            if(k == j)
                b[k] = 1;
            else
                b[k] = 0;
        }
        gj((double*) arr, (double*) b, n);
        for(int k = 0; k < n; k++)
            inv[k][j] = b[k];
        // print((double*) b, 1,n);
    }
    // print((double*) arr, n,n);
    print((double*) inv, n, n);	
	return 0;
}