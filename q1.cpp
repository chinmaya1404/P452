#include "utility.cpp"
int main(){
    int n = 6;
    double arr1[n][n], b1[n], arr[n][n] = {{1,-1,4,0,2,9},{0,5,-2,7,8,4},{1,0,5,7,3,-2},{6,-1,2,3,0,8},{-4,2,0,5,-5,3},{0,7,-1,5,4,-2}}, b[n] = {19,2,13,-7,-9,2};
    
    cout<<"Gauss Jordan"<<endl;
    assign((double*) arr1, (double*) arr, n, n);
    assign((double*) b1, (double*) b, n, 1);
    //print((double*) b1, n, 1);
    gj((double*) arr1, (double*) b1, n);
    cout<<"x = "<<endl;
    print((double*) b1, n, 1);

    cout<<"LU Decomposition"<<endl;
    assign((double*) arr1, (double*) arr, n, n);
    assign((double*) b1, (double*) b, n, 1);
    double x[n];
    //print((double*) b1, n, 1);
    LU_decomposition((double*) arr1, (double*) b1, n);
    forward_backward((double*) arr1, (double*) b1, (double*) x, n);
    cout<<"x = "<<endl;
    print((double*) x, n, 1);

    return 0;
}
/*
Gauss Jordan
x =
-1.76182
0.896228
4.05193
-1.61713
2.04191
0.151832 

LU Decomposition
x =
-1.76182
0.896228
4.05193
-1.61713
2.04191
0.151832
*/