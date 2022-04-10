#include "utility.cpp"
int main(){
    int n = 6;
    double arr1[n][n], b1[n], arr[n][n] = {{2,-3,0,0,0,0},{-1,4,-1,0,-1,0},{0,-1,4,0,0,-1},{0,0,0,2,-3,0},{0,-1,0,-1,4,-1},{0,0,-1,0,-1,4}}, b[n] = {-5/3,2/3,3,-4/3,-1/3,5/3};
    
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

    FILE *file;

    file = fopen("q2.txt","w");
    cout<<"Jacobi"<<endl;
    assign((double*) arr1, (double*) arr, n, n);
    assign((double*) b1, (double*) b, n, 1);
    double x1[n] = {1,1,1,1,1,1};
    //cout<<Isdiagonallydominant((double*)arr1, 6)<<endl;
    Jacobi((double*) arr1, (double*) b1, (double*) x1, n, 1e-4, file);
    cout<<"x = "<<endl;
    print((double*) x1, n, 1);
    fclose(file);

    
    cout<<"Jacobi"<<endl;
    assign((double*) arr1, (double*) arr, n, n);
    // print((double*) arr, n,n);
    double c[n], inv[n][n];
    for(int j = 0; j < n; j++){
        assign((double*) arr, (double*) arr1, n, n);
        for(int k = 0; k < n; k++){
            if(k == j)
                c[k] = 1;
            else
                c[k] = 0;
        }
        double x2[n] = {1,1,1,1,1,1};
        file = fopen("q2_J.txt","w");
        //cout<<Isdiagonallydominant((double*)arr1, 6)<<endl;
        Jacobi((double*) arr1, (double*) c, (double*) x2, n, 1e-4, file);
        fclose(file);
        for(int k = 0; k < n; k++)
            inv[k][j] = x2[k];
        // print((double*) b, 1,n);
    }
    // print((double*) arr, n,n);
    cout<<"Inv(A) = "<<endl;
    print((double*) inv, n, n);
    

    
    cout<<"Gauss Seidel"<<endl;
    assign((double*) arr1, (double*) arr, n, n);
    // print((double*) arr, n,n);
    for(int j = 0; j < n; j++){
        assign((double*) arr, (double*) arr1, n, n);
        for(int k = 0; k < n; k++){
            if(k == j)
                c[k] = 1;
            else
                c[k] = 0;
        }
        double x2[n] = {1,1,1,1,1,1};
        file = fopen("q2_G.txt","w");
        //cout<<Isdiagonallydominant((double*)arr1, 6)<<endl;
        Gauss_Seidel((double*) arr1, (double*) c, (double*) x2, n, 1e-4, file);
        fclose(file);
        for(int k = 0; k < n; k++)
            inv[k][j] = x2[k];
        // print((double*) b, 1,n);
    }
    // print((double*) arr, n,n);
    cout<<"Inv(A) = "<<endl;
    print((double*) inv, n, n);

/*
    cout<<"Conjugate Gradient"<<endl;
    assign((double*) arr1, (double*) arr, n, n);
    // print((double*) arr, n,n);
    for(int j = 0; j < n; j++){
        assign((double*) arr, (double*) arr1, n, n);
        for(int k = 0; k < n; k++){
            if(k == j)
                c[k] = 1;
            else
                c[k] = 0;
        }
        double x2[n] = {1,2,3,4,5,6};
        //cout<<Isdiagonallydominant((double*)arr1, 6)<<endl;
        CG_method((double*) arr1, (double*) c, (double*) x2, n, 1e-8);
        for(int k = 0; k < n; k++)
            inv[k][j] = x2[k];
        // print((double*) b, 1,n);
    }
    // print((double*) arr, n,n);
    print((double*) inv, n, n);
*/  

    return 0;
}

/*
Gauss Jordan
x =
-0.194805
0.203463
0.926407
-0.376623
0.0822511
0.502165

LU Decomposition
x =
-0.194805
0.203463
0.926407
-0.376623
0.0822511
0.502165

Jacobi
x =
-0.194603
0.203584
0.926455
-0.376403
0.0823614
0.502217

Jacobi
Inv(A) = 
0.935322 0.870316 0.259957 0.208033 0.415796 0.169058
0.290175 0.580202 0.173284 0.138668 0.277158 0.112673
0.0866417 0.173205 0.320398 0.0563349 0.112605 0.108279
0.208033 0.415796 0.169058 0.935322 0.870316 0.259957
0.138668 0.277158 0.112673 0.290175 0.580202 0.173284
0.0563349 0.112605 0.108279 0.0866417 0.173205 0.320398

Gauss Seidel
Inv(A) = 
0.935193 0.870237 0.259867 0.207915 0.415682 0.168952
0.2901 0.580134 0.173217 0.138583 0.2771 0.112608 
0.0866007 0.173177 0.320367 0.0562969 0.11257 0.108245
0.207897 0.415672 0.168935 0.935166 0.870209 0.25984
0.138575 0.277095 0.1126 0.290088 0.580122 0.173204
0.0562939 0.112568 0.108242 0.0865963 0.173173 0.320362
*/