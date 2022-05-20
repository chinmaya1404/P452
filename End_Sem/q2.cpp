#include "utility.cpp"
#include <functional>

void Polynomial_Fit_Least_Square(double * x, double * y, double * a, double * cov, int n, int N){                  //n-th order polynomial
    double arr[n+1][n+1];
    for(int i = 0; i < n+1; i++){
        for(int j = 0; j < n+1; j++){
            for(int k = 0; k < N; k++){
                arr[i][j] += pow(*(x + k), i+j);
            }
        }
        for(int k = 0; k < N; k++)
            *(a + i) += (*(y + k))*(pow(*(x + k), i)); //change
    }

    double arr1[n+1][n+1];
    assign((double*) arr1, (double*) arr, n+1, n+1);
    
    double b[n+1];
    for(int j = 0; j < n+1; j++){
        assign((double*) arr1, (double*) arr, n+1, n+1);
        for(int k = 0; k < n+1; k++){
            if(k == j)
                b[k] = 1;
            else
                b[k] = 0;
        }
        gj((double*) arr1, (double*) b, n+1);
        for(int k = 0; k < n+1; k++)
           *(cov + k*(n+1) + j) = b[k];
    }

    mm((double*) cov, (double*) a, (double*) b, n+1, n+1, 1);
    assign((double*) a, (double*) b, n+1, 1);

    
}



double P (int i,double x){
    if (i==0) return 1;
    else if (i==1)return x;
    else if (i==2)return (3*pow(x,2)-1)/2;
    else if (i==3)return (5*pow(x,3) -3*x)/2;
    else if (i==4)return (35*pow(x,4) - 30*pow(x,2) +3)/8;
    else if (i==5)return (63*pow(x,5) - 70*pow(x,3) + 15*x)/8;
    else if (i==6)return (231*pow(x,6) - 315*pow(x,4) + 105*pow(x,2) - 5)/16;
    else return 0;
}

void Polynomial_Fit_Least_Square1(double * x, double * y, double * a, double * cov, int n, int N){     //n-th order polynomial
    double arr[n+1][n+1];
    for(int i = 0; i < n+1; i++){
        for(int j = 0; j < i+1; j++){
            for(int k = 0; k < N; k++){
                arr[i][j] += P(j,x[k])*P(i,x[k]);
            }
            arr[j][i] = arr[i][j];
        }
        for(int k = 0; k < N; k++)
            *(a + i) += y[k]*(i,x[k]);
    }

    double arr1[n+1][n+1];
    assign((double*) arr1, (double*) arr, n+1, n+1);
    
    // gj((double*) arr1, (double*) a, n+1);
    // print((double*) arr, n+1, n+1);

    double b[n+1];
    for(int j = 0; j < n+1; j++){
        assign((double*) arr1, (double*) arr, n+1, n+1);
        for(int k = 0; k < n+1; k++){
            if(k == j)
                b[k] = 1;
            else
                b[k] = 0;
        }
    
        gj((double*) arr1, (double*) b, n+1);

        for(int k = 0; k < n+1; k++)
           *(cov + k*(n+1) + j) = b[k];
    }
    mm((double*) cov, (double*) a, (double*) b, n+1, n+1, 1);
    assign((double*) a, (double*) b, n+1, 1);
}


int main(){
    int N = 26, n = 6;
    double arr[N][2], x[N], y[N], cov[n+1][n+1], r2,a1[n+1], invcov[n+1][n+1];
    r2 = 0;
    import((double*) arr, "esem4fit.txt", N, 2);
    for(int i = 0; i < N; i++){
        x[i] = arr[i][0];
        y[i] = arr[i][1];
    }

        r2 = 0;
        Polynomial_Fit_Least_Square1((double*) x, (double*) y, (double*) a1, (double*) cov, n, N);
        // print((double*) a1, n+1, 1);
        // print((double*) cov, n+1, n+1);
        cout<<"a0   = "<<a1[0]<<endl;
        cout<<"a1   = "<<a1[1]<<endl;
        cout<<"a2   = "<<a1[2]<<endl;
        cout<<"a3   = "<<a1[3]<<endl;
        cout<<"a4   = "<<a1[4]<<endl;
        cout<<"a5   = "<<a1[5]<<endl;
        cout<<"a6   = "<<a1[6]<<endl;

        cout<<"Covarience Matrix :"<<endl;
        print((double*) cov, n+1, n+1);
        


    return 0;
}

/*
For order n = 3
a0   = 0.0736485
a1   = 0.00362402
a2   = 0.00854266
a3   = 0.0114262
Covarience Matrix :
0.0387265 -4.03063e-18 -0.00662436 -8.86256e-19
-4.03063e-18 0.109801 6.89458e-19 -0.0253752
-0.00662436 6.89458e-19 0.165609 1.51598e-19
-8.86256e-19 -0.0253752 1.51598e-19 0.217253

For order n = 4
a0   = 0.0696578
a1   = 0.00362402
a2   = -0.0120826
a3   = 0.0114262
a4   = 0.110492
Covarience Matrix :
0.0390716 -3.83362e-18 -0.00484104 -9.47981e-19 -0.00955349
-3.83362e-18 0.109801 1.70767e-18 -0.0253752 -5.45471e-18
-0.00484104 1.70767e-18 0.174826 -1.67416e-19 -0.0493757
-9.47981e-19 -0.0253752 -1.67416e-19 0.217253 1.70901e-18
-0.00955349 -5.45471e-18 -0.0493757 1.70901e-18 0.264513

For order n = 5
a0   = 0.0696578
a1   = 0.00430169
a2   = -0.0120826
a3   = 0.0130837
a4   = 0.110492
a5   = -0.00672697
Covarience Matrix :
0.0390716 -4.01674e-18 -0.00484104 -1.3959e-18 -0.00955349 1.81786e-18
-4.01674e-18 0.112951 1.76572e-18 -0.0176709 -5.5933e-18 -0.0312674
-0.00484104 1.76572e-18 0.174826 -2.54264e-20 -0.0493757 -5.76257e-19
-1.3959e-18 -0.0176709 -2.54264e-20 0.236097 1.37002e-18 -0.076478
-0.00955349 -5.5933e-18 -0.0493757 1.37002e-18 0.264513 1.37574e-18
1.81786e-18 -0.0312674 -5.76257e-19 -0.076478 1.37574e-18 0.310381

For order n = 6
a0   = 0.070032
a1   = 0.00430169
a2   = -0.0101667
a3   = 0.0130837
a4   = 0.114119
a5   = -0.00672697
a6   = -0.0123846
Covarience Matrix :
0.0393985 -4.13184e-18 -0.00316679 -1.39483e-18 -0.00638461 1.85562e-18 -0.0108227
-4.13184e-18 0.112951 1.17639e-18 -0.0176709 -6.70874e-18 -0.0312674 3.80956e-18
-0.00316679 1.17639e-18 0.183398 -1.99706e-20 -0.03315 -3.82867e-19 -0.0554156
-1.39483e-18 -0.0176709 -1.99706e-20 0.236097 1.38035e-18 -0.076478 -3.52677e-20
-0.00638461 -6.70874e-18 -0.03315 1.38035e-18 0.295223 1.74177e-18 -0.104886
1.85562e-18 -0.0312674 -3.82867e-19 -0.076478 1.74177e-18 0.310381 -1.25011e-18
-0.0108227 3.80956e-18 -0.0554156 -3.52677e-20 -0.104886 -1.25011e-18 0.358217


*/