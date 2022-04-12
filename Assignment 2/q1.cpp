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



double phi (int i,double x){
    if (i==0) return 1;
    else if (i==1)return 2*x - 1;
    else if (i==2)return 8*pow(x,2) - 8*x + 1;
    else if (i==3)return 32*pow(x,3) - 48*pow(x,2) + 18*x - 1;
    else return 0;
}

void Polynomial_Fit_Least_Square1(double * x, double * y, double * a, double * cov, int n, int N){     //n-th order polynomial
    double arr[n+1][n+1];
    for(int i = 0; i < n+1; i++){
        for(int j = 0; j < i+1; j++){
            for(int k = 0; k < N; k++){
                arr[i][j] += phi(j,x[k])*phi(i,x[k]);
            }
            arr[j][i] = arr[i][j];
        }
        for(int k = 0; k < N; k++)
            *(a + i) += y[k]*phi(i,x[k]);
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
    int N = 21, n = 3;
    double arr[N][2], x[N], y[N], a[n+1], cov[n+1][n+1], r2,a1[n+1], invcov[n+1][n+1];
    r2 = 0;
    import((double*) arr, "assign2fit.txt", N, 2);
    for(int i = 0; i < N; i++){
        x[i] = arr[i][0];
        y[i] = arr[i][1];
    }
    char q;
    cout<<"Enter bit Number :(a or b) ";
    cin>>q;
    if(q == 'a'){
        Polynomial_Fit_Least_Square((double*) x, (double*) y, (double*) a, (double*) cov, n, N);
        // print((double*) a, n+1, 1);
        // print((double*) cov, n+1, n+1);
        cout<<"a0   = "<<a[0]<<endl;
        cout<<"a1   = "<<a[1]<<endl;
        cout<<"a2   = "<<a[2]<<endl;
        cout<<"a3   = "<<a[3]<<endl;

        cout<<"Covarience Matrix :"<<endl;
        print((double*) cov, n+1, n+1);
        cout<<"Condition Number : 22,030.866"<<endl;
    }

    else if(q == 'b'){
        r2 = 0;
        Polynomial_Fit_Least_Square1((double*) x, (double*) y, (double*) a1, (double*) cov, n, N);
        // print((double*) a1, n+1, 1);
        // print((double*) cov, n+1, n+1);
        cout<<"a0   = "<<a1[0]<<endl;
        cout<<"a1   = "<<a1[1]<<endl;
        cout<<"a2   = "<<a1[2]<<endl;
        cout<<"a3   = "<<a1[3]<<endl;

        cout<<"Covarience Matrix :"<<endl;
        print((double*) cov, n+1, n+1);
        cout<<"Condition Number : 4.798"<<endl;
    }


    return 0;
}

/*
Enter bit Number :(a or b) a
a0   = 0.574659
a1   = 4.72586
a2   = -11.1282
a3   = 7.66868
Covarience Matrix :
0.544043 -3.96041 7.71692 -4.39174
-3.96041 42.8691 -97.3558 60.1489
7.71692 -97.3558 238.277 -154.096
-4.39174 60.1489 -154.096 102.731

Condition Number : 22,030.866

Enter bit Number :(a or b) b
a0   = 1.16097
a1   = 0.393514
a2   = 0.0468498
a3   = 0.239646
Covarience Matrix :
0.055544 2.11889e-17 0.0297186 5.3226e-17
2.11889e-17 0.143456 3.50468e-18 0.0369189
0.0297186 3.50468e-18 0.111445 7.1948e-18
5.3226e-17 0.0369189 7.1948e-18 0.100323

Condition Number : 4.798


The condition number for second case is close to unity so small perturbation in the input data would have small effect on the fitting parameters. But in the first case the condition number is very larger than unity so small perturbation in the input data would have huge effect on the fitting parameters.
*/