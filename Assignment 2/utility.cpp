#include<iostream>
#include<cmath>
#include<complex>
#include<iomanip>
#include<fstream>
using namespace std;
typedef complex<double> dcomp;

//PRINT

void print(double * arr,int n, int m){				
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			cout<<*(arr + i * m + j)<<" ";	
		}
		cout<<endl;	
	}
cout<<endl;	
}

//IMPORT

void import(double *arr,string filename, int m, int n){			
    ifstream file(filename);					
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
 			file >> *(arr+i*n+j);
	file.close();	
}

//ROW_COUNT

int row_count(string file_name){
	int count = 0;
    string line;
    /* Creating input filestream */ 
    ifstream file(file_name);
    while (getline(file, line))
        count++;
	file.close();
return count;	
}

//SWAP

void swap(double * x, double * y){
    double c = *x;
    *x = *y;
    *y = c;
}

//ASSIGN

void assign(double *arr1, double *arr2, int m, int n){
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++)
            *(arr1 + i*n +j) = *(arr2 + i*n + j);
    }
}

//TRANSPOSE

void transpose(double* Mat, double* Mat_T, int m, int n){
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++)
            *(Mat_T + i*m + j) = *(Mat + j*n + i);
    }
}

//MATRIX MULTIPLICATION

void mm(double *arr1, double *arr2, double *arr, int m, int l, int n){		//arr1[m*l], arr2[l*n], arr[m*n]
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			double sum = 0;
			for(int k = 0; k < l; k++)
				sum += (*(arr1 + i * l + k)) * ( *(arr2 + k * n + j));
			*(arr + i * n + j) = sum;
		}
	}
}

//COMPLEX MATRIX MULTIPLICATION

void cmm(dcomp *arr1, dcomp *arr2, dcomp *arr, int m, int l, int n){		//arr1[m*l], arr2[l*n], arr[m*n]
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			dcomp sum = {0,0};
			for(int k = 0; k < l; k++){
				sum += (*(arr1 + i * l + k)) * ( *(arr2 + k * n + j));
			}
			*(arr + i * n + j) = sum;
		}
	}
}

//MATRIX ADDITION

void ma(double * Mat1, double* Mat2, double* Mat, int m, int n){
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++)
            *(Mat + i*n + j) = *(Mat1 + i*n + j) + *(Mat2 + i*n + j);
    }
}

//MATRIX SUBTRACTION

void ms(double * Mat1, double* Mat2, double* Mat, int m, int n){
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++)
            *(Mat + i*n + j) = *(Mat1 + i*n + j) - *(Mat2 + i*n + j);
    }
}

//NORM OF VECTOR

double norm(double * x, int n){
    double sum = 0;
    for(int i = 0; i < n; i++){
        sum += pow(*(x + i),2);
    }
    return sqrt(sum);
}

//NORM OF MATRIX

// double matnorm(double * arr, int n){
//     double sum = 0;
//     for(int i = 0; i < n; i++){
//         for(int j = 0; j < n; j++)
//             sum += pow(*(arr + i*n + j),2);
//     }
//     return sqrt(sum);
// }

double matnorm(double * arr, int n){
    double sum = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++)
            sum += abs(*(arr + i*n + j));
    }
    return sqrt(sum);
}

//PARTIAL PIVOT

int pp(double *arr,double * b ,int i, int n){			//partial pivot 
	int l;
	if(*(arr+i*n+i) == 0){								//if the pivot is zero then swap row with nonzero pivot 
		for(int j = i+1; j<n; j++){						//search for nonzero pivot after the rows of that pivot
			if(*(arr+j*n+i) != 0){
				for(int k = 0; k<n;k++){				//swap the rows
					swap((arr+i*n+k), (arr+j*n+k));				
				}
				swap((b+i), (b+j));
				l = j;
				j = n;									//end the process to search
			}
			else{
				for(int k = i+1; (*(arr+j*n+k) == 0) && k<n; k++){		//check across the row to know whether no solution or infinite solution
					if(k == (n-1))
						exit(0);	
				}				
			}
		}
		if(i == n-1)
			exit(0);
	}
	return l;
}

//GAUSS JORDAN

void gj(double *arr, double* b, int n){
    for(int i = 0; i < n; i++){
        pp((double*) arr, (double*) b, i, n);
        double pivot = *(arr+i*n+i);
        for(int j = 0; j < n; j++){				//divide the ith row by pivot
            *(arr+i*n+j) = *(arr+i*n+j)/pivot;
        }
        *(b+i) = *(b+i)/pivot;					
        
        for(int k = 0; k < n; k++){				//make the elements of the ith column to be zero by row-k <-> row-k - A[i][k]*row-i
            if(k == i){
                continue;						//avoid the i-th row
            }
            double q = *(arr+k*n+i);						
            for(int l =0; l<n; l++){	
                *(arr+k*n+l) = *(arr+k*n+l) - q*(*(arr+i*n+l));
            }
            *(b+k) = *(b+k) - q*(*(b+i));
        }
    }
}

void invgj(double* arr, double * inv, int n){
    double arr1[n][n], b[n];
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
            *(inv + k*n + j) = b[k];
    }
}

/*
    int n = 3;
    double inv[n][n], arr1[n][n], arr[n][n] = {{1,2,4},{2,0,1},{2,7,3}}, b[n] = {1,3,2};
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
*/


//LU DECOMPOSITION

void LU_decomposition(double * arr,double * b, int n){			
	for(int i = 0; i < n; i++){				
		pp((double*) arr,(double*) b, i , n);			
		for(int j = 0; j < n; j++){
			if(i == j){
				double t = *(arr+i*n+i);
				double sum = 0;
				for(int k = 0; k < i; k++)
					sum += *(arr+ i*n + k)* *(arr+ k*n + i); 
				*(arr + i*n+i) = *(arr + i*n+i) - sum;
				if(*(arr+i*n+i) == 0){					//for Ujj == 0
					int l = pp((double*) arr,(double*) b, i , n);
					*(arr+l*n+i) = t;
				}
			}
			else if(i > j){					// for L matrix
				double sum =0;
				for(int k = 0; k < j; k++)
					sum += *(arr+ i*n + k)* *(arr+ k*n + j); 
				*(arr+ i*n + j) = (*(arr+ i*n + j) - sum)/(*(arr + j*n + j));
			}
			else if(i < j){					// for U matrix
				double sum =0;
				for(int k = 0; k < i; k++)
					sum += *(arr+ i*n + k)* *(arr+ k*n + j); 
				*(arr+ i*n + j) = *(arr+ i*n + j) - sum;
			}
		}
	}
}

//FORWARD_BACKWARD

void forward_backward(double * arr, double* b, double* x, int n){		
	double y[n];
	for(int i = 0; i < n; i++){
		double sum = 0;
		for(int j = 0; j < i; j++){
			sum += *(arr+ i*n+j) * y[j];
		}
		y[i] = *(b+i) - sum;
	}
	for(int i = n-1; i >= 0 ; i--){
		double sum = 0;
		for(int j = i+1; j < n; j++){
			sum += *(arr+i*n + j)* *(x+j);
		}
		*(x+i) = (y[i] - sum)/(*(arr+ i*n+i));	
	}
}

/*
    int n = 3;
    double inv[n][n], arr1[n][n], arr[n][n] = {{1,2,4},{2,0,1},{2,7,3}}, b[n] = {1,3,2}, x[n];
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
        LU_decomposition((double*) arr, (double*) b, n);
        forward_backward((double*) arr, (double*) b, (double*) x, n);
        for(int k = 0; k < n; k++)
            inv[k][j] = x[k];
        // print((double*) b, 1,n);
    }
    // print((double*) arr, n,n);
    print((double*) inv, n, n);
*/

//IS SYMMETRIC

bool Is_symmetric(double * arr, int n){
	for(int i = 0; i < n; i++){
		for(int j = i + 1; j < n; j++){
			if(*(arr + i * n + j) == *(arr + j * n + i))			//Check
				int doSomething();
			else
			return 0;												//Break case
		}
	}
	return 1;
}

//IS ANTISYMMETRIC

bool Is_antisymmetric(double * arr, int n){
	for(int i = 0; i < n; i++){
		for(int j = i + 1; j < n; j++){
			if(*(arr + i * n + j) == -*(arr + j * n + i))			//Check
				int doSomething();					
			else
			return 0;												//Break case
		}
	}
	return 1;
}

//CHOLESKY DECOMPOSITION

void cholesky(double * arr, double *L, int n){
	for (int i = 0; i < n; i++) {
    	for (int j = 0; j <= i; j++) {
        	float sum = 0;
        	for (int k = 0; k < j; k++)
            	sum += *(L + i*n + k) * *(L + j *n + k);

        	if (i == j)
            	*( L + i*n + j) = abs(sqrt(*(arr + i*n +i) - sum));								//Fill diagonal terms
        	else
            	*( L + i*n + j) = (1.0 / *( L + j*n + j) * (*(arr + i*n + j) - sum));			//Fill off-diagonal terms
    	}
	}
}

//POSITIVE DEFINITE

bool Is_positive_definite(double * arr, int n){
	double L[n][n];
	cholesky((double*) arr, (double*) L, n);
	double L_T[n][n], arr1[n][n];
	transpose((double*) L, (double*) L_T, n, n);												//Find L^T
	mm((double*) L, (double*) L_T, (double*) arr1, n, n, n);
	for(int i = 0; i < n; i++){																	//Check A == L * L^T
		for(int j = 0; j < n; j++){
			if(*(arr + i*n + j) == arr1[i][j])
				continue;
			else
				return 0;
		}
	}
	return 1;
}

/*
	int n = 3;
    double L[n][n], arr[n][n] = {{4,12,-16},{3,37,-43},{-16,-43,98}};
	print((double*) arr, n, n);
    cholesky((double*) arr, (double*) L, n);
	print((double*) L, n, n);
	cout<<Is_positive_definite((double*) arr, n)<<endl;
*/

//GAUSS SEIDEL

void Gauss_Seidel(double * Mat, double * b, double * x, int n, double e){
    double sum1, error = 10;
    int count = 0;
    for(int k = 1; (k < 10000)&&(error > e); k++){
        count++;
        sum1 = 0;
        //cout<<k<<endl;
        for(int i = 0; i < n; i++){
            double sum = *(b + i);
            for(int j = 0; j < n; j++){
                if(i == j)
                    continue;
                sum -= *(Mat + i*n + j)* (*(x + j));
            }
            double c = *(x + i);
            *(x + i) = 1/(*(Mat + i*n + i)) * sum;
            sum1 += pow(*(x + i) - c,2);
        }
        error = sqrt(sum1);
    }
    //cout<<count<<endl;
}

/*
    double arr[3][3] = {{20,-4,2},{-4,14,2},{2,2,5}};
    double b[3] = {2,1,3};
    double x[3] = {1,1,1};
    cout<<Isdiagonallydominant((double*)arr, 3)<<endl;
    Gauss_Seidel((double*) arr, (double*) b, (double*) x, 3, 1e-4);
    print((double*) x, 3, 1);
*/

//IS STRICTLY DIAGONALLY DOMINANT

bool Isdiagonallydominant(double * Mat, int n){
    double sum;
    for(int i = 0; i < n; i++){
        sum = 0;
        for(int j = 0; j < n; j ++){
            if(i ==j)
                continue;
            else
                sum += abs(*(Mat + i *n + j));
        }
        if(sum > abs(*(Mat + i*n + i)))
            return 0;
    }
    return 1;
}

//JACOBI METHOD

void Jacobi( double * Mat, double *b, double *x, int n, double e){
    double y[n], sum1, error = 10;
    int count = 0;
    for(int k = 1; (k < 10000)&&(error > e); k++){
        count++;
        sum1 = 0;
        for(int i = 0; i < n; i++){
            double sum = *(b + i);
            for(int j = 0; j < n; j++){
                if(i == j)
                    continue;
                sum -= *(Mat + i*n + j)* (*(x + j));
            }
            y[i] = 1/(*(Mat + i*n + i)) * sum;
        }
        ms((double*) x, (double*) y, (double*) x, 1, n);
        error = norm((double *) x, n);
        assign((double*) x, (double*) y, 1, n);
    }
    // cout<<count<<endl;
}

/*
    double arr[3][3] = {{20,-4,2},{-4,14,2},{2,2,5}};
    double b[3] = {2,1,3};
    double x[3] = {1,1,1};
    cout<<Isdiagonallydominant((double*)arr, 3)<<endl;
    Jacobi((double*) arr, (double*) b, (double*) x, 3, 1e-4);
    print((double*) x, 3, 1);
*/

//INNER PRODUCT OF VECTORS

double Inner_Product(double * x, double * y, int n){
    double sum = 0;
    for(int i = 0; i < n; i++)
        sum += *(x + i)* *(y + i);
    return sum;
}

//CONJUGATE GRADIENT METHOD

void CG_method (double* Mat, double* b, double*x, int n, double e){

    double mul[n], p[n], r[n], r1, alpha, beta;

    mm((double*) Mat, (double*) x, (double*) mul, n, n, 1);         //Ax_0
    ms((double*) b, (double*) mul,(double*) r, 1, n);               //r_0 = Ax_0 - b
    if(norm((double*) r, n) < e)                                    //Check convergence
        return;
    assign((double*) p, (double*) r, 1, n);                         //p_0 = r_0
    
    for(int k = 1; k <= n; k++){
        r1 = Inner_Product((double*) r,(double*)r,n);
        mm((double*) Mat, (double*) p, (double*) mul, n, n, 1);     //mul = Ap_0
        alpha = r1/Inner_Product((double*) p, (double*) mul, n);    //Find alpha
        
        for(int i = 0; i < n; i++){
            *(x + i) = *(x + i) + alpha * p[i];                     //x(k+1) = x(k) + alpha * p(k)
            r[i] = r[i] - alpha * mul[i];                           //r(k+1) = r(k) - alpha * A * p_k
            if(norm((double*) r, n) < e)                            //Check convergence
                return;   
        }

        beta = Inner_Product((double*) r, (double*) r, n)/r1;       //Find beta
        for(int i = 0; i < n; i++)
            p[i] = r[i] + beta * p[i];                              //p(k+1) = r(k+1) + beta * p(k)
    }

}

/*
    int n = 3;
    double arr[n][n] = {{-2,-4,2},{-4,-4,2},{2,2,5}};
    double b[n] = {2,1,3};
    double x[n] = {1,1,1};
    CG_method((double*) arr,(double*) b, (double*) x, 3, 1e-4);
    print((double*) x, 1,3);
*/
//POWER METHOD

void Power_Method(double* Mat, double * ev, double* v, int n, double e){
    double x[n], Mat1[n][n], Mat2[n][n], Mat3[n][n];
    double y[n], x1[n];
    for(int l = 0; l < 5; l++){
        for(int i = 0; i < n; i++)                                          //Initial guess of eigenvector
            x[i] = i;

        assign((double*) y, (double*) x, 1, n);
        double lambda = 1, lambda1, error = 1;
        for(int k = 1; (k < 100)&&( error > e); k++ ){
            mm((double*) Mat, (double*) x, (double*) x1, n, n, 1);
            assign((double*) x, (double*) x1, 1, n);
            lambda1 = lambda;
            lambda = Inner_Product((double*) x, (double*) y, n)/Inner_Product((double*) y, (double*) y, n);
            assign((double*) y, (double*) x, 1, n);
            error = abs(lambda1 - lambda);                                  //Check error
        }
        *(ev + l) = lambda; 
                                                                            //Eigenvalue
        double c = norm((double*) x, n);
        for(int i = 0; i < n; i++){
            x[i] = x[i]/c;                                                  //Normalize v
            *(v + l*n + i) = x[i];                                          //Eigenvector
        }
        mm((double*) x, (double*) x, (double*) Mat1, n , 1, n);             //U = x x^T
        // print((double*) Mat1, n,n);
        transpose((double*) Mat1, (double*) Mat2, n, n);                    //Find U^T
        mm((double*) Mat1, (double*) Mat2, (double*) Mat3, n , n, n);       //Find UU^T
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++)
                *(Mat + i*n + j) -= lambda * Mat3[i][j];                    //Update Mat
        }
    }
}

/*
    int n = 3;
    double arr[n][n] = {{-2,-4,2},{-4,-4,2},{2,2,5}};
    double b[n] = {2,1,3};
    double x[n] = {1,1,1};
    double v[n][n];
    double ev[1][n];
    Power_Method((double*) arr,(double*) ev, (double*) v, 3, 1e-4);
    print((double*) arr, 3,3);
    print((double*) v, 3,3);
    for(int i = 0; i < n; i++){
        cout<<ev[0][i]<<" -> [";
        for(int j = 0; j < n; j++)
            cout<<v[j][i]<<" ";
        cout<<"]"<<endl;
    }
*/
//FIND MAXIMUM OFFDIAGONAL ELEMENT

void Find_Maximum_Offdiagonal(double * Mat, int * index, int n){
    double max = abs(*(Mat + n));
    *(index) = 1;
    *(index + 1) = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < i; j++){
            if(i == j)
                continue;
            else{
                if(abs(*(Mat + i*n + j)) > abs(max)){
                    max = *(Mat + i*n + j);
                    *(index) = i;
                    *(index + 1) = j;
                    //cout<<max<<endl;
                }
            }
        }
    }
}

//CONSTRUCT GIVENS ROTATION MATRIX

void Construct_S(double * s, int* index, double  theta, int n){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i == j){
                if((*index == i)||(*(index + 1) == i))
                    *(s + i*n + j) = cos(theta);
                else
                    *(s + i*n + j) = 1;
            }
            else{
                if((i == *(index))&&(j == *(index + 1)))
                    *(s + i*n + j) = -sin(theta);
                else if((i == *(index + 1))&&(j == *(index)))
                    *(s + i*n + j) = sin(theta);
                else
                    *(s + i*n + j) = 0;
            }
        }
    }
}

//SIMILARITY TRANSFORMATION

void Similarity_Transformation(double* Mat, double* s, int n){
    double s_T[n][n];                                                               
    transpose((double*) s, (double*) s_T, n, n);                                        //Find s^T
    double B[n][n], Mat1[n][n];
    mm((double*) s_T, (double*) Mat, (double*) Mat1, n, n, n);
    mm((double*) Mat1, (double*) s, (double*) Mat, n, n, n);
}

//JACOBI EIGENVALUE METHOD

void Jacobi_Eigenvalue_Method(double * Mat, double * ev, int n, double e){
    
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i == j)
                *(ev + i*n + j) = 1;
            else
                *(ev + i*n + j) = 0;
        }
    }
    double s[n][n], ev1[n][n];
    double error = 1;
    for(int l = 0; (l < 10)&&(error > e); l++){
        double sum = 0;
        int index[2];
        Find_Maximum_Offdiagonal((double*) Mat, (int*) index, n);                       //Find Maximum A_ij
        
        //Find theta
        double b = (*(Mat + index[0]*n + index[0]) - *(Mat + index[1]*n +index[1]))/ *(Mat + index[0]*n + index[1]);
        double t = (-b + sqrt(b*b + 4))/2;
        double theta = atan(t);                                                         //Find theta
        
        Construct_S((double*) s, (int*) index, theta, n);                               //Construct S matrix
        mm((double*) ev, (double*) s, (double*) ev1, n, n, n);
        assign((double*) ev, (double*) ev1, n, n);
        Similarity_Transformation((double*) Mat, (double*) s, n);
        // print((double*) Mat, n,n);
        sum = 0;
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(i == j)
                    continue;
                else
                    sum += pow(*(Mat + i*n + j), 2);
            }
        }
        error = sqrt(sum);
    }
}

/*
int n = 3;
    double arr[n][n] = {{-2,-4,2},{-4,-4,2},{2,2,5}};
    double b[n] = {2,1,3};
    double x[n] = {1,1,1};
    double v[n][n];
    double ev[1][n];
    Jacobi_Eigenvalue_Method((double*) arr,(double*) v, 3, 1e-4);
    print((double*) arr, 3,3);
    print((double*) v, 3,3);
    for(int i = 0; i < n; i++){
        cout<<arr[i][i]<<" -> [";
        for(int j = 0; j < n; j++)
            cout<<v[j][i]<<" ";
        cout<<"]"<<endl;
    }
*/

//I_MID-POINT

float I_mid_point(float (*f)(float),float a, float b, int n){
	float h = (b-a)/n;						//Divide into intervals
	float tot = 0;
	for(int i = 0; i < n; i++ ){					//Sum over the midpoint values
		tot += f((a + h/2) + h*i);
	}
return h*tot;
}

//I_TRAPEZOIDAL

float I_trapezoidal(float (*f)(float),float a, float b, int n){
	float h = (b-a)/n;						//Divide into intervals
	float tot = f(a);
	for(int i = 1; i < n; i++ ){
		tot += 2*f(a + i*h);					//sum of middle terms
	}
	tot += f(b);
	return (h/2)*tot;
}

//I_SIMPSON

float I_simpson(float (*f)(float),float a, float b, int n){
	float h = (b-a)/n;						//Divide into intervals
	float tot = f(a);
	for(int i = 1; i < n; i++ ){
		if((i%2) == 1){
			tot += 4*f(a + i*h);				//*4 for i odd
		}
		else if((i%2) == 0){
			tot += 2*f(a + i*h);				//*2 for i even
		}
	}
	tot += f(b);
return (h/3)*tot;
}

bool k_delta(int i, int j){
    if(i == j)
        return 1;
    else
        return 0;
}

//JACK KNIFE

void Jack_Knife(double (*f)(double), double * arr, int n, double * mu, double * sigma, double * mu_f, double * sigma_f){
    double sum, fsum;
    double arr1[n], farr[n];                                //Construct the Jack Knife Samples
    for(int i = 0; i < n; i++){
        sum = 0;
        for(int j = 0; j < n; j++){
            if(j == i)
                continue;
            else
                sum += arr[j];
        }
        arr1[i] = sum/(n-1);
        farr[i] = f(arr1[i]);
    }

    sum = 0, fsum = 0;                                      //Find Jack Knife Mean
    for(int i = 0; i < n; i++){
        sum += arr1[i];
        fsum += farr[i];
    }
    *mu = sum/n;
    *mu_f = fsum/n;

    sum = 0, fsum = 0;                                      //Find Jack Knife Standard Error
    for(int i = 0; i < n; i++){
        sum += pow(arr1[i] - *mu, 2);
        fsum += pow(farr[i] - *mu_f, 2);
    }
    *sigma = sqrt(sum*(n-1)/n);
    *sigma_f = sqrt(fsum*(n-1)/n);
}

//LINEAR REGRESSION

void Linear_Regression(double * x, double * y, double * sigma, int N, double * a, double * b,double * s_a, double * s_b, double * cov, double * r2, double * chi2){
    double S = 0, S_x = 0, S_y = 0, S_xx = 0, S_yy = 0, S_xy = 0;
    *chi2 = 0;
    for(int i = 0; i < N; i++){
        S += 1/pow(*(sigma + i),2);
        S_x += *(x + i)/pow(*(sigma + i),2);
        S_y += *(y + i)/pow(*(sigma + i),2);
        S_xx += pow(*(x + i),2)/pow(*(sigma + i),2);
        S_yy += pow(*(y + i),2)/pow(*(sigma + i),2);
        S_xy += (*(x + i))*(*(y + i))/pow(*(sigma + i),2);
    }
    double delta = S * S_xx - pow(S_x, 2);
    *a = (S_xx * S_y - S_xy * S_x)/delta;
    *b = (S * S_xy - S_x * S_y)/delta;
    *s_a = sqrt(S_xx/delta);
    *s_b = sqrt(S/delta);
    *cov = -S_x/delta;
    *r2 = pow(S_xy, 2)/(S_xx * S_yy);
    for(int i = 0; i < N; i++){
        *chi2 += pow((*(y + i) - *a - *(x + i)*(*b))/(*(sigma + i)), 2);
    }
    /*
    FILE *file;
	file = fopen("chisqurefit.txt","w");
		fprintf(file,"%s	%d\n", "dof     ", N - 2);		//write in file
		fprintf(file,"%s	%lf\n", "Chi^2      ", chi2);
        fprintf(file,"%s	%lf\n", "Chi^2/dof      ", chi2/(N-2));
        fprintf(file,"%s	%lf     %s      %lf\n", "a      ", a, "±", s_a);
        fprintf(file,"%s	%lf     %s      %lf\n", "b      ", b, "±", s_b);
        fprintf(file,"%s	%lf\n", "r^2        ",r2);
        fprintf(file,"%s	%lf\n", "cov()      ", cov);
	fclose(file);
    */
}

/*
    int N = 16;
    double arr[N][3], x[N], y[N], sigma[N];
    import((double*) arr, "v1.txt", N, 3);
    for(int i = 0; i < N; i++){
        x[i] = arr[i][0];
        y[i] = log(arr[i][1]);
        sigma[i] = arr[i][2]/arr[i][1];
    }
    double a, b, s_a, s_b, cov,r2,chi2;
    Linear_Regression((double*) x, (double*) y, (double*) sigma, N, &a, &b, &s_a, &s_b, &cov, &r2, &chi2);
    a = exp(a);
    b = -b;
    s_a = a*s_a;
    cout<<"a    = "<<a<<" ± "<<s_a<<endl;
    cout<<"b    = "<<b<<" ± "<<s_b<<endl;
    cout<<"𝜒2/ν = "<<chi2/(N - 2)<<endl;
    cout<<"r2   = "<<r2<<endl;
    cout<<"cov(a,b) = "<<cov<<endl;
*/

//POLYNOMIAL FIT

void Polynomial_Fit(double * x, double * y, double * sigma, double * a, double * cov, int n, int N, double* chi2){                  //n-th order polynomial
    double arr[n+1][n+1];
    for(int i = 0; i < n+1; i++){
        for(int j = 0; j < n+1; j++){
            for(int k = 0; k < N; k++){
                arr[i][j] += pow(*(x + k), i+j)/pow(*(sigma + k), 2);
            }
        }
        for(int k = 0; k < N; k++)
            *(a + i) += (*(y + k))*(pow(*(x + k), i))/pow(*(sigma + k), 2); //change
    }
    
    double arr1[n+1][n+1];
    assign((double*) arr1, (double*) arr, n+1, n+1);
    gj((double*) arr1, (double*) a, n+1);

    double b[n+1], xx[n+1];
    for(int j = 0; j < n+1; j++){
        //assign((double*) arr, (double*) arr1, n+1, n+1);
        for(int k = 0; k < n+1; k++){
            if(k == j)
                b[k] = 1;
            else
                b[k] = 0;
        }
        LU_decomposition((double*) arr, (double*) b, n+1);
        forward_backward((double*) arr, (double*) b, (double*) xx, n+1);
        for(int k = 0; k < n+1; k++)
           *(cov + k*(n+1) + j) = xx[k];
    }
    
    //print((double*) cov, 2,2);

    for(int i = 0; i < N; i++){
        double sum = 0;
        for(int j = 0; j < n+1; j++){
            sum += *(a + j)*pow(*(x + i), j);
        }
        *chi2 += pow((*(y + i) - sum)/(*(sigma + i)), 2);
    }
}

/*
    int N = 16, n = 1;
    double arr[N][3], x[N], y[N], sigma[N];
    import((double*) arr, "v1.txt", N, 3);
    for(int i = 0; i < N; i++){
        x[i] = arr[i][0];
        y[i] = log(arr[i][1]);
        sigma[i] = arr[i][2]/arr[i][1];
    }
    double a[n+1], cov[n+1][n+1],chi2;
    Polynomial_Fit((double*) x, (double*) y, (double*) sigma, (double*) a, (double*) cov, n, N, &chi2);
    a[0] = exp(a[0]);
    a[1] = -a[1];
    print((double*) cov, 2,2);
    double s_a = a[0]*sqrt(cov[0][0]);
    double s_b = sqrt(cov[1][1]);
    cout<<"a    = "<<a[0]<<" ± "<<s_a<<endl;
    cout<<"b    = "<<a[1]<<" ± "<<s_b<<endl;
    cout<<"𝜒2/ν = "<<chi2/(N - n-1)<<endl;
    //cout<<"r2   = "<<r2<<endl;
    cout<<"cov(a,b) = "<<cov[0][1]<<endl;
*/

double fun(int x, int y){
    double m = 0.2;
    return (1/2)*(k_delta((x+1)%20, y) - k_delta((x-1)%20, y) + 2*k_delta( x, y)) + pow(m,2)*k_delta(x, y);
}

//DFT TRANSFORM

void DFT_Transform(dcomp (*fun)(dcomp), dcomp x_0, double h, dcomp * F, int N){
	dcomp i = {0, 1}, o = {0, 0};
	dcomp w = exp(-i*2.0*M_PI/double(N));
	dcomp f[N];
	for(int j = 0; j < N; j++)
		f[j] = fun(x_0 + j * h);
	dcomp W[N][N];
	for(int j = 0; j < N; j++){
		for(int k = 0; k < N; k++)
			W[j][k] = pow( w, j+k);
	}
	cmm((dcomp*) W, (dcomp*) f, (dcomp*) F, N, N, 1);
}

//DFT INVERSE

void DFT_Inverse(dcomp (*Fun)(dcomp), dcomp x_0, double h, dcomp * f, int N){
	dcomp i = {0, 1}, o = {0, 0};
	dcomp w = exp(i*2.0*M_PI/double(N));
	dcomp F[N];
	for(int j = 0; j < N; j++)
		F[j] = Fun(x_0 + j * h);
	dcomp W[N][N];
	for(int j = 0; j < N; j++){
		for(int k = 0; k < N; k++)
			W[j][k] = pow( w, j+k)/double(N);
	}
	cmm((dcomp*) W, (dcomp*) F, (dcomp*) f, N, N, 1);
}
