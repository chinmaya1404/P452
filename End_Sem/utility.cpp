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


//MONTE_CARLO

float Monte_Carlo(float (*f)(float), float a, float b, int n, float  &s){
	float x;
	float tot = 0;
	float tot_1 = 0;
	for(long int i = 0; i < n; i++){
		x = a + (b - a)*((float) rand()/RAND_MAX); 
		tot += f(x);					    //f_n
		tot_1 += pow(f(x),2);				//sigma^2
	}
	float f_n = ((b - a)/n) * tot;			//find f_n
	float s_f2 = tot_1/n - pow(tot/n , 2);	//find sigma^2
	s = sqrt(s_f2);						    //find sigma
	return f_n;
}

//EULER_EXPLICIT

double euler_explicit(double (*f)(double, double), double x_0, double y_0, double xmax, double dx, FILE *file){
	double x = x_0/dx, y = y_0;				    //converting x to x/dx 
	while(x <= xmax/dx){
		fprintf(file,"%lf	%lf\n",x*dx,y);		//print x and y
		y += dx*f(x*dx, y);
		x += 1;
	}
return y;
}

//RK_4

double RK_4(double (*f)(double, double), double (*g)(double, double), double x_0, double y_0, double z_0, double h, string filename, double xmin, double xmax){
	
	FILE *file;
	file = fopen(filename.c_str(), "w");
	double k1y, k2y, k3y, k4y, k1z, k2z, k3z, k4z;
	double x = x_0, y = y_0, z = z_0;
	float h1 = -abs(h);						            //decreasing x
	for(x = x_0/h; x > xmin/h; x -= 1){				    //changing x to x/h to avoid floating increment
		
									                    //find k1, k2, k3. k4
		k1y = h1*f(z, x*h);
		k1z = h1*g(z, x*h);

		k2y = h1*f(z+k1z/2, x*h+h1/2);
		k2z = h1*g(z+k1z/2, x*h+h1/2);		
		
		k3y = h1*f(z+k2z/2, x*h+h1/2);
		k3z = h1*g(z+k2z/2, x*h+h1/2);

		k4y = h1*f(z+k3z/2, x*h+h1);
		k4z = h1*g(z+k3z/2, x*h+h1);

		y += (k1y + 2*k2y + 2*k3y + k4y)/6;
		z += (k1z + 2*k2z + 2*k3z + k4z)/6;

		fprintf(file,"%lf	%lf\n",x*h,y);	
	}

	x = x_0, y = y_0, z = z_0;				        //increasing x
	h = abs(h);
	for(x = x_0/h; x <= xmax/h; x += 1){			//changing x to x/h to avoid floating increment
		
								                    //find k1, k2, k3. k4
		k1y = h*f(z, x*h);
		k1z = h*g(z, x*h);

		k2y = h*f(z+k1z/2, x*h+h/2);
		k2z = h*g(z+k1z/2, x*h+h/2);		
		
		k3y = h*f(z+k2z/2, x*h+h/2);
		k3z = h*g(z+k2z/2, x*h+h/2);

		k4y = h*f(z+k3z/2, x*h+h);
		k4z = h*g(z+k3z/2, x*h+h);

		y += (k1y + 2*k2y + 2*k3y + k4y)/6;
		z += (k1z + 2*k2z + 2*k3z + k4z)/6;
		fprintf(file,"%lf	%lf\n",x*h,y);
	}

	fclose(file);
return y;
}

//INTERPOLATION

double interpolation(double a, double b, double ya, double yb, double y_1){
	double t = a + (b -a)/(yb - ya)*(y_1 - ya);
return t;
}

//CHECK SHOOTING

void check_shooting(double (*f)(double, double), double (*g)(double, double), double x_0, double y_0,double x_1, double y_1, double slope, double dx, string filename){
	double y = RK_4(f,g,x_0,y_0,slope,dx,filename,x_0,x_1);
	if(abs(y - y_1) < 0.001){					//error accepted 0.001
		cout<<"The solution"<<endl;
	}
	else if(y < y_1){						    //Under shoot
		cout<<"Under Shoot"<<endl;
	}
	else if(y > y_1){						    //Over shoot
		cout<<"Over Shoot"<<endl;
	}
}

//Shooting Method

void shooting(double (*f)(double, double), double (*g)(double, double), double x_0, double y_0, double x_1, double y_1, double lb, double ub, double dx, string filename){
										                    //find yl and yu
	double yl = RK_4(f,g,x_0,y_0,lb,dx,filename,x_0,x_1);
	double yu = RK_4(f,g,x_0,y_0,ub,dx,filename,x_0,x_1);

	double yc = yl;
	if(abs(y_1 - yl) <= 0.0001){					    	//for any boundary value as solution
		RK_4(f,g,x_0,y_0,yl,dx,filename,x_0,x_1);
	}
	else if(abs(y_1 - yu) <= 0.0001){
		RK_4(f,g,x_0,y_0,yu,dx,filename,x_0,x_1);
	}
	while(abs(yc - y_1) >= 1e-4 ){				    	    //loop for interpolation
			double c;
			if((yl < y_1)&&(yu > y_1)){
				c = interpolation(lb, ub, yl, yu, y_1);
				yc = RK_4(f,g,x_0,y_0,c,dx,filename,x_0,x_1);
				if(yc >= y_1){
					yu = yc;
				}
				else{
					yl = yc;
				}
			}
			else if ((yl > y_1)&&(yu < y_1)){
				c = interpolation(ub, lb, yu, yl, y_1);
				yc = RK_4(f,g,x_0,y_0,c,dx,filename,x_0,x_1);
				if(yc >= y_1){
					yl = yc;
				}
				else{
					yu = yc;
				}
			}
		
	}
}

//Diffusion Equation
void Forward_Diffusion(double (*g)(double), double (*a)(double), double (*b)(double), double D, double dx, double dt, double L, double t, FILE* file){
    double alpha = D*dt/(dx*dx);
    int n = int(L/dx);
    //cout<<t/dt<<endl;
    double v[n+1];
    for(int i = 0; i <= n; i++){
        if(i == 0)
            v[i] = a(0);
        else if(i == n)
            v[i] = b(0);
        else
            v[i] = g(i*dx);
    }
    
    double u[n+1];

    assign((double*) u, (double*) v, 1, n+1);
    //print((double*) u,1,n+1);
    

    for(int i = 0; i*dt <= t; i++){
        for(int j = 0; j <= n; j++){
            //fprintf(file,"%lf %lf %lf\n",i*dt, j*dx, u[j]);
            if(j == 0)
                v[j] = a(i*dt);
            else if(j == n)
                v[j] = b(i*dt);
            else
                v[j] = alpha*(u[j-1]+u[j+1])+(1-2*alpha)*u[j];
        }
        if(i == 0){
            for(int j = 0; j <= n; j++)
            fprintf(file,"%lf %lf\n",j*dx, u[j]);
        }
        else if(i == 10){
            for(int j = 0; j <= n; j++)
            fprintf(file,"%lf %lf\n",j*dx, u[j]);
        }
        else if(i == 20){
            for(int j = 0; j <= n; j++)
            fprintf(file,"%lf %lf\n",j*dx, u[j]);
        }
        else if(i == 50){
            for(int j = 0; j <= n; j++)
            fprintf(file,"%lf %lf\n",j*dx, u[j]);
        }
        else if(i == 100){
            for(int j = 0; j <= n; j++)
            fprintf(file,"%lf %lf\n",j*dx, u[j]);
        }
        else if(i == 200){
            for(int j = 0; j <= n; j++)
            fprintf(file,"%lf %lf\n",j*dx, u[j]);
        }
        else if(i == 500){
            for(int j = 0; j <= n; j++)
            fprintf(file,"%lf %lf\n",j*dx, u[j]);
        }
        assign((double*) u, (double*) v, 1, n+1);
    }


/*
    //cout<<n<<endl;
    for(int i = 0; i <= n; i++)
        v[i] = g(i*dx);
    
    for(int i = 0; i <= n; i++){
        for(int j = 0; j <= n; j++){
            if(i == j)
                A[i][j] = 1-2*alpha;
            else if((i == j+1)||(i == j-1))
                A[i][j] = alpha;
            else
                A[i][j] = 0;
        }
    }
    //print((double*) A,n+1,n+1);
    
    for(int i = 0; i*dt <= t; i++){
        double u[n+1];
        mm((double*) A, (double*) v, (double*) u, n+1, n+1, 1);
        assign((double*) v, (double*) u, 1, n+1);
        for(int j = 0; j <= n; j++)
            fprintf(file,"%lf	",v[j]);
        fprintf(file,"\n");
    }
*/
}

void Backward_Diffusion(double (*g)(double), double (*a)(double), double (*b)(double), double D, double dx, double dt, double L, double t, FILE* file){
    double alpha = dt/(dx*dx);
    int n = int(L/dx);
    double v[n+1], A[n-1][n-1], u[n-1];

    for(int i = 0; i < n-1; i++)
        u[i] = v[i+1];

    //cout<<n<<endl;
    for(int i = 0; i <= n; i++){
        if(i == 0)
            v[i] = a(0);
        else if(i == n)
            v[i] = b(0);
        else
            v[i] = g(i*dx);
    }

    for(int i = 0; i < n-1; i++)
        u[i] = v[i+1];
    
    for(int i = 0; i < n-1; i++){
        for(int j = 0; j < n-1; j++){
            if(i == j)
                A[i][j] = 1+2*alpha;
            else if((i == j+1)||(i == j-1))
                A[i][j] = -alpha;
            else
                A[i][j] = 0;
        }
    }
    //print((double*) A,n+1,n+1);

    double B[n-1][n-1], A1[n-1][n-1], x[n-1], b1[n-1];
    assign((double*) A1, (double*) A, n-1, n-1);
    // print((double*) arr, n,n);
    for(int j = 0; j < n-1; j++){
        for(int k = 0; k < n-1; k++)
            x[k] = 1;
        assign((double*) A, (double*) A1, n-1, n-1);
        for(int k = 0; k < n-1; k++){
            if(k == j)
                b1[k] = 1;
            else
                b1[k] = 0;
        }
        CG_method((double*) A,(double*) b1, (double*) x, n-1, 1e-4);
        for(int k = 0; k < n-1; k++)
            B[k][j] = x[k];
        // print((double*) b, 1,n);
    }
    // print((double*) arr, n,n);
    //print((double*) B, n+1, n+1);

    for(int i = 0; i*dt <= t; i++){
        double w[n-1], z[n-1];
        for(int j = 0; j < n-1; j++){
            if(j == 0)
                w[j] = alpha*a((i+1)*dt);
            else if(j == n-2)
                w[j] = alpha*b((i+1)*dt);
            else
                w[j] = 0;
        }

        ma((double*) u, (double*) w, (double*) u, n-1, 1);

        mm((double*) B, (double*) u, (double*) z, n-1, n-1, 1);

        assign((double*) u, (double*) z, 1, n-1);
        for(int j = 0; j <= n; j++){
            if(j == 0)
                fprintf(file,"%lf %lf, %lf\n",j*dx, i*dt,a(i*dt));
            else if(j == n)
                fprintf(file,"%lf %lf, %lf\n",j*dx, i*dt,b(i*dt));
            else
                fprintf(file,"%lf %lf, %lf\n",j*dx, i*dt,u[j]);
        }
        //fprintf(file,"\n");
    }
}

//RANDOM_INT
int random_int(int min, int max){
    return min + rand()%(max - min + 1); 					
}


//RANDOM_DOUBLE

double random_double(double min, double max){
    return min + (max - min)*((double) rand()/RAND_MAX);
}





/*
    srand(time(NULL));
	FILE *file1, *file2;
	string file = "random.txt";					        //file "random.txt" for storing the positions inside the random_walk function
	int N;
	file1 = fopen("random_walk.txt","w");				//file "random_walk.txt" to store the step number and final position related data
	fprintf(file1,"%s	%s		%s		%s		%s		%s\n", "N","sqrt(N)" ,"R", "R_rms", "x_av", "y_av");
	for(int i = 0, N = 250; i < 5; i++, N *= 2){	    //loop for different number of steps
		float R_tot = 0,x_tot = 0, y_tot = 0, x_2_tot = 0, y_2_tot = 0;
		for(int i = 0; i < 100; i++){				    //loop for doning the random walk for constant N for 100 times
			double x = 0, y = 0;
			random_walk(&x, &y,N,1,file);			    //call the function
			R_tot += sqrt(x*x +y*y);
			x_tot += x;
			y_tot += y;
			x_2_tot += x*x;
			y_2_tot += y*y; 
		}
		float R_mean = R_tot/100;				            //calculate R_mean
		float R_rms = sqrt((x_2_tot/100) + (y_2_tot/100));	//calculate R_rms
		float x_av = x_tot/100;					            //calculate x_av
		float y_av = y_tot/100;					            //calculate y_av
		fprintf(file1,"%d	%lf	%lf	%lf	%lf	%lf\n", N, sqrt(N), R_mean, R_rms, x_av, y_av);		//store data in "random_walk.txt"
		
	}
	
	fclose(file1);
	remove(file.c_str());						            //remove the file "random.txt"
	for(int i = 0, N = 250; i < 5; i++, N *= 2){			//loop for five N values
		for(int i = 1; i <= 5; i++){				        //loop for storing positions for 5 walks
			double x = 0, y = 0;
			std::string file_name = "random_walk_" + std::to_string(N)+"_" + std::to_string(i) + ".txt";
			file2 = fopen(file_name.c_str(),"w");
			random_walk(&x, &y,N,1, file_name);		        //call the function
		}
	}
	fclose(file1);
*/

