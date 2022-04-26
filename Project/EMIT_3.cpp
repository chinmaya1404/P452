#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>
#include<string>
#include<complex>
using namespace std;
typedef complex<double> dcomp;

void print_r(double * arr,int n, int m){				
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			cout<<*(arr + i * m + j)<<" ";	//print rounded value upto 3 decimal place
		}
		cout<<endl;	
	}
cout<<endl;	
}

//---------------------------------------------------SWAP------------------------------------------------------------------

void swap(dcomp * a, dcomp* b){
	dcomp c;
	c = *a;
	*a = *b;
	*b = c;
}

//----------------------------------------------------PRINT-----------------------------------------------------------

void print(dcomp * arr,int n, int m){				
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			cout<<*(arr + i * m + j)<<" ";
		}
		cout<<endl;	
	}
cout<<endl;	
}

//-------------------------------------------------ROW_COUNT----------------------------------------------------------
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


//----------------------------------------------------IMPORT------------------------------------------------------------

void import(dcomp *arr,string filename, int m, int n){			
    ifstream file(filename);					
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
 			file >> *(arr+i*n+j);
	file.close();	
}

//ASSIGN

void assign(dcomp *arr1, dcomp *arr2, int m, int n){
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++)
            *(arr1 + i*n +j) = *(arr2 + i*n + j);
    }
}

//----------------------------------------------MATRIX ADDITION------------------------------------------------------

void ma(dcomp *arr1, dcomp *arr2, dcomp *arr, int m, int n){		//arr1[m*l], arr2[l*n], arr[m*n]
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			*(arr + i * n + j) = *(arr1 + i * n + j) + *(arr2 + i * n + j);
		}
	}
}


//----------------------------------------------MATRIX MULTIPLICATION------------------------------------------------------

void mm(dcomp *arr1, dcomp *arr2, dcomp *arr, int m, int l, int n){		//arr1[m*l], arr2[l*n], arr[m*n]
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

int pp(dcomp *arr,dcomp * b ,int i, int n){			    //partial pivot 
	int l;
	if(abs(*(arr+i*n+i)) == 0){								//if the pivot is zero then swap row with nonzero pivot 
		for(int j = i+1; j<n; j++){						//search for nonzero pivot after the rows of that pivot
			if(abs(*(arr+j*n+i)) != 0){
				for(int k = 0; k<n;k++){				//swap the rows
					swap((arr+i*n+k), (arr+j*n+k));				
				}
				swap((b+i), (b+j));
				l = j;
				j = n;									//end the process to search
			}
			else{
				for(int k = i+1; (abs(*(arr+j*n+k)) == 0) && k<n; k++){		//check across the row to know whether no solution or infinite solution
					if(k == (n-1))
						exit(0);	
				}				
			}
		}
		if(i == n-1)
			exit(0);		//edit
	}
	return l;
}

//GAUSS JORDAN

void gj(dcomp *arr, dcomp* b, int n){
    for(int i = 0; i < n; i++){
        pp((dcomp*) arr, (dcomp*) b, i, n);
        dcomp pivot = *(arr+i*n+i);
        for(int j = 0; j < n; j++){				//divide the ith row by pivot
            *(arr+i*n+j) = *(arr+i*n+j)/pivot;
        }
        *(b+i) = *(b+i)/pivot;					
        
        for(int k = 0; k < n; k++){				//make the elements of the ith column to be zero by row-k <-> row-k - A[i][k]*row-i
            if(k == i){
                continue;						//avoid the i-th row
            }
            dcomp q = *(arr+k*n+i);						
            for(int l =0; l<n; l++){	
                *(arr+k*n+l) = *(arr+k*n+l) - q*(*(arr+i*n+l));
            }
            *(b+k) = *(b+k) - q*(*(b+i));
        }
    }
}

//LU DECOMPOSITION

void LU_decomposition(dcomp * arr,dcomp * b, int n){			
	for(int i = 0; i < n; i++){				
		pp((dcomp*) arr,(dcomp*) b, i , n);			
		for(int j = 0; j < n; j++){
			if(i == j){
				dcomp t = *(arr+i*n+i);
				dcomp sum = {0, 0};
				for(int k = 0; k < i; k++)
					sum += *(arr+ i*n + k)* *(arr+ k*n + i); 
				*(arr + i*n+i) = *(arr + i*n+i) - sum;
				if(abs(*(arr+i*n+i)) == 0){					//for Ujj == 0
					int l = pp((dcomp*) arr,(dcomp*) b, i , n);
					*(arr+l*n+i) = t;
				}
			}
			else if(i > j){					// for L matrix
				dcomp sum = {0, 0};
				for(int k = 0; k < j; k++)
					sum += *(arr+ i*n + k)* *(arr+ k*n + j); 
				*(arr+ i*n + j) = (*(arr+ i*n + j) - sum)/(*(arr + j*n + j));
			}
			else if(i < j){					// for U matrix
				dcomp sum = {0, 0};
				for(int k = 0; k < i; k++)
					sum += *(arr+ i*n + k)* *(arr+ k*n + j); 
				*(arr+ i*n + j) = *(arr+ i*n + j) - sum;
			}
		}
	}
}

//FORWARD_BACKWARD

void forward_backward(dcomp * arr, dcomp* b, dcomp* x, int n){		
	dcomp y[n];
	for(int i = 0; i < n; i++){
		dcomp sum = 0;
		for(int j = 0; j < i; j++){
			sum += *(arr+ i*n+j) * y[j];
		}
		y[i] = *(b+i) - sum;
	}
	for(int i = n-1; i >= 0 ; i--){
		dcomp sum = {0,0};
		for(int j = i+1; j < n; j++){
			sum += *(arr+i*n + j)* *(x+j);
		}
		*(x+i) = (y[i] - sum)/(*(arr+ i*n+i));	
	}
}

//RK_4

double RK_4(double (*f)(double, double), double (*g)(double, double), double x_0, double y_0, double z_0, double h, string filename, double xmin, double xmax){
	
	FILE *file;
	file = fopen(filename.c_str(), "w");
	double k1y, k2y, k3y, k4y, k1z, k2z, k3z, k4z;
	double x = x_0, y = y_0, z = z_0;
	float h1 = -abs(h);						            //decreasing x
	for(x = x_0/h; x >= xmin/h; x -= 1){				//changing x to x/h to avoid floating increment
		
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



void euler(dcomp *arr, dcomp *x, dcomp *b1, double h, int n, FILE *file){
    dcomp x1[n];
    double t = 0;
    for(int i = 0; i < 1e6; i++){
        fprintf(file,"%lf	",t);
        for(int j = 0; j < n; j++)
            fprintf(file,"%lf	",real(*(x + j)));	
        fprintf(file,"	\n");
        mm((dcomp*) arr, (dcomp*) x, (dcomp*) x1, n, n, 1);
        ma((dcomp*) x1, (dcomp*) b1, (dcomp*) x1, n, 1);
        for(int j = 0; j < n; j++)
            x1[j] *= h;
        ma((dcomp*) x1, (dcomp*) x, (dcomp*) x1, n, 1);
        assign((dcomp*) x, (dcomp*) x1, n, 1);
        t += h;
    }
}




int main(){
	int n = 8;
    dcomp i = {0, 1}, o = {0, 0}, I = {1, 0};
    dcomp x[8];
    
    dcomp Omega_p = {1,0}, Omega_c = {0,0}, Gamma_21 = {6,0}, Gamma_32 = {0.1, 0}, Gamma_31 = {0.1,0},  Delta_c = {0,0}, Delta_p = {0, 0};
    dcomp Gamma_2 = Gamma_31 + Gamma_32;
    dcomp Gamma_3 = Gamma_21 + Gamma_31 + Gamma_32;
    
    
	//dcomp arr[2][2] = {{{1,-2},{7,0}},{{-1,3},{-4,1}}};
    
    //double mat[200][2];
    //print((dcomp*)arr1,5,5);
    FILE *file;

	file = fopen("data_Level_3_Real.txt","w");
    

    
    
    double h = 0.1;
    for(int k = 0; k < 200/h; k++){
        dcomp Delta_p = {double(k*h)-100,0};
        dcomp delta = Delta_p + Delta_c;
        dcomp arr[n][n] = {{-Gamma_31,-i*Omega_p/2.0,o, i*Omega_p/2.0, Gamma_21 - Gamma_31, o, o, o},{-i*Omega_p/2.0, -i*Delta_p - Gamma_21/2.0, -i*Omega_c/2.0, o, i*Omega_p/2.0,o,o,o},{o, -i*Omega_c/2.0, i*delta - Gamma_2/2.0, o, o, i*Omega_p/2.0,o, o}, { i*Omega_p/2.0, o, o, i*Delta_p - Gamma_21/2.0, -i*Omega_p/2.0, o, i*Omega_c/2.0, o}, {-Gamma_32, i*Omega_p/2.0, o, -i*Omega_p/2.0, -Gamma_32 - Gamma_21, -i*Omega_c/2.0 , o, i*Omega_c/2.0}, {-i*Omega_c/2.0, o, i*Omega_p/2.0, o, -i*Omega_c, -i*delta + i*Delta_p - Gamma_3/2.0, o, o},{ o, o, o, i*Omega_c/2.0, o, o, -i*delta - Gamma_2/2.0, -i*Omega_p/2.0}, {i*Omega_c/2.0, o, o, o, i*Omega_c, o, -i*Omega_p/2.0, i*Delta_p}};

        dcomp b[n] = {-Gamma_31,o,o,o, -Gamma_32, -i*Omega_c/2.0,o,i*Omega_c/2.0};
        //dcomp b1[8] = {Gamma_31,o,o,o, Gamma_32, i*Omega_c/2.0,o,-i*Omega_c/2.0};
		gj((dcomp*)arr, (dcomp*) b, n);
            // LU_decomposition((dcomp*)arr, (dcomp*) b, 8);
            // forward_backward((dcomp*) arr, (dcomp*) b, (dcomp*) x, 8);
        //mat[k][0] = double(k-100);
        //mat[k][1] = real(b[2]);
		fprintf(file,"%lf	",double(k*h-100));
        for(int j = 0; j < n; j++)
            fprintf(file,"%lf	",real(b[j]));
			fprintf(file,"%lf	",real(I-b[0]-b[4]));
        fprintf(file,"	\n");
        // fprintf(file,"%lf	%lf\n",double(k*h-100),imag(b[6]));
        

    }
    
    fclose(file);
    return 0;
}
