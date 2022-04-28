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
			exit(0);
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

//EULER METHOD

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
    int n = 15;
    dcomp i = {0, 1}, o = {0, 0}, I = {1, 0};

    
    dcomp Omega_p = {1,0}, Omega_c = {10,0}, Omega_s = {5,0}, Gamma_21 = {6.1,0}, Gamma_32 = {0.68, 0}, Gamma_43 = {0.1,0}, Delta_s = {0,0},  Delta_c = {0,0}, Delta_p = {0, 0};
    
	//dcomp arr[2][2] = {{{1,-2},{7,0}},{{-1,3},{-4,1}}};
    
    //double mat[200][2];
    //print((dcomp*)arr1,5,5);
    FILE *file;

    file = fopen("data_Level_4_Real_1_10_5.txt","w");
    
    


    double h = 0.1;
    for(int k = 0; k < 200/h; k++){
        Delta_p = {double(k*h)-100,0};
        dcomp delta = Delta_p + Delta_c;
        dcomp arr[n][n] = {{o, -i*Omega_p/2.0,o,o,i*Omega_p/2.0, Gamma_21,o,o,o,o,o,o,o,o,o},{-i*Omega_p/2.0, -i*Delta_p-Gamma_21/2.0, -i*Omega_c/2.0,o,o,i*Omega_p/2.0,o,o,o,o,o,o,o,o,o},{o, -i*Omega_c/2.0, -i*(Delta_p+Delta_c) -Gamma_32/2.0, -i*Omega_s/2.0,o,o,i*Omega_p/2.0,o,o,o,o,o,o,o,o}, {o,o,-i*Omega_s/2.0, -i*(Delta_p + Delta_c + Delta_s) -Gamma_43/2.0,o,o,o,i*Omega_p/2.0,o,o,o,o,o,o,o}, {i*Omega_p/2.0, o,o,o,i*Delta_p-Gamma_21/2.0, -i*Omega_p/2.0,o,o,i*Omega_c/2.0,o,o,o,o,o,o}, {o,i*Omega_p/2.0, o, o, -i*Omega_p/2.0, -Gamma_21, -i*Omega_c/2.0, o,o, i*Omega_c/2.0, Gamma_32,o,o,o,o}, {o,o, -i*Omega_s/2.0, o, o, -i*Omega_c/2.0, -i*Delta_c-(Gamma_32+Gamma_21)/2.0, -i*Omega_s/2.0, o,o,i*Omega_c/2.0, o,o,o,o}, {o,o,o,i*Omega_p/2.0, o,o, -i*Omega_s/2.0, -i*(Delta_c+Delta_s)/2.0 -(Gamma_43+Gamma_21)/2.0, o,o,o,i*Omega_c/2.0, o,o,o}, {o,o,o,o,i*Omega_c/2.0, o,o,o, i*(Delta_p+Delta_c)/2.0 -Gamma_32/2.0, -i*Omega_p/2.0, o,o, i*Omega_s/2.0,o,o}, {o,o,o,o,o,i*Omega_c/2.0, o,o, -i*Omega_p/2.0, i*Delta_c-(Gamma_32+Gamma_21)/2.0, -i*Omega_c/2.0,o,o,i*Omega_s/2.0,o}, {-Gamma_43,o,o,o,o,-Gamma_43,i*Omega_c/2.0, o,o,-i*Omega_c/2.0, -Gamma_32-Gamma_43, -i*Omega_s/2.0, o,o,i*Omega_s/2.0}, {-i*Omega_s/2.0, o,o,o, o, -i*Omega_s/2.0, o, i*Omega_c/2.0, o,o,-i*Omega_s, -i*Delta_s - (Gamma_43+Gamma_32)/2.0, o,o,o}, {o,o,o,o,o,o,o,o,i*Omega_s/2.0, o,o,o,i*(Delta_p+Delta_c+Delta_s)-Gamma_43/2.0, -i*Omega_p/2.0, o}, {o,o,o,o,o,o,o,o,o,i*Omega_s/2.0, o,o,-i*Omega_p/2.0, i*(Omega_c+Omega_s)-(Gamma_43+Gamma_21)/2.0, -i*Omega_c/2.0}, {i*Omega_s/2.0,o,o,o,o,i*Omega_s/2.0,o,o,o,o,i*Omega_s,o,o,-i*Omega_c/2.0, i*Omega_s-(Gamma_43+Gamma_32)/2.0}};

        dcomp b[n] = {o,o,o,o,o,o,o,o,o,o,-Gamma_43,-i*Omega_s/2.0,o,o,i*Omega_s/2.0};
        //dcomp b1[15] = {o,o,o,o,o,o,o,o,o,o,Gamma_43,i*Omega_s/2.0,o,o,-i*Omega_s/2.0};
            gj((dcomp*)arr, (dcomp*) b, n);
        //mat[k][0] = double(k-100);
        //mat[k][1] = real(b[2]);

        fprintf(file,"%lf	",double(k*h-100));
        for(int j = 0; j < n; j++)
            fprintf(file,"%lf	",real(b[j]));
			fprintf(file,"%lf	",real(I-b[0]-b[5]-b[10]));
        fprintf(file,"	\n");
        
        // fprintf(file,"%lf	%lf\n",double(k*h-100),imag(b[4]));
        

    }
    fclose(file);
    

    //print_r((double*)mat,200,2);
    //print((dcomp*)arr,5,5);
    //mm((dcomp*) arr, (dcomp*) b, (dcomp*) x, 5,5,1);
    //print((dcomp*)x,5,1);
    //dcomp inv[5][5];


    
    return 0;
}
/*
    dcomp Omega_p = {10,0}, Omega_c = {5,0}, Gamma_21 = {6,0}, Gamma_32 = {0.001, 0}, Gamma_31 = {0.01,0}, delta = {1,0}, Delta_p = {0,0};
    dcomp Gamma_2 = Gamma_31 + Gamma_32;
    dcomp Gamma_3 = Gamma_21 + Gamma_31 + Gamma_32;
    dcomp i = {0, 1}, o = {0, 0};
    dcomp arr[9][9] = {{o,-i*Omega_p/2.0,o, i*Omega_p/2.0, Gamma_21, o, o, o, Gamma_31},{-i*Omega_p/2.0, -i*Delta_p - Gamma_21/2.0, -i*Omega_c/2.0,o, i*Omega_p/2.0,o,o,o,o},{o, -i*Omega_c/2.0, -i*delta - Gamma_2/2.0, o, o, i*Omega_p/2.0,o, o, o}, { o, -i*Omega_p/2.0, o, -i*Omega_p/2.0, -Gamma_21, i*Omega_c/2.0, o, i*Omega_c/2.0, Gamma_32}, {o, o, i*Omega_p/2.0, o, -i*Omega_c/2.0, i*Delta_p-Gamma_3/2.0+i*delta, o, o, i*Omega_c/2.0}, {o, o, o, o, o, i*Omega_c/2.0, o, -i*Omega_c/2.0, -Gamma_2}, {i*Omega_p/2.0, o, o, i*Delta_p-Gamma_21/2.0, -i*Delta_p/2.0, o, i*Omega_c/2.0, o, o}, { o, o, o, i*Omega_c/2.0, o, o, i*delta-Gamma_2/2.0,-i*Omega_p/2.0, o},{ o, o, o, o, i*Omega_c/2.0, o, -i*Omega_p/2.0, -i*Delta_p-Gamma_3/2.0-i*delta, -i*Omega_c/2.0}};
	//dcomp arr[2][2] = {{{1,-2},{7,0}},{{-1,3},{-4,1}}};
    dcomp arr1[9][9];
    for(int i = 0; i<9;i++){
        for(int j =0; j<9;j++)
            arr1[i][j] = arr[i][j];
    }
    
    
    print((dcomp*)arr1,9,9);
    dcomp b[9] = {o,o,o,o,o,o,o,o,o};
    dcomp x[9];
    for(int i=0;i<9;i++){
        pp((dcomp*)arr, (dcomp*) b, i, 9);
        gj((dcomp*)arr, (dcomp*) b, i, 9);
        print((dcomp*)arr,9,9);
    }
    print((dcomp*)x,9,1);
    //print((dcomp*)arr,5,5);
    //mm((dcomp*) arr1, (dcomp*) b, (dcomp*) x, 5,5,1);
    //print((dcomp*)x,5,1);
    //dcomp inv[5][5];
    
    return 0;
}
*/