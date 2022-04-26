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
    int n = 3;
    dcomp i = {0, 1}, o = {0, 0}, I = {1, 0};
    dcomp Omega = {1,0}, Gamma = {6,0};

    
    FILE *file;

	file = fopen("data_Level_2_Imag.txt.txt","w");
    

    double h = 0.1;
    for(int k = 0; k < 200/h; k++){
        dcomp Delta = {double(k*h)-100,0};
        dcomp arr[n][n] = {{-Gamma, i*Omega/2.0, -i*Omega/2.0}, {i*Omega, i*Delta - Gamma/2.0, o}, {-i*Omega, o, -i*Delta - Gamma/2.0}};
        dcomp b[n] = {-Gamma, i*Omega/2.0, -i*Omega/2.0};
        gj((dcomp*)arr, (dcomp*) b, n);
        fprintf(file,"%lf	",double(k*h-100));
        for(int j = 0; j < n; j++)
            fprintf(file,"%lf	",imag(b[j]));
			fprintf(file,"%lf	",imag(I-b[0]));
        fprintf(file,"	\n");
    }
    fclose(file);
    
    
    return 0;
}
