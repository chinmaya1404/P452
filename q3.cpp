#include "utility.cpp"

const int n = 20;
const int N = n*n;
const double mu = 0.04;

double OA(int i, double *arr){
    int x = i%n;
	int y = i/n;
	double sum = 0;
	
	if (x == 0) sum += 0.5 * *(arr+y*n+n-1);          //For x component
	else sum += 0.5 * *(arr+y*n+x-1);

	if (x == n-1) sum += 0.5* *(arr+y*n);
	else sum += 0.5 * *(arr+y*n+x+1);

	if ( y== 0) sum += 0.5 * *(arr+n*n-n+x);          //For y component
	else sum += 0.5 * *(arr+y*n-n+x);

	if(y == n-1) sum += 0.5 * *(arr+x);
	else sum += 0.5 * *(arr+y*n+n+x); 

	return sum + (mu-2) * *(arr+i);
}

void CG_method_fun( double * b, double * x, double e, FILE* file){

    double mul[N], p[N], r[N], r1, alpha, beta;

    for(int i = 0; i < N; i++)
        mul[i] = OA(i, x);
    ms((double*) b, (double*) mul,(double*) r, 1, N);               //r_0 = Ax_0 - b
    if(norm((double*) r, N) < e)                                    //Check convergence
        return;
    assign((double*) p, (double*) r, 1, N);                         //p_0 = r_0
    for(int k = 1; k <= N; k++){
        r1 = Inner_Product((double*) r,(double*)r, N);
        for(int j = 0; j < N; j++)
            mul[j] = OA(j, p);
        //print((double*) mul, N, 1);
        alpha = r1/Inner_Product((double*) p, (double*) mul, N);    //Find alpha
        for(int i = 0; i < N; i++){
            *(x + i) = *(x + i) + alpha * p[i];                     //x(k+1) = x(k) + alpha * p(k)
            r[i] = r[i] - alpha * mul[i];                           //r(k+1) = r(k) - alpha * A * p_k
            if(norm((double*) r, N) < e)                            //Check convergence
                return;   
        }

        beta = Inner_Product((double*) r, (double*) r, N)/r1;       //Find beta
        for(int i = 0; i < N; i++)
            p[i] = r[i] + beta * p[i];                              //p(k+1) = r(k+1) + beta * p(k)
            //cout<<norm((double*) r, N)<<endl;
        fprintf(file,"%d   %lf\n", k,norm((double*) r, N));
    }
}

int main(){
    FILE* file;
    double b[N], x1[N], x2[N];

    for(int i = 0; i < N; i++)  x1[i] = 1;
    for(int i = 0; i < N; i++){
        if(i == 0)  b[i] = 1;
        else    b[i] = 0;
    }
    

    file = fopen("q3_1.txt","w");
    CG_method_fun((double*) b, (double*) x1, 1e-6, file);
    fclose(file);

    for(int i = 0; i < N; i++)  x2[i] = 1;
    for(int i = 0; i < N; i++){
        if(i == 1)  b[i] = 1;
        else    b[i] = 0;
    }
    
    file = fopen("q3_2.txt","w");
    CG_method_fun((double*) b, (double*) x2, 1e-6, file);
    fclose(file);



    
    return 0;
}