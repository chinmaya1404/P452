#include "utility.cpp"

double phi1(double x){
    return 1;
}

double phi2(double x){
    return 0;
}

double mnorm(double * arr, int n){
    double sum = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++)
            sum += pow(*(arr + i*n + j),2);
    }
    return sum;
}

int main(){
    int n = 40;
    double e = 1e-4;
    double dx = 1/double(n),dt = 1/double(n);
    double arr[n+1][n+1];
    
    for(int i = 0; i <= n; i++)
        arr[i][0] = phi1(i*dx);

    for(int i = 0; i <= n; i++)
        arr[n][i] = phi2(i*dx);

    for(int j = 1; j < n; j++)
        arr[0][j] = phi2(j*dt);

    for(int j = 1; j < n; j++)
        arr[n][j] = phi2(j*dt);
    
    for(int i = 1; i < n; i++){
        for(int j = 1; j < n; j++)
            arr[i][j] = 0;
    }

    
    double a2 = 0, a1 = mnorm((double*) arr, n+1);
    
    while(abs(a1 - a2) > e){
        a2 = a1;
        for(int i = 1; i < n; i++){
            for(int j = 1; j < n; j++)
                arr[i][j] = (arr[i+1][j] + arr[i-1][j] + arr[i][j+1] + arr[i][j-1])/4;
        }
        a1 = mnorm((double*) arr, n+1);
        
    }

    FILE* file;
    file = fopen("Laplace.txt", "w");
    for(int i = 0; i <= n; i++){
        for(int j = 0; j < n; j++)
        fprintf(file,"%lf %lf, %lf\n",i*dx, j*dt,arr[i][j]);
    }
    
    fclose(file);


    return 0;
}