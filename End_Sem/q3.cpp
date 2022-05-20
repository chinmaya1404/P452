#include "utility.cpp"
#include<cmath>
double g(double x){
    return 20*abs(sin(M_PI*x));
}
double b(double x){
    return 0;
}
double a(double x){
    return 0;
}

int main(){
    
	FILE* file;
    file = fopen("wt.txt","w");
	Forward_Diffusion(g, a, b, 1, 0.1, 0.0008, 2, 4, file);
    fclose(file);
    return 0;
}

/*
The temperature of the rod is constant at the boundary but the temperature is changing at other points. At the middle point x = 1, initially the temperature is 0, but with time it's temperature is increasing. There are two peaks at points x = 0.5 and x = 1.5 for nt < 100 and after that the middle point is the only peak.

Over time the temperature of the peaks are decreasing and reaches a saturation state.
*/