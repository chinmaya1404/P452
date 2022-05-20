#include "utility.cpp"
#include<cmath>

double f(double x){
    return 1/sqrt(1+x*x);
}

double GQ_Integrator(double (*f)(double), double x1, double x2, double * roots, double * weights, int n, double e){
	double sum = 0;
	double m = (x2-x1)/2;
	double c = (x2+x1)/2;
	for(int i = 0; i < n; i++)
		sum += m*weights[i]*f(m*roots[i]+c);
	return sum;
}

double Integrator(double (*f)(double),double x1,double x2,double * roots, double * weights, int n, double e){
	double c = (x1+x2)/2;
	double I0 = GQ_Integrator(f,x1,x2,roots,weights, n, e);
	double I1 = GQ_Integrator(f,x1,c,roots,weights, n, e);
	double I2 = GQ_Integrator(f,c,x2,roots,weights, n, e);
	if (abs(I0-I1-I2) < e)
		return I1+I2;
	else
		return Integrator(f,x1,c,roots,weights,n,e) + Integrator(f,c,x2,roots,weights,n,e);
}

int main(){
    int n;
    double x[n], w[n];
    double e = 1e-9;
    
    n = 4;
    cout<<"n = "<<n<<", I = ";
    double x4[n] = {0.861136311, 0.339981043,-0.339981043,-0.861136311};
    double w4[n] = {0.347854845, 0.652145154, 0.652145154, 0.347854845};
    cout<<Integrator(f,-1,1,x4,w4, n, e)<<endl;

    n = 5;
    cout<<"n = "<<n<<", I = ";
    double x5[n] = {0.906179845, 0.538469310,0.0,-0.538469310,-0.906179845};
    double w5[n] = {0.236926885, 0.478628670, 0.568888889, 0.478628670, 0.236926885};
    cout<<Integrator(f,-1,1,x5,w5, n, e)<<endl;

    n = 6;
    cout<<"n = "<<n<<", I = ";
    double x6[n] = {0.932469514, 0.661209386, 0.238619186, -0.238619186,-0.661209386,-0.932469514};
    double w6[n] = {0.171324492,0.360761573,0.467913934,0.467913934,0.360761573,0.171324492};
    cout<<Integrator(f,-1,1,x6,w6, n, e)<<endl;

    return 0;
}

/*
n = 4, I = 1.76275
n = 5, I = 1.76275
n = 6, I = 1.76275
*/