#include "utility.cpp"


double f(double x, double y, double z){
	return z;
}



double es1(double x, double y, double z){
    return -4*M_PI*M_PI*y;
}

double gs(double x, double y, double z){
    return -M_PI*M_PI*y;
}

void Normalization(double * y, int n){
	double tot = 0;
	tot = abs(*y);
	for(int i = 1; i < n; i++)
		tot += 2* abs(*(y + i));
	tot += abs(*(y + n));
	tot = tot/(2*double(n));
	//cout<<tot<<endl;
	for(int i = 0; i <= n; i++)
		*(y + i) = *(y + i)/tot;
}

int main()
{
	
	double h = 1e-3;
	int n = 1/h;
	double lb = 1, ub = 1.505, slope;
	double arr[n+1][2], y[n+1], x[n+1];
	FILE* file;
	double tot = 0;

	
	string str1 = "q2.txt";
	shooting(f, gs, 0, 0, 1,0,lb, ub, h,str1);		
	import((double*) arr, str1, n+1,2);
	for(int i = 0; i <= n; i++){
		x[i] = arr[i][0];
		y[i] = arr[i][1];
	}
	Normalization((double*) y, n);
	tot = y[0];
	for(int i = 1; i < n; i++)
		tot += 2*y[i];
	tot += y[n];
	tot = tot*h/2;
	cout<<tot<<endl;
		file = fopen("gs.txt", "w");
		for(int i = 0; i <= n; i++)
			fprintf(file,"%lf	%lf\n",x[i],y[i]);
		fclose(file);
	
	
	string str2 = "q2.txt";
	shooting(f, es1, 0, 0, 1,0,lb, ub, h,str2);		
	import((double*) arr, str2, n+1,2);
	for(int i = 0; i <= n; i++){
		x[i] = arr[i][0];
		y[i] = arr[i][1];
	}
	Normalization((double*) y, n);
	
	file = fopen("es1.txt", "w");
	for(int i = 0; i <= n; i++)
		fprintf(file,"%lf	%lf\n",x[i],y[i]);
	fclose(file);
	
return 0;

}

