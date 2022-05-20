#include "utility.cpp"

//RANDOM NUMBER GENERATOR

double randgen(double x_0, int a, int m){
    return double(int(a*x_0)%m);
}



int main(){
    double z = 2.2, theta, a = 572, m = 16381, min = 0, max = 2*M_PI, sum = 0, x = 0, y = 0;
    int n = 500, N = 200;

    x = 2.2;
    
    FILE* file;

    file = fopen("random_walk.txt","w");				//file "random_walk.txt" to store the step number and final position related data
	fprintf(file,"%s	%s		%s		%s\n", "N","sqrt(N)" ,"R", "R_rms");
		float R_tot = 0,x_tot = 0, y_tot = 0, x_2_tot = 0, y_2_tot = 0;
		for(int j = 0; j < 500; j++){				//loop for doning the random walk for constant N for 100 times
            z = double(j/10+1);
			double x = 0, y = 0;
			//random_walk(&x, &y,N,1,file);			//call the function
            for(int i = 0; i < N; i++){
            z = randgen(z,a,m);
            theta = min + (max - min)*z/(m-1);
            x += cos(theta);					        //find x coordinate
            y += sin(theta); 					        //find y coordinate
            
            }
			R_tot += sqrt(x*x +y*y);
			x_tot += x;
			y_tot += y;
			x_2_tot += x*x;
			y_2_tot += y*y; 
		}
		float R_mean = R_tot/500;				//calculate R_mean
		float R_rms = sqrt((x_2_tot/500) + (y_2_tot/500));	//calculate R_rms
		fprintf(file,"%d	%lf	%lf	%lf\n", N, sqrt(N), R_mean, R_rms);		//store data in "random_walk.txt"
        cout<<"N = "<<N<<" , sqrt(N) = "<<sqrt(N)<<", R_mean = "<<R_mean<<", R_rms = "<<R_rms<<endl;
		
	
	fclose(file);

    return 0;
}

/*
N = 200 , sqrt(N) = 14.1421, R_mean = 12.1976, R_rms = 13.4652
*/