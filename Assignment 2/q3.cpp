#include "utility.cpp"

double randgen(double x_0, int a, int m){
    return double(int(a*x_0)%m);
}

double c1(double x, double y, double z){
    return x*x + y*y;
}

double c2(double x, double y, double z){
    return z*z + y*y;
}

int main(){
    double x1, y1, z1, x = 2.2, y = 4.3, z = 3.1, a = 572, m = 16381, sum = 0;
    int n = 1e8;
    
    int count = 0;
    for(int i = 0; i < n; i++){
        x = randgen(x,a,m);
        y = randgen(y,a,m);
        z = randgen(z,a,m);
        
        x1 = -1 + 2*x/(m-1);
        y1 = -1 + 2*y/(m-1);
        z1 = -1 + 2*z/(m-1);
        //cout<<x1<<", "<<y1<<", "<<z1<<endl;
        if((c1(x1,y1,z1) <= 1)&&(c2(x1,y1,z1) <= 1))
            count++;
    }
    double v = pow(2,3);        //volume of cube of side 2
    cout<<"Volume of Steinmetz solid = "<<double(count)/n*v<<endl;
    return 0;
}

// Volume of Steinmetz solid = 5.36264