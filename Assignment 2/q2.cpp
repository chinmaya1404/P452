#include "utility.cpp"

double randgen(double x_0, int a, int m){
    return double(int(a*x_0)%m);
}

double f(double x){
    return sqrt(1-x*x);
}

double g(double x, double y){
    return x*x + y*y;
}

int main(){
    double x1, y1, x = 2.2, y = 4.3, a1 = 65, m1 = 1021, a2 = 572, m2 = 16381, sum = 0;
    int n = 1e8;
    
    cout<<"Throwing Points"<<endl;
    for(int i = 0; i < n; i++){
        x = randgen(x,a1,m1);
        // cout<<x/(m1-1)<<endl;
        x1 = x/(m1-1);
        sum += f(x1);
    }
    cout<<"Set1: ℼ/4 = "<<sum/n<<endl;

    x = 2.2;
    sum = 0;
    for(int i = 0; i < n; i++){
        x = randgen(x,a2,m2);
        // cout<<x/(m1-1)<<endl;
        x1 = x/(m2-1);
        sum += f(x1);
    }
    cout<<"Set2: ℼ/4 = "<<sum/n<<endl;
    cout<<"solving the integral by Monte Carlo"<<endl;
    x = 2.2;
    y = 4.3;
    int count = 0;
    for(int i = 0; i < n; i++){
        x = randgen(x,a1,m1);
        y = randgen(y,a1,m1);
        // cout<<x/(m1-1)<<endl;
        x1 = x/(m1-1);
        y1 = y/(m1-1);
        if(g(x1,y1) <= 1)
            count++;
    }
    cout<<"Set1: ℼ/4 = "<<double(count)/n<<endl;

    x = 2.2;
    y = 4.3;
    count = 0;
    for(int i = 0; i < n; i++){
        x = randgen(x,a2,m2);
        y = randgen(y,a2,m2);
        // cout<<x/(m1-1)<<endl;
        x1 = x/(m2-1);
        y1 = y/(m2-1);
        if(g(x1,y1) <= 1)
            count++;
    }
    cout<<"Set2: ℼ/4 = "<<double(count)/n<<endl;
    return 0;
}

/*
Throwing Points
Set1: ℼ/4 = 0.784899
Set2: ℼ/4 = 0.785367
solving the integral by Monte Carlo
Set1: ℼ/4 = 0.788235
Set2: ℼ/4 = 0.785531
*/