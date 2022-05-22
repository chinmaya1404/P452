#include "utility.cpp"

double f(double x){
    return exp(-x*x);
}

double sample(double x, double alpha){
    return alpha*exp(-x);
}
int main(){
    int n = 1e8;
    cout<<"Without Importance Sampling"<<endl;
    double sum = 0;
    for(int i = 0; i < n; i++){
        double x = (double) rand()/RAND_MAX;
        sum += f(x);
    }
    cout<<"I = "<<sum/n<<endl;

    cout<<"With Importance Sampling"<<endl;
    sum = 0;
    for(int i = 0; i < n; i++){
        double x = (double) rand()/RAND_MAX;
        double alpha = exp(1)/(exp(1)-1);
        double y = -log(1-x/alpha);
        sum += f(y)/sample(y,alpha);
    }
    cout<<"I = "<<sum/n<<endl;
    return 0;
}

/*
Without Importance Sampling
I = 0.74682
With Importance Sampling
I = 0.746823
*/