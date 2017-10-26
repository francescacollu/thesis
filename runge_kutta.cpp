//RUNGE KUTTA FUNCTION FOR EDO CALCULUS

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

double F(double t, double y)
{
    return y-pow(t,2)+1;
}

int main()
{
    double h, t, w;
    h = 0.5;
    t = 0.;
    w = 0.5;
    
    for(double i=0.; i<2; i = i+h)
    {
        double k1, k2, k3, k4;
        k1 = h * F(t, w);
        k2 = h * F(t+0.5*h, w+0.5*k1);
        k3 = h * F(t+0.5*h, w+0.5*k2);
        k4 = h * F(t+h, w+k3);

        w = w + 1./6. * (k1 + 2*k2 + 2*k3 + k4);
        
        t = t + h;
        
        cout << t << "\t" << "w" << i << "\t" << w << endl;
    }
    return 0;
}
