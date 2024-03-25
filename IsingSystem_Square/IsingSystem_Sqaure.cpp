#include <iostream>
#include "IsingSystem_Square.hpp"

using namespace std;

int main()
{
    double start = 0.05; 
    double end = 4.0; 
    double step = 0.05; 

    vector<double> Temperature;
   
     for (double i = start; i < end; i += step) 
    {
         Temperature.push_back(i); 
    }
     vector<double> beta(Temperature.size());
    
    for(int i = 0; i<Temperature.size(); i++)
        {
            beta[i] = 1.0/Temperature[i];
        }

    const vector<int> system_size = { 5,5 };

    IsingSystem_Square model(system_size,beta);
    
    model.exact();
    model.print_exact();

    return 0;

};