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

    // const vector<double> beta = {0.1, 1 , 4};

    const vector<int> system_size = { 2,2 };

    IsingSystem_Square model(system_size,beta);
    
    model.exact();
    model.print_exact();
    // model.set_state_by_code(1);
    // model.weight_unnormalized(1);
    // model.weight_unnormalized(1,0);
    // for (int i=0;i<16;i++){
    // cout << "for configuration site_idx" << i <<"energy weights: "<< model._exact_energy_q(1,i) << "\n";
    // }
    return 0;

};