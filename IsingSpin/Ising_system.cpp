#include <iostream>
#include "Ising_system.hpp"

int main()
{
    using namespace std;
    int n_spins;
    int rep_state;
    cout << "Enter the total spin sites: ";
    cin >> n_spins;
    cout << "Enter the spin state:";
    cin >> rep_state;
    IsingSystem Ising(n_spins);
    Ising.set_state_by_code(rep_state);
    cout<< "The total magnetization: "<< Ising.eval_mz() << "\n";
    cout<< "The total energy(1D): "<< Ising.eval_energy_1D() << "\n";
    return 0;
};