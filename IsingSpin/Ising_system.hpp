#ifndef Ising_system_hpp
#define Ising_system_hpp

#include <iostream>
#include <vector>
#include"spin.hpp"
#include <cmath>

class IsingSystem{
private:
    const double J;
    const int n_spins;
    const long long maxrep_state;
    std::vector<IsingSpin> spin;

public:
   IsingSystem(const int n_spins_spec):J(-1.0),n_spins(n_spins_spec),
   maxrep_state(static_cast<long long>(std::pow(2,n_spins))-1){
   spin.resize(n_spins);};
   virtual ~IsingSystem(){};

   double _J() const {return J;};
   int _n_spins() const {return n_spins;};
   long long _maxrep_state() const { return maxrep_state;};
  
   int _sz(const int site_idx) const {return spin[site_idx]._sz();};
   void set_up_spin(const int site_idx) {spin[site_idx].set_up();};
   void set_dw_spin(const int site_idx) {spin[site_idx].set_down();};
   void set_spin(const int site_idx, int s_spec)
       {  spin[site_idx].set_sz(s_spec);};
   void flip_spin(const int site_idx) { spin[site_idx].flip();};

   void set_state_by_code(long long rep_state)
   {
    for (int i = 1;i<=_n_spins();i++) 
    {
        if (rep_state>=std::pow(2,_n_spins()-i))
        {
            set_up_spin(i-1);
            rep_state -= std::pow(2,_n_spins()-i);
        }
        else
        {
            set_dw_spin(i-1);
        }
    } 
   }
   
   double eval_mz() const
   {
    double M = 0;
    for (int i = 1;i<=_n_spins();i++)
    {
        M += spin[i-1]._sz(); 
    }
    return M;

   }

   double eval_energy_1D() const
   { 
    double energy = 0;
    for (int i = 1;i<_n_spins();i++) 
   {
    energy += spin[i-1]._sz() * spin[i]._sz();
   }
    energy += spin[_n_spins() - 1]._sz() * spin[0]._sz();
    energy *= J;
    return energy;

   }




};

#endif