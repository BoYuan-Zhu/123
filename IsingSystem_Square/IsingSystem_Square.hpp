#ifndef IsingSpin_Square_hpp
#define IsingSpin_Square_hpp

#include <iostream>
#include <vector>
#include "Ising_system.hpp"

class  IsingSystem_Square:public IsingSystem
{
private:
   std::vector<int> system_size;
   void setup_NN()
   {
    for (int site_idx = 0; site_idx < n_spins; site_idx++)
    {
        std::vector<int> r = lattice_coordinate(site_idx);
        spin[site_idx].set_NN(0,site_index(shift_pos_x(r)));
        spin[site_idx].set_NN(1,site_index(shift_pos_y(r)));
        spin[site_idx].set_NN(2,site_index(shift_neg_x(r)));
        spin[site_idx].set_NN(3,site_index(shift_neg_y(r)));
    }
   }
public:
     IsingSystem_Square(const std::vector<int> system_size_spec):
       IsingSystem(system_size_spec[0]*system_size_spec[1]),
       system_size(system_size_spec){setup_NN();};
       
    ~ IsingSystem_Square(){};
    

    int site_index(const std::vector<int> lattice_coordinate) const
    {
        if (lattice_coordinate[0]<=system_size[0]-1 && lattice_coordinate[1]<=system_size[1])
        {
            return system_size[0]*(lattice_coordinate[1])+lattice_coordinate[0];
        }
        else
        {
            std::cout << "lattice_coordinate exceeds the size of the lattice.\n";
            std::cout << "Assigned site_index to 0.";
            return 0;    
        }
    };

    std::vector<int> lattice_coordinate(int site_idx) const
    {
        if (site_idx<=system_size[0]*system_size[1])
        {
            return {site_idx%system_size[0],int(site_idx/system_size[0])};
        }
        else
        {
           std::cout << "lattice_coordinate exceeds the size of the lattice.\n";
           std::cout << "Assigned lattice_index to {0,0}.";
           return {0,0};   
        }
    };

    std::vector<int> shift_pos_x(const std::vector<int> r_spec)const
    {
        std::vector<int> r(r_spec);
        r[0] = (r[0] + 1) % system_size[0];
        return r;
    };

    std::vector<int> shift_pos_y(const std::vector<int> r_spec)const
    {
        std::vector<int> r(r_spec);
        r[1] = (r[1] + 1) % system_size[1];
        return r;
    };

    std::vector<int> shift_neg_x(const std::vector<int> r_spec)const
    {
        std::vector<int> r(r_spec);
        r[0] = (r[0] - 1 + system_size[0]) % system_size[0];
        return r;
    };

    std::vector<int> shift_neg_y(const std::vector<int> r_spec)const
    {
        std::vector<int> r(r_spec);
        r[1] = (r[1] - 1 + system_size[1]) % system_size[1];
        return r;
    };

    int NN(const int site_idx, const int bond_idx) const
    {
        return spin[site_idx]._NN(bond_idx);
    };

    double eval_energy() const
   { 
    double energy = 0;
    for (int i = 0;i<_n_spins();i++) 
   {
      for (int j = 0; j<4; j++ )
          {energy += spin[i]._sz() * spin[spin[i]._NN(j)]._sz();};
   }
   
    energy *= J/2;
    return energy;

   };

    
};

#endif