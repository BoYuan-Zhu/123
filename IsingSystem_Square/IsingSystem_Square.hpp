#ifndef IsingSpin_Square_hpp
#define IsingSpin_Square_hpp

#include <iostream>
#include <vector>
#include "Ising_system.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>

class  IsingSystem_Square:public IsingSystem
{
private:
   std::vector<int> system_size;
   std::vector<double> internal_E;
   std::vector<double> internal_E_sq;
   std::vector<double> M_sq;
   std::vector<double> Z;
   std::vector<double> C;
   const double ground_state_energy;
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
     IsingSystem_Square(const std::vector<int> system_size_spec,const std::vector<double> beta_spec):
       IsingSystem(system_size_spec[0]*system_size_spec[1],beta_spec),
       ground_state_energy(-4*system_size_spec[0]*system_size_spec[1]/2),
       system_size(system_size_spec)
       {
       setup_NN();
       internal_E.assign(beta.size(),0);
       internal_E_sq.assign(beta.size(),0);
       M_sq.assign(beta.size(),0);
       Z.assign(beta.size(),0);
       C.assign(beta.size(),0);
       };
       
       
    ~ IsingSystem_Square(){};
    

    int site_index(const std::vector<int> lattice_coordinate) const
    {
        if (lattice_coordinate[0]<=system_size[0]-1 && lattice_coordinate[1]<=system_size[1])
        {
            return system_size[0] * lattice_coordinate[1]+lattice_coordinate[0];
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
    
    double ground_state() const
    {
        return ground_state_energy;
    };

    double weight_unnormalized(const std::size_t beta_idx)const
    {
        return exp(-beta[beta_idx]*(eval_energy()-ground_state()));
    };

    double _exact_energy_Z(const std::size_t beta_idx)const
    {
        return weight_unnormalized(beta_idx);
    };

    double _exact_energy_q(const std::size_t beta_idx)const
    {
        return _exact_energy_Z(beta_idx)*eval_energy();
    }; 
    
    double _exact_energy_q_sq(const std::size_t beta_idx)const
    {
        return _exact_energy_q(beta_idx)*eval_energy();
    }; 
   
    double _exact_magz_Z(const std::size_t beta_idx)const
    {
        return weight_unnormalized(beta_idx);
    }; 

    double _exact_magz_q_sq(const std::size_t beta_idx)const
    {
        return _exact_magz_Z(beta_idx)*eval_mz()*eval_mz();
    }; 
    
    void exactly_evaluate_given()
    {
        for(int i = 0; i<beta.size(); i++)
        {
            internal_E[i] += _exact_energy_q(i);
            internal_E_sq[i] += _exact_energy_q_sq(i);
            Z[i] += weight_unnormalized(i);
            M_sq[i] += _exact_magz_q_sq(i);
        }
        
    };

    // For state in vector form
    void exactly_evaluate(const std::vector<bool>& state)
    {
         set_state(state),
         exactly_evaluate_given();
    };
    
    // For state in integer form
    void exactly_evaluate(const long long& rep_state)
    {
        std::vector<bool> state = state_by_code(rep_state);
        exactly_evaluate(state);
    };
    
    //going through all the state
    void exact()
    {
        long long rep_state =0;
        while (rep_state <= maxrep_state)
        {
            exactly_evaluate(rep_state++);
        }
        normalize_direct();
    };

    void normalize_direct()
    {
      for(int i = 0; i<beta.size(); i++)
        {
            internal_E[i] *= 1/Z[i];
            internal_E_sq[i] *= 1/Z[i];
            M_sq[i] *= 1/Z[i];
            M_sq[i] *= 1.0/(n_spins*n_spins);

        }
        for(int i = 0; i<beta.size(); i++)
        {
            C[i] = beta[i] * beta[i] * (internal_E_sq[i]-internal_E[i]*internal_E[i]);
            C[i] *= 1.0/n_spins;
        }
        
    };
    void print_exact() const
    {
        std::cout << "Specific Heat: ";
    for (double value : C) 
    {
        std::cout << value << " ";
    }
    std::cout << "." << std::endl;
        
        std::cout << "Magnetization (Squared): ";
    for (double value : M_sq) 
    {
        std::cout << value << " ";
    }
    std::cout << "." << std::endl;
    };



};

#endif