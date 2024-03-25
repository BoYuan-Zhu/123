#include<catch2/catch_test_macros.hpp>
#include<catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>
#include "IsingSpinOnLattice.hpp"
#include "Ising_system.hpp"
#include "IsingSystem_Square.hpp"
#include <vector>

using namespace std;
TEST_CASE("IsingSpinOnLattice","[single spin]")
{
    IsingSpinOnLattice spin;
    SECTION("spin position/sublattice (initial)")
    {
        constexpr int dim = 2;
        spin.set_dim(dim);
        REQUIRE(spin._position() == vector<int>({0,0}));
        REQUIRE(spin._NN() == std::vector<int>({-1}));
        REQUIRE(spin._NN(0) == -1);
    };
};


TEST_CASE("IsingSystem","[examples of 10 spins]")
{
    const int n_spin = 10;
    std::vector<double> beta = {0.1,1.0};
    IsingSystem spin(n_spin,beta);
    
    SECTION("basics")
    {
      REQUIRE(spin._sz(0) == 1);
      
      spin.set_up_spin(0);
      REQUIRE(spin._sz(0)== 1);
      
      spin.set_dw_spin(0);
      REQUIRE(spin._sz(0)== -1);
      
      spin.set_spin(0,1);
      REQUIRE(spin._sz(0)== 1);
     
      spin.set_spin(0,-1);
      REQUIRE(spin._sz(0)==-1);
    
      spin.flip_spin(0);
      REQUIRE(spin._sz(0)==1);

      spin.flip_spin(0);
      spin.flip_spin(0);
      REQUIRE(spin._sz(0)==1);
  
      REQUIRE(spin._beta()[0] == 0.1); 
      REQUIRE(spin._beta()[1] == 1.0); 
    };
      
    SECTION("spin state code #7")
    {
        spin.set_state_by_code(7);
        REQUIRE(spin.eval_mz() == -4);
        REQUIRE(spin.eval_energy_1D() == -6);
    };

    SECTION("spin state code #77")
    {
        spin.set_state_by_code(77);
        REQUIRE(spin.eval_mz() == -2);
        REQUIRE(spin.eval_energy_1D() == 2);
    };
    
    SECTION("spin state code #777")
    {
        spin.set_state_by_code(777);
        REQUIRE(spin.eval_mz() == -2);
        REQUIRE(spin.eval_energy_1D() == -2);
    };
    
    IsingSystem model(n_spin,beta);
    SECTION("spin position initialization")
    {
        constexpr int dim = 2;
        model.set_dim(dim);
        for (int i = 0;i<n_spin; i++)
        {
            REQUIRE(model._spin_position(i) == std::vector<int>({0,0}));
            REQUIRE(model._spin_NN(i) == std::vector<int>({-1}));
            REQUIRE(model._spin_NN(i,0) == -1);
        }
    };
};

TEST_CASE("IsingSystem_Square","[examples of 6 x 6 spins]")
{
    const std::vector<int> system_size = {6,6};
    const std::vector<double> beta = {1.0};
    IsingSystem_Square model(system_size,beta);

    SECTION("basics")
    {
        REQUIRE(model._n_spins() == 36);
        REQUIRE_THAT(model._J(),Catch::Matchers::WithinULP(-1.0,4));
        REQUIRE(model.ground_state() == -4*6*6/2);
    };
      SECTION("site index")
    {
        const std::vector<int> lattice_coordinate = {3,4};
        REQUIRE(model.site_index(lattice_coordinate) == 27);
        REQUIRE(model.lattice_coordinate(27) == lattice_coordinate);
    };
    SECTION("neighboring site coordinates")
    {
        const std::vector<int> lattice_coordinate ={3,3};
        const std::vector<int> lattice_coordinate_pos_x={4,3};
        const std::vector<int> lattice_coordinate_pos_y={3,4};
        const std::vector<int> lattice_coordinate_neg_x={2,3};
        const std::vector<int> lattice_coordinate_neg_y={3,2};
        REQUIRE(model.shift_pos_x(lattice_coordinate)== lattice_coordinate_pos_x);
        REQUIRE(model.shift_pos_y(lattice_coordinate)== lattice_coordinate_pos_y);
        REQUIRE(model.shift_neg_x(lattice_coordinate)== lattice_coordinate_neg_x);
        REQUIRE(model.shift_neg_y(lattice_coordinate)== lattice_coordinate_neg_y);


    };

    SECTION("neighboring site coordinates2")
    {
        const std::vector<int> lattice_coordinate ={0,5};
        const std::vector<int> lattice_coordinate_pos_x={1,5};
        const std::vector<int> lattice_coordinate_pos_y={0,0};
        const std::vector<int> lattice_coordinate_neg_x={5,5};
        const std::vector<int> lattice_coordinate_neg_y={0,4};
        REQUIRE(model.shift_pos_x(lattice_coordinate)== lattice_coordinate_pos_x);
        REQUIRE(model.shift_pos_y(lattice_coordinate)== lattice_coordinate_pos_y);
        REQUIRE(model.shift_neg_x(lattice_coordinate)== lattice_coordinate_neg_x);
        REQUIRE(model.shift_neg_y(lattice_coordinate)== lattice_coordinate_neg_y);
    };

    SECTION("connectivity")
    {
        constexpr int i = 21;
        REQUIRE(model.NN(i,0) == 22);
        REQUIRE(model.NN(i,1) == 27);
        REQUIRE(model.NN(i,2) == 20);
        REQUIRE(model.NN(i,3) == 15);
    };

    SECTION("pi state : magnetization and energy")
    {
        std::vector<bool>
          state({1,1,0,1,1,1,0,0,1,1,1,0,1,1,1,1,0,1,0,0,0,0,0,0,1,1,0,1,0,1,1,1,0,0,0,0});
        model.set_state(state);
        REQUIRE(model.eval_mz() == 2);
        REQUIRE_THAT(model.eval_energy(),Catch::Matchers::WithinULP(-4.0,4));
    };
};




// TEST_CASE("IsingSystem_Square","[tests for exact counting]")
// {
//     const vector<int> system_size = { 6, 6 };
//     std::vector<double> beta = { 0.1, 1.0, 2.0 };
//     IsingSystem_Square model(system_size, beta);
//     std::vector<bool>
//           pi_state({1,1,0,1,1,1,0,0,1,1,1,0,1,1,1,1,0,1,0,0,0,0,0,0,1,1,0,1,0,1,1,1,0,0,0,0});
//     constexpr double energy = -4.0;
//     constexpr int magz = 2;
//     constexpr double w0= 1.491824697641270;//exp(-1 *0.1*(-4.0))
//     constexpr double w1 = 54.59815003314424;// exp(-1*1.0 *(-4.0))
//     constexpr double w2=2980.957987041728;// exp(-1*2.0 *(-4.0))
//     constexpr double w[3] = {w0,w1,w2};

//     SECTION("'pi' state :M, E, and Boltzmann weight") 
//     {
//         model.set_state(pi_state);
//         REQUIRE(model.eval_mz() == magz);
//         REQUIRE_THAT(model.eval_energy(), Catch::Matchers::WithinULP(energy, 4));
//         for(std::size_t beta_idx=0; beta_idx<beta.size(); beta_idx++)
//         {
//             REQUIRE_THAT(model.weight_unnormalized(beta[beta_idx]), Catch::Matchers::WithinULP(w[beta_idx], 4));
//         };
//     };

//     SECTION("'pi' state :single term in the whole sum")
//     {
//         model.exactly_evaluate(pi_state);
//         for(std::size_t beta_idx=0;beta_idx<beta.size(); beta_idx++)
//         {
//             REQUIRE_THAT(model._exact_energy_Z(beta_idx), Catch::Matchers::WithinULP(w[beta_idx], 4));
//             REQUIRE_THAT(model._exact_energy_q(beta_idx), Catch::Matchers::WithinULP(energy * w[beta_idx], 4));
//             REQUIRE_THAT(model._exact_energy_q_sq(beta_idx), Catch::Matchers::WithinULP(energy * energy * w[beta_idx],4));
//             REQUIRE_THAT(model._exact_magz_Z(beta_idx), Catch::Matchers::WithinULP(w[beta_idx], 4));
//             REQUIRE_THAT(model._exact_magz_q_sq(beta_idx), Catch::Matchers::WithinULP(magz*magz*w[beta_idx],4));
//         }
//     };
// };
