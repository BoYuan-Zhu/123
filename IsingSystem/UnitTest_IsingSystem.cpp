#include<catch2/catch_test_macros.hpp>
#include<catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>
#include "Ising_system.hpp"

TEST_CASE("IsingSpin","[1 spin site in 1D]")
{
    int n_spins = 1;
    IsingSystem spin(n_spins);
    SECTION("spin state(initial)"){
      REQUIRE(spin._sz(0) == 1);
    }
    SECTION("set spin state as up(1)"){
        spin.set_up_spin(0);
        REQUIRE(spin._sz(0)== 1);
    }
    SECTION("set spin state as down(1)"){
        spin.set_dw_spin(0);
        REQUIRE(spin._sz(0)== -1);
    }
    SECTION("set spin state as up(2)"){
        spin.set_spin(0,1);
        REQUIRE(spin._sz(0)== 1);
    }
    SECTION("set spin state as down(2)"){
        spin.set_spin(0,-1);
        REQUIRE(spin._sz(0)==-1);
    }
    SECTION("spin flip once"){
        spin.flip_spin(0);
        REQUIRE(spin._sz(0)==-1);
    }
    SECTION("spin flip twice"){
        spin.flip_spin(0);
        spin.flip_spin(0);
        REQUIRE(spin._sz(0)==1);
    }
};

TEST_CASE("IsingSpin","[5 spin sites in 1D]")
{
    int n_spins = 5;
    IsingSystem spins(n_spins);
    SECTION("Initial spin states = 0 ")
    {
        spins.set_state_by_code(0);
        REQUIRE(spins.eval_mz() == -5);
        REQUIRE(spins.eval_energy_1D() == -5);
    }
};

TEST_CASE("IsingSpin","[10 spin sites in 1D]")
{
    int n_spins = 10;
    IsingSystem spins(n_spins);
    SECTION("Initial spin states = #7 ")
    {
        spins.set_state_by_code(7);
        REQUIRE(spins.eval_mz() == -4);
        REQUIRE(spins.eval_energy_1D() == -6);
    }
    SECTION("Initial spin states = #77 ")
    {
        spins.set_state_by_code(77);
        REQUIRE(spins.eval_mz() == -2);
        REQUIRE(spins.eval_energy_1D() == 2);
    }
    SECTION("Initial spin states = #777 ")
    {
        spins.set_state_by_code(777);
        REQUIRE(spins.eval_mz() == -2);
        REQUIRE(spins.eval_energy_1D() == -2);
    }
};