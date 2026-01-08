#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm> // for std::shuffle
#include <random> // for std::mt19937 and std::random_device
#include <mpi.h>
#include <HYPRE.h>
#include <HYPRE_struct_ls.h>
#include "ex.h"
#include "setup.h"
#include "utilities.h"
#include "geometry.h"   //defines mesh
#include "material.h"   //defines material properties
#include "physics.h"
#include "radsolver.h"



//void solve_radiation_groups(const Mesh& mesh, State& state);
//int setupsolver() ;

// RadSolve.hpp
class RadSolve : public IRadSolver {
public:
    RadSolve(int nx, int ny, int nz, Pars& pars);

    ~RadSolve();
    //const std::vector<std::vector<double>>& getTemperatureField() const;
    //void setTemperature(int i, int j, double value);
    //void readMesh(const std::string& filename);
    //void setupGrid();
    double gradenergy(const Mesh& mesh, State& state, Pars &pars);
    double larsendelimiter(const Mesh &mesh, State &state, Pars &pars, double opact, int i, int j, int k, int ifreq,int ord=2);
    double divergence(const Mesh &mesh, State &state, Pars &pars, int i, int j, int k, int ifreq);

    //interface methods
    void solveRadiationTransport(const Mesh&, State&, Pars&, double) override;
    void UpdateBEmission(const Mesh&, State&, Pars&) override;
    void UpdateRadFlux(const Mesh&, State&, Pars&) override;

    void apply_milne_boundary_conditions(Mesh& mesh, State& state, Pars &pars);
    void apply_reflect_boundary_conditions(Mesh& mesh, State& state, Pars &pars);
    int updatestate(Pars &pars, Mesh &mesh, State &state);

    std::random_device rd;


private:
    int nx, ny, nz; // Grid dimensions 

    // HYPRE grid, stencil, matrix, vectors
    HYPRE_StructGrid grid;
    HYPRE_StructStencil stencil;
    HYPRE_StructMatrix A;
    HYPRE_StructVector b, x;
    //HYPRE_StructSolver solver;
    HYPRE_StructSolver krylov_solver;  //MKG June 2025

    std::vector<int> ordered; // Vector to hold ordered indices for shuffling frequency groups
    std::vector<std::vector<double>> T; // Temperature field
    
    //std::vector<std::pair<double, double>> nodes;

    std::vector<std::vector<std::vector<double>>> grad_energy; // Gradient of energy field  // [group][cell][3] 
    std::vector<std::vector<double>> diff_coeff; // Diffusion coefficient field  // [group][cell] corrected eg using Larsen corrector
    std::vector<std::vector<double>> Bag; // Emission spectrum per group  // [group][cell]
    std::vector<std::vector<double>> ddelr; // the divergence of the (energy gradient multiplied by the corrected diffusion coefficient)  // [group][cell]
};
