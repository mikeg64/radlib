// radiation_step.cpp
#include "../include/radiation_step.h"


void run_radiation_step(IRadSolver& solver,
                        Mesh& mesh,
                        State& state,
                        Pars& pars)
{
    solver.solveRadiationTransport(mesh, state, pars, pars.time);
    solver.UpdateBEmission(mesh, state, pars);
    solver.UpdateRadFlux(mesh, state, pars);
}