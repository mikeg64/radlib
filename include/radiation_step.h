// radiation_step.hpp
#pragma once

#include "geometry.h"
#include "physics.h"
#include "setup.h"
#include "radsolver.h"


void run_radiation_step(IRadSolver& solver,
                        Mesh& mesh,
                        State& state,
                        Pars& pars);