// IRadSolver.hpp
#pragma once
#include "geometry.h"
#include "physics.h"
#include "setup.h"

class IRadSolver {
public:
    virtual ~IRadSolver() = default;

    virtual void solveRadiationTransport(const Mesh&, State&, Pars&, double) = 0;
    virtual void UpdateBEmission(const Mesh&, State&, Pars&) = 0;
    virtual void UpdateRadFlux(const Mesh&, State&, Pars&) = 0;
};

