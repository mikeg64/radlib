// MockRadSolver.hpp
#pragma once
#include <gmock/gmock.h>
#include "../../../include/solver.h"


class MockRadSolver : public IRadSolver {
public:
    MOCK_METHOD(void, solveRadiationTransport,
                (const Mesh&, State&, Pars&, double), (override));

    MOCK_METHOD(void, UpdateBEmission,
                (const Mesh&, State&, Pars&), (override));

    MOCK_METHOD(void, UpdateRadFlux,
                (const Mesh&, State&, Pars&), (override));
};