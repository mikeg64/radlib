#include <gtest/gtest.h>
#include <gmock/gmock.h>



#include "../../../include/radiation_step.h"
#include "../../../include/solver.h"
#include "../../../include/physics.h"
#include "../../../include/geometry.h"
#include "../../../include/setup.h"
#include "mockradsolver.h"

using ::testing::_;
using ::testing::Exactly;
using ::testing::Return;
using ::testing::Invoke;




TEST(RadiationIntegration, PipeTwoMaterials_MockedHeating)
{
    Mesh mesh(200, 20, 0.1, 0.1);
    Pars pars;
    pars.dt = 1e-6;
    pars.time = 0.0;

    Materials materials = initialize_materials(mesh);
    State state = initialize_physics(mesh, materials, pars);

    MockRadSolver mockSolver;

    // Fake heating: +1 K per call
    EXPECT_CALL(mockSolver, solveRadiationTransport(_, _, _, _))
        .Times(100)
        .WillRepeatedly(::testing::Invoke(
            [](const Mesh&, State& s, Pars&, double){
                for (auto& T : s.temperature)
                    T += 1.0;
            }
        ));

    EXPECT_CALL(mockSolver, UpdateBEmission(_, _, _))
        .Times(100);

    EXPECT_CALL(mockSolver, UpdateRadFlux(_, _, _))
        .Times(100);

    // Run the loop using dependency injection
    for (int step = 0; step < 100; step++) {
        run_radiation_step(mockSolver, mesh, state, pars);
        pars.time += pars.dt;
    }

    // Check that heating happened
    EXPECT_NEAR(state.temperature[0], 100.0, 1e-12);
}