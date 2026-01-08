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



TEST(RadiationIntegration, PipeTwoMaterials_MockedSolver)
{
    Pars pars;
    pars.tini = 2000.0;
    pars.nx = 200;
    pars.ny = 20;
    pars.nz = 1;
    pars.dt = 1e-6;

    int nsteps=10;
    int nint=2; //interval to write vtk files

    Mesh mesh(200, 20, 0.1, 0.1);
    Materials materials = initialize_materials(mesh);

    MaterialProperties copper{0.15, 385.0};
    materials.add_material(2, copper);

    // Assign materials along pipe
    for (int i = 0; i < mesh.num_cells; i++) {
        int ix = i % mesh.nx;
        if (mesh.cells[i].in_pipe) {
            mesh.cells[i].material_id = (ix < mesh.nx/2 ? 1 : 2);
        }
    }

    State state = initialize_physics(mesh, materials, pars);
    state.temperature[10] = 500.0;

    // Create mock solver
    MockRadSolver mockSolver;

    // EXPECT solver calls
    EXPECT_CALL(mockSolver, solveRadiationTransport(_, _, _, _))
        .Times(Exactly(nsteps));

    EXPECT_CALL(mockSolver, UpdateBEmission(_, _, _))
        .Times(Exactly(nsteps));

    EXPECT_CALL(mockSolver, UpdateRadFlux(_, _, _))
        .Times(Exactly(nsteps));

    // Define behavior for solveRadiationTransport to simulate heating
    //- Every time the simulation calls solveRadiationTransport
    //- The mock will pretend to solve radiation transport
    //- And will heat the material by +1 K per call

    EXPECT_CALL(mockSolver, solveRadiationTransport(_, _, _, _))
    .WillRepeatedly(Invoke([](const Mesh&, State& s, Pars&, double){
        // Fake heating
        for (auto& T : s.temperature) T += 1.0;
    }));

    // Run the time loop (but using mocks)
    for (int step = 0; step < nsteps; step++) {
        run_radiation_step(mockSolver, mesh, state, pars);

        //mockSolver.solveRadiationTransport(mesh, state, pars, pars.time);
        //mockSolver.UpdateBEmission(mesh, state, pars);
        //mockSolver.UpdateRadFlux(mesh, state, pars);
        pars.time += pars.dt;
    }

    EXPECT_EQ(pars.time, nsteps * pars.dt);
}


