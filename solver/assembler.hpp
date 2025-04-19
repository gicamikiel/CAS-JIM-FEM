#ifndef SOLVER_ASSEMBLER
#define SOLVER_ASSEMBLER

#include <assembler.hpp>
#include <Kokkos_Core.hpp>

using team_policy = Kokkos::TeamPolicy<>;
using member_type = Kokkos::TeamPolicy<>::member_type;

using ScratchViewType = Kokkos::View<double*[4],
                            Kokkos::DefaultExecutionSpace::scratch_memory_space,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

template<typename ElmDef>
void assembleK(Kokkos::View<double*> kappa, Kokkos::View<std::size_t*> nodeTags, Kokkos::View<double*[2]> nodeCoords, 
    std::size_t numElm, Kokkos::View<double**> conn, Kokkos::View<double**> stiffness) {

    ElmDef e;
    
    int scratchSize = ScratchViewType::shmem_size(e.numQuad);
    Kokkos::parallel_for("elements", team_policy(numElm, Kokkos::AUTO).set_scratch_size(0, Kokkos::PerTeam(scratchSize)),
                            KOKKOS_LAMBDA(const member_type &teamMember) {
        const std::size_t elmIdx = teamMember.league_rank();
        ScratchViewType derivatives(e.numQuad);

        double elmKappa = kappa(elmIdx);

        Kokkos::parallel_for("jacobianEval", Kokkos::TeamThreadRange(teamMember, e.jacAddendEvalsN()), [&] (const int i) {
            std::size_t node;
            std::size_t quadPt;
            std::size_t derivative;
            e.unwrapAddendEvalsN(i, node, quadPt, derivative);
            double xi;
            double eta;
            double wt;
            e.quadPtInfo(quadPt, xi, eta, wt);

            double addend = 0;
            std::size_t coordIdx = conn(elmIdx, node)-1;
            switch(derivative) {
                case 0: addend = e.jacAddendDelxDelXi(node, xi, eta, nodeCoords(coordIdx, 0)); break; // takes xa, a in ad-bc
                case 1: addend = e.jacAddendDelyDelEta(node, xi, eta, nodeCoords(coordIdx, 1)); break; // takes ya, d in ad-bc
                case 2: addend = e.jacAddendDelxDelEta(node, xi, eta, nodeCoords(coordIdx, 0)); break; // takes xa, b in ad-bc
                case 3: addend = e.jacAddendDelyDelXi(node, xi, eta, nodeCoords(coordIdx, 1)); break; // takes ya, c in ad-bc
            }
            Kokkos::atomic_add(&derivatives(quadPt, derivative), addend);
        });

        teamMember.team_barrier();

        Kokkos::parallel_for("stiffnessAsy", Kokkos::TeamThreadRange(teamMember, e.stiffnessEvalsN()), [&] (const int i) {
            std::size_t quadPt;
            std::size_t localRow;
            std::size_t localCol;
            e.unwrapStiffnessN(i, quadPt, localRow, localCol);
            double xi;
            double eta;
            double wt;
            e.quadPtInfo(quadPt, xi, eta, wt);

            double jacobian = derivatives(quadPt,0)*derivatives(quadPt,1)-derivatives(quadPt,2)*derivatives(quadPt,3);
            double eval = wt*elmKappa*e.stiffnessIntegrand(localRow, localCol, xi, eta, jacobian);
            
            std::size_t globalRow = conn(elmIdx, localRow)-1;
            std::size_t globalCol = conn(elmIdx, localCol)-1;
            Kokkos::atomic_add(&stiffness(globalRow, globalCol), eval);
            if(globalRow != globalCol) {
                Kokkos::atomic_add(&stiffness(globalCol, globalRow), eval);
            }
        });
    });
}

template<typename ElmDef>
void assembleF(Kokkos::View<double*> elmForces, Kokkos::View<std::size_t*> nodeTags, Kokkos::View<double*[2]> nodeCoords, 
    std::size_t numElm, Kokkos::View<double**> conn, Kokkos::View<double*> forcing) {

    ElmDef e;
    
    int scratchSize = ScratchViewType::shmem_size(e.numQuad);
    Kokkos::parallel_for("elements", team_policy(numElm, Kokkos::AUTO).set_scratch_size(0, Kokkos::PerTeam(scratchSize)),
                            KOKKOS_LAMBDA(const member_type &teamMember) {
        const std::size_t elmIdx = teamMember.league_rank();
        ScratchViewType derivatives(e.numQuad);

        double elmForce = elmForces(elmIdx);

        Kokkos::parallel_for("jacobianEval", Kokkos::TeamThreadRange(teamMember, e.jacAddendEvalsN()), [&] (const int i) {
            std::size_t node;
            std::size_t quadPt;
            std::size_t derivative;
            e.unwrapAddendEvalsN(i, node, quadPt, derivative);
            double xi;
            double eta;
            double wt;
            e.quadPtInfo(quadPt, xi, eta, wt);

            double addend = 0;
            std::size_t coordIdx = conn(elmIdx, node)-1;
            switch(derivative) {
                case 0: addend = e.jacAddendDelxDelXi(node, xi, eta, nodeCoords(coordIdx, 0)); break; // takes xa, a in ad-bc
                case 1: addend = e.jacAddendDelyDelEta(node, xi, eta, nodeCoords(coordIdx, 1)); break; // takes ya, d in ad-bc
                case 2: addend = e.jacAddendDelxDelEta(node, xi, eta, nodeCoords(coordIdx, 0)); break; // takes xa, b in ad-bc
                case 3: addend = e.jacAddendDelyDelXi(node, xi, eta, nodeCoords(coordIdx, 1)); break; // takes ya, c in ad-bc
            }
            Kokkos::atomic_add(&derivatives(quadPt, derivative), addend);
        });

        teamMember.team_barrier();

        Kokkos::parallel_for("fAsy", Kokkos::TeamThreadRange(teamMember, e.forcingTotalN()), [&] (const int i) {
            std::size_t node;
            std::size_t quadPt;
            e.unwrapForcingN(i, node, quadPt);

            double xi;
            double eta;
            double wt;
            e.quadPtInfo(i, xi, eta, wt);

            double jacobian = derivatives(quadPt,0)*derivatives(quadPt,1)-derivatives(quadPt,2)*derivatives(quadPt,3);
            double eval = elmForce*e.shape(node, xi, eta)*wt;
            
            std::size_t globalIdx = conn(elmIdx, node)-1;
            Kokkos::atomic_add(&forcing(globalIdx), eval);
        });
    });
}

#endif