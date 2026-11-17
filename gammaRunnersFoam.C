/*---------------------------------------------------------------------------*\
    gammaRunnersFoam
    Solver adjunto 2D para optimizacion topologica conceptual de runners
    de inyeccion. Laminar, incompresible, no-newtoniano termico.

    Base: OpenFOAM-14 applications/legacy/incompressible/
          adjointShapeOptimisationFoam

    Objetivo: J = -sum_out int (p + 0.5|U|^2) U.n dGamma   (minimizar)
    Primal BCs: p fija en inlet, U fija en outlets, T=Twall en paredes.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "viscosityModel.H"
#include "incompressibleMomentumTransportModels.H"
#include "simpleControl.H"
#include "pressureReference.H"
#include "findRefCell.H"
#include "constrainPressure.H"
#include "constrainHbyA.H"
#include "adjustPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"

#include "fvcDdt.H"
#include "fvcGrad.H"
#include "fvcFlux.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

#include <vector>
#include <algorithm>
#include <cmath>

#include "OFstream.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline void zeroCells
(
    VolField<Type>& vf,
    const labelList& cells
)
{
    forAll(cells, i)
    {
        vf[cells[i]] = Zero;
    }
}

template<class Type>
inline void setCells
(
    VolField<Type>& vf,
    const labelList& cells,
    const Type& value
)

{
    forAll(cells, i)
    {
        vf[cells[i]] = value;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"

    #include "createFields.H"
    #include "initOptimization.H"
    #include "initContinuityErrs.H"

    turbulence->validate();
    
    Info<< "\nStarting time loop\n" << endl;
    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        turbulence->predict();

        #include "UEqn.H"
        #include "TEqn.H"
        
        #include "TaEqn.H"
        #include "UaEqn.H"
  
        #include "sens.H"
        #include "valFunctions.H"
        #include "runSolverSDG.H"
        #include "lastUpdateEta.H"
        #include "penaltiesFunc.H"

        viscosity->correct();
        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
