{
    volVectorField adjointTransposeConvection(fvc::grad(Ua) & U);

    fvVectorMatrix UaEqn
    (
        fvm::div(-phi, Ua)
      - adjointTransposeConvection
      + turbulence->divDevSigma(Ua)     // <- misma nuEff, gradiente de Ua
      + fvm::Sp(alpha, Ua)
    );

    UaEqn.relax();
    solve(UaEqn == -fvc::grad(pa));

    volScalarField rAUa(1.0/UaEqn.A());
    volVectorField HbyAa("HbyAa", rAUa*UaEqn.H());
    surfaceScalarField phiHbyAa("phiHbyAa", fvc::flux(HbyAa));

    // En las salidas primales: [dJ]_p = -Un  =>  phia = -phi
    phiHbyAa.boundaryFieldRef()[outletPatchID] == -phi.boundaryField()[outletPatchID];

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix paEqn(fvm::laplacian(rAUa, pa) == fvc::div(phiHbyAa));
        paEqn.setReference(paRefCell, paRefValue);
        paEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phia = phiHbyAa - paEqn.flux();
        }
    }

    //#include "adjointContinuityErrs.H"

    pa.relax();
    Ua = HbyAa - rAUa*fvc::grad(pa);
    Ua.correctBoundaryConditions();
}