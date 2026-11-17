{
    fvVectorMatrix UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevSigma(U)      // <- laplaciano + dev2(T(grad U))
      + fvm::Sp(alpha, U)
    );

    UEqn.relax();
    solve(UEqn == -fvc::grad(p));

    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

    adjustPhi(phiHbyA, U, p);
    constrainPressure(p, U, phiHbyA, rAU);

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(fvm::laplacian(rAU, p) == fvc::div(phiHbyA));
        pEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );
        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "continuityErrs.H"

    p.relax();
    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
}