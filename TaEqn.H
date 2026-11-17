{
    const volSymmTensorField D (symm(fvc::grad(U)));
    const volSymmTensorField Da(symm(fvc::grad(Ua)));

    // Unico termino de acoplamiento retenido:  2*(dnu/dT)*(D:Da)
    volScalarField thermalSource
    (
        IOobject("thermalSource", runTime.name(), mesh),
        2.0*dnuDT*(D && Da)
    );

    fvScalarMatrix TaEqn
    (
        fvm::div(-phi, Ta)          // adveccion invertida
      - fvm::laplacian(DT, Ta)
     ==
      - thermalSource
    );

    TaEqn.relax();
    TaEqn.solve();
    Ta.correctBoundaryConditions();  // Ta = 0 en todo el contorno
}