{
    // Adveccion-difusion pura. SIN disipacion viscosa (no hay Cp).
    fvScalarMatrix TEqn
    (
        fvm::div(phi, T)
      - fvm::laplacian(DT, T)
    );

    TEqn.relax();
    TEqn.solve();
    T.correctBoundaryConditions();
}