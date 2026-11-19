fvScalarMatrix sensFilterEqn_DP
(
    fvm::laplacian(fSens)
  - fvm::Sp(1.0/sqr(R), fSens)
  ==
  - npSens / (R*R)
);
sensFilterEqn_DP.solve();

fvScalarMatrix sensFilterEqn_vf
(
    fvm::laplacian(fSensV)
  - fvm::Sp(1.0/sqr(R), fSensV)
  ==
  - npSensV / (R*R)
);
sensFilterEqn_vf.solve();

fSens.correctBoundaryConditions();
fSensV.correctBoundaryConditions(); 