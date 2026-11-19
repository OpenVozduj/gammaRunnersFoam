// Adjoint sensitivities
diffAlpha = -alphaMax * (1.0 + qAlpha) * qAlpha 
            / Foam::pow(eta + qAlpha, 2);

sens = -(U&Ua) * diffAlpha;

// Normalize sensitivities
scalar normS = gMax(mag(sens.primitiveField()));
normS = max(normS, 1e-8);
npSens = sens / normS;

npSensV.primitiveFieldRef() = 1.0;
npSensV.correctBoundaryConditions();

#include "diffHTProj.H"

#include "sensHFilter.H"
