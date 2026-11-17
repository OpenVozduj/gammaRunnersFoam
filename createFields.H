Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field pa\n" << endl;
volScalarField pa
(
    IOobject
    (
        "pa",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Ua\n" << endl;
volVectorField Ua
(
    IOobject
    (
        "Ua",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Ta\n" << endl;
volScalarField Ta
(
    IOobject
    (
        "Ta",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading/calculating face flux field phia\n" << endl;
surfaceScalarField phia
(
    IOobject
    (
        "phia",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Ua) & mesh.Sf()
);

label paRefCell = 0;
scalar paRefValue = 0.0;
setRefCell
(
    pa,
    simple.dict(),
    paRefCell,
    paRefValue
);
mesh.schemes().setFluxRequired(pa.name());

volScalarField eta
(
    IOobject
    (
        "eta",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

pressureReference pressureReference(p, simple.dict());
mesh.schemes().setFluxRequired(p.name());

autoPtr<viscosityModel> viscosity(viscosityModel::New(mesh));
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, viscosity)
);

IOdictionary configOptimization
(
    IOobject
    (
        "configOptimization",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

volScalarField nuEff
(
    IOobject("nuEff", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("0", dimArea/dimTime, Zero),
    calculatedFvPatchScalarField::typeName
);

dimensionedScalar zeroSensitivity("0", dimVelocity*dimVelocity, 0.0);

dimensionedScalar etaMin("etaMin", dimless, 0.0);
dimensionedScalar etaMax("etaMax", dimless, 1.0);

dimensionedScalar alphaMax(viscosity->lookup("alphaMax"));
dimensionedScalar DtFluid(viscosity->lookup("DtFluid"));
dimensionedScalar DtSolid(viscosity->lookup("DtSolid"));

const scalar qAlpha = readScalar(configOptimization.lookup("qAlphaStart"));
const scalar qD = readScalar(configOptimization.lookup("qD"));

volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless/dimTime, 0.0)
);

volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength*dimVelocity, 0.0)
);

alpha = alphaMax - alphaMax*eta*(1+qAlpha)/(eta+qAlpha);
DT = DtSolid + (DtFluid - DtSolid)*eta*(1+qD)/(eta+qD);

volScalarField dnuDT
(
    IOobject("dnuDT", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("0", dimArea/dimTime/dimTemperature, Zero),
    calculatedFvPatchScalarField::typeName
);

#include "createFvModels.H"
#include "createFvConstraints.H"