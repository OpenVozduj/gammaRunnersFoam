const labelList& inletCells = mesh.boundary()["inlet"].faceCells();
const labelList& outletCells = mesh.boundary()["outlet"].faceCells();

const scalar tarVF = readScalar(configOptimization.lookup("targetFluidVF"));

volScalarField diffAlpha
(
    IOobject
    (
        "diffAlpha",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless/dimTime, 0.0)
);

diffAlpha = -alphaMax * (1.0 + qAlpha) * qAlpha 
            / Foam::pow(eta + qAlpha, 2);

volScalarField sens
(
    IOobject
    (
        "sens",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    (U & Ua) * diffAlpha
);

volScalarField npSens
(
    IOobject
    (
        "npSens",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    sens
);

volScalarField npSensV
(
    IOobject
    (
        "npSensV",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 1.0)
);

volScalarField etaFil
(
    IOobject
    (
        "etaFil",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    eta,
    "zeroGradient"
);

int opt = 0;
scalar beta = readScalar(configOptimization.lookup("beta"));
scalar eta0 = tarVF;

scalar Jpower = 0.0;
scalar Jout   = 0.0;
