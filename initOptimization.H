const labelList& inletCells = mesh.boundary()["inlet"].faceCells();
const labelList& outletCells = mesh.boundary()["outlet"].faceCells();

label inletPatchID  = -1;
label outletPatchID = -1;

forAll(mesh.boundary(), patchi)
{
    const word& pname = mesh.boundary()[patchi].name();

    if (pname == "inlet")  inletPatchID  = patchi;
    if (pname == "outlet") outletPatchID = patchi;
}

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
    -(U & Ua) * diffAlpha
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

volScalarField fSens
(
    IOobject
    (
        "fSens",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    npSens,
    "zeroGradient"
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

volScalarField fSensV
(
    IOobject
    (
        "fSensV",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    npSensV,
    "zeroGradient"
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
const label flowIter = readLabel(configOptimization.lookup("flowIter"));

scalar radiusAbs = readScalar(configOptimization.lookup("radiusAbs"));
dimensionedScalar R("R", dimLength, radiusAbs);

scalar Jpower = 0.0;
scalar Jout   = 0.0;
