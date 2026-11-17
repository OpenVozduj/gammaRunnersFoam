// ---------------------------------------------------------------------------
// Objetivo de potencia y chequeo independiente de disipacion volumetrica
//
//   J = -sum_out int (p + 0.5|U|^2) U.n dGamma     [m^5/s^3] cinematico
//
// Nota: las operaciones sobre patch fields devuelven scalar plano (sin
// dimensiones). El numero resultante esta en m^5/s^3; multiplicar por rho
// para obtener Watts.
// ---------------------------------------------------------------------------

Jout = 0.0;

forAll(outletCells, i)
{
    const label pi = outletCells[i];

    Jout += sum
    (
        phi.boundaryField()[pi]
       *(p.boundaryField()[pi] + 0.5*magSqr(U.boundaryField()[pi]))
    );                                       // <- sin .value()
}

Jpower = -Jout;


// Chequeo independiente: disipacion volumetrica int(2*nu*|D|^2) dV
const volScalarField diss
(
    IOobject
    (
        "diss",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    2.0*nuEff*magSqr(symm(fvc::grad(U)))
);

scalar Diss = 0.0;
{
    // DimensionedField IS-A Field -> se puede enlazar a scalarField&
    const scalarField& V  = mesh.V();
    const scalarField& df = diss.primitiveField();

    forAll(df, celli)
    {
        Diss += V[celli]*df[celli];
    }
    // reduce(Diss, sumOp<scalar>());   // solo si algun dia corres en paralelo
}

Info<< "OBJECTIVE  J = " << Jpower << " m^5/s^3"
    << "   |  sum_out(p_t*Q) = " << Jout
    << "   |  int(2*nu*|D|^2)dV = " << Diss
    << endl;

// Metric Fluid Volume fraction
scalar totalVol = gSum(mesh.V());
scalar fluidVol = gSum(eta.primitiveField() * mesh.V());
scalar currentVF = fluidVol / totalVol;
Info<< "Fluid Vf: Initial = " << currentVF << endl; 
