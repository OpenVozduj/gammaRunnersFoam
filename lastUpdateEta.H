scalar globalTotalVol = gSum(mesh.V());
scalar targetVol = tarVF * globalTotalVol;
scalar tol = 0.01 * targetVol; // Tolerancia del 1% sobre el volumen objetivo

const label iterHP = readLabel(configOptimization.lookup("iterHP"));
scalar eta0L = 1e-4;
scalar eta0U = 1.0 - 1e-4;
scalar eta0 = 0.5; // Valor inicial para la bisección

auto projectedVolTanh = [&](scalar eta0_val) -> scalar
{
    scalar tanh_beta_eta0 = Foam::tanh(beta * eta0_val);
    scalar tanh_beta_1_minus_eta0 = Foam::tanh(beta * (1.0 - eta0_val));
    scalar denom = tanh_beta_1_minus_eta0 + tanh_beta_eta0;
    denom = Foam::max(denom, SMALL);

    scalar volProjLocal = 0.0;
    const scalarField& etaField = etaFil.primitiveField();
    const scalarField& V = mesh.V();
    
    forAll(etaField, celli)
    {
        scalar x = etaField[celli];
        scalar tanh_arg = beta * (x - eta0_val);
        scalar etaProj = (Foam::tanh(tanh_arg) + tanh_beta_eta0) / denom;
        etaProj = Foam::max(0.0, Foam::min(1.0, etaProj));
        volProjLocal += etaProj * V[celli];
    }
    
    // ¡CRÍTICO PARA PARALELO! Reducir el volumen local a global en todos los procesadores
    reduce(volProjLocal, sumOp<scalar>());
    
    return volProjLocal;
};

scalar vol = projectedVolTanh(eta0);
for (label iter = 0; iter < iterHP; ++iter)
{
    if (mag(vol - targetVol) < tol) break;
    
    if (vol > targetVol)
        eta0L = eta0;  
    else
        eta0U = eta0; 
        
    eta0 = 0.5 * (eta0L + eta0U);
    vol = projectedVolTanh(eta0);
}

// Leer factor de relajación del diccionario (si no existe, usa 1.0 = sin relajación)

scalar relaxFactor = readScalar(configOptimization.lookup("relaxFactorHT"));

scalar tanh_beta_eta0 = Foam::tanh(beta * eta0);
scalar tanh_beta_1_minus_eta0 = Foam::tanh(beta * (1.0 - eta0));
scalar denominator = tanh_beta_1_minus_eta0 + tanh_beta_eta0;
denominator = Foam::max(denominator, SMALL);

scalarField& eta_I = eta.primitiveFieldRef();
const scalarField& etaFilI = etaFil.primitiveField();

forAll(eta_I, celli)
{
    scalar x = etaFilI[celli];
    scalar tanh_arg = beta * (x - eta0);
    scalar etaProj = (Foam::tanh(tanh_arg) + tanh_beta_eta0) / denominator;
    etaProj = Foam::max(0.0, Foam::min(1.0, etaProj));

    // Aplicar sub-relajación (Under-relaxation)
    scalar etaNew = (1.0 - relaxFactor) * eta_I[celli] + relaxFactor * etaProj;

    // Aplicar thresholding (corte duro) DESPUÉS de la relajación 
    // para asegurar que el campo físico final sea binario (blanco/negro)
    if (etaNew > 0.999) etaNew = 1.0;
    if (etaNew < 0.001) etaNew = 0.0;
    
    // Clamp final de seguridad
    eta_I[celli] = Foam::max(0.0, Foam::min(1.0, etaNew));
}

setCells(eta, inletCells, scalar(1.0));
setCells(eta, outletCells, scalar(1.0));

scalar projVol = gSum(eta.primitiveField() * mesh.V());
scalar projVF = projVol / totalVol;

Info<< "Fluid vf: Final = " << projVF << endl; 