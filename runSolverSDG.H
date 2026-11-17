const label nLocalCells  = mesh.nCells();
const label nGlobalCells = returnReduce(nLocalCells, sumOp<label>());
const label nConstraints = 1;

List<scalarField> allEta   (Pstream::nProcs());
List<scalarField> allSens (Pstream::nProcs());
List<scalarField> allSensV (Pstream::nProcs());

allEta   [Pstream::myProcNo()] = eta.primitiveField();
allSens [Pstream::myProcNo()] = npSens.primitiveField();
allSensV [Pstream::myProcNo()] = npSensV.primitiveField();

Pstream::gatherList(allEta);
Pstream::gatherList(allSens);
Pstream::gatherList(allSensV);

std::vector<double> xval  (nGlobalCells);
std::vector<double> dfdx  (nGlobalCells);
std::vector<double> gx    (nConstraints);
std::vector<double> dgdx  (nGlobalCells * nConstraints, 0.0);
std::vector<double> xmin  (nGlobalCells, 0.0001); 
std::vector<double> xmax  (nGlobalCells, 0.9999);


if (Pstream::master())
{
    label globalIdx = 0;
    for (label proci = 0; proci < Pstream::nProcs(); ++proci)
    {
        const scalarField& eProc = allEta[proci];
        const scalarField& pProc = allSens[proci];
        const scalarField& vProc = allSensV[proci];

        for (label i = 0; i < eProc.size(); ++i)
        {
            xval[globalIdx] = eProc[i];
            dfdx[globalIdx] = pProc[i];
            dgdx[globalIdx * nConstraints + 0] = vProc[i];
            ++globalIdx;
        }
    }

    // Calcular restricciones (Asegúrate de que currentDP, tarDP, etc. estén definidos antes)
    gx[0] = currentVF - tarVF;

    // Debug info
    scalar sumDfdx = 0.0, sumDgdx0 = 0.0;
    for (label i = 0; i < nGlobalCells; ++i)
    {
        sumDfdx += dfdx[i];
        sumDgdx0 += dgdx[i * nConstraints + 0];
    }
    Info<< "DEBUG SDG INPUT (GLOBAL):" << nl
        << "  dfdx sum=" << sumDfdx << nl
        << "  dgdx[0] sum=" << sumDgdx0 << nl
        << "  gx[0]=" << gx[0] << nl << endl;
}


// Parámetros de control para SDG (Ajustables)
scalar penV = readScalar(configOptimization.lookup("penalizationVF"));   // Factor de penalización para Volumen
scalar stepSize  = readScalar(configOptimization.lookup("stepSize"));   // Tamaño de paso base (Learning rate)
scalar moveLimit = readScalar(configOptimization.lookup("moveLimit"));  // Límite máximo de cambio de eta por iteración (5%)

if (Pstream::master())
{
    for (label i = 0; i < nGlobalCells; ++i)
    {
        double sens_total = dfdx[i]; 

        if (gx[0] > 0.0)
        {
            sens_total += 2.0 * penV * gx[0] * dgdx[i * nConstraints + 0];
        }

        // Dirección de descenso (restamos el gradiente)
        double step = -stepSize * sens_total;

        // Aplicar límite de movimiento (Move Limit) para evitar oscilaciones
        if (step > moveLimit)  step = moveLimit;
        if (step < -moveLimit) step = -moveLimit;

        // Actualizar variable de diseño
        double x_new = xval[i] + step;

        // Clamp a los límites físicos [0.0001, 0.9999]
        xval[i] = std::max(0.0001, std::min(0.9999, x_new));
    }

    Info<< "SDG (global): nVars=" << nGlobalCells
        << " min=" << *std::min_element(xval.begin(), xval.end())
        << " max=" << *std::max_element(xval.begin(), xval.end()) << nl;
}


scalarField xvalLocal(nLocalCells);
if (Pstream::nProcs() > 1)
{
    if (Pstream::master())
    {
        label offset = nLocalCells;
        for (label proci = 1; proci < Pstream::nProcs(); ++proci)
        {
            label nProcCells = allEta[proci].size();
            scalarField procData(nProcCells);
            for (label i = 0; i < nProcCells; ++i)
            {
                procData[i] = xval[offset + i];
            }
            OPstream toProc(Pstream::commsTypes::blocking, proci);
            toProc << procData;
            offset += nProcCells;
        }
        for (label i = 0; i < nLocalCells; ++i)
        {
            xvalLocal[i] = xval[i];
        }
    }
    else
    {
        IPstream fromMaster(Pstream::commsTypes::blocking, Pstream::masterNo());
        fromMaster >> xvalLocal;
    }
}
else
{
    for (label i = 0; i < nLocalCells; ++i)
    {
        xvalLocal[i] = xval[i];
    }
}

// Actualizar el campo eta local
scalarField& etaI = eta.primitiveFieldRef();
for (label i = 0; i < nLocalCells; ++i)
{
    etaI[i] = xvalLocal[i];
}

setCells(eta, inletCells, scalar(1.0));
setCells(eta, outletCells, scalar(1.0));
eta.correctBoundaryConditions();