// diagnostico.H  -- incluir despues de valFunctions.H
{
    const scalarField& V = mesh.V();
    const scalar tz   = mesh.bounds().span().z();
    const scalar Vavg = gSum(V)/returnReduce(mesh.nCells(), sumOp<label>());
    const scalar dx   = Foam::sqrt(Vavg/tz);

    const scalar nuAvg = gAverage(nuEff.primitiveField());
    const scalar nuHi  = gMax(nuEff.primitiveField());
    const scalar aHi   = gMax(alpha.primitiveField());

    Info<< "=== DIAGNOSTICO ===" << nl
        << "  dx                    = " << dx*1e3 << " mm" << nl
        << "  nuEff min/avg/max     = " << gMin(nuEff.primitiveField())
        << " / " << nuAvg << " / " << nuHi << " m2/s" << nl
        << "  alpha max             = " << aHi << " 1/s" << nl
        << "  delta = sqrt(nu/alph) = " << Foam::sqrt(nuAvg/max(aHi,SMALL))*1e3
        << " mm   (necesita << dx)" << nl
        << "  Pi = alpha*dx^2/nu    = " << aHi*sqr(dx)/max(nuAvg,SMALL)
        << "      (necesita >> 1)" << nl
        << "  |U| max               = " << gMax(mag(U.primitiveField()))  << " m/s" << nl
        << "  |Ua| max              = " << gMax(mag(Ua.primitiveField())) << " m/s" << nl
        << "  p inlet / outlet      = " << gAverage(p.boundaryField()[inletPatchID])
        << " / " << gAverage(p.boundaryField()[outletPatchID]) << nl
        << "  Q inlet / outlet      = " << gSum(phi.boundaryField()[inletPatchID])
        << " / " << gSum(phi.boundaryField()[outletPatchID]) << nl
        << "  phia/phi outlet       = "
        << gSum(phia.boundaryField()[outletPatchID])
          /max(gSum(phi.boundaryField()[outletPatchID]),SMALL)
        << "   (debe ser -1.0)" << nl
        << endl;
}