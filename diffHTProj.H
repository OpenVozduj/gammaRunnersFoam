// Derivative Hyperbolic Tangent projection
if (opt > 1) // Asegúrate de que esta condición coincida con el inicio de tu proyección
{
    // 1. Pre-calcular constantes (Seguro en paralelo si eta0 es globalmente consistente)
    scalar tanh_beta_eta0 = Foam::tanh(beta * eta0);
    scalar tanh_beta_1_minus_eta0 = Foam::tanh(beta * (1.0 - eta0));
    scalar denominator = tanh_beta_1_minus_eta0 + tanh_beta_eta0;
    denominator = Foam::max(denominator, SMALL);

    // 2. Referencias directas a los campos internos (Evita crear un volScalarField temporal)
    scalarField& sensI = npSens.primitiveFieldRef();
    scalarField& sensVI = npSensV.primitiveFieldRef();
    const scalarField& etaFilI = etaFil.primitiveField();

    // 3. Bucle local y aplicación de la regla de la cadena
    forAll(etaFilI, celli)
    {
        scalar x = etaFilI[celli];
        scalar tanh_arg = beta * (x - eta0);
        scalar tanh_val = Foam::tanh(tanh_arg);
        
        // Derivada de la función tanh proyectada
        scalar dHdx = beta * (1.0 - tanh_val * tanh_val) / denominator;
        
        // Evitar derivadas exactamente cero cuando beta es muy grande
        dHdx = Foam::max(dHdx, SMALL); 

        // Aplicar regla de la cadena directamente a las sensibilidades
        sensI[celli] *= dHdx;
        sensVI[celli] *= dHdx;
    }
}