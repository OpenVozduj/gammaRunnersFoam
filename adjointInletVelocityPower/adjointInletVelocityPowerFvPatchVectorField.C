/*---------------------------------------------------------------------------*\
    adjointRunnerFoam - BC adjunta de velocidad en la entrada primal
\*---------------------------------------------------------------------------*/

#include "adjointInletVelocityPowerFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointInletVelocityPowerFvPatchVectorField::
adjointInletVelocityPowerFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}


Foam::adjointInletVelocityPowerFvPatchVectorField::
adjointInletVelocityPowerFvPatchVectorField
(
    const adjointInletVelocityPowerFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::adjointInletVelocityPowerFvPatchVectorField::
adjointInletVelocityPowerFvPatchVectorField
(
    const adjointInletVelocityPowerFvPatchVectorField& pivpvf,
    const DimensionedField<vector, fvMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointInletVelocityPowerFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi");

    const fvsPatchField<scalar>& phiap =
        patch().lookupPatchField<surfaceScalarField, scalar>("phia");

    const fvPatchField<vector>& Uap =
        patch().lookupPatchField<volVectorField, vector>("Ua");

    const fvPatchField<scalar>& nuP =
        patch().lookupPatchField<volScalarField, scalar>("nuEff");

    const scalarField& magSf = patch().magSf();
    const vectorField  nf(patch().nf());

    const scalarField Un (phip /magSf);   // U.n  (negativo en la entrada)
    const scalarField Uan(phiap/magSf);   // Ua.n (fijado por continuidad adjunta)
    const scalarField nuD(nuP*patch().deltaCoeffs());

    // Tangencial del vecino interior
    const vectorField UaN(Uap.patchInternalField());
    const vectorField UaNt(UaN - (UaN & nf)*nf);

    // Robin con [dJ]_t = 0:  (Un + nu/d)*Ua_t = (nu/d)*Ua_t^N
    // Si el denominador se anula -> zeroGradient (fallback seguro)
    vectorField Uat(UaNt);

    forAll(Un, i)
    {
        const scalar den = Un[i] + nuD[i];

        if (mag(den) > SMALL)
        {
            Uat[i] = nuD[i]*UaNt[i]/den;
        }
    }

    operator==(Uat + Uan*nf);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::adjointInletVelocityPowerFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointInletVelocityPowerFvPatchVectorField
    );
}

// ************************************************************************* //