/*---------------------------------------------------------------------------*\
    adjointRunnerFoam - BC adjunta de presion en la entrada primal
\*---------------------------------------------------------------------------*/

#include "adjointInletPressurePowerFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointInletPressurePowerFvPatchScalarField::
adjointInletPressurePowerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    viscousCoeff_(dict.lookupOrDefault<scalar>("viscousCoeff", 0.0))
{}


Foam::adjointInletPressurePowerFvPatchScalarField::
adjointInletPressurePowerFvPatchScalarField
(
    const adjointInletPressurePowerFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    viscousCoeff_(ptf.viscousCoeff_)
{}


Foam::adjointInletPressurePowerFvPatchScalarField::
adjointInletPressurePowerFvPatchScalarField
(
    const adjointInletPressurePowerFvPatchScalarField& ptfsf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptfsf, iF),
    viscousCoeff_(ptfsf.viscousCoeff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointInletPressurePowerFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi");

    const fvsPatchField<scalar>& phiap =
        patch().lookupPatchField<surfaceScalarField, scalar>("phia");

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U");

    const fvPatchField<vector>& Uap =
        patch().lookupPatchField<volVectorField, vector>("Ua");

    const scalarField& magSf = patch().magSf();

    const scalarField Un (phip /magSf);   // U.n   (negativo en entrada)
    const scalarField Uan(phiap/magSf);   // Ua.n  (fijado por continuidad adjunta)

    // Termino viscoso opcional: viscousCoeff*nu*d(Ua.n)/dn
    scalarField visc(Un.size(), Zero);

    if (viscousCoeff_ > 0.0)
    {
        const fvPatchField<scalar>& nuP =
            patch().lookupPatchField<volScalarField, scalar>("nuEff");

        const scalarField UanN(Uap.patchInternalField() & patch().nf());

        visc = viscousCoeff_*nuP*patch().deltaCoeffs()*(Uan - UanN);
    }

    // [dJ]_n = 0 en la entrada (el objetivo vive en las salidas)
    operator==(Un*Uan + (Up & Uap) + visc);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::adjointInletPressurePowerFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "viscousCoeff", viscousCoeff_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        adjointInletPressurePowerFvPatchScalarField
    );
}

// ************************************************************************* //