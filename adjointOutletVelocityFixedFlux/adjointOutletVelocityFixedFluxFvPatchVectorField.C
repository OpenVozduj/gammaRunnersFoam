/*---------------------------------------------------------------------------*\
    adjointRunnerFoam - BC adjunta de velocidad en las salidas primales
\*---------------------------------------------------------------------------*/

#include "adjointOutletVelocityFixedFluxFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointOutletVelocityFixedFluxFvPatchVectorField::
adjointOutletVelocityFixedFluxFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}


Foam::adjointOutletVelocityFixedFluxFvPatchVectorField::
adjointOutletVelocityFixedFluxFvPatchVectorField
(
    const adjointOutletVelocityFixedFluxFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::adjointOutletVelocityFixedFluxFvPatchVectorField::
adjointOutletVelocityFixedFluxFvPatchVectorField
(
    const adjointOutletVelocityFixedFluxFvPatchVectorField& pivpvf,
    const DimensionedField<vector, fvMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointOutletVelocityFixedFluxFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi");

    const fvPatchField<vector>& Uap =
        patch().lookupPatchField<volVectorField, vector>("Ua");

    const vectorField nf(patch().nf());
    const scalarField Un(phip/patch().magSf());

    // Tangencial: zeroGradient (delta(U) = 0 en la salida primal)
    const vectorField UaN(Uap.patchInternalField());
    const vectorField UatN(UaN - (UaN & nf)*nf);

    // Normal: Ua.n = -U.n   (porque [dJ]_p = -Un)
    operator==(UatN - Un*nf);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::adjointOutletVelocityFixedFluxFvPatchVectorField::write(Ostream& os) const
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
        adjointOutletVelocityFixedFluxFvPatchVectorField
    );
}

// ************************************************************************* //