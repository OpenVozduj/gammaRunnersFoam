/*---------------------------------------------------------------------------*\
    adjointRunnerFoam
    BC adjunta de velocidad en las SALIDAS primales (velocidad prescrita).

    Como delta(U) = 0 alli, el objetivo solo aporta via delta(p):
        [dJ]_p = -Un   =>   Ua.n = -U.n   =>   phia = -phi
    Tangencial: zeroGradient.
\*---------------------------------------------------------------------------*/

#ifndef adjointOutletVelocityFixedFluxFvPatchVectorField_H
#define adjointOutletVelocityFixedFluxFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class adjointOutletVelocityFixedFluxFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class adjointOutletVelocityFixedFluxFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("adjointOutletVelocityFixedFlux");


    // Constructors

        //- Construct from patch, internal field and dictionary
        adjointOutletVelocityFixedFluxFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        adjointOutletVelocityFixedFluxFvPatchVectorField
        (
            const adjointOutletVelocityFixedFluxFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, fvMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        adjointOutletVelocityFixedFluxFvPatchVectorField
        (
            const adjointOutletVelocityFixedFluxFvPatchVectorField&,
            const DimensionedField<vector, fvMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, fvMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new adjointOutletVelocityFixedFluxFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using fixedValueFvPatchVectorField::operator=;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //