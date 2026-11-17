/*---------------------------------------------------------------------------*\
    adjointRunnerFoam
    BC adjunta de velocidad en la ENTRADA primal (presion prescrita).

    Tangencial: (Un + nu/d)*Ua_t = (nu/d)*Ua_t^N     ([dJ]_t = 0 aqui)
    Normal:     Ua.n = phia/|Sf|                     (continuidad adjunta)

    Para fundido nu/d >> |Un|  ->  practicamente zeroGradient.
\*---------------------------------------------------------------------------*/

#ifndef adjointInletVelocityPowerFvPatchVectorField_H
#define adjointInletVelocityPowerFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class adjointInletVelocityPowerFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class adjointInletVelocityPowerFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("adjointInletVelocityPower");


    // Constructors

        //- Construct from patch, internal field and dictionary
        adjointInletVelocityPowerFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        adjointInletVelocityPowerFvPatchVectorField
        (
            const adjointInletVelocityPowerFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, fvMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        adjointInletVelocityPowerFvPatchVectorField
        (
            const adjointInletVelocityPowerFvPatchVectorField&,
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
                new adjointInletVelocityPowerFvPatchVectorField(*this, iF)
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