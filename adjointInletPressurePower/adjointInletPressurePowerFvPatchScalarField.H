/*---------------------------------------------------------------------------*\
 *  H   H   55555                                                              
 *  H   H   5                                                                  
 *  HHHHH   55555                                                             
 *  H   H       5                                                              
 *  H   H   55555 
 *
 *  Solver: hsTopOptBeta
 *  Autors: 
\*---------------------------------------------------------------------------*/

#ifndef adjointInletPressurePowerFvPatchScalarField_H
#define adjointInletPressurePowerFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class adjointInletPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class adjointInletPressurePowerFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    scalar viscousCoeff_; 
    
public:

    //- Runtime type information
    TypeName("adjointInletPressurePower");


    // Constructors

        //- Construct from patch, internal field and dictionary
        adjointInletPressurePowerFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given adjointInletPressurePowerFvPatchScalarField
        //  onto a new patch
        adjointInletPressurePowerFvPatchScalarField
        (
            const adjointInletPressurePowerFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        adjointInletPressurePowerFvPatchScalarField
        (
            const adjointInletPressurePowerFvPatchScalarField&,
            const DimensionedField<scalar, fvMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, fvMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new adjointInletPressurePowerFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
