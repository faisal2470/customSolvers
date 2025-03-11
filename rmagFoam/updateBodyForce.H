#include "fvCFD.H"
#include "constants.H"
#include "updateMagneticProperties.H" // Update chi for each species based on the updated temperature

if (!mesh.foundObject<uniformDimensionedVectorField>("g"))
{
    FatalErrorInFunction << "Gravity field 'g' not found in constant/" << exit(FatalError);
}

const uniformDimensionedVectorField& g = mesh.lookupObject<uniformDimensionedVectorField>("g");
const volVectorField& gradB2 = mesh.lookupObject<volVectorField>("gradB2");
forAll(Y, i)
{
    const word &speciesName = Y[i].name();

    // Retrieve the chi and bodyForce fields for each species
    volScalarField &chi = const_cast<volScalarField&>(
        mesh.lookupObject<volScalarField>("chi_" + speciesName)
    );

    volVectorField &bodyForce = const_cast<volVectorField&>(
        mesh.lookupObject<volVectorField>("bodyForce_" + speciesName)
    );

    // Update bodyForce based on the updated chi and grad(B^2)
    // Info << "Updating bodyForce for species " << speciesName << endl;
    bodyForce = (0.5 / mu0) * chi * gradB2 + g;

}

