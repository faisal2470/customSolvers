#include "fvCFD.H"
#include "constants.H"
#include "updateMagneticProperties.H" // Update chi for each species based on the updated temperature

Foam::vector g(0, -9.81, 0);
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
    forAll(bodyForce, cellI){
        bodyForce[cellI] = (0.5 / mu0.value() * chi[cellI] * gradB2[cellI]) + g;
    }

}

