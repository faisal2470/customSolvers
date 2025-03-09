#include "fvCFD.H"
#include "constants.H"

Foam::vector g(0, -9.81, 0);
const volVectorField& gradB2 = mesh.lookupObject<volVectorField>("gradB2");
thermo.correct();




// Update chi and magneticForce for each species based on the updated temperature
forAll(Y, i)
{
    const word &speciesName = Y[i].name();

    // Retrieve magnetic properties from the dictionary
    const dictionary& speciesDict = magneticProperties.subDict("speciesProperties");
    const dictionary& speciesSubDict = speciesDict.subDict(speciesName);

    scalar Si = speciesSubDict.lookupOrDefault<scalar>("Si", 1.0);
    dimensionedScalar mi(
        "mi_" + speciesName, dimMolarMass,
        speciesSubDict.lookupOrDefault<scalar>("mi", 1.0e-3)
    );

    // Retrieve the chi and magneticForce fields for each species
    volScalarField &chi = const_cast<volScalarField&>(
        mesh.lookupObject<volScalarField>("chi_" + speciesName)
    );
    volVectorField &magneticForce = const_cast<volVectorField&>(
        mesh.lookupObject<volVectorField>("magneticForce_" + speciesName)
    );

    // Update chi based on the current temperature
    //const scalar T_min = SMALL; // Prevent division by zero

    
    forAll(chi, cellI)
    {
        scalar T_cell = thermo.T()[cellI];
           
        chi[cellI] = (NA.value() * pow(gL.value(), 2) * pow(muB.value(), 2) * Si * (Si + 1) * mu0.value())
                     / (3 * k.value() * T_cell * mi.value());
       //Info<<"T_cell = " <<T[cellI]<<" for chi = "<<chi[cellI]<<endl;
       
       
    }
    
    // Update magneticForce based on the updated chi and grad(B^2)
    forAll(magneticForce, cellI)
{
    //const scalar y = mesh.C()[cellI].y(); // Get the y-coordinate of the cell center

    // Check if gradB2 has meaningful values
         

    // Compute the magnetic force
    magneticForce[cellI] = (0.5 / mu0.value() * chi[cellI] * gradB2[cellI]) + g;

    // Debug the magneticForce value

}
    
}
