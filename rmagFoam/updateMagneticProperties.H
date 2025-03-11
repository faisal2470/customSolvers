#include "fvCFD.H"
#include "constants.H"


thermo.correct();


// Update chi for each species based on the updated temperature
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

    // Retrieve the chi field for each species
    volScalarField &chi = const_cast<volScalarField&>(
        mesh.lookupObject<volScalarField>("chi_" + speciesName)
    );

    forAll(chi, cellI){
        chi[cellI] = (NA.value() * pow(gL.value(), 2) * pow(muB.value(), 2) * Si * (Si + 1) * mu0.value())
                     / (3 * k.value() * thermo.T()[cellI] * mi.value());
    }

}
