#include "fvCFD.H"
#include "constants.H"


thermo.correct();


// Update chi for each species based on the updated temperature
forAll(Y, i){
    const word &speciesName = Y[i].name();

    // Retrieve magnetic properties from the dictionary
    const dictionary& speciesDict = magneticProperties.subDict("speciesProperties");
    const dictionary& speciesSubDict = speciesDict.subDict(speciesName);

    scalar Si = speciesSubDict.lookupOrDefault<scalar>("Si", 1.0);

    // Retrieve the chi field for each species
    volScalarField &chi = const_cast<volScalarField&>(
        mesh.lookupObject<volScalarField>("chi_" + speciesName)
    );

    // Define molecular weight as a dimensioned scalar
    const dimensionedScalar W_i("W_i", dimensionSet(1, 0, 0, 0, -1, 0, 0), thermo.composition().W(i));

    // Update chi based on the current temperature
    chi = (1000 * NA * pow(gL, 2) * pow(muB, 2) * Si * (Si + 1) * mu0) / (3 * k * thermo.T() * W_i);

}
