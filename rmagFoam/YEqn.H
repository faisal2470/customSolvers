#include "dimensionedTypes.H"

// Define the dimensioned types for the mean molar mass and scaling factor
const dimensionSet dimAmount(0, 0, 0, 0, 0, 1, 0);   // Moles
dimensionSet dimScalingFactor(0, -2, 2, 0, 0, 0, 0); // [M^0, L^-2, T^2, Theta^0, I^0, N^0, J^0]

// Create the mean molar mass field
volScalarField meanMolarMass(
    IOobject(
        "meanMolarMass",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0));

// Compute the mean molar mass
forAll(Y, i){
    meanMolarMass += Y[i]/thermo.composition().W(i);
}

// Invert the mean molar mass
meanMolarMass = 1/meanMolarMass;

tmp<fv::convectionScheme<scalar>> mvConvection(
    fv::convectionScheme<scalar>::New(
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

// Precompute total weighted body force
volVectorField totalWeightedForce(
    IOobject(
        "totalWeightedForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedVector("zero", dimLength/dimTime/dimTime, vector::zero)
);

forAll(Y, j){
    const volVectorField &bodyForce_Yj = mesh.lookupObject<volVectorField>("bodyForce_" + Y[j].name());
    totalWeightedForce += Y[j] * bodyForce_Yj;
}

// Compute the reaction rate and heat release
reaction->correct();
Qdot = reaction->Qdot();
volScalarField Yt(0.0 * Y[0]);

forAll(Y, i){
    if (i != inertIndex && composition.active(i)){

        volScalarField &Yi = Y[i];
        const word &speciesName = Y[i].name();

        // Initialize the summation term for species i
        const volVectorField &bodyForce_Yi = mesh.lookupObject<volVectorField>("bodyForce_" + Yi.name());
        volVectorField summationTerm(
            IOobject(
                "summationTerm",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh,
            dimensionedVector("zero", dimLength / dimTime / dimTime, vector::zero) // Initialize to zero
        );

        // Access predefined body force field as const reference
        const volVectorField &bodyForce_Yi = mesh.lookupObject<volVectorField>("bodyForce_" + Yi.name());

        // Compute the summation term for species i
        summationTerm = (scalar(1) - Yi) * bodyForce_Yi - (totalWeightedForce - Yi * bodyForce_Yi);

        // Initialize the scaling factor for species i
        volScalarField scalingFactor(
            IOobject(
                "scalingFactor_" + Yi.name(),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("scalingFactor", dimScalingFactor, 0.0) // Initialize with zero
        );
        scalingFactor = thermo.composition().W(i) * rho / (meanMolarMass * p);

        // Initialize the magnetic diffusion fkux for species i
        volVectorField magdiff_Yi(
            IOobject(
                "magdiff" + Yi.name(),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh,
            dimensionedVector("zero", dimless / dimLength, vector::zero) // Explicit dimensions as inverse length (1/m)
        );

        // Assign values to magdiff_Yi based on the scalingFactor, Yi, and summationTerm
        magdiff_Yi = scalingFactor * Yi * summationTerm;

        volVectorField &diffFlux_Yi = const_cast<volVectorField &>(mesh.lookupObject<volVectorField>("diffFlux_" + speciesName));

        // Compute the modified diffusion flux for species i
        diffFlux_Yi = -(((turbulence->alphaEff()* fvc::grad(Yi)) / rho) - ((turbulence->alphaEff() * magdiff_Yi) / rho));

        fvScalarMatrix YiEqn(
            fvm::ddt(rho, Yi) + mvConvection->fvmDiv(phi, Yi) - fvm::laplacian(turbulence->alphaEff(), Yi) + fvc::div(turbulence->alphaEff() * magdiff_Yi) ==
            reaction->R(Yi) + fvOptions(rho, Yi));

        YiEqn.relax();

        fvOptions.constrain(YiEqn);

        YiEqn.solve(mesh.solver("Yi"));

        fvOptions.correct(Yi);

        Yi.max(0.0);
        Yt += Yi;
    }
}

Y[inertIndex] = scalar(1) - Yt;
Y[inertIndex].max(0.0);
