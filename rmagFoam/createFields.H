#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}


volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& p = thermo.p();

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> reaction
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
);

#include "createMRF.H"
#include "createFvOptions.H"

//Magnetic field block
Info << "Reading magnetic Properties\n" << endl;

IOdictionary magneticProperties(
    IOobject(
        "magneticProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

scalar B0_value = magneticProperties.lookupOrDefault<scalar>("B0", 0.52);                   // Default to 0.52 T
scalar k1_value = magneticProperties.lookupOrDefault<scalar>("k1", 120.0);                  // Default to 120 m^-1
vector gradDir = magneticProperties.lookupOrDefault<vector>("gradDir", vector(0, 1, 0));    // Default to (0, 1, 0)

// Assign dimensions [M, L, T, K, N, A, J]
dimensionedScalar B0("B0", dimensionSet(1, 0, -2, 0, 0, -1, 0), B0_value); // [T]
dimensionedScalar k1("k1", dimensionSet(0, -1, 0, 0, 0, 0, 0), k1_value);  // [m^-1]

// Normalize the gradient direction
scalar gradMag = mag(gradDir);
if (gradMag < SMALL) gradMag = 1.0;     // Avoid division by zero
gradDir /= gradMag;

Info << "Configured B0: " << B0 << " T" << endl;
Info << "Configured k1: " << k1 << " mm^-1" << endl;
Info << "Configured gradDir: " << gradDir << endl;

// Compute projection r = mesh.C() & gradDir
volScalarField r(
    IOobject(
        "r",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh.C() & gradDir
);

// Create the magnetic field B and initialize it to zero
Info << "\nCreating magnetic field B" << endl;
const dimensionSet dimMagneticField(1, 0, -2, 0, 0, -1, 0);     // Tesla
volVectorField B(
    IOobject(
        "B",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector("zero", dimMagneticField, vector::zero) // Initialize with zero values
);

// Compute B internal field
// Info << "Computing internal magnetic field\n" << endl;
B = B0 * Foam::exp(-k1 * r) * gradDir;

// Compute B boundary field
// Info << "Computing boundary magnetic field\n" << endl;
forAll(B.boundaryField(), patchI)
{
    scalarField rPatch = B.boundaryField()[patchI].patch().Cf() & gradDir; // Project the face centers onto the gradient direction
    B.boundaryFieldRef()[patchI] = B0.value() * Foam::exp(-k1.value() * rPatch) * gradDir; // Compute the magnetic field on the patch

}

// Create the gradient of B^2 field
Info << "Creating gradient of B^2 field\n" << endl;
const dimensionSet dimGradB2(2, -1, -4, 0, -2, 0, 0);

volVectorField gradB2(
    IOobject(
        "gradB2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimGradB2, vector::zero) // Set all components to zero
);

// Compute the gradient of B^2 internal field
// Info << "Computing internal gradient of B^2 field\n" << endl;
gradB2 = -2.0 * k1 * B0 * B0 * Foam::exp(-2.0 * k1 * r) * gradDir;

// Compute the gradient of B^2 boundary field
// Info << "Computing boundary gradient of B^2 field\n" << endl;
forAll(gradB2.boundaryField(), patchI)
{
    scalarField rPatch = gradB2.boundaryField()[patchI].patch().Cf() & gradDir;                         // Project the face centers onto the gradient direction
    gradB2.boundaryFieldRef()[patchI] = -2.0 * (k1 * B0 * B0).value() * Foam::exp(-2.0 * k1.value() * rPatch) * gradDir;  // Compute the gradient of B^2 on the patch

}


// Create dimensioned scalars for physical constants

// Avogadro's number (NA) [1/mol]
const dimensionSet dimInvMole(0, 0, 0, 0, -1, 0, 0);                        // Dimensions for mol^-1
dimensionedScalar NA("NA", dimInvMole, 6.022e23);

// Lande's g-factor (dimensionless, so no units are required)
dimensionedScalar gL("gL", dimless, 2.0); // dimensionless

// Bohr magneton (muB) [J/T]
const dimensionSet dimBohrMagneton = dimEnergy / (dimMagneticField);        // [M L^2 T^-2 A^-1]
dimensionedScalar muB("muB", dimBohrMagneton, 9.274e-24);

// Boltzmann constant (k) [J/K]
dimensionedScalar k("k", dimEnergy / dimTemperature, 1.38e-23);

// Magnetic permeability of free space (mu0) [H/m]
const dimensionSet dimPermeability = dimForce / (dimCurrent * dimCurrent);  // [M L T^-2 A^-2]
dimensionedScalar mu0("mu0", dimPermeability, 4 * M_PI * 1e-7);

// Loop over species and create individual chi and bodyForce fields
forAll(Y, i)
{
    const word &speciesName = Y[i].name();

    // Retrieve the sub-dictionary for the specific species

        new volScalarField(
        IOobject(
            "chi_" + speciesName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            true // registerObject = true
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    // Create and register bodyForce field without assigning to a variable
    new volVectorField(
        IOobject(
            "bodyForce_" + speciesName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            true // registerObject = true
        ),
        mesh,
        dimensionedVector("zero", dimAcceleration, vector::zero)
    );

        new volVectorField(
        IOobject(
            "diffFlux_" + speciesName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            true // registerObject = true
        ),
        mesh,
        dimensionedVector("zero", dimVelocity, vector::zero)
    );

}
