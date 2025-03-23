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
const dimensionSet dimrho(1, -3, 0, 0, 0, 0, 0);                           // [kg m^-3]

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
const dimensionSet dimGradB2(2, -1, -4, 0, 0, -2, 0);

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
