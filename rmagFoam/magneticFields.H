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

// Read magnetic properties
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

int nFields = magneticProperties.lookupOrDefault<int>("nFields", 1);                        // Default to 1 if not found
Info << "Number of fields: " << nFields << endl;

// Define a small dimensioned offset for >= and <= conditions
dimensionedScalar smallOffset("smallOffset", dimensionSet(0, 1, 0, 0, 0, 0, 0), 1e-10);     // [m], adjust if needed

for (int fieldI = 1; fieldI <= nFields; fieldI++)
{
    word fieldName = "field" + Foam::name(fieldI);
    if (!magneticProperties.isDict(fieldName))
    {
        FatalErrorIn("magneticFields.H")
            << "magneticProperties dictionary specifies nFields = " << nFields
            << ", but sub-dictionary " << fieldName << " is not defined!" << endl
            << "Please define all fields from field1 to field" << nFields << "." << exit(FatalError);
    }

    Info << "\nDefining " << fieldI << endl;
    const dictionary &fieldDict = magneticProperties.subDict(fieldName);

    word fieldType = fieldDict.lookupOrDefault<word>("type", "exp"); // Default to exp if not found
    vector gradDir = fieldDict.lookupOrDefault<vector>("gradDir", vector(0, 1, 0)); // Default to (0, 1, 0)
    dimensionedScalar r0("r0", dimensionSet(0, 1, 0, 0, 0, 0, 0), fieldDict.lookupOrDefault<scalar>("r0", 0.0));
    dimensionedScalar rf("rf", dimensionSet(0, 1, 0, 0, 0, 0, 0), fieldDict.lookupOrDefault<scalar>("rf", 0.01));
    bool inv = fieldDict.lookupOrDefault<int>("inv", 0) == 1;

    if (r0.value() > rf.value())
    {
        FatalErrorIn("magneticFields.H")
            << "Field " << fieldI << " has r0 = " << r0 << " > rf = " << rf << endl
            << "Please ensure r0 <= rf for all fields." << exit(FatalError);
    }

    scalar gradMag = mag(gradDir);
    if (gradMag < SMALL) gradMag = 1.0; // Avoid division by zero
    gradDir /= gradMag;

    Info << "\tField type: " << fieldType << endl;
    Info << "\tGradient direction: " << gradDir << endl;
    Info << "\tInitial position: " << r0.value() << " m" << endl;
    Info << "\tFinal position: " << rf.value() << " m" << endl;
    Info << "\tInverse: " << inv << endl;

    dimensionedScalar B0("B0", dimensionSet(1, 0, -2, 0, 0, -1, 0), 0.52);
    dimensionedScalar k("k", dimensionSet(0, -1, 0, 0, 0, 0, 0), 120.0);

    if (fieldType == "exp")
    {
        const dictionary &expDict = fieldDict.subDict("exp");
        B0 = dimensionedScalar("B0", dimMagneticField, expDict.lookupOrDefault<scalar>("B0", 0.52));
        k = dimensionedScalar("k", dimensionSet(0, -1, 0, 0, 0, 0, 0), expDict.lookupOrDefault<scalar>("k", 120.0));
        Info << "\tB0: " << B0.value() << " T" << endl;
        Info << "\tk: " << k.value() << " m^-1" << endl;

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

        // Compute B contribution for this field
        volVectorField B_field = (
            inv
            ? pos(r - r0 + smallOffset) * pos(rf - r + smallOffset) * B0 * Foam::exp(-k * (rf - r))
            : pos(r - r0 + smallOffset) * pos(rf - r + smallOffset) * B0 * Foam::exp(-k * (r - r0))
        ) * gradDir;
        B += B_field;

        // Compute gradB2 contribution for this field
        volVectorField gradB2_field = (
            inv
            ? 2.0 * pos(r - r0 + smallOffset) * pos(rf - r + smallOffset) * k * B0 * B0 * Foam::exp(-2.0 * k * (rf - r))
            : -2.0 * pos(r - r0 + smallOffset) * pos(rf - r + smallOffset) * k * B0 * B0 * Foam::exp(-2.0 * k * (r - r0))
        ) * gradDir;
        gradB2 += gradB2_field;

    }
    else
    {
        FatalErrorIn("magneticFields") << "Invalid type: " << fieldType << endl;
        exit(FatalError);
    }
}

