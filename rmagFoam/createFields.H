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


IOdictionary magneticProperties(
    IOobject(
        "magneticProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));

//Magnetic field block
Info << "Reading magnetic flux density B\n"
     << endl;

scalar B0 = magneticProperties.lookupOrDefault<scalar>("B0", 0.52);
scalar k1 = magneticProperties.lookupOrDefault<scalar>("k1", 120.0);// Default to 100 mm^-1

const dimensionSet dimMagneticField(1, 0, -2, 0, 0, -1, 0); // Tesla
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
Info << "Configured B0: " << B0 << " T" << endl;
Info << "Configured k1: " << k1 << " mm^-1" << endl;

forAll(B, cellI)
{
    const scalar y = mesh.C()[cellI].y(); // Get the y-coordinate of each cell's center

    if (y < 0)
    {
        B[cellI] = vector(0, B0, 0); // Set B to 0.52 if y < 0
    }
    else
    {
        B[cellI] = vector(0, B0 * Foam::exp(-k1 * y), 0); // Use the exponential expression if y >= 0
    }

}

forAll(B.boundaryField(), patchI)
{
    fvPatchField<vector>& patchField = const_cast<fvPatchField<vector>&>(B.boundaryField()[patchI]); // Cast away const for modification
    const pointField& faceCenters = patchField.patch().Cf(); // Get face centers of the patch

    forAll(patchField, faceI)
    {
        const scalar y = faceCenters[faceI].y(); // y-coordinate of the face center

        if (y < 0)
        {
            patchField[faceI] = vector(0, B0, 0);
        }
        else
        {
            patchField[faceI] = vector(0, B0 * Foam::exp(-k1 * y), 0);
        }
    }
}




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

forAll(gradB2, cellI)
{
    const scalar y = mesh.C()[cellI].y(); // Get y-coordinate of cell center

    if (y <= 0) // For y <= 0, ensure zero gradient
    {
        gradB2[cellI] = vector(0, 0, 0);
    }
    else
    {
        scalar By = B0 * Foam::exp(-k1 * y); // Example B_y = 0.52 * exp(-100 * y)
        scalar dBy_dy = -k1 * By;             // Derivative of B_y
        gradB2[cellI] = vector(0, 2.0 * By * dBy_dy, 0); // Analytical gradient
    }
}

forAll(gradB2.boundaryField(), patchI)
{
    fvPatchField<vector>& patchField = const_cast<fvPatchField<vector>&>(gradB2.boundaryField()[patchI]); // Cast away const for modification
    const pointField& faceCenters = patchField.patch().Cf(); // Get face centers of the patch

    forAll(patchField, faceI)
    {
        const scalar y = faceCenters[faceI].y(); // y-coordinate of the face center

        if (y < 0)
        {
            patchField[faceI] = vector(0, 0, 0);
        }
        else
        {
            scalar By = B0 * Foam::exp(-k1 * y); // Example B_y = 0.52 * exp(-100 * y)
            scalar dBy_dy = -k1 * By;           // Derivative of B_y
            patchField[faceI] = vector(0, 2.0 * By * dBy_dy, 0);
        }
    }
}



const dimensionSet dimInvMole(0, 0, 0, 0, -1, 0, 0); // Dimensions for mol^-1

// Avogadro's number (NA) [1/mol]
dimensionedScalar NA("NA", dimInvMole, 6.022e23);

// Lande's g-factor (dimensionless, so no units are required)
dimensionedScalar gL("gL", dimless, 2.0); // dimensionless

// Bohr magneton (muB) [J/T] = [M L^2 T^-2 A^-1]
const dimensionSet dimBohrMagneton = dimEnergy / (dimMagneticField); // [M L^2 T^-2 A^-1]
dimensionedScalar muB("muB", dimBohrMagneton, 9.274e-24);

// Boltzmann constant (k) [J/K]
dimensionedScalar k("k", dimEnergy / dimTemperature, 1.38e-23);

// Magnetic permeability of free space (mu0) [H/m] = [M L T^-2 A^-2]
const dimensionSet dimPermeability = dimForce / (dimCurrent * dimCurrent);
dimensionedScalar mu0("mu0", dimPermeability, 4 * M_PI * 1e-7);

// Loop over species and create individual chi and magneticForce fields

const dimensionSet dimMolarMass(1, 0, 0, 0, 0, -1, 0); // kg/mol
const dimensionSet dimVolumePerMass( -1, 3, 0, 0, 0, 0, 0);


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

    // Create and register magneticForce field without assigning to a variable
    new volVectorField(
        IOobject(
            "magneticForce_" + speciesName,
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

    // Store magneticForce in the object registry

}
