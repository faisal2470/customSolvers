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
#include "magneticFields.H"

// Create the gravity field
Info << "Creating gravity field g\n" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);



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
        dimensionedScalar("zero", dimless / dimrho, 0.0)
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
