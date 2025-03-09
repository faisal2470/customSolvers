#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<fluidReactionThermo> pThermo(fluidReactionThermo::New(mesh));
fluidReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

dimensionedScalar pRef("pRef", dimPressure, 1.01325e5);
Info << "Reference pressure pRef: " << pRef << nl;

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    pRef/((thermo.Cp() - thermo.Cv())*thermo.T())
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

pressureReference pressureReference(p, pimple.dict(), false);

mesh.setFluxRequired(p.name());


Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidReactionThermophysicalTransportModel>
thermophysicalTransport
(
    fluidReactionThermophysicalTransportModel::New
    (
        turbulence(),
        thermo
    )
);

Info<< "Creating reaction model\n" << endl;
autoPtr<combustionModel> reaction(combustionModel::New(thermo, turbulence()));


/* Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
); */

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

dimensionedScalar initialMass = fvc::domainIntegrate(rho);


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

// forAll(Y, i)
// {
//     ROP_[i].field() = -(reaction->R(Y[i])()).source()/mesh.V();
// }

PtrList<volScalarField> ROP_(Y.size());
forAll(Y, i)
{
    const word RRname = "ROP." + Y[i].name();
    ROP_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                RRname,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );    
}
    
#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
