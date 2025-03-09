// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);
dimensionSet dimrho(1, -3, 0, 0, 0, 0, 0);
#include "fvcGrad.H"
#include "fvMatrix.H"

volVectorField magneticForceMom(
    IOobject(
        "magneticForceMom",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector("zero", dimForce / dimVolume / dimrho, vector::zero) // Initialize to zero
);

forAll(Y, i)
{
    volScalarField &Yi = Y[i];

    // Check if magneticForce_Yi is available before accessing it
    const word magneticForceName = "magneticForce_" + Yi.name();
    if (mesh.foundObject<volVectorField>(magneticForceName))
    {
        const volVectorField &magneticForce_Yi = mesh.lookupObject<volVectorField>(magneticForceName);

        // Accumulate the magnetic force for each species, scaled by its fraction Yi
        magneticForceMom += Yi * magneticForce_Yi;
    }
    else
    {
        Info << "Warning: " << magneticForceName << " not found in objectRegistry. Skipping." << endl;
    }
}


tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevTau(U)
 ==
    fvModels.source(rho, U) + rho * magneticForceMom
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvConstraints.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvConstraints.constrain(U);
    K = 0.5*magSqr(U);
}
