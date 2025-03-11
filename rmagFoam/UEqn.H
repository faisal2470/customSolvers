// Solve the Momentum equation
MRF.correctBoundaryVelocity(U);
dimensionSet dimrho(1, -3, 0, 0, 0, 0, 0);
#include "fvcGrad.H"
#include "fvMatrix.H"

volVectorField bodyForceMom(
    IOobject(
        "bodyForceMom",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector("zero", dimForce / dimVolume / dimrho, vector::zero) // Initialize to zero
);

forAll(Y, i){
    volScalarField &Yi = Y[i];

    // Check if bodyForce_Yi is available before accessing it
    const word bodyForceName = "bodyForce_" + Yi.name();
    if (mesh.foundObject<volVectorField>(bodyForceName))
    {
        const volVectorField &bodyForce_Yi = mesh.lookupObject<volVectorField>(bodyForceName);

        // Accumulate the magnetic force for each species, scaled by its fraction Yi
        bodyForceMom += Yi * bodyForce_Yi;
    }
    else
    {
        Info << "Warning: " << bodyForceName << " not found in objectRegistry. Skipping." << endl;
    }
}

tmp<fvVectorMatrix> tUEqn(
    fvm::ddt(rho, U) + fvm::div(phi, U) + MRF.DDt(rho, U) + turbulence->divDevRhoReff(U) ==
    fvOptions(rho, U) + rho * bodyForceMom);
fvVectorMatrix &UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
    K = 0.5 * magSqr(U);
}
