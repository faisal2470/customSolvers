#include "dimensionedTypes.H"
const dimensionSet dimAmount(0, 0, 0, 0, 0, 1, 0);   // Moles
dimensionSet dimScalingFactor(0, -2, 2, 0, 0, 0, 0); // [M^0, L^-2, T^2, Theta^0, I^0, N^0, J^0]
volScalarField meanMolarMass(
    IOobject(
        "meanMolarMass",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0));

 forAll(Y, i)
    {
       meanMolarMass += Y[i]/thermo.composition().W(i);
    }

    meanMolarMass = 1/meanMolarMass;

tmp<fv::convectionScheme<scalar>> mvConvection(
    fv::convectionScheme<scalar>::New(
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")));

{
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0 * Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField &Yi = Y[i];

            const word &speciesName = Y[i].name();
            // Access predefined magnetic force field as const reference
            const volVectorField &magneticForce_Yi = mesh.lookupObject<volVectorField>("magneticForce_" + Yi.name());
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

            forAll(Y, j)
            {
                volScalarField &Yj = Y[j];

                if (i != j)
                {
                    // Use const references for chi_Yj and magneticForce_Yj, and remove chi_Yj if unused
                    const volVectorField &magneticForce_Yj = mesh.lookupObject<volVectorField>("magneticForce_" + Y[j].name());

                    // Add Y_i * Y_j * (f_i - f_j) to the summation term
                    summationTerm += Yj * (magneticForce_Yi - magneticForce_Yj);
                }
            }

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
            // Compute the modified diffusion velocity for species i
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
}
