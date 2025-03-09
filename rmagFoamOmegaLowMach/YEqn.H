#include "dimensionedTypes.H"
const dimensionSet dimAmount(0, 0, 0, 0, 0, 1, 0);   // Moles
dimensionSet dimScalingFactor(0, -2, 2, 0, 0, 0, 0); // [M^0, L^-2, T^2, Theta^0, I^0, N^0, J^0]

tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

volScalarField meanMolarMass(
    IOobject(
        "meanMolarMass",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimMass/dimAmount, 0.0));
   
volScalarField temp(
    IOobject(
        "temp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimAmount/dimMass, 0.0));
    
forAll(Y, i)
{   const word &speciesName = Y[i].name();
    
    const dictionary& speciesDict = magneticProperties.subDict("speciesProperties");
    const dictionary& speciesSubDict = speciesDict.subDict(speciesName);

    dimensionedScalar mi(
        "mi_" + speciesName, dimMolarMass,
        speciesSubDict.lookup<scalar>("mi")
    );

    temp += Y[i]/ mi;
}

meanMolarMass = 1/temp;

reaction->correct();

forAll(Y, i)
{
    if (composition.solve(i))
    {
        volScalarField& Yi = Y[i];
        
        
        const word &speciesName = Y[i].name();
            // Access predefined magnetic force field as const reference
            const volVectorField &magneticForce_Yi = mesh.lookupObject<volVectorField>("magneticForce_" + Yi.name());
            const dictionary& speciesDict = magneticProperties.subDict("speciesProperties");
            const dictionary& speciesSubDict = speciesDict.subDict(speciesName);
        	dimensionedScalar mi(
        	"mi_" + speciesName, dimMolarMass,
        	speciesSubDict.lookup<scalar>("mi") );
            
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
            scalingFactor = mi * rho / (meanMolarMass * p);
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

            volVectorField &diffflux_Yi = const_cast<volVectorField &>(mesh.lookupObject<volVectorField>("diffflux_" + speciesName));

            diffflux_Yi = -(((thermophysicalTransport->DEff(Yi)* fvc::grad(Yi)) / rho) - ((thermophysicalTransport->DEff(Yi) * magdiff_Yi) / rho));
               
         

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho, Yi)
          + mvConvection->fvmDiv(phi, Yi)
          + thermophysicalTransport->divj(Yi)
          + fvc::div(thermophysicalTransport->DEff(Yi) * magdiff_Yi)
         ==
            reaction->R(Yi)
          + fvModels.source(rho, Yi)
        );

        YiEqn.relax();

        fvConstraints.constrain(YiEqn);

        YiEqn.solve("Yi");

        fvConstraints.constrain(Yi);
    }
}

composition.normalise();

forAll(Y, i)
{
    ROP_[i].field() = -(reaction->R(Y[i])()).source()/mesh.V();
}
