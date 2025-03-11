
volScalarField &he = thermo.he();

volScalarField bodyForceWork(
    IOobject(
        "bodyForceWork",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimEnergy / dimVolume / dimTime / dimrho, 0.0));

forAll(Y, j)
{
  volScalarField &Yj = Y[j];

  // Access predefined fields for magnetic susceptibility and body force

  const volVectorField &bodyForce_Yj = mesh.lookupObject<volVectorField>("bodyForce_" + Yj.name());
  const volVectorField &diffFlux_Yj = mesh.lookupObject<volVectorField>("diffFlux_" + Yj.name());

  // Accumulate the magnetic work contribution: rho * Yi * (bodyForce_Yi · U)
  bodyForceWork += (bodyForce_Yj & diffFlux_Yj);
}

fvScalarMatrix EEqn(
    fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he) + fvc::ddt(rho, K) + fvc::div(phi, K) + (he.name() == "e" ? fvc::div(fvc::absolute(phi / fvc::interpolate(rho), U), p, "div(phiv,p)") : -dpdt) - fvm::laplacian(turbulence->alphaEff(), he) ==
    Qdot + fvOptions(rho, he) + rho * bodyForceWork);

EEqn.relax();

fvOptions.constrain(EEqn);

EEqn.solve();

fvOptions.correct(he);

thermo.correct();

Info << "min/max(T) = "
      << min(T).value() << ", " << max(T).value() << endl;
