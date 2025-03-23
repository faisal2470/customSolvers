{
    volScalarField& he = thermo.he();
    dimensionSet dimrho(1, -3, 0, 0, 0, 0, 0);
     volScalarField magneticWork(
      IOobject(
          "magneticWork",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE),
      mesh,
      dimensionedScalar("zero", dimEnergy / dimVolume / dimTime / dimrho, 0.0));

  forAll(Y, j)
  {
    volScalarField &Yj = Y[j];

    // Access predefined fields for magnetic susceptibility and magnetic force

    const volVectorField &magneticForce_Yj = mesh.lookupObject<volVectorField>("magneticForce_" + Yj.name());
    const volVectorField &diffflux_Yj = mesh.lookupObject<volVectorField>("diffflux_" + Yj.name());

    // Accumulate the magnetic work contribution: rho * Yi * (magneticForce_Yi · U)
    magneticWork += (magneticForce_Yj & diffflux_Yj);
  }
	
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
/*       + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? mvConvection->fvcDiv(fvc::absolute(phi, rho, U), p/rho)
          : -dpdt
        ) */
      + thermophysicalTransport->divq(he)
     ==
        reaction->Qdot()
      + fvModels.source(rho, he)
      + rho * magneticWork
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
