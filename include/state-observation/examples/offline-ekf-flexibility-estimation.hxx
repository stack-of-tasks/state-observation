stateObservation::DiscreteTimeArray offlineEKFFlexibilityEstimation(
            const stateObservation::DiscreteTimeArray & y,
            const stateObservation::DiscreteTimeArray & u,
            const Matrix & xh0,
            double dt)
{

        ///Sizes of the states for the state, the measurement, and the input vector
    const unsigned stateSize=18;
    const unsigned measurementSize=6;
    const unsigned inputSize=15;

    flexibilityEstimation::FixedContactEKFFlexEstimatorIMU estimator;

    ///
    estimator.setFlexibilityGuess(xh0);

    ///the array of the state estimations over time
    stateObservation::DiscreteTimeArray xh;
    xh.pushBack(xh0,y.getFirstTime()-1);

    ///the reconstruction of the state
    for (int i=y.getFirstTime();i<=y.getLastTime();++i)
    {
        std::cout << i << std::endl;

        ///introduction of the measurement
        estimator.setMeasurement(y[i]);

        estimator.setMeasurementInput(u[i-1]);

        ///get the estimation and give it to the array
        xh.pushBack(estimator.getFlexibility());

    }

    return xh;
}


stateObservation::DiscreteTimeArray offlineEKFFlexibilityEstimation(
            const stateObservation::DiscreteTimeArray & y,
            const Matrix & xh0,
            double dt)
{
    const unsigned inputSize=15;

    ///initialization of a zero input
    stateObservation::DiscreteTimeArray u;
    for (int k=y.getFirstTime()-1; k<=y.getLastTime(); ++k)
    {
        u.pushBack(Vector::Zero(inputSize,1),k);
    }

    return offlineEKFFlexibilityEstimation (y, u, xh0, dt);
}

