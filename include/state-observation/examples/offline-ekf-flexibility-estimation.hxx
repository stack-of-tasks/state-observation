stateObservation::DiscreteTimeArray offlineEKFFlexibilityEstimation(
            const stateObservation::DiscreteTimeArray & y,
            const stateObservation::DiscreteTimeArray & u,
            const Matrix & xh0,
            unsigned numberOfContacts,
            const std::vector<Vector3> & contactsPositions,
            double dt,
            DiscreteTimeArray * ino,
            DiscreteTimeArray * premea)
{

        ///Sizes of the states for the state, the measurement, and the input vector
    const unsigned stateSize=18;
    const unsigned measurementSize=6;
    const unsigned inputSize=15;

    flexibilityEstimation::FixedContactEKFFlexEstimatorIMU estimator;


    Matrix R(estimator.getEKF().getRmatrixIdentity());
    R=R*1.e-16;

    Matrix Q(estimator.getEKF().getQmatrixIdentity());
    Q=Q*1.e-4;
    Q.block(6,6,3,3)=Matrix3::Identity()*1.e-2;
    Q.block(15,15,3,3)=Matrix3::Identity()*1.e-2;


    estimator.setFlexibilityCovariance(Q);
    estimator.setMeasurementNoiseCovariance(R);
    estimator.setProcessNoiseCovariance(Q);


    estimator.setContactsNumber(numberOfContacts);

    for (unsigned i = 0; i<numberOfContacts ; ++i)
    {
        estimator.setContactPosition(i,contactsPositions[i]);
    }

    ///
    estimator.setFlexibilityGuess(xh0);

    ///the array of the state estimations over time
    stateObservation::DiscreteTimeArray xh;
    xh.setValue(xh0,y.getFirstTime()-1);

    ///the reconstruction of the state
    for (int i=y.getFirstTime();i<=y.getLastTime();++i)
    {
        std::cout << i << std::endl;

        ///introduction of the measurement
        estimator.setMeasurement(y[i]);

        estimator.setMeasurementInput(u[i]);

        ///get the estimation and give it to the array
        xh.pushBack(estimator.getFlexibilityVector());

        if (ino != 0)
        {
            ino->setValue(estimator.getInovation(),i);
        }

        if (premea != 0)
        {
            premea->setValue(estimator.getPredictedMeaurement(),i);
        }

    }

    return xh;
}


stateObservation::DiscreteTimeArray offlineEKFFlexibilityEstimation(
            const stateObservation::DiscreteTimeArray & y,
            const Matrix & xh0,
            unsigned numberOfContacts,
            const std::vector<Vector3> & contactsPositions,
            double dt)
{
    const unsigned inputSize=15;

    ///initialization of a zero input
    stateObservation::DiscreteTimeArray u;
    for (int k=y.getFirstTime()-1; k<=y.getLastTime(); ++k)
    {
        u.setValue(Vector::Zero(inputSize,1),k);
    }

    return offlineEKFFlexibilityEstimation
                        (y, u, xh0, numberOfContacts, contactsPositions, dt);
}

