DiscreteTimeArray imuAttitudeTrajectoryReconstruction
    (
    const DiscreteTimeArray & y,
    const DiscreteTimeArray & u,
    const Vector & xh0,
    const Matrix & p,
    const Matrix & q,
    const Matrix & r,
    double dt)
{
    ///Sizes of the states for the state, the measurement, and the input vector
    const unsigned stateSize=18;
    const unsigned measurementSize=6;
    const unsigned inputSize=6;

    ///initialization of the extended Kalman filter
    ExtendedKalmanFilter filter(stateSize, measurementSize, inputSize, false);

    ///initalization of the functor
    IMUDynamicalSystem imuFunctor;
    imuFunctor.setSamplingPeriod(dt);
    filter.setFunctor(& imuFunctor);

    ///the initalization of the estimation of the initial state
    filter.setState(xh0,y.getFirstTime()-1);

    ///computation and initialization of the covariance matrix of the initial state
    filter.setStateCovariance(p);

    ///set initial input
    filter.setInput(u[y.getFirstTime()-1],y.getFirstTime()-1);

    ///The covariance matrix of the process noise and the measurement noise
    /// for the extended Kalman filter
    filter.setR(r);
    filter.setQ(q);

    ///set the derivation step for the finite difference method
    Vector dx=filter.stateVectorConstant(1)*1e-8;

    ///the array of the state estimations over time
    DiscreteTimeArray xh;
    xh.setValue(xh0,y.getFirstTime()-1);

    ///the reconstruction of the state
    for (int i=y.getFirstTime();i<=y.getLastTime();++i)
    {
        ///introduction of the measurement
        filter.setMeasurement(y[i],i);

        ///introduction of the input
        if (i<y.getLastTime())
            filter.setInput(u[i],i);

        ///get the jacobians by finite differences and provide them to the Kalman filter
        Matrix a=filter.getAMatrixFD(dx);
        Matrix c= filter.getCMatrixFD(dx);
        filter.setA(a);
        filter.setC(c);

        ///get the estimation and give it to the array
        Vector xhk=filter.getEstimateState(i);

        ///regulate the part of orientation vector in the state vector
        xhk.segment(kine::ori,3)=kine::regulateOrientationVector(xhk.segment(9,3));

        ///give the new value of the state to the kalman filter.
        ///This step is usually unnecessary, unless we modify the value of the state esimation
        ///which is the case here.
        filter.setState(xhk,i);

        xh.setValue(xhk,i);
    }

    return xh;
}

DiscreteTimeArray imuAttitudeTrajectoryReconstruction(
    const DiscreteTimeArray & y,
    const Vector & xh0,
    const Matrix & p,
    const Matrix & q,
    const Matrix & r,
    double dt)
{
    const unsigned inputSize=6;

    ///initialization of a zero input
    DiscreteTimeArray u;
    for (int k=y.getFirstTime()-1; k<y.getLastTime(); ++k)
    {
        u.setValue(Vector::Zero(inputSize,1),k);
    }

    return imuAttitudeTrajectoryReconstruction (y, u, xh0, p, q, r, dt);
}




