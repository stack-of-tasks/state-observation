#include <iostream>
#include <fstream>

#include <state-observation/sensors-simulation/accelerometer-gyrometer.hpp>
#include <state-observation/noise/gaussian-white-noise.hpp>
#include <state-observation/dynamical-system/imu-dynamical-system.hpp>
#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>
#include <state-observation/observer/extended-kalman-filter.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>

using namespace stateObservation;

int test()
{
    /// The number of samples
    const unsigned kmax=3000;

    ///sampling period
    const double dt=1e-3;

    ///Sizes of the states for the state, the measurement, and the input vector
    const unsigned stateSize=18;
    const unsigned measurementSize=6;
    const unsigned inputSize=6;

    ///The array containing all the states, the measurements and the inputs
    DiscreteTimeArray x;
    DiscreteTimeArray y;
    DiscreteTimeArray u;

    ///The covariance matrix of the process noise and the measurement noise
    Matrix q;
    Matrix r;

    {
        ///simulation of the signal
        /// the IMU dynamical system functor
        IMUDynamicalSystem imu;

        ///The process noise initialization
        Matrix q1=Matrix::Identity(stateSize,stateSize)*0.01;
        GaussianWhiteNoise processNoise(imu.getStateSize());
        processNoise.setStandardDeviation(q1);
        imu.setProcessNoise( & processNoise );
        q=q1*q1.transpose();

        ///The measurement noise initialization
        Matrix r1=Matrix::Identity(measurementSize,measurementSize)*0.01;
        GaussianWhiteNoise MeasurementNoise(imu.getMeasurementSize());
        MeasurementNoise.setStandardDeviation(r1);
        imu.setMeasurementNoise( & MeasurementNoise );
        r=r1*r1.transpose();

        ///the simulator initalization
        DynamicalSystemSimulator sim;
        sim.setDynamicsFunctor(&imu);

        ///initialization of the state vector
        Vector x0=Vector::Zero(stateSize,1);
        sim.setState(x0,0);

        ///construction of the input
        /// the input is constant over 10 time samples
        for (int i=0;i<kmax/10;++i)
        {
            Vector uk=Vector::Zero(imu.getInputSize(),1);

            uk[0]=0.4 * sin(M_PI/10*i);
            uk[1]=0.6 * sin(M_PI/12*i);
            uk[2]=0.2 * sin(M_PI/5*i);

            uk[3]=10  * sin(M_PI/12*i);
            uk[4]=0.07  * sin(M_PI/15*i);
            uk[5]=0.05 * sin(M_PI/5*i);

            ///filling the 10 time samples of the constant input
            for (int j=0;j<10;++j)
            {
                u.pushBack(uk,i*10+j);
            }

            ///give the input to the simulator
            ///we only need to give one value and the
            ///simulator takes automatically the appropriate value
            sim.setInput(uk,10*i);

        }

        ///set the sampling perdiod to the functor
        imu.setSamplingPeriod(dt);

        ///launched the simulation to the time kmax+1
        sim.simulateDynamicsTo(kmax+1);

        ///extract the array of measurements and states
        y = sim.getMeasurementArray(1,kmax);
        x = sim.getStateArray(1,kmax);
    }

    ///initialization of the extended Kalman filter
    ExtendedKalmanFilter filter(stateSize, measurementSize, measurementSize, false);

    ///initalization of the functor
    IMUDynamicalSystem imuFunctor;
    imuFunctor.setSamplingPeriod(dt);
    filter.setFunctor(& imuFunctor);

    ///the initalization of a random estimation of the initial state
    Vector xh0=filter.stateVectorRandom()*3.14;
    xh0[9]=3.14;
    filter.setState(xh0,0);

    ///computation and initialization of the covariance matrix of the initial state
    Matrix p=filter.getPmatrixZero();
    for (unsigned i=0;i<filter.getStateSize();++i)
    {
        p(i,i)=xh0[i];
    }
    p=p*p.transpose();
    filter.setStateCovariance(p);

    ///set initial input
    filter.setInput(u[0],0);

    ///set the covariance matrices for the extended Kalman filter
    filter.setR(r);
    filter.setQ(q);

    ///set the derivation step for the finite difference method
    Vector dx=filter.stateVectorConstant(1)*1e-8;

    ///the array of the state estimations over time
    DiscreteTimeArray xh;
    xh.pushBack(xh0,0);

    ///file of output
    std::ofstream f;
    f.open("trajectory.dat");

    ///the reconstruction of the state
    for (int i=y.getFirstTime();i<=y.getLastTime();++i)
    {
        ///introduction of the measurement
        filter.setMeasurement(y[i],i);
        if (i<x.getLastTime())
            filter.setInput(u[i],i);

        ///get the jacobians by finite differences and provide them to the Kalman filter
        Matrix a=filter.getAMatrixFD(dx);
        Matrix c= filter.getCMatrixFD(dx);
        filter.setA(a);
        filter.setC(c);

        ///get the estimation and give it to the array
        Vector xhk=filter.getEstimateState(i);

        ///regulate the part of orientation vector in the state vector
        xhk.segment(9,3)=tools::regulateOrientationVector(xhk.segment(9,3));

        ///give the new value of the state to the kalman filter.
        ///This step is usually unnecessary, unless we modify the value of the state esimation
        ///which is the case here.
        filter.setState(xhk,i);


        xh.pushBack(xhk,i);


        {///display part, useless
            Vector3 g;
            {
                Matrix3 R;
                Vector3 orientationV=Vector(x[i]).segment(9,3);
                double angle=orientationV.norm();
                if (angle > cst::epsilonAngle)
                    R = AngleAxis(angle, orientationV/angle).toRotationMatrix();
                else
                    R = Matrix3::Identity();
                g=R.transpose()*Vector3::UnitZ();
                g.normalize();
            }

            Vector3 gh;
            {
                Matrix3 Rh;
                Vector3 orientationV=Vector(xh[i]).segment(9,3);
                double angle=orientationV.norm();
                if (angle > cst::epsilonAngle)
                    Rh = AngleAxis(angle, orientationV/angle).toRotationMatrix();
                else
                    Rh = Matrix3::Identity();
                gh=Rh.transpose()*Vector3::UnitZ();
                gh.normalize();
            }


            f << i<< " \t "<< acos(double(g.transpose()*gh)) * 180 / M_PI << " \t\t\t "
            << g.transpose() << " \t\t\t " << gh.transpose() << std::endl;
        }
    }
    return 0;
}

int main()
{


    return test();

}
