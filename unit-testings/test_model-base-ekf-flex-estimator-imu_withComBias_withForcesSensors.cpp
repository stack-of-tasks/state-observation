#include <iostream>
#include <fstream>

//#include <state-observation/noise/gaussian-white-noise.hpp>
//#include <state-observation/examples/offline-ekf-flexibility-estimation.hpp>
//#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>
//#include <state-observation/tools/miscellaneous-algorithms.hpp>

#include <state-observation/flexibility-estimation/model-base-ekf-flex-estimator-imu.hpp>

#include <time.h>
#include <iostream>


using namespace stateObservation;

int test()
{
    std::cout << "Starting" << std::endl;

    // For measurement vector
    bool withUnmodeledForces_ = false;
    bool withForceSensors_=true;
    bool withAbsolutePose_ = false;

    // For state vector
    bool withComBias_=false;

    // Time
    const double dt=5e-3;
    const unsigned kinit=3;
    const unsigned kmax=2144;

    // Fix sizes
    const unsigned measurementSizeBase=6;
    const unsigned stateSize=35;

    /// Definitions of input vectors
     // Measurement
     IndexedMatrixArray y;
     std::cout << "Loading measurements file" << std::endl;
     y.getFromFile("source_measurement.dat",1,30);

     // Input
     IndexedMatrixArray u;
     std::cout << "Loading input file" << std::endl;
     u.getFromFile("source_input.dat",1,66);

     // Number of support contacts
     IndexedMatrixArray nbSupport;
     std::cout << "Loading the number of supports file" << std::endl;
     nbSupport.getFromFile("source_nbSupport.dat",1,1);

    /// Definition of ouptut vectors 
     // State: what we want
     IndexedMatrixArray x_output;
     // Measurement
     IndexedMatrixArray y_output;
     // Input
     IndexedMatrixArray u_output;

    std::cout << "Creating estimator" <<std::endl;

    stateObservation::flexibilityEstimation::ModelBaseEKFFlexEstimatorIMU est;
    est.setSamplingPeriod(dt);

    // Model
    est.setContactModel(1);
    est.setRobotMass(56.8);
    est.setKfe(40000*Matrix3::Identity());
    est.setKte(600*Matrix3::Identity());
    est.setKfv(600*Matrix3::Identity());
    est.setKtv(60*Matrix3::Identity());

    // Config
    est.setWithUnmodeledMeasurements(withUnmodeledForces_);
    est.setWithForcesMeasurements(withForceSensors_);
    est.setWithAbsolutePos(withAbsolutePose_);
    est.setWithComBias(withComBias_);

    std::cout << "Setting covariances" << std::endl;

    // Measurement noise covariance
    stateObservation::Matrix R_; R_.resize(measurementSizeBase,measurementSizeBase); R_.setIdentity();
    R_.block(0,0,3,3)*=1.e-3;
    R_.block(3,3,3,3)*=1.e-6;
    est.setMeasurementNoiseCovariance(R_);
    est.setUnmodeledForceVariance(1e-13);
    est.setForceVariance(1.e-8);
    est.setAbsolutePosVariance(1e-4);

    // Process noise covariance
    stateObservation::Matrix Q_; Q_.resize(stateSize,stateSize); Q_.setIdentity();
    Q_.block(0,0,12,12)*=1.e-8;
    Q_.block(12,12,12,12)*=1.e-4;
    Q_.block(24,24,6,6)*=1.e-2;
    Q_.block(30,30,2,2)*=1.e-15;
    Q_.block(32,32,3,3)*=1.e-8;
    est.setProcessNoiseCovariance(Q_);

//    stateObservation::Vector inputInit, measurementInit;
//    measurementInit.resize(est.getMeasurementSize()); measurementInit.setZero();
//    inputInit.resize(est.getInputSize()); inputInit.setZero();
//    for(unsigned i=0;i<=kinit;++i)
//    {
//        est.setMeasurementInput(inputInit);
//        est.setMeasurement(measurementInit);
//        flexibility = est.getFlexibilityVector();
//    }

    // Temporary variables
    stateObservation::Vector input, measurement;
    Vector x; x.resize(stateSize);
    unsigned contactNbr, inputSize, measurementSize;

    std::cout << "Beginning reconstruction "<<std::endl;

    for (unsigned k=kinit;k<kmax;++k)
    {
        std::cout << k << std::endl;

        if(nbSupport[k](0)!=est.getContactsNumber())
        {
            contactNbr = nbSupport[k](0);
            est.setContactsNumber(contactNbr);
        }

        measurementSize = est.getMeasurementSize();
        measurement.resize(measurementSize); measurement.setZero();
        measurement = (y[k].block(0,0,1,measurementSize)).transpose();
        est.setMeasurement(measurement);

        inputSize = est.getInputSize();
        input.resize(inputSize);
        input=(u[k].block(0,0,1,inputSize)).transpose();
        est.setMeasurementInput(input);

        x = est.getFlexibilityVector();

        x_output.setValue(x,k);
        y_output.setValue(y[k],k);
        u_output.setValue(u[k],k);
    }

    std::cout << "Completed "<<std::endl;

    x_output.writeInFile("state.dat");
    y_output.writeInFile("measurement.dat");
    u_output.writeInFile("input.dat");

    return 1;
}

int main()
{
    return test();
}

