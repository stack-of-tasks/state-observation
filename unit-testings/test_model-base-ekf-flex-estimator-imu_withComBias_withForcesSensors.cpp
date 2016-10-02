#include <iostream>
#include <fstream>

#include <state-observation/flexibility-estimation/model-base-ekf-flex-estimator-imu.hpp>
#include <state-observation/flexibility-estimation/imu-elastic-local-frame-dynamical-system.hpp>

#include <time.h>

using namespace stateObservation;

Vector computeZmp (unsigned footNumber, IndexedMatrixArray& forces, IndexedMatrixArray& sensorPositions)
{
    double fnormal = 0;
    double sumZmpx = 0;
    double sumZmpy = 0;
    double sumZmpz = 0;
    Vector zmp; zmp.setZero(); zmp.resize (3);

    for (unsigned int i=0; i<footNumber; ++i)
    {
        const Vector& f = forces [i];
        // Check that force is of dimension 6
        if (f.size () != 6)
        {
            zmp.fill (0.);
            return zmp;
        }

        const Matrix& M = sensorPositions[i];
        double fx = M (0,0) * f (0) + M(0,1) * f (1) + M (0,2) * f (2);
        double fy = M (1,0) * f (0) + M(1,1) * f (1) + M (1,2) * f (2);
        double fz = M (2,0) * f (0) + M(2,1) * f (1) + M (2,2) * f (2);

        if (fz > 0)
        {
            double Mx = M (0,0)*f(3) + M (0,1)*f(4) + M (0,2)*f(5) + M (1,3)*fz - M (2,3)*fy;
            double My = M (1,0)*f(3) + M (1,1)*f(4) + M (1,2)*f(5) + M (2,3)*fx - M (0,3)*fz;

            fnormal += fz;
            sumZmpx -= My;
            sumZmpy += Mx;
            sumZmpz += fz * M (2,3);
        }
    }

    if (fnormal != 0)
    {
        zmp (0) = sumZmpx / fnormal;
        zmp (1) = sumZmpy / fnormal;
        zmp (2) = sumZmpz / fnormal;
    }
    else
    {
        zmp.fill (0.);
    }

    return zmp;
}

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
    const unsigned kinit=3500; // 10000; //
    const unsigned kmax=5700; // 11500; //

    // Fix sizes
    const unsigned measurementSizeBase=6;
    const unsigned inputSizeBase=42;
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

     // Zmp from sensors
     IndexedMatrixArray sensorZmp_output;
     // Zmp from filter
     IndexedMatrixArray filteredZmp_output;

    std::cout << "Creating estimator" <<std::endl;

    stateObservation::flexibilityEstimation::ModelBaseEKFFlexEstimatorIMU est;
    est.setSamplingPeriod(dt);

    // Model
    est.setContactModel(1);
    est.setRobotMass(56.8);
    est.setKfe(40000*Matrix3::Identity());
    est.setKte(350*Matrix3::Identity());
    est.setKfv(600*Matrix3::Identity());
    est.setKtv(10*Matrix3::Identity());

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
    est.setForceVariance(1.e-6);
    est.setAbsolutePosVariance(1e-4);

    // Process noise covariance
    stateObservation::Matrix Q_; Q_.resize(stateSize,stateSize); Q_.setIdentity();
    Q_.block(0,0,12,12)*=1.e-8;
    Q_.block(12,12,12,12)*=1.e-4;
    Q_.block(24,24,6,6)*=1.e-2;
    Q_.block(30,30,2,2)*=1.e-11;
    Q_.block(32,32,3,3)*=1.e-8;
    est.setProcessNoiseCovariance(Q_);

    // Temporary variables
    stateObservation::Vector input, measurement;
    Vector x; x.resize(stateSize);
    unsigned contactNbr, inputSize, measurementSize;

    stateObservation::IndexedMatrixArray measurementForces, filteredForces, inputFeetPositions;
    stateObservation::Vector sensorZmp, filteredZmp;

    std::cout << "Beginning reconstruction "<<std::endl;

    for (unsigned k=kinit;k<kmax;++k)
    {
        std::cout << "\n" << k << std::endl;

        if(nbSupport[k](0)!=est.getContactsNumber())
        {
            contactNbr = nbSupport[k](0);
            est.setContactsNumber(contactNbr);
        }

        measurementSize = est.getMeasurementSize();
        measurement.resize(measurementSize);
        measurement = (y[k].block(0,0,1,measurementSize)).transpose();
        est.setMeasurement(measurement);

        inputSize = est.getInputSize();
        input.resize(inputSize);
        input = (u[k+1].block(0,0,1,inputSize)).transpose();
        est.setMeasurementInput(input);

        x = est.getFlexibilityVector();

        // Compute Zmp from sensor and filter
        for(unsigned i=0;i<contactNbr;++i)
        {
            measurementForces.setValue(y[k].block<1,6>(0,measurementSizeBase+withUnmodeledForces_*6+i*6).transpose(),i);
            filteredForces.setValue(x.segment<6>(stateObservation::flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::state::fc+i*6),i);
            inputFeetPositions.setValue(kine::vector6ToHomogeneousMatrix(u[k].block<1,6>(0,inputSizeBase+i*12)).transpose(),i);
        }
        sensorZmp=computeZmp(contactNbr, measurementForces, inputFeetPositions);
        filteredZmp=computeZmp(contactNbr, filteredForces, inputFeetPositions);

        x_output.setValue(x,k);
        y_output.setValue(y[k],k);
        u_output.setValue(u[k],k);
        sensorZmp_output.setValue(sensorZmp,k);
        filteredZmp_output.setValue(filteredZmp,k);
    }

    std::cout << "Completed "<<std::endl;

    x_output.writeInFile("state.dat");
    y_output.writeInFile("measurement.dat");
    u_output.writeInFile("input.dat");
    sensorZmp_output.writeInFile("sensorZmp.dat");
    filteredZmp_output.writeInFile("filteredZmp.dat");

    return 1;
}

int main()
{
    return test();
}

