#include <iostream>
#include <fstream>

//#include <state-observation/noise/gaussian-white-noise.hpp>
//#include <state-observation/examples/offline-ekf-flexibility-estimation.hpp>
//#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>
//#include <state-observation/tools/miscellaneous-algorithms.hpp>

#include <state-observation/flexibility-estimation/model-base-ekf-flex-estimator-imu.hpp>

#include <time.h>


using namespace stateObservation;

int test()
{
    std::cout << "Debut du test" << std::endl;

  /// sampling period
    const double dt=5e-3;

  /// Measurement noise covariance
    Matrix Cov;
    Cov.resize(6,6);
    double unitCov = 1e0;
    Cov <<  unitCov,0,0,0,0,0,
            0,unitCov,0,0,0,0,
            0,0,unitCov,0,0,0,
            0,0,0,unitCov,0,0,
            0,0,0,0,unitCov,0,
            0,0,0,0,0,unitCov;

  /// Initializations
    // Dimensions things
    const unsigned kmax=14012;
    const unsigned measurementSize=6;
    const unsigned inputSize=54;
    const unsigned stateSize=18;
    unsigned contactNbr = 2;
    // State initialization => not used here because it is set in model-base-ekf-flex-estimator-imu
    Vector x0=Vector::Zero(stateSize,1);
    x0 <<   4.19416e-06,
            2.39345e-07,
            -0.00713916,
            -0.000638939,
            -0.0184693,
            5.30652e-06,
            -2.65572e-05,
            -3.50845e-05,
            -0.000156354,
            -0.000233892,
            -0.00268997,
            -7.22751e-05,
            2.11239e-05,
            0.000119309,
            -0.000310061,
            -0.000189487,
            -0.00423664,
            0.000263969;
     // Input initialization
     Vector u0=Vector::Zero(inputSize-6*contactNbr,1);
     u0 <<  0.0135672,
             0.001536,
             0.80771,
             -2.50425e-06,
             -1.03787e-08,
             5.4317e-08,
             -2.50434e-06,
             -1.03944e-08,
             5.45321e-08,
             48.1348,
             46.9498,
             1.76068,
             -0.0863332,
             -0.59487,
             -0.0402246,
             0,
             0,
             0,
             0,
             0,
             0,
             0,
             0,
             0,
             0,
             0,
             0,
             -0.098,
             -1.21619e-10,
             1.1174,
             3.06752e-22,
             -1.06094e-20,
             7.75345e-22,
             -2.84609e-06,
             -1.18496e-08,
             -4.52691e-18,
             2.95535e-20,
             -1.0346e-18,
             7.58731e-20,
             -0.000284609,
             -1.18496e-06,
             -4.52691e-16;

   /// Definitions of test vectors
     // State: what we want
     DiscreteTimeArray x;
     // Measurement
     DiscreteTimeArray y;
     y.getFromFile("source_measurement.dat",1,measurementSize);
     // Input
     DiscreteTimeArray u;
     u.getFromFile("source_input.dat",1,inputSize);

    stateObservation::flexibilityEstimation::ModelBaseEKFFlexEstimatorIMU est;
    est.setSamplingPeriod(dt);
    est.setInput(u0);
    est.setMeasurementInput(u0);

    est.setMeasurementNoiseCovariance(Cov);

    est.setContactsNumber(contactNbr);
    est.setInputSize(inputSize);

    Vector flexibility;
    flexibility.resize(18);

    for (int k=2;k<kmax;++k)
    {
        est.setMeasurement(y[k].transpose());
        est.setMeasurementInput(u[k].transpose());

        flexibility = est.getFlexibilityVector();
        x.setValue(flexibility,k);
    }

    x.writeInFile("state.dat");
    y.writeInFile("measurement.dat");
    u.writeInFile("input.dat");

    std::cout << "Fin du test" << std::endl;

}

int main()
{

    return test();

}





