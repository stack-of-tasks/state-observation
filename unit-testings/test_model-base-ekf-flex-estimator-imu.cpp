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
    double unitCov = 1e-2;
    Cov <<  unitCov,0,0,0,0,0,
            0,unitCov,0,0,0,0,
            0,0,unitCov,0,0,0,
            0,0,0,unitCov,0,0,
            0,0,0,0,unitCov,0,
            0,0,0,0,0,unitCov;

  /// Initializations
    // Dimensions things
    const unsigned kinit=0;
    const unsigned kmax=1400;
    const unsigned measurementSize=6;
    const unsigned inputSize=54;
    const unsigned stateSize=18;
    unsigned contactNbr = 2;
    // State initialization => not used here because it is set in model-base-ekf-flex-estimator-imu
//    Vector x0=Vector::Zero(stateSize,1);
//    x0 <<   4.19416e-06,
//            2.39345e-07,
//            -0.00713916,
//            -0.000638939,
//            -0.0184693,
//            5.30652e-06,
//            -2.65572e-05,
//            -3.50845e-05,
//            -0.000156354,
//            -0.000233892,
//            -0.00268997,
//            -7.22751e-05,
//            2.11239e-05,
//            0.000119309,
//            -0.000310061,
//            -0.000189487,
//            -0.00423664,
//            0.000263969;
     // Input initialization
     Vector u0=Vector::Zero(inputSize-6*contactNbr,1);
     u0 <<  0.0135673,
             0.001536,
             0.80771,
             -2.63605e-06,
             -1.09258e-08,
             5.71759e-08,
             2.71345,
             0.3072,
             161.542,
             48.1348,
             46.9498,
             1.76068,
             -0.0863332,
             -0.594871,
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
             -6.23712e-11,
             1.1174,
             1.58984e-22,
             -5.43636e-21,
             3.9598e-22,
             -2.99589e-06,
             -1.24742e-08,
             -4.7647e-18,
             3.17968e-20,
             -1.08727e-18,
             7.91959e-20,
             -0.000299589,
             -1.24742e-06,
             -4.7647e-16;

   /// Definitions of input vectors
     // Measurement
     DiscreteTimeArray y;
     y.getFromFile("source_measurement.dat",1,measurementSize);
     // Input
     DiscreteTimeArray u;
     u.getFromFile("source_input.dat",1,inputSize);

   /// Definition of ouptut vectors
     // State: what we want
     DiscreteTimeArray x_output;
     // Measurement
     DiscreteTimeArray y_output;
     // Input
     DiscreteTimeArray u_output;



    stateObservation::flexibilityEstimation::ModelBaseEKFFlexEstimatorIMU est;
    est.setSamplingPeriod(dt);
    est.setInput(u0);
    est.setMeasurementInput(u0);
 //   est.getEKF().setState(x0,0);

    est.setMeasurementNoiseCovariance(Cov);

    est.setContactsNumber(contactNbr);
    est.setInputSize(inputSize);

//    u0.resize(inputSize);
//    u0 <<  0.013567,
//            0.001536,
//            0.80771,
//            -4.85723e-16,
//            4.33681e-18,
//            -1.08843e-14,
//            0,
//            0,
//            0,
//            48.1348,
//            46.9498,
//            1.76068,
//            -0.0863332,
//            -0.594861,
//            -0.0402246,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            -0.0980005,
//            -1.25418e-09,
//            1.1174,
//            2.25124e-24,
//            1.21594e-20,
//            -8.70945e-22,
//            4.82571e-18,
//            9.83976e-18,
//            -1.08802e-14,
//            -4.54866e-24,
//            -2.45645e-20,
//            1.75948e-21,
//            1.97163e-15,
//            5.58275e-21,
//            -2.46868e-20,
//            0.00949046,
//            -0.095,
//            1.98197e-07,
//            -2.06795e-24,
//            -7.44034e-16,
//            -1.73252e-24,
//            0.00949046,
//            0.095,
//            1.98197e-07,
//            -1.39829e-24,
//            -4.00152e-16,
//            -7.3383e-25;
//    est.setInput(u0);
//    est.setMeasurementInput(u0);

    Vector flexibility;
    flexibility.resize(18);

    for (int k=kinit+2;k<kmax;++k)
    {
        est.setMeasurement(y[k].transpose());
        est.setMeasurementInput(u[k].transpose());

        flexibility = est.getFlexibilityVector();
        x_output.setValue(flexibility,k);
        y_output.setValue(y[k],k);
        u_output.setValue(u[k],k);
    }

    x_output.writeInFile("state.dat");
    y_output.writeInFile("measurement.dat");
    u_output.writeInFile("input.dat");

    std::cout << "Fin du test" << std::endl;

}

int main()
{

    return test();

}





