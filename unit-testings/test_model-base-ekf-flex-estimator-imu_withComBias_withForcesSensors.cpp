#include <iostream>
#include <fstream>

//#include <state-observation/noise/gaussian-white-noise.hpp>
//#include <state-observation/examples/offline-ekf-flexibility-estimation.hpp>
//#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>
//#include <state-observation/tools/miscellaneous-algorithms.hpp>

#include <state-observation/flexibility-estimation/model-base-ekf-flex-estimator-imu.hpp>

#include <time.h>


using namespace stateObservation;

timespec diff(const timespec & start, const timespec & end)
{
        timespec temp;
        if ((end.tv_nsec-start.tv_nsec)<0) {
                temp.tv_sec = end.tv_sec-start.tv_sec-1;
                temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
        } else {
                temp.tv_sec = end.tv_sec-start.tv_sec;
                temp.tv_nsec = end.tv_nsec-start.tv_nsec;
        }
        return temp;
}

int test()
{
    std::cout << "Starting" << std::endl;

    bool withComBias_=false;
    bool withForceSensors_=false;

    // Dimensions
    const double dt=5e-3;
    const unsigned kinit=8;
    const unsigned kmax=1400;
    unsigned contactNbr = 2;
    const unsigned inputSize=42+contactNbr*6;
    const unsigned measurementSizeBase=6;
    const unsigned measurementSize=measurementSizeBase+withForceSensors_*contactNbr*6;
    const unsigned stateSizeBase_=18;
    const unsigned stateSize=stateSizeBase_+(int)withComBias_*2;

    // State initialization => not used here because it is set in model-base-ekf-flex-estimator-imu

  /// Initializations
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

    Vector m0=Vector::Zero(measurementSize);
    m0 <<  -0.00191194,
           -0.00417779,
           9.80011,
           -0.000159848,
           0.00118677,
           -1.01118e-05;

    Vector x0=Vector::Zero(stateSize);
    x0 <<  0.00579374,
           4.77747e-05,
           -0.00503099,
           2.8051e-05,
           0.00842914,
           -0.000329279,
           -0.377671,
           -0.00429404,
           -0.167778,
           -0.00281867,
           0.288473,
           0.000692284,
           3.10531,
           0.0816574,
           0.986816,
           0.0726735,
           -4.03091,
           0.228881;

    stateObservation::Matrix P0_; P0_.resize(stateSize,stateSize);
    P0_ <<  2.53574e-08, 9.8097e-12, 3.09505e-09, -5.82127e-11, 1.43065e-08, 1.4134e-10, -7.11114e-07, -7.0412e-10, -5.4998e-08, -4.24379e-10, 7.07223e-07, -8.74669e-10, -3.09802e-05, 4.00308e-08, -2.09651e-06, 5.02917e-08, 3.0132e-05, -2.45227e-08,
           9.8097e-12, 3.03951e-08, 5.02406e-11, -1.59497e-08, 2.49844e-11, -2.17655e-09, -2.97152e-10, -8.64265e-07, 2.81677e-11, -8.36905e-07, 3.50777e-10, 1.09398e-07, -3.27145e-08, -4.53405e-05, -1.77178e-09, -4.302e-05, 3.93155e-08, 3.51997e-06,
           3.09505e-09, 5.02406e-11, 1.5331e-08, 5.83244e-12, 2.82322e-09, 1.1099e-11, -7.72645e-08, -7.52285e-11, -3.90394e-08, -3.80282e-11, 7.85662e-08, 3.12856e-10, -1.54327e-05, -3.28099e-07, -6.25128e-06, -3.1672e-07, 1.56903e-05, 1.09632e-07,
           -5.82127e-11, -1.59497e-08, 5.83244e-12, 1.01798e-07, -1.06038e-10, -4.08813e-09, 3.65534e-09, 8.21972e-07, -1.27679e-10, 7.54921e-07, -3.7944e-09, -1.67646e-08, -4.09666e-08, 6.45776e-06, -3.75853e-10, 1.84885e-06, 5.83376e-08, 4.14927e-06,
           1.43065e-08, 2.49844e-11, 2.82322e-09, -1.06038e-10, 1.03951e-07, -1.18105e-11, -7.54383e-07, -1.42177e-09, -8.7564e-08, -9.27792e-11, 7.33019e-07, -1.47306e-09, -9.68387e-06, -5.1679e-09, -1.84061e-06, 6.19559e-09, 6.84397e-06, -5.22238e-08,
           1.4134e-10, -2.17655e-09, 1.1099e-11, -4.08813e-09, -1.18105e-11, 6.75919e-08, -6.68363e-09, 2.96585e-08, -3.55745e-10, 2.98775e-08, 5.28615e-09, -5.90598e-07, -2.58879e-07, 1.11262e-05, -1.52355e-08, 1.0784e-05, 2.03524e-07, -2.30801e-05,
           -7.11114e-07, -2.97152e-10, -7.72645e-08, 3.65534e-09, -7.54383e-07, -6.68363e-09, 4.13448e-05, 2.09017e-08, 2.99499e-06, 4.94288e-09, -4.10041e-05, 6.25423e-08, 0.000204982, -9.13586e-07, 1.6441e-05, -1.01823e-06, -9.80051e-05, 8.77898e-06,
           -7.0412e-10, -8.64265e-07, -7.52285e-11, 8.21972e-07, -1.42177e-09, 2.96585e-08, 2.09017e-08, 5.22679e-05, -5.7584e-09, 5.0475e-05, -2.14392e-08, -5.7461e-06, 2.55252e-06, -0.000124243, 4.84995e-07, -0.000282105, -1.83147e-06, 0.00038937,
           -5.4998e-08, 2.81677e-11, -3.90394e-08, -1.27679e-10, -8.7564e-08, -3.55745e-10, 2.99499e-06, -5.7584e-09, 1.93935e-06, -9.36749e-09, -2.98692e-06, 1.69319e-09, 4.22052e-05, 4.78036e-07, 1.15707e-05, 4.87627e-07, -3.48383e-05, 4.69863e-07,
           -4.24379e-10, -8.36905e-07, -3.80282e-11, 7.54921e-07, -9.27792e-11, 2.98775e-08, 4.94288e-09, 5.0475e-05, -9.36749e-09, 4.90788e-05, -7.12519e-09, -5.69618e-06, 2.49196e-06, -0.000102352, 4.679e-07, -0.000254325, -1.83484e-06, 0.000381884,
           7.07223e-07, 3.50777e-10, 7.85662e-08, -3.7944e-09, 7.33019e-07, 5.28615e-09, -4.10041e-05, -2.14392e-08, -2.98692e-06, -7.12519e-09, 4.09618e-05, -3.50854e-08, -0.000216523, 7.09046e-07, -1.81174e-05, 8.2821e-07, 0.000110377, -9.06836e-06,
           -8.74669e-10, 1.09398e-07, 3.12856e-10, -1.67646e-08, -1.47306e-09, -5.90598e-07, 6.25423e-08, -5.7461e-06, 1.69319e-09, -5.69618e-06, -3.50854e-08, 8.88517e-06, -1.19927e-06, -7.61875e-05, -3.68304e-07, -5.81568e-05, 1.45156e-06, 2.79828e-05,
           -3.09802e-05, -3.27145e-08, -1.54327e-05, -4.09666e-08, -9.68387e-06, -2.58879e-07, 0.000204982, 2.55252e-06, 4.22052e-05, 2.49196e-06, -0.000216523, -1.19927e-06, 0.242923, -0.000215012, 0.0153966, -0.000286973, -0.251132, -0.00110584,
           4.00308e-08, -4.53405e-05, -3.28099e-07, 6.45776e-06, -5.1679e-09, 1.11262e-05, -9.13586e-07, -0.000124243, 4.78036e-07, -0.000102352, 7.09046e-07, -7.61875e-05, -0.000215012, 0.369984, -0.000122324, 0.372829, 7.30365e-05, -0.0780296,
           -2.09651e-06, -1.77178e-09, -6.25128e-06, -3.75853e-10, -1.84061e-06, -1.52355e-08, 1.6441e-05, 4.84995e-07, 1.15707e-05, 4.679e-07, -1.81174e-05, -3.68304e-07, 0.0153966, -0.000122324, 0.0082466, -0.000142214, -0.0157927, -2.33657e-05,
           5.02917e-08, -4.302e-05, -3.1672e-07, 1.84885e-06, 6.19559e-09, 1.0784e-05, -1.01823e-06, -0.000282105, 4.87627e-07, -0.000254325, 8.2821e-07, -5.81568e-05, -0.000286973, 0.372829, -0.000142214, 0.376433, 0.000144442, -0.0794679,
           3.0132e-05, 3.93155e-08, 1.56903e-05, 5.83376e-08, 6.84397e-06, 2.03524e-07, -9.80051e-05, -1.83147e-06, -3.48383e-05, -1.83484e-06, 0.000110377, 1.45156e-06, -0.251132, 7.30365e-05, -0.0157927, 0.000144442, 0.260155, 0.00124888,
           -2.45227e-08, 3.51997e-06, 1.09632e-07, 4.14927e-06, -5.22238e-08, -2.30801e-05, 8.77898e-06, 0.00038937, 4.69863e-07, 0.000381884, -9.06836e-06, 2.79828e-05, -0.00110584, -0.0780296, -2.33657e-05, -0.0794679, 0.00124888, 0.0384546;

    stateObservation::flexibilityEstimation::ModelBaseEKFFlexEstimatorIMU est;
    est.setSamplingPeriod(dt);

    // Measurement noise covariance
    stateObservation::Matrix R_; R_.resize(measurementSizeBase,measurementSizeBase); R_.setIdentity();
    R_.block(0,0,3,3)*=1.e-2;
    R_.block(3,3,3,3)*=1.e-2;
    est.setMeasurementNoiseCovariance(R_);   
    est.setWithForcesMeasurements(withForceSensors_);
    est.setForceVariance(1.e-4);

    // Process noise covariance
    est.setWithComBias(withComBias_);
    stateObservation::Matrix Q_; Q_.resize(stateSize,stateSize); Q_.setIdentity();
    Q_.block(0,0,12,12)*=1.e-8;
    Q_.block(12,12,6,6)*=1.e-4;
    if(withComBias_) Q_.block(18,18,2,2)*=1e-13;
    est.setProcessNoiseCovariance(Q_);

    // Flexibility noise covariance
    est.setFlexibilityCovariance(P0_);

    // estimator state
    est.setInput(u0.block(0,0,inputSize-6*contactNbr,1));
//    est.setMeasurementInput(u0.block(0,0,inputSize-6*contactNbr,1));
//    est.setMeasurement(m0.block(0,0,measurementSizeBase,1));
    est.setFlexibilityGuess(x0);

    // Set contacts number
    est.setContactsNumber(contactNbr);
    est.setContactModel(stateObservation::flexibilityEstimation::
                ModelBaseEKFFlexEstimatorIMU::contactModel::elasticContact);

    // Set stifness and damping
    est.setKfe(40000*Matrix3::Identity());
    est.setKte(600*Matrix3::Identity());
    est.setKfv(600*Matrix3::Identity());
    est.setKtv(60*Matrix3::Identity());

   /// Definitions of input vectors
     // Measurement
    IndexedMatrixArray y;
    std::cout << "Loading measurements file" << std::endl;
    y.getFromFile("source_measurement.dat",1,measurementSize);
     // Input
    IndexedMatrixArray u;
     std::cout << "Loading input file" << std::endl;
    u.getFromFile("source_input.dat",1,inputSize);
     //state
    IndexedMatrixArray xRef;
      std::cout << "Loading reference state file" << std::endl;
    xRef.getFromFile("source_state.dat",stateSize,1);

   /// Definition of ouptut vectors
     // State: what we want
    IndexedMatrixArray x_output;
     // Measurement
    IndexedMatrixArray y_output;
     // Input
    IndexedMatrixArray u_output;
    IndexedMatrixArray deltax_output;

    Vector flexibility;
    flexibility.resize(stateSize);
    Vector xdifference(flexibility);

    timespec time1, time2, time3;
    IndexedMatrixArray computationTime_output;
    double computationTime_moy=0;
    Vector computeTime;
    computeTime.resize(1);

    double norm=0;

    std::cout << "Beginning reconstruction "<<std::endl;
    for (unsigned k=kinit+2;k<kmax;++k)
    {
        est.setMeasurement(y[k].transpose());
        est.setMeasurementInput(u[k].transpose());

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);

        flexibility = est.getFlexibilityVector();

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time3);
        computeTime[0]=(double)diff(time2,time3).tv_nsec-(double)diff(time1,time2).tv_nsec;

        xdifference =flexibility-xRef[k];
        norm += xdifference.squaredNorm();

        if(k==10) std::cout << y[k].transpose() << std::endl;

        x_output.setValue(flexibility,k);
        y_output.setValue(y[k],k);
        u_output.setValue(u[k],k);
        deltax_output.setValue(xdifference,k);

        computeTime[0]=est.getComputeFlexibilityTime();
        computationTime_output.setValue(computeTime,k);
        computationTime_moy+=computeTime[0];
    }

    std::cout << "Completed "<<std::endl;

    computeTime[0]=computationTime_moy/(kmax-kinit-10);
    computationTime_output.setValue(computeTime,kmax);
    computationTime_output.writeInFile("computationTime.dat");

    x_output.writeInFile("state.dat");
    y_output.writeInFile("measurement.dat");
    u_output.writeInFile("input.dat");

    std::cout << "Mean computation time " << computeTime[0] <<std::endl;

    std::cout << "Mean quadratic error " << norm/(kmax-kinit-2)<<std::endl;

    if (norm/(kmax-kinit-2)>1e-04)
    {
      std::cout << "Failed : error is too big !!"<< std::endl <<"The end" << std::endl;
      return 1;
    }
#ifdef NDEBUG
    if (computeTime[0]>2e5)
     {
      std::cout << "Failed : Computation time is too long !!"<< std::endl <<"The end" << std::endl;
      return 2;
    }
#endif

    std::cout << "Succeed !!"<< std::endl <<"The end" << std::endl;
    return 1;
}

int main()
{
    return test();
}





