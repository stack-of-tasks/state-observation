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

    bool withComBias_=true;
    bool withForceSensors_=true;

    // Dimensions
    const unsigned kinit=0;
    const unsigned kmax=1848;
    const unsigned inputSize=54;
    const unsigned measurementSizeBase=6;
    const unsigned measurementSize=measurementSizeBase+withForceSensors_*12;
    const unsigned stateSizeBase_=18;
    const unsigned stateSize=stateSizeBase_+(int)withComBias_*2;
    unsigned contactNbr = 2;
    // State initialization => not used here because it is set in model-base-ekf-flex-estimator-imu

  /// sampling period
    const double dt=5e-3;

  /// Initializations
     // Input initialization
    Vector u0=Vector::Zero(inputSize);//-6*contactNbr,1);
    u0 <<   0.0142039,
            0.00165921,
            0.80771,
            0.127323,
            0.0246383,
            3.19163e-05,
            12.7323,
            2.46383,
            0.00319163,
            45.2594,
            44.0532,
            1.73905,
            -0.00185306,
            -0.375901,
            -0.0996394,
            0,
            0,
            0,
            0,
            0,
            0,
            -4.53201,
            2.90831,
            0.0296007,
            -113.171,
            584.828,
            0.788792,
            -0.0972751,
            0.000140574,
            1.11742,
            -3.41081e-16,
            2.71888e-16,
            -3.80512e-17,
            0.144989,
            0.0281149,
            0.00319886,
            -6.82163e-14,
            5.43775e-14,
            -7.61025e-15,
            14.4989,
            2.81149,
            0.319886,
            0.00949046,
            0.095,
            1.98197e-07,
            0,
            0,
            0,
            0.00949605,
            -0.095,
            1.98197e-07,
            0,
            0,
            0;

//    u0 <<   0.0135672,
//            0.001536,
//            0.80771,
//            -2.50425e-06,
//            -1.03787e-08,
//            5.4317e-08,
//            -2.50434e-06,
//            -1.03944e-08,
//            5.45321e-08,
//            48.1348,
//            46.9498,
//            1.76068,
//            -0.0863332,
//            -0.59487,
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
//            -0.098,
//            -1.21619e-10,
//            1.1174,
//            3.06752e-22,
//            -1.06094e-20,
//            7.75345e-22,
//            -2.84609e-06,
//            -1.18496e-08,
//            -4.52691e-18,
//            2.95535e-20,
//            -1.0346e-18,
//            7.58731e-20,
//            -0.000284609,
//            -1.18496e-06,
//            -4.52691e-16,
//            0.00949046,
//            0.095,
//            1.98197e-07,
//            0,
//            0,
//            0,
//            0.00949605,
//            -0.095,
//            1.98197e-07,
//            0,
//            0,
//            0;

    Vector m0=Vector::Zero(measurementSize);
    m0 <<       -0.0330502,
               -0.169031,
               9.91812,
               0.0137655,
               0.0797922,
               0.000778988,
                   6.15302,
                   -8.44315,
                   245.826,
                   0.749795,
                   2.59329,
                   0.140388,
                   5.59534,
                   7.49882,
                   227.461,
                   0.24084,
                   2.74922,
                   -0.120347;

//    m0 << -7.03234,
//          0.0805895,
//          13.5925,
//          0.000333733,
//          -0.157283,
//          -0.00480441,
//          45.1262,
//          -21.367,
//          361.344,
//          1.12135,
//          -14.5562,
//          1.89125,
//          44.6005,
//          21.7871,
//          352.85,
//          -1.00715,
//          -14.5158,
//          -1.72017;

    stateObservation::flexibilityEstimation::ModelBaseEKFFlexEstimatorIMU est;
    est.setSamplingPeriod(dt);

    est.setInput(u0.block(0,0,inputSize-6*contactNbr,1));
    est.setMeasurementInput(u0.block(0,0,inputSize-6*contactNbr,1));
    est.setMeasurement(m0.block(0,0,measurementSizeBase,1));

    // Measurement noise covariance
    stateObservation::Matrix R_; R_.resize(measurementSizeBase,measurementSizeBase); R_.setIdentity();
    R_.block(0,0,3,3)*=1.e-3;
    R_.block(3,3,3,3)*=1.e-6;
    est.setMeasurementNoiseCovariance(R_);   
    est.setForceVariance(1.e-4);
    est.setWithForcesMeasurements(withForceSensors_);

    // Process noise covariance
    est.setWithComBias(withComBias_);
    stateObservation::Matrix Q_; Q_.resize(stateSize,stateSize); Q_.setIdentity();
    Q_.block(0,0,12,12)*=1.e-8;
    Q_.block(12,12,6,6)*=1.e-4;
    if(withComBias_) Q_.block(18,18,2,2)*=1e-13;
    est.setProcessNoiseCovariance(Q_);

    est.setContactsNumber(contactNbr);
    est.setContactModel(stateObservation::flexibilityEstimation::
                ModelBaseEKFFlexEstimatorIMU::contactModel::elasticContact);

//    est.setInput(u0);
//    est.setMeasurementInput(u0);
//    est.setMeasurement(m0);

    est.setKfe(40000*Matrix3::Identity());
    est.setKte(350*Matrix3::Identity());
    est.setKfv(600*Matrix3::Identity());
    est.setKtv(10*Matrix3::Identity());

   /// Definitions of input vectors
     // Measurement
    IndexedMatrixArray y;
    std::cout << "Loading measurements file" << std::endl;
    y.getFromFile("inputFiles/source_measurement.dat",1,measurementSize);
     // Input
    IndexedMatrixArray u;
     std::cout << "Loading input file" << std::endl;
    u.getFromFile("inputFiles/source_input.dat",1,inputSize);
      //state
    IndexedMatrixArray xRef;
      std::cout << "Loading reference state file" << std::endl;
    xRef.getFromFile("inputFiles/source_state.dat",stateSize,1);

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

        x_output.setValue(flexibility,k);
        y_output.setValue(y[k],k);
        u_output.setValue(u[k],k);
        deltax_output.setValue(xdifference,k);

        computeTime[0]=est.getComputeFlexibilityTime();
        computationTime_output.setValue(computeTime,k);
        computationTime_moy+=computeTime[0];
    }

    std::cout << "Completed "<<std::endl;

    computeTime[0]=computationTime_moy/(kmax-kinit-2);
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





