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
    // Dimensions
    const unsigned kinit=0;
    const unsigned kmax=1400;
    const unsigned measurementSize=6;
    const unsigned inputSize=54;
    const unsigned stateSize=18;
    unsigned contactNbr = 2;
    // State initialization => not used here because it is set in model-base-ekf-flex-estimator-imu

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

    stateObservation::flexibilityEstimation::ModelBaseEKFFlexEstimatorIMU est;
    est.setSamplingPeriod(dt);
    est.setInput(u0);
    est.setMeasurementInput(u0);

    est.setMass(56.8679920);
    stateObservation::Vector3 v1; v1.setOnes();
    v1=1000*v1;
    est.setAngularAccelerationLimit(v1);
    stateObservation::Vector3 v2; v2.setOnes();
    v2=1000*v2;
    est.setLinearAccelerationLimit(v2);
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


    est.setMeasurementNoiseCovariance(Cov);

    est.setContactsNumber(contactNbr);


    Vector flexibility;
    flexibility.resize(18);
    Vector xdifference(flexibility);

    timespec time1, time2, time3;
    IndexedMatrixArray computationTime_output;
    double computationTime_moy=0;
    Vector computeTime;
    computeTime.resize(1);

    double norm=0;

    est.setContactModel(stateObservation::flexibilityEstimation::
                ModelBaseEKFFlexEstimatorIMU::contactModel::elasticContact);

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





