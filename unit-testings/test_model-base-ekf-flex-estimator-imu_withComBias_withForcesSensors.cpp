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

    // For state vector
    bool withComBias_=false;

    // For measurement vector
    bool withUnmodeledForces_ = false;
    bool withForceSensors_=true;
    bool withAbsolutePose_ = false;

    // time
    const double dt=5e-3;
    const unsigned kinit=3;
    const unsigned kmax=1100;

    // fix sizes
    const unsigned inputSizeBase = 42;
    const unsigned measurementSizeBase=6;
    const unsigned stateSize=35;

    // variable sizes
    unsigned contactNbr = 2;
    unsigned inputSize=inputSizeBase+contactNbr*12;
    unsigned measurementSize=measurementSizeBase+withUnmodeledForces_*6+withAbsolutePose_*6+withForceSensors_*contactNbr*6;

    // State initialization => not used here because it is set in model-base-ekf-flex-estimator-imu

    /// Definitions of input vectors
     // Measurement
     IndexedMatrixArray y;
     std::cout << "Loading measurements file" << std::endl;
     y.getFromFile("source_measurement.dat",1,30);

     // Input
     IndexedMatrixArray u;
     std::cout << "Loading input file" << std::endl;
     u.getFromFile("source_input.dat",1,inputSize);

     // Input
     IndexedMatrixArray nbSupport;
     std::cout << "Loading the number of supports file" << std::endl;
     nbSupport.getFromFile("source_nbSupport.dat",1,1);

//      //state
//     IndexedMatrixArray xRef;
//       std::cout << "Loading reference state file" << std::endl;
//     xRef.getFromFile("source_state.dat",18,1);


    /// Definition of ouptut vectors
     // State: what we want
     IndexedMatrixArray x_output;
     // State covariance
     IndexedMatrixArray xCov_output;
     // Measurement
     IndexedMatrixArray y_output;
     // Input
     IndexedMatrixArray u_output;
     IndexedMatrixArray deltax_output;

     // Inovation
     IndexedMatrixArray i_output;
     // Inovation
     IndexedMatrixArray iInt_output;
     // Predicted sensors
     IndexedMatrixArray ps_output;
     // Simulated densors ss_output;
     IndexedMatrixArray ss_output;
     // Predicted sensors
     IndexedMatrixArray xPredicted_output;

     Vector flexibility;
     flexibility.resize(stateSize);
     Vector xdifference(flexibility);

     Vector inov;
     inov.resize(stateSize);
     Vector inovInt;
     inovInt.resize(stateSize);

     Vector simulatedMeasurements;
     simulatedMeasurements.resize(measurementSize);

     Vector predictedMeasurements;
     predictedMeasurements.resize(measurementSize);

     Vector predictedState;
     predictedState.resize(stateSize);

     Vector flexCovariance;
     flexCovariance.resize(stateSize);

     timespec time1, time2, time3;
     IndexedMatrixArray computationTime_output;
     double computationTime_moy=0;
     Vector computeTime;
     computeTime.resize(1);

    /// Initializations
    stateObservation::flexibilityEstimation::ModelBaseEKFFlexEstimatorIMU est;
    est.setSamplingPeriod(dt);

    std::cout << "setting covariances" << std::endl;

    // Measurement noise covariance
    est.setWithUnmodeledMeasurements(withUnmodeledForces_);
    est.setWithForcesMeasurements(withForceSensors_);
    est.setWithAbsolutePos(withAbsolutePose_);

    stateObservation::Matrix R_; R_.resize(measurementSizeBase,measurementSizeBase); R_.setIdentity();
    R_.block(0,0,3,3)*=1.e-3;
    R_.block(3,3,3,3)*=1.e-6;
    est.setMeasurementNoiseCovariance(R_);
    est.setUnmodeledForceVariance(1e-13);
    est.setForceVariance(1.e-4);
    est.setAbsolutePosVariance(1e-4);

    // Process noise covariance
    est.setWithComBias(withComBias_);
    stateObservation::Matrix Q_; Q_.resize(stateSize,stateSize); Q_.setIdentity();
    Q_.block(0,0,12,12)*=1.e-8;
    Q_.block(12,12,6,6)*=1.e-4;
    if(withComBias_) Q_.block(18,18,2,2)*=0;//1.e-13;
    est.setProcessNoiseCovariance(Q_);

//    // Flexibility noise covariance
//    stateObservation::Matrix P0_(stateSize,stateSize); P0_.setZero();
//    for(int i=0;i<stateSize;++i) P0_(i,i)=stateCovariance[kinit+2](i);
//    est.setFlexibilityCovariance(P0_);

    // Estimator state
    est.setInput(u[kinit+2].block(0,0,1,est.getInputSize()).transpose());
//    est.setFlexibilityGuess(xRef[kinit+2].block(0,0,est.getStateSize(),1));

    // Set contacts number
    est.setContactsNumber(contactNbr);
    est.setContactModel(stateObservation::flexibilityEstimation::
                ModelBaseEKFFlexEstimatorIMU::contactModel::elasticContact);

    est.setRobotMass(56.8);//48.6);//
//    stateObservation::Vector3 v1; v1.setOnes();
//    v1=100*v1;
//    est.setAngularAccelerationLimit(v1);
//    stateObservation::Vector3 v2; v2.setOnes();
//    v2=10*v2;
//    est.setLinearAccelerationLimit(v2);

    // Set stifness and damping
    est.setKfe(40000*Matrix3::Identity());
    est.setKte(350*Matrix3::Identity());
    est.setKfv(600*Matrix3::Identity());
    est.setKtv(10*Matrix3::Identity());

    double normState=0;
    Vector errorsum=Vector::Zero(est.getEKF().getStateSize());
    stateObservation::Vector input, measurement;

    std::cout << "Beginning reconstruction "<<std::endl;
    for (unsigned k=kinit+3;k<kmax;++k)
    {
        std::cout << k << std::endl;

        contactNbr = nbSupport[k](0);
        est.setContactsNumber(contactNbr);

        inputSize=inputSizeBase+contactNbr*12;
        input.resize(inputSize);
        input=(u[k].block(0,0,1,inputSize)).transpose();
        est.setMeasurementInput(input);

        measurementSize=measurementSizeBase+withUnmodeledForces_*6+withAbsolutePose_*6+withForceSensors_*contactNbr*6;
        measurement.resize(measurementSize);
        measurement.segment(0,6) = (y[k].block(0,0,1,6)).transpose();
        int i = 6;
        if(withUnmodeledForces_)
        {
            measurement.segment<6>(i) = (y[k].block(0,0,1,6)).transpose();
            i+=6;
        }
        if(withForceSensors_)
        {
            measurement.segment(i,6*contactNbr) = (y[k].block(0,12,1,6*contactNbr)).transpose();
            i+=6*contactNbr;
        }
        if(withAbsolutePose_)
        {
            measurement.segment(i,6) = (y[k].block(0,24,1,6)).transpose();
        }
        est.setMeasurement(measurement);

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);

        flexibility = est.getFlexibilityVector();

        for(int i=0;i<stateSize;++i){
            flexCovariance[i]=est.getFlexibilityCovariance()(i,i);
        }

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time3);
        computeTime[0]=(double)diff(time2,time3).tv_nsec-(double)diff(time1,time2).tv_nsec;

//        xdifference =flexibility-xRef[k];
//        normState += xdifference.squaredNorm();
//        errorsum += xdifference.cwiseProduct(xdifference);

        inov=est.getInovation();
        inovInt=inov*dt;
        simulatedMeasurements=est.getSimulatedMeasurement();
        predictedMeasurements=est.getLastPredictedMeasurement();
        predictedState=est.getPrediction();

        x_output.setValue(flexibility,k);
        xCov_output.setValue(flexCovariance,k);
        y_output.setValue(y[k],k);
        u_output.setValue(u[k],k);
        deltax_output.setValue(xdifference,k);

        i_output.setValue(inov,k);
        iInt_output.setValue(inovInt,k);
        ps_output.setValue(predictedMeasurements,k);
        ss_output.setValue(simulatedMeasurements,k);
        xPredicted_output.setValue(predictedState,k);

        computeTime[0]=est.getComputeFlexibilityTime();
        computationTime_output.setValue(computeTime,k);
        computationTime_moy+=computeTime[0];
    }

    std::cout << "Completed "<<std::endl;

    computeTime[0]=computationTime_moy/(kmax-kinit-3);
    computationTime_output.setValue(computeTime,kmax);
    computationTime_output.writeInFile("computationTime.dat");

    x_output.writeInFile("state.dat");
    xCov_output.writeInFile("stateCov.dat");
    y_output.writeInFile("measurement.dat");
    u_output.writeInFile("input.dat");

    i_output.writeInFile("inovation.dat");
    iInt_output.writeInFile("inovationInt.dat");
    ps_output.writeInFile("predictedSensors.dat");
    ss_output.writeInFile("simulatedSensors.dat");
    xPredicted_output.writeInFile("statePredicted.dat");

    std::cout << "Mean computation time " << computeTime[0] <<std::endl;

//    std::cout << "flexibility mean quadratic error " << normState/(kmax-kinit-3)<<std::endl;
//    errorsum = errorsum/(kmax-kinit-3);

//    Vector error(6);
//
//    error(0)= sqrt((errorsum(kine::pos) + errorsum(kine::pos+1) + errorsum(kine::pos+2))/(kmax-kinit-3));
//    error(1) = sqrt((errorsum(kine::linVel) + errorsum(kine::linVel+1) + errorsum(kine::linVel+2))/(kmax-kinit-3));
//    error(2) = sqrt((errorsum(kine::linAcc) + errorsum(kine::linAcc+1) + errorsum(kine::linAcc+2))/(kmax-kinit-3));
//    error(3) = sqrt((errorsum(kine::ori) + errorsum(kine::ori+1) + errorsum(kine::ori+2))/(kmax-kinit-3));
//    error(4) = sqrt((errorsum(kine::angVel) + errorsum(kine::angVel+1) + errorsum(kine::angVel+2))/(kmax-kinit-3));
//    error(5) = sqrt((errorsum(kine::angAcc) + errorsum(kine::angAcc+1) + errorsum(kine::angAcc+2))/(kmax-kinit-3));
//
//    double syntherror = 666*error(0)+10*error(3)+(10*error(1)+1*error(4))+(1*error(2)+1*error(5));
//    std::cout << "synthError=" << syntherror << std::endl;
//    std::cout << "Mean error " << error.transpose() <<std::endl;

//    if (normState/(kmax-kinit-3)>1e-04)
//    {
//      std::cout << "Failed : error is too big !!"<< std::endl <<"The end" << std::endl;
//      return 1;
//    }
//#ifdef NDEBUG
//    if (computeTime[0]>2e5)
//     {
//      std::cout << "Failed : Computation time is too long !!"<< std::endl <<"The end" << std::endl;
//      return 2;
//    }
//#endif

    std::cout << "Succeed !!"<< std::endl <<"The end" << std::endl;
    return 1;
}

int main()
{
    return test();
}

