 /*
 * dynamical-system.cpp
 *
 *  Created on: 19 mai 2014
 *      Author: alexis
 */

#include <state-observation/flexibility-estimation/imu-elastic-local-frame-dynamical-system.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>

#include <iostream>


namespace stateObservation
{
namespace flexibilityEstimation
{
    using namespace stateObservation;



   IMUElasticLocalFrameDynamicalSystem::
    	IMUElasticLocalFrameDynamicalSystem(double dt):
        processNoise_(0x0), dt_(dt),orientationVector_(Vector3::Zero()),
        curRotation_(Matrix3::Identity()),
        measurementSize_(measurementSizeBase_),
        robotMassInv_(1/hrp2::m),
        robotMass_(hrp2::m)
    {
#ifdef STATEOBSERVATION_VERBOUS_CONSTRUCTORS
       // std::cout<<std::endl<<"IMUElasticLocalFrameDynamicalSystem Constructor"<<std::endl;
#endif //STATEOBSERVATION_VERBOUS_CONSTRUCTOR

      calculationState <<  -1,
                           -1,
                           -1;

      sensor_.setMatrixMode(true);

      kcurrent_=-1;
      disp=true;

    }

    IMUElasticLocalFrameDynamicalSystem::
                    ~IMUElasticLocalFrameDynamicalSystem()
    {
        //dtor
    }

    void IMUElasticLocalFrameDynamicalSystem::setRobotMass(double m)
    {
      robotMass_ = m;
      robotMassInv_=1/m;
    }


    void IMUElasticLocalFrameDynamicalSystem::getForcesAndMoments
                              (const IndexedMatrixArray& contactPosArray,
                               const IndexedMatrixArray& contactOriArray,
                               const Vector3& position, const Vector3& linVelocity,
                               const Vector3& oriVector, const Matrix3& orientation,
                               const Vector3& angVel,
                               Vector3& forces, Vector3& moments)
    {
      unsigned nbContacts(getContactsNumber());
      fc_.resize(nbContacts*3);
      tc_.resize(nbContacts*3);
      forces.setZero();
      moments.setZero();
      Matrix3 Rci;
      Matrix3 Rcit;
      Vector3 contactPos;

      Vector3 forcei;
      Vector3 momenti;
      Vector3 globalContactPos;



      if (disp)
      {
        std::cout << "Forces and moments ============="<<std::endl;

        std::cout << "oriVector "<<oriVector.transpose()<<std::endl;
        std::cout << "angVel "<< angVel.transpose()<<std::endl;
        std::cout << "position "<<position.transpose()<<std::endl;
        std::cout << "linVelocity "<<linVelocity.transpose()<<std::endl;
      }

      for (int i = 0; i<nbContacts ; ++i)
      {
        if (disp)
          std::cout << "contactNumber ============== "<<i <<std::endl;

        contactPos = contactPosArray[i];

        Rci = contactOriArray[i];
        Rcit= Rci.transpose();

        globalContactPos = position ;
        globalContactPos.noalias() += orientation*contactPos ;

        forcei.noalias() = - Rci*Kfe_*Rcit*(globalContactPos-contactPos);
        forcei.noalias() += - Rci*Kfv_*Rcit*(kine::skewSymmetric(angVel)*orientation*contactPos+linVelocity);

        fc_.segment<3>(3*i)= forcei;

        forces += forcei;

        momenti.noalias() = -Rci*Kte_*Rcit*oriVector;
        momenti.noalias() += -Rci*Ktv_*Rcit*angVel;
        momenti.noalias() += kine::skewSymmetric(globalContactPos)*forcei;

        tc_.segment<3>(3*i)= momenti;

        moments += momenti;

        if (disp)
        {

          std::cout << "contactPos "<<contactPos.transpose()<<std::endl;
          std::cout << "globalContactPos "<<globalContactPos.transpose()<<std::endl;
          std::cout << "globalContactPos-contactPos "<<(globalContactPos-contactPos).transpose()<<std::endl;
          std::cout << "omega^orientation*contactPos+v "<<(kine::skewSymmetric(angVel)*orientation*contactPos+linVelocity).transpose()<<std::endl;
          std::cout << "RKpR'*dsup_i"<<(- Rci*Kfe_*Rcit*(globalContactPos-contactPos)).transpose() <<std::endl;
          std::cout << "RKdR'*dotsup_i"<<((- Rci*Kfv_*Rcit*(kine::skewSymmetric(angVel)*orientation*contactPos+linVelocity))).transpose() <<std::endl;


        //std::cout << "orientation "<<std::endl<<orientation<<std::endl;

          std::cout << "forcei "<<forcei.transpose()<<std::endl;
          std::cout << "momenti "<<momenti.transpose()<<std::endl;
        }


      }

      if (disp)
      {
        std::cout << "Total ============== "<<std::endl;

        std::cout << "forces "<<forces.transpose()<<std::endl;
        std::cout << "moments "<<moments.transpose()<<std::endl;

        std::cout << "End ============== "<<std::endl;
      }
    }


    void IMUElasticLocalFrameDynamicalSystem::test()
    {
        Vector x, u;
        x.resize(18);
        u.resize(54);
//        x << 1.95675e-08,-2.22298e-10,-0.0109575,1.54011e-12,1.35566e-10,3.81192e-23,1.41111e-06,-1.60311e-08,-0.0489989,6.56636e-12,-3.14893e-09,-9.62102e-11,-5.76107e-05,6.54491e-07,9.47828,-8.92673e-08,-8.41481e-06,-7.82905e-09;
//        u << 0.0135673,0.001536,0.80771,-2.63605e-06,-1.09258e-08,5.71759e-08,2.71345,0.3072,161.542,1.55342e-318,-6.53681e-66,-6.39612e-66,-6.54978e-66,-6.53171e-66,-6.55957e-66,0,0,0,0,0,0,0,0,0,0,0,0,-0.098,-6.23712e-11,1.1174,1.56933e-22,-5.40778e-21,3.86235e-22,-2.99589e-06,-1.24742e-08,-4.76329e-18,3.13865e-20,-1.08156e-18,7.72471e-20,-0.000299589,-1.24742e-06,-4.76329e-16,0.00949046,-0.095,1.98197e-07,8.75907e-23,-0,7.17001e-22,0.00949046,0.095,1.98197e-07,6.29302e-23,1.01427e-16,-6.61701e-22;

        x << 0,0,-0.0109575,0,0,0,0,0,-0.0489989,0,0,0,0,0,0,0,0,0;
        u << 0.0135673,0.001536,0.80771,0,0,0,2.71345,0.3072,161.542,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.098,0,1.1174,0,0,0,0,0,0,0,0,0,-0.000299589,0,0,0.00949046,-0.095,0,0,-0,0,0.00949046,0.095,0,0,0,0;

//        std::cout << "AccAng" << computeAccelerationAngular(x,u,1) << std::endl;
//        std::cout << "AccLin" << computeAccelerationLinear(x,u,1) << std::endl;
    }


    Vector3 IMUElasticLocalFrameDynamicalSystem::computeAccelerations
       (const Vector3& positionCom, const Vector3& velocityCom,
        const Vector3& accelerationCom, const Vector3& AngMomentum,
        const Vector3& dotAngMomentum,
        const Matrix3& Inertia, const Matrix3& dotInertia,
        const IndexedMatrixArray& contactPosV,
        const IndexedMatrixArray& contactOriV,
        const Vector3& position, const Vector3& linVelocity, Vector3& linearAcceleration,
        const Vector3 &oriVector ,const Matrix3& orientation,
        const Vector3& angularVel, Vector3& angularAcceleration)
    {

        if (disp)
        {
          std::cout << "positionCom "<<positionCom.transpose()<<std::endl;
          std::cout << "velocityCom "<< velocityCom.transpose()<<std::endl;
          std::cout << "accelerationCom "<<accelerationCom.transpose()<<std::endl;
        }

        Vector3 fc;
        Vector3 tc;

        Matrix3 skewV(kine::skewSymmetric(angularVel));
        Matrix3 skewV2(kine::skewSymmetric2(angularVel));
        Matrix3 skewVR(skewV * orientation);
        Matrix3 skewV2R(skewV2 * orientation);

        Matrix3 orientationT=orientation.transpose();

        getForcesAndMoments (contactPosV, contactOriV,
                          position, linVelocity, oriVector, orientation,
                             angularVel, fc, tc);

        Vector3 vf (robotMassInv_*fc);
        vf.noalias() -= orientation*accelerationCom;
        vf.noalias() -= 2*skewVR*velocityCom;
        vf.noalias() -= skewV2R*positionCom;
        vf.noalias() -= cst::gravity;

        Vector3 vt (tc);
        vt.noalias() -= skewVR * (Inertia* (orientationT * angularVel));
        vt.noalias() -= orientation * (dotInertia * (orientationT * angularVel)) ;
        vt.noalias() -= orientation * dotAngMomentum;
        vt.noalias() -= skewVR * AngMomentum;
        vt.noalias() -= robotMass_* (kine::skewSymmetric(position) *
                            (skewV2R * positionCom + 2*(skewVR * velocityCom) + orientation * accelerationCom ));
        vt.noalias() -= robotMass_* kine::skewSymmetric(orientation * positionCom + position) * cst::gravity;


        angularAcceleration = (orientation*Inertia*orientationT + robotMass_*kine::skewSymmetric2(orientation * positionCom)).inverse()
                                *(vt - robotMass_*kine::skewSymmetric(orientation * positionCom + position)*vf);


        linearAcceleration = vf;
        linearAcceleration += kine::skewSymmetric(orientation*positionCom)*angularAcceleration;

        if (disp)
        {

          std::cout << "robotMassInv_*fc " << robotMassInv_*fc.transpose()<<std::endl;
          std::cout << "orientation*accelerationCom " << (orientation*accelerationCom).transpose()<<std::endl;
          std::cout << "2*skewVR*velocityCom " << (2*skewVR*velocityCom).transpose()<<std::endl;
          std::cout << "skewV2R*positionCom " << (skewV2R*positionCom).transpose()<<std::endl;
          std::cout << "cst::gravity " << (cst::gravity).transpose()<<std::endl;
          std::cout << "vf "<<vf.transpose()<<std::endl;

          std::cout << "w^ R I RT w "  << (skewVR * (Inertia* (orientationT * angularVel))).transpose()<<std::endl;
          std::cout << " R dotI RT w " << (orientation * (dotInertia * (orientationT * angularVel)) ).transpose()<<std::endl;
          std::cout << " R * dotAngMomentum " << (orientation * dotAngMomentum ).transpose()<<std::endl;
          std::cout << " m * (t)^(w^w^Rc+2w^Rcdot + Rcddot) " << (robotMass_* (kine::skewSymmetric(position) *
                        (skewV2R * positionCom + 2*(skewVR * velocityCom) + orientation * accelerationCom )) ).transpose()<<std::endl;
          std::cout <<  "m*c ^g " <<(robotMass_* kine::skewSymmetric(orientation * positionCom + position) * cst::gravity).transpose()<<std::endl;
          std::cout << "vt "<<vt.transpose()<<std::endl;

          std::cout << "(R*I*RT + m*[R * cl]x)^-1 " <<std::endl << (orientation*Inertia*orientationT + robotMass_*kine::skewSymmetric2(orientation * positionCom)).inverse()<<std::endl;

          std::cout << "c "<<(orientation * positionCom + position).transpose()<<std::endl;
          std::cout << "c^vf "<<(kine::skewSymmetric(orientation * positionCom + position)*vf).transpose()<<std::endl;
          std::cout << "mc^vf "<<(robotMass_*kine::skewSymmetric(orientation * positionCom + position)*vf).transpose()<<std::endl;
          std::cout << "vt-mc^vf "<<(vt - robotMass_*kine::skewSymmetric(orientation * positionCom + position)*vf).transpose()<<std::endl;
          std::cout << "R*cl "<<(orientation*positionCom).transpose()<<std::endl;
          std::cout << "(R*cl)^omegadot "<<(kine::skewSymmetric(orientation*positionCom)*angularAcceleration).transpose()<<std::endl;
          std::cout << "linearAcceleration "<<linearAcceleration.transpose()<<std::endl;
          std::cout << "angularAcceleration "<<angularAcceleration.transpose()<<std::endl;

        }
    }

    void IMUElasticLocalFrameDynamicalSystem::iterateDynamics
             (const Vector3& positionCom, const Vector3& velocityCom,
              const Vector3& accelerationCom, const Vector3& AngMomentum,
              const Vector3& dotAngMomentum,
              const Matrix3& inertia, const Matrix3& dotInertia,
              const IndexedMatrixArray& contactPos,
              const IndexedMatrixArray& contactOri,
              Vector3& position, Vector3& linVelocity, Vector3& linearAcceleration,
              Vector3 &oriVector, Vector3& angularVel, Vector3& angularAcceleration
              )
    {

        Matrix3 orientationFlex(computeRotation_(oriVector));
        //compute new acceleration with the current flex position and velocity and input
        computeAccelerations (positionCom, velocityCom,
        accelerationCom, AngMomentum, dotAngMomentum,
        inertia, dotInertia,  contactPosV_, contactOriV_, position, linVelocity, linearAcceleration,
                       oriVector, orientationFlex, angularVel, angularAcceleration);

        //integrate kinematics with the last acceleration
        integrateKinematics(position, linVelocity, linearAcceleration,
            orientationFlex,angularVel, angularAcceleration, dt_);

        AngleAxis orientationAA(orientationFlex);
        oriVector=orientationAA.angle()*orientationAA.axis();
    }

    Vector IMUElasticLocalFrameDynamicalSystem::stateDynamics
        (const Vector& x, const Vector& u, unsigned k)
    {
        if (kcurrent_!= k )
        {
          std::cout <<std::endl << "k "<<k<<std::endl;
          kcurrent_ = k;
          disp =true;
        }
        else
          disp =false;

        assertStateVector_(x);

        Vector3 positionFlex(x.segment(kine::pos,3));
        Vector3 velocityFlex(x.segment(kine::linVel,3));
        Vector3 accelerationFlex(x.segment(kine::linAcc,3));
        Vector3 orientationFlexV(x.segment(kine::ori,3));
        Vector3 angularVelocityFlex(x.segment(kine::angVel,3));
        Vector3 angularAccelerationFlex(x.segment(kine::angAcc,3));


        const Matrix3 inertia(kine::computeInertiaTensor(u.segment<6>(input::inertia)));
        const Matrix3 dotInertia(kine::computeInertiaTensor(u.segment<6>(input::dotInertia)));

        unsigned nbContacts(getContactsNumber());
        Vector contact(u.tail(6*getContactsNumber()));
        for (int i = 0; i<nbContacts ; ++i)
        {
          contactPosV_.setValue(contact.segment<3>(6*i),i);
          contactOriV_.setValue(computeRotation_(contact.segment<3>(6*i+3)),i);
        }

        Vector3 positionCom(u.segment<3>(input::posCom));
        Vector3 velocityCom(u.segment<3>(input::velCom));
        Vector3 accelerationCom(u.segment<3>(input::accCom));
        Vector3 AngMomentum(u.segment<3>(input::angMoment));
        Vector3 dotAngMomentum(u.segment<3>(input::dotAngMoment));

         if (false)
        {
        std::cout << "accelerationFlex "<<accelerationFlex.transpose()<<std::endl;
        std::cout << "angularAccelerationFlex "<<angularAccelerationFlex.transpose()<<std::endl;
        }



        iterateDynamics (positionCom, velocityCom,
                          accelerationCom, AngMomentum, dotAngMomentum,
                      inertia, dotInertia,  contactPosV_, contactOriV_, positionFlex, velocityFlex, accelerationFlex,
                       orientationFlexV, angularVelocityFlex, angularAccelerationFlex);



//        std::cout << "accelerationFlex" << accelerationFlex.transpose() << std::endl;
//        std::cout << "angularAccelerationFlec" << angularAccelerationFlex.transpose() << std::endl;

        //x_{k+1}
         if (false)
        {
        std::cout << "accelerationFlex "<<accelerationFlex.transpose()<<std::endl;
        std::cout << "angularAccelerationFlex "<<angularAccelerationFlex.transpose()<<std::endl;
        std::cout << "velocityFlex "<<velocityFlex.transpose()<<std::endl;
        std::cout << "angularVelocityFlex "<<angularVelocityFlex.transpose()<<std::endl;
        std::cout << "positionFlex "<<positionFlex.transpose()<<std::endl;
        std::cout << "orientationFlexV "<<orientationFlexV.transpose()<<std::endl;
        }

        Vector xk1(x);

        xk1.segment(kine::pos,3) = positionFlex;
        xk1.segment(kine::linVel,3) = velocityFlex;
        xk1.segment(kine::linAcc,3) = accelerationFlex;



        xk1.segment(kine::ori,3) =  orientationFlexV;
        xk1.segment(kine::angVel,3) = angularVelocityFlex;
        xk1.segment(kine::angAcc,3) = angularAccelerationFlex;

         if (processNoise_!=0x0)
            return processNoise_->addNoise(xk1);
        else
            return xk1;
    }

    inline Matrix3 IMUElasticLocalFrameDynamicalSystem::computeRotation_
                (const Vector3 & x)
    {
        if (orientationVector_!=x)
        {
            orientationVector_ = x;
            curRotation_ = kine::rotationVectorToAngleAxis(x).toRotationMatrix();
        }

        return curRotation_;
    }

    Vector IMUElasticLocalFrameDynamicalSystem::measureDynamics
                (const Vector& x, const Vector& u, unsigned k)
    {
        assertStateVector_(x);

        Vector3 positionFlex(x.segment(kine::pos,3));
        Vector3 velocityFlex(x.segment(kine::linVel,3));
        Vector3 accelerationFlex(x.segment(kine::linAcc,3));
        Vector3 orientationFlexV(x.segment(kine::ori,3));
        Vector3 angularVelocityFlex(x.segment(kine::angVel,3));
        Vector3 angularAccelerationFlex(x.segment(kine::angAcc,3));

        Matrix3 rFlex (computeRotation_(orientationFlexV));

        assertInputVector_(u);

        Vector3 positionControl(u.segment(input::posIMU,3));
        Vector3 velocityControl(u.segment(input::linVelIMU,3));
        Vector3 accelerationControl(u.segment(input::linAccIMU,3));
        Vector3 orientationControlV(u.segment(input::oriIMU,3));
        Vector3 angularVelocityControl(u.segment(input::angVelIMU,3));

        Matrix3 rControl(computeRotation_(orientationControlV));
        Matrix3 r = rFlex * rControl;

//        std::cout << "Mesure x " << x.transpose() << std::endl;
//        std::cout << "Mesure u " << u.transpose() << std::endl;

        // Translation sensor dynamic
        Vector3 acceleration;

        acceleration << 0,0,0;

        acceleration += 2*kine::skewSymmetric(angularVelocityFlex) * rFlex * velocityControl;

//        std::cout << "acceleration " << acceleration.transpose() << std::endl;

        acceleration += accelerationFlex + rFlex * accelerationControl;

//        std::cout << "acceleration " << acceleration.transpose() << std::endl;

        acceleration += (
                    kine::skewSymmetric(angularAccelerationFlex)
                    + tools::square(kine::skewSymmetric(angularVelocityFlex))
                )*rFlex * positionControl;

//        std::cout << "acceleration " << acceleration.transpose() << std::endl;

        // Rotation sensor dynamic
        Vector3 angularVelocity( angularVelocityFlex + rFlex * angularVelocityControl);

        // Set sensor state before measurement
        Vector v(Vector::Zero(15,1));
        v.head<9>() = Eigen::Map<Eigen::Matrix<double, 9, 1> >(&rFlex(0,0));
        v.segment<3>(9)=acceleration;
        v.tail<3>()=angularVelocity;
        sensor_.setState(v,k);

        // Measurement
        Vector y (Matrix::Zero(measurementSize_,1));
        y.head(sensor_.getMeasurementSize()) = sensor_.getMeasurements();

//        std::cout << "===> v state sensor  " << v.transpose() << std::endl;
//        std::cout << "===> y  " << y.transpose() << std::endl;

        return y;
    }

    void IMUElasticLocalFrameDynamicalSystem::setProcessNoise(NoiseBase * n)
    {
        processNoise_=n;
    }

    void IMUElasticLocalFrameDynamicalSystem::resetProcessNoise()
    {
        processNoise_=0x0;
    }

    void IMUElasticLocalFrameDynamicalSystem::setMeasurementNoise
                ( NoiseBase * n)
    {
        sensor_.setNoise(n);
    }

    void IMUElasticLocalFrameDynamicalSystem::resetMeasurementNoise()
    {
        sensor_.resetNoise();
    }

    void IMUElasticLocalFrameDynamicalSystem::setSamplingPeriod(double dt)
    {
        dt_=dt;
    }

    unsigned IMUElasticLocalFrameDynamicalSystem::getStateSize()
    {
        return stateSize_;
    }

    unsigned IMUElasticLocalFrameDynamicalSystem::getInputSize()
    {
        return inputSize_;
    }

    void IMUElasticLocalFrameDynamicalSystem::setInputSize(unsigned i)
    {
        inputSize_=i;
    }

    unsigned IMUElasticLocalFrameDynamicalSystem::getMeasurementSize()
    {
        return measurementSize_;
    }

    NoiseBase * IMUElasticLocalFrameDynamicalSystem::getProcessNoise() const
    {
        return processNoise_;
    }

    NoiseBase * IMUElasticLocalFrameDynamicalSystem::
                                                    getMeasurementNoise() const
    {
        return sensor_.getNoise();
    }

    void IMUElasticLocalFrameDynamicalSystem::setContactsNumber(unsigned i)
    {
        nbContacts_=i;
    }

    void IMUElasticLocalFrameDynamicalSystem::setContactPosition
                                        (unsigned i, const Vector3 & position)
    {
        BOOST_ASSERT( i< contactPositions_.size() &&
                    "ERROR: The index of contact is out of range.");

        contactPositions_[i] = position;
    }

    Vector3 IMUElasticLocalFrameDynamicalSystem::getContactPosition(unsigned i)
    {
        return contactPositions_[i];
    }


    inline unsigned IMUElasticLocalFrameDynamicalSystem::getContactsNumber(void) const
    {
        return  nbContacts_;
    }

    void IMUElasticLocalFrameDynamicalSystem::setKfe(const Matrix3 & m)
    {
        Kfe_=m;
    }

    void IMUElasticLocalFrameDynamicalSystem::setKfv(const Matrix3 & m)
    {
        Kfv_=m;
    }

    void IMUElasticLocalFrameDynamicalSystem::setKte(const Matrix3 & m)
    {
        Kte_=m;
    }

    void IMUElasticLocalFrameDynamicalSystem::setKtv(const Matrix3 & m)
    {
        Ktv_=m;
    }
}
}


