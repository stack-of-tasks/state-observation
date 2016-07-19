 /*
 * dynamical-system.cpp
 *
 *  Created on: 19 mai 2014
 *      Author: alexis
 */

#include <state-observation/flexibility-estimation/imu-elastic-local-frame-dynamical-system.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>
#include <stdexcept>

namespace stateObservation
{
namespace flexibilityEstimation
{
    using namespace stateObservation;

   IMUElasticLocalFrameDynamicalSystem::
    	IMUElasticLocalFrameDynamicalSystem(double dt):
        processNoise_(0x0), dt_(dt),
        robotMass_(hrp2::m),
        robotMassInv_(1/hrp2::m),
        measurementSize_(measurementSizeBase_),
        withForceMeasurements_(false), withComBias_(false), withAbsolutePos_(false),
        withUnmodeledMeasurements_(false)
    {
#ifdef STATEOBSERVATION_VERBOUS_CONSTRUCTORS
       // std::cout<<std::endl<<"IMUElasticLocalFrameDynamicalSystem Constructor"<<std::endl;

#endif //STATEOBSERVATION_VERBOUS_CONSTRUCTOR
      Kfe_=40000*Matrix3::Identity();
      Kte_=600*Matrix3::Identity();
      Kfv_=600*Matrix3::Identity();
      Ktv_=60*Matrix3::Identity();

      sensor_.setMatrixMode(true);
      sensor_.concatenateWithInput(6);
      contactModel_=contactModel::none;

      nbContacts_=0;
      inputSize_=42;

      kcurrent_=-1;

      fc_.resize(2*3); fc_.setZero();
      tc_.resize(2*3); tc_.setZero();

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

    double IMUElasticLocalFrameDynamicalSystem::getRobotMass() const
    {
      return robotMass_;
    }


   void IMUElasticLocalFrameDynamicalSystem::setContactModel(unsigned nb)
   {
       contactModel_=nb;
   }

    void IMUElasticLocalFrameDynamicalSystem::computeElastPendulumForcesAndMoments
                              (const IndexedMatrixArray& PeArray,
                               const IndexedMatrixArray& PrArray,
                               const Vector3& position, const Vector3& linVelocity,
                               const Vector3& oriVector, const Matrix3& orientation,
                               const Vector3& angVel,
                               Vector6& forces, Vector6& moments)
    {
        unsigned nbContacts(getContactsNumber());
        fc_.resize(2*3); fc_.setZero();
        tc_.resize(2*3); tc_.setZero();
        forces.setZero();
        moments.setZero();

        Vector3 forcei;
        Vector3 momenti;
        Vector3 globalContactPos;

        Vector3 contactOriUnitVector;

        double stringLength;
        double modifiedStringLength;

        for (unsigned i = 0; i<nbContacts ; ++i)
        {
            globalContactPos = position ;
            globalContactPos.noalias() += orientation*PrArray[i] ;

            stringLength=(PrArray[i]-PeArray[i]).norm();
            modifiedStringLength=(globalContactPos-PeArray[i]).norm();
            contactOriUnitVector= (PeArray[i]-globalContactPos)/modifiedStringLength;

            forcei = (modifiedStringLength-stringLength)*Kfe_*contactOriUnitVector;
            fc_.segment<3>(3*i)= forcei;
            forces.segment<3>(0) += forcei;

            momenti.noalias() = kine::skewSymmetric(globalContactPos)*forcei;
            tc_.segment<3>(3*i)= momenti;
            moments.segment<3>(0) += momenti;
        }

        moments.segment<3>(0).noalias() += - Kte_*oriVector;
        moments.segment<3>(0).noalias() += - Ktv_*angVel;
        forces.segment<3>(0).noalias() += -Kfv_*linVelocity;
    }

    void IMUElasticLocalFrameDynamicalSystem::computeElastPendulumForcesAndMoments1
                              (const IndexedMatrixArray& PeArray,
                               const IndexedMatrixArray& PrArray,
                               const Vector3& position, const Vector3& linVelocity,
                               const Vector3& oriVector, const Matrix3& orientation,
                               const Vector3& angVel,
                               Vector6& forces, Vector6& moments)
    {
        unsigned nbContacts(getContactsNumber());
        fc_.resize(2*3); fc_.setZero();
        tc_.resize(2*3); tc_.setZero();
        forces.setZero();
        moments.setZero();

        stateObservation::Vector3 forcei;
        stateObservation::Vector3 momenti;
        stateObservation::Vector3 u, du;

        double stringLength;
        double modifiedStringLength;
        double lengthRate;

        for (unsigned i = 0; i<nbContacts ; ++i)
        {
            u = position;
            u.noalias() += orientation*PrArray[i] ;
            u.noalias() -= PeArray[i];

            du = linVelocity;
            du.noalias() += kine::skewSymmetric(angVel)*orientation*PrArray[i] ;

            stringLength=(PrArray[i]-PeArray[i]).norm();
            modifiedStringLength=u.norm();

            lengthRate=(modifiedStringLength-stringLength)/modifiedStringLength;

            forcei = -lengthRate*(Kfe_*u+Kfv_*du);
            fc_.segment<3>(3*i)= forcei;
            forces.segment<3>(0) += forcei;

            momenti.noalias() = kine::skewSymmetric(u)*forcei;
            tc_.segment<3>(3*i)= momenti;
            moments.segment<3>(0) += momenti;
        }

        moments.segment<3>(0).noalias() += - Kte_*oriVector;
        moments.segment<3>(0).noalias() += - Ktv_*angVel;
    }

    void IMUElasticLocalFrameDynamicalSystem::computeElastPendulumForcesAndMoments2
                              (const IndexedMatrixArray& PeArray,
                               const IndexedMatrixArray& PrArray,
                               const Vector3& position, const Vector3& linVelocity,
                               const Vector3& oriVector, const Matrix3& orientation,
                               const Vector3& angVel,
                               Vector6& forces, Vector6& moments)
    {
        unsigned nbContacts(getContactsNumber());
        fc_.resize(2*3); fc_.setZero();
        tc_.resize(2*3); tc_.setZero();
        forces.setZero();
        moments.setZero();

        stateObservation::Vector3 forcei;
        stateObservation::Vector3 momenti;
        stateObservation::Vector3 u, du;

        double stringLength;
        double modifiedStringLength;
        double lengthRate;
        unsigned lengthRateVel;

        for (unsigned i = 0; i<nbContacts ; ++i)
        {
            u = position;
            u.noalias() += orientation*PrArray[i] ;
            u.noalias() -= PeArray[i];

            du = linVelocity;
            du.noalias() += kine::skewSymmetric(angVel)*orientation*PrArray[i] ;

            stringLength=(PrArray[i]-PeArray[i]).norm();
            modifiedStringLength=u.norm();

            lengthRate=(modifiedStringLength-stringLength)/modifiedStringLength;
            lengthRateVel=(stringLength/std::pow(modifiedStringLength,3))*du.sum();

            forcei = -Kfe_*lengthRate*u-Kfv_*(lengthRateVel*u+lengthRate*du);
            fc_.segment<3>(3*i)= forcei;
            forces.segment<3>(0) += forcei;

            momenti.noalias() = kine::skewSymmetric(u)*forcei;
            tc_.segment<3>(3*i)= momenti;
            moments.segment<3>(0) += momenti;
        }

        moments.segment<3>(0).noalias() += - Kte_*oriVector;
        moments.segment<3>(0).noalias() += - Ktv_*angVel;
    }


    Vector IMUElasticLocalFrameDynamicalSystem::getForcesAndMoments()
    {
        unsigned nbContacts(getContactsNumber());
        Vector x(6*nbContacts);
        x.setZero();

        for (unsigned int i=0; i<nbContacts; ++i)
        {
          x.segment<3>(6*i) = fc_.segment<3>(3*i);
          x.segment<3>(6*i+3) = tc_.segment<3>(3*i);
        }

        return x;
    }

    Vector IMUElasticLocalFrameDynamicalSystem::getForcesAndMoments(const Vector& x, const Vector& u)
    {
        computeForcesAndMoments(x,u);
        return getForcesAndMoments();
    }

    void IMUElasticLocalFrameDynamicalSystem::computeElastContactForcesAndMoments
                              (const IndexedMatrixArray& contactPosArray,
                               const IndexedMatrixArray& contactOriArray,
                               const IndexedMatrixArray& contactVelArray,
                               const IndexedMatrixArray& contactAngVelArray,
                               const Vector3& position, const Vector3& linVelocity,
                               const Vector3& oriVector, const Matrix3& orientation,
                               const Vector3& angVel,
                               Vector6& fc, Vector6& tc)
    {
        unsigned nbContacts(getContactsNumber());
        fc_.resize(2*3); fc_.setZero();
        tc_.resize(2*3); tc_.setZero();

      for (unsigned i = 0; i<nbContacts ; ++i)
      {
        op_.contactPos = contactPosArray[i];
        op_.contactVel = contactVelArray[i];

        op_.Rci = computeRotation_(contactOriArray[i],i+2);
        op_.Rcit.noalias()= op_.Rci.transpose();

        op_.RciContactPos.noalias()= orientation*op_.contactPos;

        op_.globalContactPos = position;
        op_.globalContactPos.noalias() += op_.RciContactPos ;

        op_.forcei.noalias() = - op_.Rci*Kfe_*op_.Rcit*(op_.globalContactPos-op_.contactPos);
        op_.forcei.noalias() += - op_.Rci*Kfv_*op_.Rcit*(kine::skewSymmetric(angVel)*op_.RciContactPos
                                  +linVelocity + (orientation+stateObservation::Matrix3::Identity())*op_.contactVel);

        fc_.segment<3>(3*i)= op_.forcei;

        op_.momenti.noalias() = -op_.Rci*Kte_*op_.Rcit*oriVector;
        op_.momenti.noalias() += -op_.Rci*Ktv_*op_.Rcit*angVel;

        tc_.segment<3>(3*i)= op_.momenti + kine::skewSymmetric(op_.globalContactPos)*op_.forcei;
      }

      fc=fc_;
      tc=tc_;
    }


    inline void IMUElasticLocalFrameDynamicalSystem::computeForcesAndMoments
                          (const IndexedMatrixArray& contactpos,
                           const IndexedMatrixArray& contactori,
                           const IndexedMatrixArray& contactvel,
                           const IndexedMatrixArray& contactangvel,
                           const Vector3& position, const Vector3& linVelocity,
                           const Vector3& oriVector, const Matrix3& orientation,
                           const Vector3& angVel,
                           Vector6& fc, Vector6& tc)
    {
        switch(contactModel_){
        case contactModel::elasticContact : computeElastContactForcesAndMoments(contactpos, contactori, contactvel, contactangvel, position, linVelocity, oriVector, orientation, angVel, fc, tc);
                 break;

        case contactModel::pendulum : computeElastPendulumForcesAndMoments(contactpos, contactori, position, linVelocity, oriVector, orientation, angVel, fc, tc);
                 break;

        case contactModel::pendulum1 : computeElastPendulumForcesAndMoments1(contactpos, contactori, position, linVelocity, oriVector, orientation, angVel, fc, tc);
                 break;

        case contactModel::pendulum2 : computeElastPendulumForcesAndMoments2(contactpos, contactori, position, linVelocity, oriVector, orientation, angVel, fc, tc);
                 break;

        default: throw std::invalid_argument("IMUElasticLocalFrameDynamicalSystem: the contact model is incorrectly set.");
        }
    }

    inline void IMUElasticLocalFrameDynamicalSystem::computeForcesAndMoments (const Vector& x,const Vector& u)
    {

        op_.positionFlex=x.segment(state::pos,3);
        op_.velocityFlex=x.segment(state::linVel,3);
        op_.orientationFlexV=x.segment(state::ori,3);
        op_.angularVelocityFlex=x.segment(state::angVel,3);
        op_.rFlex = computeRotation_(op_.orientationFlexV,0);

        unsigned nbContacts(getContactsNumber());

        for (unsigned i = 0; i<nbContacts ; ++i)
        {
            // Positions
            op_.contactPosV.setValue(u.segment<3>(input::contacts + 12*i),i);
            op_.contactOriV.setValue(u.segment<3>(input::contacts +12*i+3),i);

            // Velocities
            op_.contactVelArray.setValue(u.segment<3>(input::contacts + 12*i +6),i);
            op_.contactAngVelArray.setValue(u.segment<3>(input::contacts +12*i +9),i);
        }

        computeForcesAndMoments (op_.contactPosV, op_.contactOriV, op_.contactVelArray, op_.contactAngVelArray,
                          op_.positionFlex, op_.velocityFlex, op_.orientationFlexV, op_.rFlex,
                             op_.angularVelocityFlex, fc_, tc_);
    }

    void IMUElasticLocalFrameDynamicalSystem::computeAccelerations
       (const Vector3& positionCom, const Vector3& velocityCom,
        const Vector3& accelerationCom, const Vector3& AngMomentum,
        const Vector3& dotAngMomentum,
        const Matrix3& Inertia, const Matrix3& dotInertia,
        const IndexedMatrixArray& contactPosV,
        const IndexedMatrixArray& contactOriV,
        const Vector3& position, const Vector3& linVelocity, Vector3& linearAcceleration,
        const Vector3 &oriVector ,const Matrix3& orientation,
        const Vector3& angularVel, Vector3& angularAcceleration,
        Vector6& fc, Vector6& tc,
        Vector3 & fm, Vector3& tm)
    {

        kine::skewSymmetric(angularVel,op_.skewV);
        kine::skewSymmetric2(angularVel,op_.skewV2);
        op_.skewVR.noalias()=op_.skewV * orientation;
        op_.skewV2R.noalias()=op_.skewV2 * orientation;

        op_.rFlexT.noalias()=orientation.transpose();

        op_.f=orientation*fm;
        op_.t=orientation*tm+kine::skewSymmetric(position)*orientation*fm;

        for (unsigned i = 0; i<getContactsNumber() ; ++i)
        {
            op_.globalContactPos.noalias() = orientation*contactPosV[i] + position ;
            op_.f+= fc.segment<3>(i*3);
            op_.t+= tc.segment<3>(i*3)+kine::skewSymmetric(op_.globalContactPos)*fc.segment<3>(i*3);
        }    

        op_.wx2Rc.noalias()=op_.skewV2R*positionCom;
        op_._2wxRv.noalias()=2*op_.skewVR*velocityCom;
        op_.Ra.noalias() = orientation*accelerationCom;
        op_.Rc.noalias() = orientation * positionCom;
        op_.Rcp.noalias() = op_.Rc+position;
        op_.RIRT.noalias() = orientation*Inertia*op_.rFlexT;

        op_.vf.noalias() =robotMassInv_*op_.f;

        op_.vf.noalias() -= op_.Ra;
        op_.vf.noalias() -= op_._2wxRv;
        op_.vf.noalias() -= op_.wx2Rc;
        op_.vf.noalias() -= cst::gravity;

        op_.vt =op_.t;
        op_.vt.noalias() -= op_.skewV * (op_.RIRT * angularVel);
        op_.vt.noalias() -= orientation * (dotInertia * (op_.rFlexT * angularVel)+dotAngMomentum) ;
        op_.vt.noalias() -= op_.skewVR * AngMomentum;
        op_.vt.noalias() -= robotMass_* (kine::skewSymmetric(position) *
                            (op_.wx2Rc+ op_._2wxRv + op_.Ra ));
        op_.vt.noalias() -= robotMass_* kine::skewSymmetric(op_.Rcp) * cst::gravity;

        op_.orinertia.noalias() = op_.RIRT + robotMass_*kine::skewSymmetric2(op_.Rc);
        op_.invinertia.compute(op_.orinertia);

        angularAcceleration.noalias() = op_.vt - robotMass_*op_.Rcp.cross(op_.vf);

        op_.invinertia.matrixL().solveInPlace(angularAcceleration);
        op_.invinertia.matrixL().transpose().solveInPlace(angularAcceleration);

        //angularAcceleration = ( op_.RIRT + robotMass_*kine::skewSymmetric2(op_.Rc)).llt().solve(
        //                        (op_.vt - robotMass_*kine::skewSymmetric(op_.Rcp)*op_.vf));

        linearAcceleration = op_.vf;
        linearAcceleration += kine::skewSymmetric(op_.Rc)*angularAcceleration;
    }

    void IMUElasticLocalFrameDynamicalSystem::iterateDynamicsEuler
             (const Vector3& positionCom, const Vector3& velocityCom,
              const Vector3& accelerationCom, const Vector3& AngMomentum,
              const Vector3& dotAngMomentum,
              const Matrix3& inertia, const Matrix3& dotInertia,
              const IndexedMatrixArray& contactPos,
              const IndexedMatrixArray& contactOri,
              Vector3& position, Vector3& linVelocity, Vector6& fc,
              Vector3 &oriVector, Vector3& angularVel, Vector6& tc,
              Vector3 & fm, Vector3& tm,
              double dt)
    {

        op_.rFlex = computeRotation_(oriVector,0);
        //compute new acceleration with the current flex position and velocity and input

//        computeForcesAndMoments (contactPos, contactOri,
//                          position, linVelocity, oriVector, op_.rFlex,
//                             angularVel, fc, tc);
        
        computeAccelerations (positionCom, velocityCom,
        accelerationCom, AngMomentum, dotAngMomentum,
        inertia, dotInertia,  contactPos, contactOri, position, linVelocity, op_.linearAcceleration,
                       oriVector, op_.rFlex, angularVel, op_.angularAcceleration,
                       fc, tc, fm,tm);

        //integrate kinematics with the last acceleration
        integrateKinematics(position, linVelocity, op_.linearAcceleration,
            op_.rFlex,angularVel, op_.angularAcceleration, dt);

        op_.orientationAA=op_.rFlex;
        oriVector.noalias()=op_.orientationAA.angle()*op_.orientationAA.axis();

        computeForcesAndMoments (contactPos, contactOri, op_.contactVelArray, op_.contactAngVelArray,
                          position, linVelocity, oriVector, op_.rFlex,
                             angularVel, fc, tc);
    }

    void IMUElasticLocalFrameDynamicalSystem::iterateDynamicsRK4
             (const Vector3& positionCom, const Vector3& velocityCom,
              const Vector3& accelerationCom, const Vector3& AngMomentum,
              const Vector3& dotAngMomentum,
              const Matrix3& inertia, const Matrix3& dotInertia,
              const IndexedMatrixArray& contactPos,
              const IndexedMatrixArray& contactOri,
              Vector3& position, Vector3& linVelocity, Vector6& fc,
              Vector3 &oriVector, Vector3& angularVel, Vector6& tc,
              Vector3 & fm, Vector3& tm,
              double dt)
    {

        Matrix3 orientationFlex(computeRotation_(oriVector,0));

        Vector3 linAcc1;
        Vector3 angAcc1;


        Vector3 pos2;
        Vector3 oriv2;
        AngleAxis oriaa2;
        Matrix3 ori2;

        Vector3 linVelocity2;
        Vector3 angVelocity2;

        Vector3 linAcc2;
        Vector3 angAcc2;


        Vector3 pos3;
        Vector3 oriv3;
        AngleAxis oriaa3;
        Matrix3 ori3;

        Vector3 linVelocity3;
        Vector3 angVelocity3;

        Vector3 linAcc3;
        Vector3 angAcc3;


        Vector3 pos4;
        Vector3 oriv4;
        AngleAxis oriaa4;
        Matrix3 ori4;

        Vector3 linVelocity4;
        Vector3 angVelocity4;

        Vector3 linAcc4;
        Vector3 angAcc4;

        //////////1st/////////////    
        computeAccelerations (positionCom, velocityCom,
                        accelerationCom, AngMomentum, dotAngMomentum,
                        inertia, dotInertia,  contactPos, contactOri,
                        position, linVelocity, linAcc1, oriVector,
                        orientationFlex, angularVel, angAcc1,
                        fc, tc, fm, tm);


        //////////2nd//////////////
        pos2 = position;
        ori2 = orientationFlex;

        integrateConfiguration(pos2, linVelocity, ori2, angularVel, dt/2);

        oriaa2.fromRotationMatrix(ori2);
        oriv2 = oriaa2.angle()*oriaa2.axis();

        linVelocity2 = linVelocity + (dt/2)*linAcc1;
        angVelocity2 = angularVel + (dt/2)*angAcc1;


        computeAccelerations (positionCom, velocityCom,
                        accelerationCom, AngMomentum, dotAngMomentum,
                        inertia, dotInertia,  contactPos, contactOri,
                        pos2, linVelocity2, linAcc2, oriv2,
                        ori2, angVelocity2, angAcc2,
                        fc, tc, fm, tm);

        ////////////3rd/////////////

        pos3 = position;
        ori3 = orientationFlex;

        integrateConfiguration( pos3, linVelocity2, ori3, angVelocity2, dt/2);

        oriaa3.fromRotationMatrix(ori3);
        oriv3 = oriaa3.angle()*oriaa3.axis();

        linVelocity3 = linVelocity + (dt/2)*linAcc2;
        angVelocity3 = angularVel + (dt/2)*angAcc2;

        computeAccelerations (positionCom, velocityCom,
                        accelerationCom, AngMomentum, dotAngMomentum,
                        inertia, dotInertia,  contactPos, contactOri,
                        pos3, linVelocity3, linAcc3, oriv3,
                        ori3, angVelocity3, angAcc3,
                        tc, fc, fm, tm);


        ////////////4th/////////////
        pos4 = position;
        ori4 = orientationFlex;

        integrateConfiguration( pos4, linVelocity3, ori4, angVelocity3, dt);

        oriaa4.fromRotationMatrix(ori4);
        oriv4 = oriaa4.angle()*oriaa4.axis();

        linVelocity4 = linVelocity + (dt)*linAcc3;
        angVelocity4 = angularVel + (dt)*angAcc3;

        computeAccelerations (positionCom, velocityCom,
                        accelerationCom, AngMomentum, dotAngMomentum,
                        inertia, dotInertia,  contactPos, contactOri,
                        pos4, linVelocity4, linAcc4, oriv4,
                        ori4, angVelocity4, angAcc4,
                        fc, tc, fm, tm);

        /////////////////////////////

        integrateConfiguration( position, linVelocity, orientationFlex, angularVel, dt/6);
        integrateConfiguration( position, linVelocity2, orientationFlex, angVelocity2, dt/3);
        integrateConfiguration( position, linVelocity3, orientationFlex, angVelocity3, dt/3);
        integrateConfiguration( position, linVelocity4, orientationFlex, angVelocity4, dt/6);

        linVelocity+= (dt/6)*(linAcc1+2*linAcc2+2*linAcc3+linAcc4);
        angularVel+= (dt/6)*(angAcc1+2*angAcc2+2*angAcc3+angAcc4);

        AngleAxis orientationAA(orientationFlex);
        oriVector=orientationAA.angle()*orientationAA.axis();

        // Getting forces
        op_.rFlex = computeRotation_(oriVector,0);
        computeForcesAndMoments (contactPos, contactOri, op_.contactVelArray, op_.contactAngVelArray,
                          position, linVelocity, oriVector, op_.rFlex,
                             angularVel, fc, tc);
    }

    Vector IMUElasticLocalFrameDynamicalSystem::stateDynamics
        (const Vector& x, const Vector& u, unsigned )
    {
        assertStateVector_(x);
        assertInputVector_(u);

        xk_=x;
        uk_=u;

        op_.positionFlex=x.segment(state::pos,3);
        op_.orientationFlexV=x.segment(state::ori,3);
        op_.velocityFlex=x.segment(state::linVel,3);
        op_.angularVelocityFlex=x.segment(state::angVel,3);
        op_.fc1=x.segment(state::fc1,6);
        op_.fc2=x.segment(state::fc2,6);
        op_.positionComBias <<  x.segment(state::comBias,2),
                                0;// the bias of the com along the z axis is assumed 0.
        op_.fm=x.segment(state::unmodeledForces,3);
        op_.tm=x.segment(state::unmodeledForces+3,3);

        kine::computeInertiaTensor(u.segment<6>(input::inertia),op_.inertia);
        kine::computeInertiaTensor(u.segment<6>(input::dotInertia),op_.dotInertia);

        unsigned nbContacts(getContactsNumber());

        for (unsigned i = 0; i<nbContacts ; ++i)
        {
          // Positions
          op_.contactPosV.setValue(u.segment<3>(input::contacts + 12*i),i);
          op_.contactOriV.setValue(u.segment<3>(input::contacts +12*i+3),i);

          // Velocities
          op_.contactVelArray.setValue(u.segment<3>(input::contacts + 12*i +6),i);
          op_.contactAngVelArray.setValue(u.segment<3>(input::contacts +12*i +9),i);

        }

        op_.positionCom=u.segment<3>(input::posCom);
        if(withComBias_) op_.positionCom+=op_.positionComBias;
        op_.velocityCom=u.segment<3>(input::velCom);
        op_.accelerationCom=u.segment<3>(input::accCom);
        op_.AngMomentum=u.segment<3>(input::angMoment);
        op_.dotAngMomentum=u.segment<3>(input::dotAngMoment);

        op_.fc << op_.fc1.segment<3>(0),
                  op_.fc2.segment<3>(0);
        op_.tc << op_.fc1.segment<3>(3),
                  op_.fc2.segment<3>(3);

        int subsample=1;
        for (int i=0; i<subsample; ++i)
        {
          iterateDynamicsEuler (
                          op_.positionCom, op_.velocityCom,
                          op_.accelerationCom, op_.AngMomentum, op_.dotAngMomentum,
                          op_.inertia, op_.dotInertia,  op_.contactPosV, op_.contactOriV,
                          op_.positionFlex, op_.velocityFlex, op_.fc,
                          op_.orientationFlexV, op_.angularVelocityFlex,
                          op_.tc, op_.fm, op_.tm,
                          dt_/subsample);
        }
        //x_{k+1}

        op_.fc1 << op_.fc.segment<3>(0),
                   op_.tc.segment<3>(0);
        op_.fc2 << op_.fc.segment<3>(3),
                   op_.tc.segment<3>(3);

        xk1_=x;

        xk1_.segment<3>(state::pos) = op_.positionFlex;
        xk1_.segment<3>(state::ori) =  op_.orientationFlexV;

        xk1_.segment<3>(state::linVel) = op_.velocityFlex;
        xk1_.segment<3>(state::angVel) = op_.angularVelocityFlex;

        xk1_.segment<6>(state::fc1) = op_.fc1;
        xk1_.segment<6>(state::fc2) = op_.fc2;

        // xk1_.segment<2>(state::comBias) = op_.positionComBias.head<2>();

        xk1_.segment<3>(state::unmodeledForces) = 0.99*op_.fm;
        xk1_.segment<3>(state::unmodeledForces+3) = 0.99*op_.tm;

        if (processNoise_!=0x0)
            return processNoise_->addNoise(op_.xk1);
        else
            return xk1_;
    }


    inline Matrix3& IMUElasticLocalFrameDynamicalSystem::computeRotation_
                (const Vector3 & x, int i)
    {
        Vector3 & oriV = op_.orientationVector(i);
        Matrix3 & oriR = op_.curRotation(i);
        if (oriV!=x)
        {
            oriV = x;
            oriR = kine::rotationVectorToAngleAxis(x).toRotationMatrix();
        }

        return oriR;
    }

    Vector IMUElasticLocalFrameDynamicalSystem::measureDynamics
                (const Vector& x, const Vector& u, unsigned k)
    {
        assertStateVector_(x);
        assertInputVector_(u);

        xk_fory_=x;
        uk_fory_=u;
        op_.k_fory=k;

        op_.positionFlex=x.segment(state::pos,3);
        op_.velocityFlex=x.segment(state::linVel,3);
        op_.orientationFlexV=x.segment(state::ori,3);
        op_.angularVelocityFlex=x.segment(state::angVel,3);
        op_.fc1=x.segment(state::fc1,6);
        op_.fc2=x.segment(state::fc2,6);
        op_.fm=x.segment(state::unmodeledForces,3);
        op_.tm=x.segment(state::unmodeledForces+3,3);
        op_.drift=x.segment<3>(state::drift);

        op_.rFlex =computeRotation_(op_.orientationFlexV,0);

        kine::computeInertiaTensor(u.segment<6>(input::inertia),op_.inertia);
        kine::computeInertiaTensor(u.segment<6>(input::dotInertia),op_.dotInertia);
        unsigned nbContacts(getContactsNumber());
        for (unsigned i = 0; i<nbContacts ; ++i)
        {
            // Positions
            op_.contactPosV.setValue(u.segment<3>(input::contacts + 12*i),i);
            op_.contactOriV.setValue(u.segment<3>(input::contacts +12*i+3),i);
        }
        op_.positionCom=u.segment<3>(input::posCom);
        if(withComBias_) op_.positionCom+=op_.positionComBias;
        op_.velocityCom=u.segment<3>(input::velCom);
        op_.accelerationCom=u.segment<3>(input::accCom);
        op_.AngMomentum=u.segment<3>(input::angMoment);
        op_.dotAngMomentum=u.segment<3>(input::dotAngMoment);

        op_.positionControl=u.segment(input::posIMU,3);
        op_.velocityControl=u.segment(input::linVelIMU,3);
        op_.accelerationControl=u.segment(input::linAccIMU,3);
        op_.orientationControlV=u.segment(input::oriIMU,3);
        op_.angularVelocityControl=u.segment(input::angVelIMU,3);

        op_.rControl=computeRotation_(op_.orientationControlV,1);

        op_.rimu = op_.rFlex * op_.rControl;

        op_.fc << op_.fc1.segment<3>(0),
                  op_.fc2.segment<3>(0);
        op_.tc << op_.fc1.segment<3>(3),
                  op_.fc2.segment<3>(3);

        // Get acceleration
        computeAccelerations (op_.positionCom, op_.velocityCom,
        op_.accelerationCom, op_.AngMomentum, op_.dotAngMomentum,
        op_.inertia, op_.dotInertia,  op_.contactPosV, op_.contactOriV, op_.positionFlex, op_.velocityFlex, op_.linearAcceleration,
                       op_.orientationFlexV, op_.rFlex, op_.angularVelocityFlex, op_.angularAcceleration,
                       op_.fc, op_.tc, op_.fm,op_.tm);

        // Translation sensor dynamic
        op_.imuAcc = 2*kine::skewSymmetric(op_.angularVelocityFlex) * op_.rFlex * op_.velocityControl;
        op_.imuAcc += op_.accelerationFlex + op_.rFlex * op_.accelerationControl;
        op_.imuAcc += (
                    kine::skewSymmetric(op_.angularAccelerationFlex)
                    + tools::square(kine::skewSymmetric(op_.angularVelocityFlex))
                )*op_.rFlex * op_.positionControl;

        // Rotation sensor dynamic
        op_.imuOmega = op_.angularVelocityFlex + op_.rFlex * op_.angularVelocityControl;

        // Set sensor state before measurement

        op_.sensorState.resize(sensor_.getStateSize());
        op_.sensorState.head<9>() = Eigen::Map<Eigen::Matrix<double, 9, 1> >(&op_.rimu(0,0));
        op_.sensorState.segment<3>(9)=op_.imuAcc;
        op_.sensorState.segment<3>(12)=op_.imuOmega;

        index_=15;

        if(withUnmodeledMeasurements_)
        {
          op_.sensorState.segment<6>(index_) << op_.fm,
                                                op_.tm;
          index_+=6;
        }

        if (withForceMeasurements_)
        {

          for (unsigned int i=0; i<nbContacts_; ++i)
          {
            op_.sensorState.segment(index_,6) = x.segment(state::fc1+i*6,6);
            index_+=6;
          }
        }

        if (withAbsolutePos_)
        {
          op_.cy=cos(op_.drift(2));
          op_.sy=sin(op_.drift(2));
          op_.rdrift<< op_.cy, -op_.sy, 0,
                       op_.sy,  op_.cy, 0,
                       0,       0,      1;
          op_.rtotal.noalias() = op_.rdrift*op_.rFlex;

          op_.ptotal << op_.drift(0),op_.drift(1),0;
          op_.ptotal.noalias() += op_.rdrift*op_.positionFlex;

          op_.aatotal = op_.rtotal;
          op_.oritotal.noalias() = op_.aatotal.angle()*op_.aatotal.axis();

          op_.sensorState.segment<3>(index_) = op_.ptotal;
          op_.sensorState.segment<3>(index_+3) = op_.oritotal;
        }

        sensor_.setState(op_.sensorState,k);

        //measurements
        yk_=sensor_.getMeasurements();
        return yk_;
    }

    stateObservation::Matrix IMUElasticLocalFrameDynamicalSystem::measureDynamicsJacobian()
    {
      op_.Jy.resize(getMeasurementSize(),getStateSize());
      op_.Jy.setZero();

      op_.xdx = xk_fory_;
      op_.xk_fory = xk_fory_;
      op_.yk = yk_;

      unsigned sizeAfterComBias=getStateSize()-state::comBias-2;

      for (unsigned i=0; i<state::unmodeledForces; ++i)
      {
        op_.xdx[i]+= dx_[i];

        op_.ykdy=measureDynamics(op_.xdx,uk_fory_, op_.k_fory);
        op_.ykdy-=op_.yk;
        op_.ykdy/=dx_[i];

        op_.Jy.col(i)=op_.ykdy;
        op_.xdx[i]=op_.xk_fory[i];
      }

      if(!withUnmodeledMeasurements_)
      {
          op_.Jy.block(0,state::unmodeledForces,getMeasurementSize(),6).setZero();
      }
      else
      {
          for (unsigned i=state::unmodeledForces; i<state::unmodeledForces+6; ++i)
          {
            op_.xdx[i]+= dx_[i];

            op_.ykdy=measureDynamics(op_.xdx,uk_fory_, op_.k_fory);
            op_.ykdy-=op_.yk;
            op_.ykdy/=dx_[i];

            op_.Jy.col(i)=op_.ykdy;
            op_.xdx[i]=op_.xk_fory[i];
          }
      }

      if(!withComBias_)
      {
          op_.Jy.block(0,state::comBias,getMeasurementSize(),2).setZero();
      }
      else
      {
        for (unsigned i=state::comBias; i<state::comBias+2; ++i)
        {
          op_.xdx[i]+= dx_[i];

          op_.ykdy=measureDynamics(op_.xdx,uk_fory_, op_.k_fory);
          op_.ykdy-=op_.yk;
          op_.ykdy/=dx_[i];

          op_.Jy.col(i)=op_.ykdy;
          op_.xdx[i]=op_.xk_fory[i];
        }
      }

      if (!withAbsolutePos_)
      {
         op_.Jy.block(0,state::drift,getMeasurementSize(),3).setZero();
      }
      else
      {
        for (unsigned i=state::comBias+2; i<state::comBias+2+sizeAfterComBias; ++i)
        {
          op_.xdx[i]+= dx_[i];

          op_.ykdy=measureDynamics(op_.xdx,uk_fory_, op_.k_fory);
          op_.ykdy-=op_.yk;
          op_.ykdy/=dx_[i];

          op_.Jy.col(i)=op_.ykdy;
          op_.xdx[i]=op_.xk_fory[i];
        }
      }

      //std::cout << "JACOBIAN: "<<std::endl;
      //std::cout << op_.Jy<<std::endl;

      xk_fory_ = op_.xk_fory; //last thing to do before returning to make the prediction correct

      yk_ = op_.yk;

      return op_.Jy;
    }

    stateObservation::Matrix IMUElasticLocalFrameDynamicalSystem::stateDynamicsJacobian()
    {
      op_.Jx.resize(getStateSize(),getStateSize());
      op_.Jx.setZero();

      op_.xdx = xk_;
      op_.xk = xk_;
      op_.xk1 = xk1_;

      unsigned sizeAfterComBias=getStateSize()-state::comBias-2;
      unsigned sizeAfterUnmodeledForces=getStateSize()-state::unmodeledForces-6;

      for (unsigned i=0; i<state::unmodeledForces; ++i)
      {
        op_.xdx[i]+= dx_[i];

        op_.xk1dx=stateDynamics(op_.xdx,uk_, 0);
        op_.xk1dx-=op_.xk1;
        op_.xk1dx/=dx_[i];

        op_.Jx.col(i)=op_.xk1dx;
        op_.xdx[i]=op_.xk[i];
      }

      if(!withUnmodeledMeasurements_)
      {
          op_.Jx.block(state::unmodeledForces,state::unmodeledForces,6,6).setIdentity();
          op_.Jx.block(0,state::unmodeledForces,state::unmodeledForces,6).setZero();
          op_.Jx.block(state::unmodeledForces+6,state::unmodeledForces,sizeAfterUnmodeledForces,6).setZero();
      }
      else
      {
        for (unsigned i=state::unmodeledForces; i<state::unmodeledForces+6; ++i)
        {
          op_.xdx[i]+= dx_[i];

          op_.xk1dx=stateDynamics(op_.xdx,uk_, 0);
          op_.xk1dx-=op_.xk1;
          op_.xk1dx/=dx_[i];

          op_.Jx.col(i)=op_.xk1dx;
          op_.xdx[i]=op_.xk[i];
        }
      }

      if(!withComBias_)
      {
          op_.Jx.block(state::comBias,state::comBias,2,2).setIdentity();
          op_.Jx.block(0,state::comBias,state::comBias,2).setZero();
          op_.Jx.block(state::comBias+2,state::comBias,sizeAfterComBias,2).setZero();
      }
      else
      {
        for (unsigned i=state::comBias; i<state::comBias+2; ++i)
        {
          op_.xdx[i]+= dx_[i];

          op_.xk1dx=stateDynamics(op_.xdx,uk_, 0);
          op_.xk1dx-=op_.xk1;
          op_.xk1dx/=dx_[i];

          op_.Jx.col(i)=op_.xk1dx;
          op_.xdx[i]=op_.xk[i];
        }
      }

      if (!withAbsolutePos_)
      {
         op_.Jx.block(0,state::drift,getStateSize(),3).setZero();
         op_.Jx.block(state::drift,state::drift,3,3).setIdentity();
      }
      else
      {
          for (unsigned i=state::drift; i<state::drift+3; ++i)
        {
          op_.xdx[i]+= dx_[i];

          op_.xk1dx=stateDynamics(op_.xdx,uk_, 0);
          op_.xk1dx-=op_.xk1;
          op_.xk1dx/=dx_[i];

          op_.Jx.col(i)=op_.xk1dx;
          op_.xdx[i]=op_.xk[i];
        }
      }

      //std::cout << "JACOBIAN: "<<std::endl;
      //std::cout << op_.Jx <<std::endl;

      xk_= op_.xk; //last thing to do before returning
      xk1_ = op_.xk1;

      return op_.Jx;
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

    unsigned IMUElasticLocalFrameDynamicalSystem::getStateSize() const
    {
        return stateSize_;
    }

    unsigned IMUElasticLocalFrameDynamicalSystem::getInputSize() const
    {
        return inputSize_;
    }

    void IMUElasticLocalFrameDynamicalSystem::setInputSize(unsigned i)
    {
        inputSize_=i;
    }

    unsigned IMUElasticLocalFrameDynamicalSystem::getMeasurementSize() const
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
        inputSize_ = 42+12*i;

        updateMeasurementSize_();
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

    void IMUElasticLocalFrameDynamicalSystem::updateMeasurementSize_()
    {
      measurementSize_=measurementSizeBase_;

      if(withUnmodeledMeasurements_)
      {
        measurementSize_+=6; // For unmodeled measurements
      }

      if (withForceMeasurements_)
      {
        measurementSize_+=nbContacts_*6;
      }

      if (withAbsolutePos_)
      {
        measurementSize_+=6;
      }

      sensor_.concatenateWithInput(measurementSize_-measurementSizeBase_);

    }

    void IMUElasticLocalFrameDynamicalSystem::setWithForceMeasurements(bool b)
    {
      withForceMeasurements_=b;
      updateMeasurementSize_();
    }

    void IMUElasticLocalFrameDynamicalSystem::setWithComBias(bool b)
    {
      withComBias_=b;
    }

    void IMUElasticLocalFrameDynamicalSystem::setWithAbsolutePosition(bool b)
    {
      withAbsolutePos_=b;
      updateMeasurementSize_();

    }

    void IMUElasticLocalFrameDynamicalSystem::setWithUnmodeledMeasurements(bool b)
    {
      withUnmodeledMeasurements_=b;
      updateMeasurementSize_();

    }

    bool IMUElasticLocalFrameDynamicalSystem::getWithForceMeasurements() const
    {
      return withForceMeasurements_;
    }

    bool IMUElasticLocalFrameDynamicalSystem::getWithComBias() const
    {
      return withComBias_;
    }

    bool IMUElasticLocalFrameDynamicalSystem::getWithAbsolutePosition() const
    {
      return withAbsolutePos_;
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

    void  IMUElasticLocalFrameDynamicalSystem::setFDstep(const stateObservation::Vector & dx)
    {
      dx_ = dx;
    }
}
}


