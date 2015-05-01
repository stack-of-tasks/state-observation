/**
 * \file      miscellaneous-algorithms.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief      Gathers many kinds of algorithms
 *
 *
 *
 */


#ifndef STATEOBSERVATIONTOOLSMISCELANEOUSALGORITHMS
#define STATEOBSERVATIONTOOLSMISCELANEOUSALGORITHMS

#include <state-observation/tools/definitions.hpp>


namespace stateObservation
{
    namespace tools
    {

        ///computes the square of a value of any type
        template <class T>
        inline T square (const T & x)
        {
            return T(x*x);
        }

        ///derivates any type with finite differences
        template <class T>
        inline T derivate(const T & o1 , const T & o2 , double dt)
        {
            T o(o2-o1);
            return o*(1/dt);
        }

        template <typename T> inline
        int signum(T x)
        {
          return (T(0) < x) - (x < T(0));
        }
    }

    namespace kine
    {

        /// Puts the orientation vector norm between 0 and Pi if its
        /// get close to 2pi
        inline Vector regulateOrientationVector(const Vector3 & v )
        {
            double n2=v.squaredNorm();
            if (n2 > (3./2.) * M_PI * (3./2.) * M_PI )
            {
                double n=sqrt(n2);
                unsigned k =  unsigned(ceil((n - M_PI) / (2*M_PI))) ;
                return (v / n) * ( n - k*2*M_PI );
            }
            else
                return v;
        }

        /// Transform the rotation vector into angle axis
        inline AngleAxis rotationVectorToAngleAxis(const Vector3 & v)
        {
            double angle=v.squaredNorm();
            if (angle > cst::epsilonAngle * cst::epsilonAngle)
            {
                angle=sqrt(angle);
                return AngleAxis(angle, v/angle);
            }
            else
                return AngleAxis(0.0 , Vector3::UnitZ());
        }

        /// Tranbsform the rotation vector into rotation matrix
        inline Matrix3 rotationVectorToRotationMatrix(const Vector3 & v)
        {
                return (rotationVectorToAngleAxis(Vector3(v))).toRotationMatrix();
        }

        /// Tranbsform the rotation matrix into rotation vector
        inline Vector3 rotationMatrixToRotationVector(const Matrix3 & R)
        {
            AngleAxis a;
            a=AngleAxis(Matrix3(R));
            Vector v(3);
            v=a.angle()*a.axis();
            return v;
        }

        ///transform a 3d vector into a skew symmetric 3x3 matrix
        inline Matrix3 skewSymmetric(const Vector3 & v)
        {
            //R <<     0, -v[2],  v[1],
            //      v[2],     0, -v[0],
            //     -v[1],  v[0],     0;


            Matrix3 R;

            R(0,0)=   R(1,1) = R(2,2) = 0.;
            R(0,1)= -( R(1,0)= v[2] );
            R(2,0)= -( R(0,2)= v[1] );
            R(1,2)= -( R(2,1)= v[0] );

            return R;
        }

        ///transform a 3d vector into a squared skew symmetric 3x3 matrix
        inline Matrix3 skewSymmetric2(const Vector3 & v)
        {
            Matrix3 R( v * v.transpose());

            double n = R.trace();

            R(0,0) -= n;
            R(1,1) -= n;
            R(2,2) -= n;

            return R;
        }


        inline Matrix3 computeInertiaTensor(const Vector6 inputInertia, Matrix3& inertiaTensor)
        {

            const double & Ixx=inputInertia[0];
            const double & Iyy=inputInertia[1];
            const double & Izz=inputInertia[2];
            const double & Ixy=inputInertia[3];
            const double & Ixz=inputInertia[4];
            const double & Iyz=inputInertia[5];

            inertiaTensor   <<    Ixx, Ixy, Ixz,
                                  Ixy, Iyy, Iyz,
                                  Ixz, Iyz, Izz;

            return inertiaTensor;

        }

        ///transforms a homogeneous matrix into 6d vector (position theta mu)
        inline Vector6 homogeneousMatrixToVector6(const Matrix4 & M)
        {
            Vector6 v;
            AngleAxis a (AngleAxis(Matrix3(M.block(0,0,3,3))));

            v.head(3) = M.block(0,3,3,1);
            v.tail(3) = a.angle() * a.axis();

            return v;
        }

        ///transforms a 6d vector (position theta mu) into a homogeneous matrix
        inline Matrix4 vector6ToHomogeneousMatrix(const Vector6 & v)
        {
            Matrix4 M(Matrix4::Identity());
            M.block(0,0,3,3) = rotationVectorToAngleAxis(Vector3(v.tail(3)))
                                .toRotationMatrix();
            M.block(0,3,3,1) = v.head(3);
            return M;
        }

        ///transforms a rotation into translation given a constraint of a fixed point
        inline void fixedPointRotationToTranslation
            (const Matrix3 & R, const Vector3 & rotationVelocity,
            const Vector3 & rotationAcceleration, const Vector3 & fixedPoint,
        Vector3 & outputTranslation, Vector3 & outputLinearVelocity,
        Vector3 & outputLinearAcceleration)
        {
            Matrix3 omega_x = skewSymmetric(rotationVelocity);
            outputTranslation = fixedPoint - R * fixedPoint;
            outputLinearVelocity = -omega_x * R * fixedPoint;
            outputLinearAcceleration =-(skewSymmetric(rotationAcceleration) + tools::square(omega_x))
                                           * R * fixedPoint;
        }



        ///derivates a quaternion using finite difference to get a angular velocity vector
        inline Vector3 derivateRotationFD
        (const Quaternion & q1, const Quaternion & q2, double dt)
        {
            AngleAxis aa (q2 * q1.conjugate());

            return (aa.angle()/dt)*aa.axis();
        }

        ///derivates a quaternion using finite difference to get a angular velocity vector
        inline Vector3 derivateRotationFD
        (const Vector3 & o1, const Vector3 & o2, double dt)
        {
            Quaternion q1(rotationVectorToAngleAxis(o1));
            Quaternion q2(rotationVectorToAngleAxis(o2));

            return derivateRotationFD(q1, q2, dt);
        }

        inline Vector6 derivateHomogeneousMatrixFD
                    (const Matrix4 & m1, const Matrix4 & m2, double dt )
        {
            Vector6 out;

            Matrix3 r1 = m1.block(0,0,3,3);
            Matrix3 r2 = m2.block(0,0,3,3);

            AngleAxis aa (r2 * r1.transpose());

            out.tail(3) = (aa.angle()/dt)*aa.axis();

            out.head(3) = (m2.block(0,3,3,1) - m1.block(0,3,3,1))/dt;

            return out;
        }

        inline Vector6 derivatePoseThetaUFD
                    (const Vector6 & v1, const Vector6 & v2, double dt )
        {
            Vector6 out;

            out.tail(3) = derivateRotationFD(v1.tail(3), v2.tail(3), dt);

            out.head(3) = (v2.head(3) - v1.head(3))/dt;

            return out;
        }



        ///uses the derivation to reconstruct the velocities and accelerations given
        ///trajectories in positions and orientations only
        inline IndexedMatrixArray reconstructStateTrajectory
                (const IndexedMatrixArray & positionOrientation,
                double dt)
        {
            Vector r(Vector::Zero(18,1));

            const IndexedMatrixArray & po= positionOrientation;

            unsigned i0=positionOrientation.getFirstIndex();
            unsigned i1=positionOrientation.getLastIndex()+1;

            IndexedMatrixArray a;
            a.setValue(r,i0);
            a.resize(po.size(),r);

            for (unsigned i=i0; i<i1; ++i)
            {
                Vector poi = po[i];

                r.segment(kine::pos,3) = poi.head(3);
                r.segment(kine::ori,3) = poi.tail(3);
                a.setValue(r,i);
            }

            for (unsigned i=i0; i<i1-1; ++i)
            {
                r = a[i];

                Vector poi = po[i];
                Vector poi1 = po[i+1];

                r.segment(linVel,3)  = tools::derivate(Vector3(poi.head(3)),
                                               Vector3(poi1.head(3)), dt);
                r.segment(angVel,3) = derivateRotationFD(
                                    Quaternion(rotationVectorToAngleAxis (poi.tail(3))),
                                    Quaternion(rotationVectorToAngleAxis (poi1.tail(3))),
                                    dt);

                a.setValue(r,i);

            }

            for (unsigned i=i0; i<i1-2; ++i)
            {
                r = a[i];
                Vector r2 = a[i+1];

                r.segment(linAcc,3) = tools::derivate(Vector3(r.segment(linVel,3)),
                                              Vector3(r2.segment(linVel,3)), dt);
                r.segment(angAcc,3) = tools::derivate(Vector3(r.segment(angVel,3)),
                                         Vector3(r2.segment(angVel,3)), dt);

                a.setValue(r,i);

            }

            return a;
        }

        inline Vector invertState( const Vector & state)
        {
            Matrix3 r2 = (rotationVectorToAngleAxis( - state.segment(kine::ori,3))).
                             toRotationMatrix();//inverse
            Vector3 omega1 = state.segment(kine::angVel,3);
            Vector3 omega1dot = state.segment(kine::angAcc,3);
            Matrix3 omega1x = skewSymmetric(omega1);
            Matrix3 omega1dotx = skewSymmetric(omega1dot);

            Vector3 t1 = state.segment(kine::pos,3);
            Vector3 t1dot = state.segment(kine::linVel,3);
            Vector3 t1dotdot = state.segment(kine::linAcc,3);



            Vector state2(Vector::Zero(18,1));
            state2.segment(kine::pos,3)= - r2 * t1 ;   //t2
            state2.segment(kine::linVel,3)= r2 * ( omega1x * t1 - t1dot); //t2dot
            state2.segment(kine::linAcc,3)=  r2 * ( omega1x * (2 * t1dot - omega1x * t1)
                                            - t1dotdot + omega1dotx * t1); //t2dotdot
            state2.segment(kine::ori,3)= -state.segment(kine::ori,3);   //thetaU2
            state2.segment(kine::angVel,3)= -r2 * omega1; //omega2
            state2.segment(kine::angAcc,3)=  r2 * (omega1x * omega1 - omega1dot); //omega2dot
            return state2;

        }

        inline Matrix4 invertHomoMatrix (Matrix4 m)
        {
            Matrix4 m2(Matrix4::Identity());
            Matrix3 rt = m.block(0,0,3,3).transpose();
            m2.block(0,0,3,3) = rt;
            m2.block(0,3,3,1) = - rt * m.block(0,3,3,1);
            return m2;
        }

    }

}



#endif //STATEOBSERVATIONTOOLSMISCELANEOUSALGORITHMS
