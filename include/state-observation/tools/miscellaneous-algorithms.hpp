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
            return x*x;
        }

        ///derivates any type with finite differences
        template <class T>
        inline T derivate(const T & o1 , const T & o2 , double dt)
        {
            return (o2-o1)/dt;
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
                unsigned k =  ceil((n - M_PI) / (2*M_PI)) ;
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

        ///transform a 3d vector into a skew symmetric 3x3 matrix
        inline Matrix3 skewSymmetric(const Vector3 & v)
        {
            Matrix3 R ;
            R << 0,     -v[2],    v[1],
            v[2],   0,      -v[0],
            -v[1],   v[0],    0;

            return R;
        }

        ///transforms a homogeneous matrix into 6d vector (position theta mu)
        inline Vector6 homogeneousMatrixToVector6(const Matrix4 & M)
        {
            Vector6 v;
            AngleAxis a = AngleAxis(Matrix3(M.block(0,0,3,3)));

            v.head(3) = M.block(0,3,3,1);
            v.tail(3) = a.angle() * a.axis();

            return v;
        }

        ///transforms a 6d vector (position theta mu) into a homogeneous matrix
        inline Matrix4 vector6ToHomogeneousMatrix(const Vector6 & v)
        {
            Matrix4 M;
            M.block(0,0,3,3) = rotationVectorToAngleAxis(Vector3(v.tail(3)))
                                .toRotationMatrix();
            M.block(0,3,3,1) = v.head(3);
            M(3,3) = 1.0;

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
            AngleAxis aa (q2 * q1.inverse());

            double a=aa.angle();

            Vector3 v =  aa.axis();

            return (aa.angle()/dt)*aa.axis();
        }

        ///uses the derivation to reconstruct the velocities and accelerations given
        ///trajectories in positions and orientations only
        inline DiscreteTimeArray reconstructStateTrajectory
                (const DiscreteTimeArray & positionOrientation,
                double dt)
        {
            Vector r(Vector::Zero(18,1));

            const DiscreteTimeArray & po= positionOrientation;

            unsigned i0=positionOrientation.getFirstTime();
            unsigned i1=positionOrientation.getLastTime()+1;

            DiscreteTimeArray a;
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
