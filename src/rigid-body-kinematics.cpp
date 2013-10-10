#include <state-observation/dynamical-system/algorithm/rigid-body-kinematics.hpp>

namespace stateObservation
{
    namespace algorithm

    {
        RigidBodyKinematics::~RigidBodyKinematics()
        {
            //dtor
        }

        void RigidBodyKinematics::integrateKinematics
                (Vector3 & position, Vector3 & velocity, const Vector3 & acceleration,
                    Quaternion & orientation, Vector3 & rotationVelocityVector,
                        const Vector3 & rotationVelocityVectorRate, double dt)
        {
            position +=  dt * velocity + 0.5 * dt * dt * acceleration;
            velocity +=  dt * acceleration;

            double norm=rotationVelocityVector.norm();

            Vector3 a;

            if (norm>cst::epsilonAngle)
                a = rotationVelocityVector / norm;
            else
            {
                a = Vector3::UnitX();
                norm =0;
            }



            AngleAxis rot( norm*dt, a);

            orientation = Quaternion(rot) * orientation;

            rotationVelocityVector += dt * rotationVelocityVectorRate;

        }
    }
}
