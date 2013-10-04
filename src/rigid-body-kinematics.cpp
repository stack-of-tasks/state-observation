#include <state-observation/dynamical-system/algorithm/rigid-body-kinematics.hpp>

namespace stateObservation
{
    namespace algorithm

    {
        RigidBodyKinematics::RigidBodyKinematics()
        {
            //ctor
        }

        RigidBodyKinematics::~RigidBodyKinematics()
        {
            //dtor
        }

        void RigidBodyKinematics::integrateKinematics
                (Vector3 & position, Vector3 & velocity, Vector3 & acceleration,
                    Matrix3 & orientation, Vector3 & rotationVelocityVector,
                        Vector3 & rotationVelocityVectorRate, double dt)
        {
            position = position + dt * velocity + 0.5 * dt * dt * acceleration;
            velocity = velocity + dt * acceleration;


        }

    }
}
