#include <state-observation/dynamical-system/algorithm/rigid-body-kinematics.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>

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

        orientation = Quaternion( tools::kinematics::rotationVectorToAngleAxis
                                                    (rotationVelocityVector*dt) )
                            * orientation;

        rotationVelocityVector += dt * rotationVelocityVectorRate;

    }
}
}
