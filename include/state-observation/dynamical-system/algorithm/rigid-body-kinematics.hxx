
namespace stateObservation
{
namespace algorithm

{
    inline void RigidBodyKinematics::integrateKinematics
            (Vector3 & position, Vector3 & velocity, const Vector3 & acceleration,
             Quaternion & orientation, Vector3 & rotationVelocityVector,
             const Vector3 & rotationVelocityVectorRate, double dt)
    {
        position +=  dt * velocity + 0.5 * dt * dt * acceleration;
        velocity +=  dt * acceleration;

        orientation = Quaternion( kine::rotationVectorToAngleAxis
                                                    (rotationVelocityVector*dt) )
                            * orientation;

        rotationVelocityVector += dt * rotationVelocityVectorRate;

    }

    inline void RigidBodyKinematics::integrateKinematics
            (Vector3 & position, Vector3 & velocity, const Vector3 & acceleration,
                    Matrix3 & orientation, Vector3 & rotationVelocityVector,
                        const Vector3 & rotationVelocityVectorRate, double dt)
    {
        position +=  dt * velocity + 0.5 * dt * dt * acceleration;
        velocity +=  dt * acceleration;

        orientation = kine::rotationVectorToAngleAxis (rotationVelocityVector*dt + 0.5 * dt * dt*rotationVelocityVectorRate).toRotationMatrix()
                            * orientation;

        rotationVelocityVector += dt * rotationVelocityVectorRate;

    }

    inline void RigidBodyKinematics::integrateConfiguration(Vector3 & position,
          const Vector3 & velocity, Matrix3 & orientation,
          const Vector3 & rotationVelocity, double dt)
    {
        position +=  dt * velocity;

        orientation = kine::rotationVectorToAngleAxis (rotationVelocity*dt).toRotationMatrix()
                            * orientation;
    }
}
}
