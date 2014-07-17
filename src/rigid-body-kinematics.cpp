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
        //position +=  dt * velocity + 0.5 * dt * dt * acceleration;
        position += dt*velocity+0.5*dt*dt*acceleration;
        velocity +=  dt * acceleration;

        //debug
        AngleAxis preori  = AngleAxis(orientation);
        Vector3 previousOri = preori.axis() * preori.angle();
        //eodebug

        orientation = Quaternion( kine::rotationVectorToAngleAxis
                                                    (rotationVelocityVector*dt) )
                            * orientation;



        //debug
        AngleAxis postori  = AngleAxis(orientation);
        Vector3 postOri = postori.axis() * postori.angle();
        Vector3 preOridot = rotationVelocityVector;
        //eofdebug


        rotationVelocityVector += dt * rotationVelocityVectorRate;

       // cout << "preori" << previousOri.transpose()<<endl;
       // cout << "postOri" << postOri.transpose()<<endl;
       // cout << "derivative " << kine::derivateRotationFD(previousOri,postOri,dt).transpose();
       // cout << "omega" << rotationVelocityVector.transpose() << endl;
       // cout <<  "dt * rotationVelocityVectorRate" <<  (dt * rotationVelocityVectorRate).transpose() <<endl;
       // cout << "preOridot " << preOridot.transpose() <<endl;
       // cout << "dt integr " << dt <<endl;

    }
}
}
