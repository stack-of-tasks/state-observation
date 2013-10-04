#ifndef StATEOBSERVATIONRIGIDBODYKINEMATICS_H
#define StATEOBSERVATIONRIGIDBODYKINEMATICS_H

#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
    namespace algorithm
    {
        class RigidBodyKinematics
        {
        public:
            RigidBodyKinematics();
            virtual ~RigidBodyKinematics();

            virtual void integrateKinematics
                (Vector3 & position, Vector3 & velocity, Vector3 & acceleration,
                    Matrix3 & orientation, Vector3 & rotationVelocityVector,
                        Vector3 & rotationVelocityVectorRate, double dt);

        protected:
        private:
        };
    }
}

#endif // StATEOBSERVATIONRIGIDBODYKINEMATICS_H
