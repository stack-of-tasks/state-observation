#include <state-observation/sensors-simulation/algorithm/linear-acceleration.hpp>

namespace stateObservation
{
    namespace algorithm
    {
        Vector3 LinearAcceleration::accelerationMeasure(const Vector3 & acceleration, const Matrix3 & orientation) const
        {

            Vector3 localAcceleration (acceleration);
            localAcceleration += cst::gravity;
            return  orientation.transpose()*localAcceleration;
        }
    }
}
