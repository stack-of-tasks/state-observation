#include <state-observation/sensors-simulation/algorithm/magnetic-field.hpp>

namespace stateObservation
{
    namespace algorithm
    {
        Vector3 MagneticField::earthLocalMagneticField_ = Vector3 (0, 32, -37);

        Vector3 MagneticField::magneticFieldMeasure(const Matrix3 & orientation) const
        {
            return Vector3(orientation.transpose()*earthLocalMagneticField_);
        }
    }
}
