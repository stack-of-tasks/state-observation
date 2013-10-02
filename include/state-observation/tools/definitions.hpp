/**
 * \file     definitions.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief
 *
 * \details
 *
 *
 */



#ifndef SENSORSIMULATIONDEFINITIONSHPP
#define SENSORSIMULATIONDEFINITIONSHPP

#include <Eigen/Core>
#include <Eigen/Geometry>


namespace stateObservation
{

    typedef Eigen::VectorXd Vector;

    typedef Eigen::Vector3d Vector3;

    typedef Eigen::Vector4d Vector4;

    typedef Eigen::MatrixXd Matrix;

    typedef Eigen::Matrix3d Matrix3;

    typedef Eigen::Quaterniond Quaternion;

    typedef Eigen::AngleAxis<double> AngleAxis;

    namespace cst
    {
        const Vector gravity= 9.81 * Eigen::Vector3d::UnitZ();
    }
}

#endif //SENSORSIMULATIONDEFINITIONSHPP
