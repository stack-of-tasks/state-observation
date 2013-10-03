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

#include <boost/assert.hpp>
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


    /**
     * \class    DiscreteTimeMatrix
     * \brief    This class describes a structure composed by a matrix
     *           of a given size and a time-index parameter. It can tell also if
     *           it initialized or not.
     *
     *
     */
    class DiscreteTimeMatrix
    {
    public:
        ///Default constructor
        DiscreteTimeMatrix();

        ///A constructor with a given matrix value and a time index
        DiscreteTimeMatrix(const Matrix& v, unsigned k);

        ///Set the value of the matrix and the time sample
        inline void set(const Matrix& v,unsigned k)
        {
            k_=k;
            v_=v;
            isSet_=true;
        }

        ///Get the matrix value
        inline Matrix operator()()const
        {
            check_();
            return v_;
        }

        ///Get the time index
        inline const unsigned & getTime()const
        {
            check_();
            return k_;
        }

        ///Says whether the matrix is initialized or not
        inline const bool & isSet()const
        {
            return isSet_;
        }

        ///Switch off the initalization flag, the value is no longer accessible
        inline void reset()
        {
            k_=0;
            isSet_=false;
        }

    protected:
        ///Checks whether the matrix is set or not (assert)
        ///does nothing in release mode
        inline void check_()const
        {
            BOOST_ASSERT(isSet_ && "Error: Matrix not initialized");
        }

        ///this variable ensures the matrix is initialized,
        bool isSet_;

        unsigned k_;
        Matrix v_;
    };

    namespace cst
    {
        const Vector gravity= 9.81 * Eigen::Vector3d::UnitZ();
    }
}

#endif //SENSORSIMULATIONDEFINITIONSHPP
