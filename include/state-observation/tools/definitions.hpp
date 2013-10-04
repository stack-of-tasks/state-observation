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

#include <vector>
#include <deque>

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

    class CheckedMatrix
    {
    public:
        ///Default constructor
        CheckedMatrix();

        ///A constructor with a given matrix value and a time index
        CheckedMatrix(const Matrix& v);

        ///Set the value of the matrix and the time sample
        inline void set(const Matrix& v);

        ///Get the matrix value
        inline Matrix operator()()const;

        ///Says whether the matrix is initialized or not
        inline const bool & isSet()const;

        ///Switch off the initalization flag, the value is no longer accessible
        inline void reset();

    protected:
        ///Checks whether the matrix is set or not (assert)
        ///does nothing in release mode
        inline void check_() const;

        ///this variable ensures the matrix is initialized,
        bool isSet_;

        Matrix v_;
    };


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
        inline void set(const Matrix& v,unsigned k);

        ///Get the matrix value
        inline Matrix operator()() const;

        ///Get the time index
        inline const unsigned & getTime() const;

        ///Says whether the matrix is initialized or not
        inline const bool & isSet() const;

        ///Switch off the initalization flag, the value is no longer accessible
        inline void reset();

    protected:
        ///Checks whether the matrix is set or not (assert)
        ///does nothing in release mode
        inline void check_() const;

        ///this variable ensures the matrix is initialized,
        bool isSet_;

        unsigned k_;
        Matrix v_;
    };

    class DiscreteTimeArray
    {
        public:
        ///Default constructor
        DiscreteTimeArray();

        inline void pushBack(const Matrix& v,unsigned k);

        inline void popFront();

        inline Matrix operator[](unsigned time) const;

        inline Matrix  & operator[](unsigned time);

        void truncate(unsigned time);


        ///Get the time index
        inline unsigned getLastTime() const;

        ///Get the time index
        inline unsigned getFirstTime() const;

        inline unsigned size() const;

        ///Switch off the initalization flag, the value is no longer accessible
        inline void reset();

        std::vector<Matrix> getArray() const;

        inline bool checkIndex(unsigned k) const;

    protected:
        ///Checks whether the matrix is set or not (assert)
        ///does nothing in release mode
        inline void check_(unsigned time) const;

        ///Checks whether the matrix is set or not (assert)
        ///does nothing in release mode
        inline void check_() const;

        inline void checkNext_(unsigned time) const;

        unsigned k_;

        std::deque<Matrix> v_;

    };

    namespace cst
    {
        const Vector gravity= 9.81 * Eigen::Vector3d::UnitZ();
    }

    #include <state-observation/tools/definitions.hxx>
}

#endif //SENSORSIMULATIONDEFINITIONSHPP
