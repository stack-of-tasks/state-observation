/**
 * \file      observer-base.hpp
 * \author    Mehdi Benallegue
 * \date       2012
 * \brief      Defines the base class of a state
 *             observer.
 *             The observer is destinated to any dynamical system with a vector
 *             state representation
 *             The file describes also
 *             the used data structures and exceptions of the derivated classes
 *
 * \details
 *
 *
 */



#ifndef TEMPLATEOBSERVERBASEHPP
#define TEMPLATEOBSERVERBASEHPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <boost/static_assert.hpp>
#include <boost/assert.hpp>

namespace stateObserver
{
    namespace compileTime
    {
        /**
         * \class    DiscreteTimeMatrix
         * \brief    This class describes a structure composed by a matrix
         *           of a given size and a time-index parameter. It can tell also if
         *           it initialized or not.
         *           It is templated by:
         *           \li r : number of rows
         *           \li c : number of columns
         *           r and c must be positive.
         *
         * \details
         *
         */
        template<unsigned r,unsigned c>
        class DiscreteTimeMatrix
        {
        public:

            ///Definition of matrix type
            typedef Eigen::Matrix<double, r,c> MatrixT;

            ///Default constructor
            DiscreteTimeMatrix();

            ///A constructor with a given matrix value and a time index
            DiscreteTimeMatrix(const MatrixT& v, unsigned k);

            ///Set the value of the matrix and the time sample
            inline void set(const MatrixT& v,unsigned k);

            ///Get the matrix value
            inline MatrixT operator()()const;

            ///Get the time index
            inline const unsigned & getTime()const;

            ///Says whether the matrix is initialized or not
            inline const bool & isSet()const;

            ///Switch off the initalization flag, the value is no longer accessible
            inline void reset();

            ///Special instructions to have a static-sized eigen vector as a member
            enum { NeedsToAlign = (sizeof(MatrixT)%16)==0 };
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)

        protected:

            inline void check_()const;

            ///this variable ensures the matrix is initialized,
            bool isSet_;

            unsigned k_;
            MatrixT v_;
        };


        /**
         * \class  ObserverBase
         * \brief  The base class for observers.
         *         The observer is destinated to any dynamical system with a vector
         *         state representation. This class mostly defined an abstract
         *         interface, static constants and types. It is templated by:
         *
         *         \li n : size of the state vector
         *         \li m : size of the measurements vector
         *         \li p : size of the input vector
         *
         * \details
         *
         */


        template <unsigned n,unsigned m, unsigned p=0>
        class ObserverBase
        {
        public:

            ///stateSize is the size of the state vector
            static unsigned const stateSize=n;

            ///measureSize is the size of measurements vector
            static unsigned const measureSize=m;

            ///inputSize is the size of the input vector
            static unsigned const inputSize=p;


            ///Destructor
            virtual ~ObserverBase(){};


            ///StateVector is the type of state vector
            typedef Eigen::Matrix<double, n,1> StateVector;

            ///MeasureVector is the type of measurements vector
            typedef Eigen::Matrix<double, m,1> MeasureVector;

            ///InputVector is the type of the input vector
            typedef Eigen::Matrix<double, p,1> InputVector;


            ///Set the value of the state vector at time index k
            virtual void setState(const StateVector& x_k,unsigned k)=0;

            ///Remove all the given past values of the state
            virtual void clearState()=0;

            ///Set the value of the measurements vector at time index k
            virtual void setMeasurement(const MeasureVector& x_k,unsigned k)=0;

            ///Remove all the given past values of the measurements
            virtual void clearMeasurements()=0;

            ///Set the value of the input vector at time index k
            virtual void setInput(const InputVector& x_k,unsigned k)=0;

            ///Remove all the given past values of the inputs
            virtual void clearInputs()=0;

            ///Run the observer loop and gets the state estimation of the state at
            ///instant k
            virtual StateVector getEstimateState(unsigned k)=0;

            ///Reinitializes the whole observer
            ///default behavior is to call the three "ObserverBase::clear*" methods
            virtual void reset();

        protected:

            ///Internal (protected) typedefs the timed states.
            typedef DiscreteTimeMatrix<n,1> State;

            ///Internal (protected) typedefs the timed measurements.
            typedef DiscreteTimeMatrix<m,1> Measure;

            ///Internal (protected) typedefs the timed inputs.
            typedef DiscreteTimeMatrix<p,1> Input;
        };

#include <state-observer/template/observer-base.hxx>
    }
}

#endif
