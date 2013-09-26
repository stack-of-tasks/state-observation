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



#ifndef OBSERVERBASEHPP
#define OBSERVERBASEHPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <boost/static_assert.hpp>
#include <boost/assert.hpp>

namespace stateObserver
{
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

        ///Definition of matrix type
        typedef Eigen::MatrixXd MatrixT;

        ///Default constructor
        DiscreteTimeMatrix();

        ///A constructor with a given matrix value and a time index
        DiscreteTimeMatrix(const MatrixT& v, unsigned k);


        ///Set the value of the matrix and the time sample
        inline void set(const MatrixT& v,unsigned k)
        {
            k_=k;
            v_=v;
            isSet_=true;
        }

        ///Get the matrix value
        inline MatrixT operator()()const
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
        MatrixT v_;
    };


    /**
     * \class  ObserverBase
     * \brief  The base class for observers.
     *         The observer is destinated to any dynamical system with a vector
     *         state representation. This class mostly defined an abstract
     *         interface, static constants and types. It is templated by:
     *
     *
     */


    class ObserverBase
    {
    public:

        ///StateVector is the type of state vector
        typedef Eigen::VectorXd StateVector;

        ///MeasureVector is the type of measurements vector
        typedef Eigen::VectorXd MeasureVector;

        ///InputVector is the type of the input vector
        typedef Eigen::VectorXd InputVector;


        ///constructor
        ///      \li n : size of the state vector
        ///      \li m : size of the measurements vector
        ///      \li p : size of the input vector
        ObserverBase(unsigned n, unsigned m, unsigned p=0);

        ///default constructor (default values for n,m,p are zero)
        ObserverBase();

        ///Destructor
        virtual ~ObserverBase(){};

        /// Changes the size of the state vector
        virtual void setStateSize(unsigned n);

        /// gets the size of the state vector
        virtual unsigned getStateSize() const;

        /// Changes the size of the measurement vector
        virtual void setMeasureSize(unsigned m);

        /// gets the size of the measurement vector
        virtual unsigned getMeasureSize() const;

        /// Changes the size of the input vector
        virtual void setInputSize(unsigned p);

        /// gets the size of the input vector
        virtual unsigned getInputSize() const;


        ///Set the value of the state vector at time index k
        virtual void setState(const StateVector& x_k,unsigned k)=0;

        ///Remove all the given past values of the state
        virtual void clearStates()=0;

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

        ///Gives a vector of state vector size having duplicated "c" value
        virtual Eigen::VectorXd stateVectorConstant( double c ) const;

        ///Gives a vector of state vector size having random values
        virtual Eigen::VectorXd stateVectorRandom() const;

        ///Gives a vector of state vector size having zero values
        virtual Eigen::VectorXd stateVectorZero() const;

        ///Tells whether or not the vector has the dimensions of a state vector
        virtual bool checkStateVector(const StateVector & v ) const;


        ///Gives a vector of measurement vector size having duplicated "c" value
        virtual Eigen::VectorXd measureVectorConstant( double c ) const;

        ///Gives a vector of measurement vector size having random values
        virtual Eigen::VectorXd measureVectorRandom() const;

        ///Gives a vector of measurement vector size having zero values
        virtual Eigen::VectorXd measureVectorZero() const;

        ///Tells whether or not the vector has the dimensions of a measurement vector
        virtual bool checkMeasureVector(const MeasureVector &) const;


        ///Gives a vector of input vector size having duplicated "c" value
        virtual Eigen::VectorXd inputVectorConstant( double c ) const;

        ///Gives a vector of input vector size having random values
        virtual Eigen::VectorXd inputVectorRandom() const;

        ///Gives a vector of input vector size having zero values
        virtual Eigen::VectorXd inputVectorZero() const;

        ///Tells whether or not the vector has the dimensions of a input vector
        virtual bool checkInputVector(const InputVector &) const;

    protected:


        ///stateSize is the size of the state vector
        unsigned n_;

        ///measureSize is the size of measurements vector
        unsigned m_;

        ///inputSize is the size of the input vector
        unsigned p_;

    };

}

#endif//OBSERVERBASEHPP
