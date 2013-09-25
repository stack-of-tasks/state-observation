/**
 * \file      kalman-filter-base.hpp
 * \author    Mehdi Benallegue
 * \date       2012
 * \brief      Defines the base class of a Kalman filter
 *
 *             It mostly implements the equations of Kalman filtering
 *             It is suitablle by derivation to be used incases of Linear,
 *             linearized and extended Kalman filtering. It may be
 *             derived to unscented Kalman filtering, but non-straighforwardly
 *             because the state vector is modified.
 *
 *             x_{k+1}=f(x_k,u_k)+v_k
 *
 *             y_k=h(x_k,u_k)+w_k
 *
 * \details
 *
 *
 */



#ifndef KALMANFILTERBASEHPP
#define KALMANFILTERBASEHPP

#include <state-observer/zero-delay-observer.hpp>

namespace stateObserver
{


/**
     * \class  KalmanFilterBase
     * \brief
     *        It mostly implements the equations of Kalman filtering
     *        It is suitablle by derivation to be used incases of Linear,
     *        linearized and extended Kalman filtering. It may be
     *        derived to unscented Kalman filtering, but non-straighforwardly
     *        because the state vector is modified. This class requires
     *        to be derived to overload the update routine and the measurements
     *        simulation routine.
     *
     *             x_{k+1}=f(x_k,u_k)+v_k
     *
     *             y_k=h(x_k,u_k)+w_k
     *
     *
     *
     */
    class KalmanFilterBase: public ZeroDelayObserver
    {
    public:

        /// The type of the jacobian df/dx
        typedef Eigen::MatrixXd Amatrix;

        /// The type of the jacobian dh/dx
        typedef Eigen::MatrixXd Cmatrix;

        /// The type of the covariance matrix of the process noise v
        typedef Eigen::MatrixXd Qmatrix;

        /// The type of the covariance matrix of the measurement noise w
        typedef Eigen::MatrixXd Rmatrix;

        /// The type of the covariance matrix of the state estimation error.
        typedef Eigen::MatrixXd Pmatrix;

        /// Default constructor
        KalmanFilterBase(){}

        /// The constructor
        ///  \li n : size of the state vector
        ///  \li m : size of the measurements vector
        ///  \li p : size of the input vector
        KalmanFilterBase(unsigned n,unsigned m,unsigned p=0)
            :ZeroDelayObserver(n,m,p){}



        /// Set the value of the jacobian df/dx
        virtual void setA(const Amatrix& A);

        /// Clear the jacobian df/dx
        virtual void clearA();


        /// Set the value of the Jacobian dh/dx
        virtual void setC(const Cmatrix& C);

        /// Clear the jacobian dh/dx
        virtual void clearC();


        /// Set the measurement noise covariance matrix
        virtual void setR(const Rmatrix& R);

        /// Clear the measurement noise covariance matrix
        virtual void clearR();


        /// Set the process noise covariance matrix
        virtual void setQ(const Qmatrix& Q);

        /// Clear the process noise covariance matrix
        virtual void clearQ();


        /// Set the covariance matrix of the current time state estimation error
        virtual void setStateCovariance(const Pmatrix& P);

        /// Clear the covariace matrix of the current time state estimation
        /// error
        virtual void clearStateCovariance();

        /// Get the covariance matrix of the current time state estimation
        virtual Pmatrix getStateCovariance(unsigned k);

        /// Resets all the observer
        virtual void reset();


        /// Get a matrix having the size of the A matrix having "c" values
        Amatrix getAmatrixConstant(double c) const;

        /// Get a matrix having the size of the A matrix having random values
        Amatrix getAmatrixRandom() const;

        /// Get a matrix having the size of the A matrix having zero values
        Amatrix getAmatrixZero() const;

        /// Get an identity matrix having the size of the A matrix
        Amatrix getAmatrixIdentity() const;

        ///checks whether or not a matrix has the dimensions of the A matrix
        bool checkAmatrix(const Amatrix & ) const;


        /// Get a matrix having the size of the C matrix having "c" values
        Cmatrix getCmatrixConstant(double c) const;

        /// Get a matrix having the size of the C matrix having random values
        Cmatrix getCmatrixRandom() const;

        /// Get a matrix having the size of the C matrix having zero values
        Cmatrix getCmatrixZero() const;

        ///checks whether or not a matrix has the dimensions of the C matrix
        bool checkCmatrix(const Cmatrix &) const;


        /// Get a matrix having the size of the Q matrix having "c" values
        Qmatrix getQmatrixConstant(double c) const;

        /// Get a matrix having the size of the Q matrix having random values
        Qmatrix getQmatrixRandom() const;

        /// Get a matrix having the size of the Q matrix having zero values
        Qmatrix getQmatrixZero() const;

        /// Get an identity matrix having the size of the Q matrix
        Qmatrix getQmatrixIdentity() const;

        ///checks whether or not a matrix has the dimensions of the Q matrix
        bool checkQmatrix(const Qmatrix &) const;


        /// Get a matrix having the size of the R matrix having "c" values
        Rmatrix getRmatrixConstant(double c) const;

        /// Get a matrix having the size of the R matrix having random values
        Rmatrix getRmatrixRandom() const;

        /// Get a matrix having the size of the R matrix having zero values
        Rmatrix getRmatrixZero() const;

        /// Get an identity matrix having the size of the R matrix
        Rmatrix getRmatrixIdentity() const;

        ///checks whether or not a matrix has the dimensions of the R matrix
        bool checkRmatrix(const Rmatrix &) const;


        /// Get a matrix having the size of the P matrix having "c" values
        Pmatrix getPmatrixConstant(double c) const;

        /// Get a matrix having the size of the P matrix having random values
        Pmatrix getPmatrixRandom() const;

        /// Get a matrix having the size of the P matrix having zero values
        Pmatrix getPmatrixZero() const;

        /// Get an identity matrix having the size of the P matrix
        Pmatrix getPmatrixIdentity() const;

        ///checks whether or not a matrix has the dimensions of the P matrix
        bool checkPmatrix(const Pmatrix & ) const;

        ///changes the dimension of the state vector:
        ///resets the internal container for the state vector and
        ///the containers for the matrices A, C, Q, P
        virtual void setStateSize(unsigned n);

        ///changes the dimension of the measurement vector:
        ///resets the internal container for the measurement vectors and
        ///the containers for the matrices C, R
        virtual void setMeasureSize(unsigned m);

    protected:

        /// The type of Kalman gain matrix
        typedef Eigen::MatrixXd Kmatrix;

        /// The Kalman filter loop
        virtual StateVector oneStepEstimation_();

        /// The abstract method to overload to implement f(x,u)
        virtual StateVector prediction_(unsigned k)=0;

        /// The abstract method to overload to implement h(x,u)
        virtual MeasureVector simulateSensor_(const StateVector& x, unsigned k)=0;

        /// Containers for the jacobian matrix of the process
        DiscreteTimeMatrix a_;

        /// Containers for the jacobian matrix of the measurement
        DiscreteTimeMatrix c_;

        /// Container for the process noise covariance matrice
        DiscreteTimeMatrix q_;

        /// Container for the measurement noise covariance matrice
        DiscreteTimeMatrix r_;

        /// Container for the covariance matrix of the estimation error
        DiscreteTimeMatrix pr_;

    };




}

#endif //KALMANFILTERBASEHPP
