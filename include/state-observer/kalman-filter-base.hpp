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
     *         \li n : size of the state vector
     *         \li m : size of the measurements vector
     *         \li p : size of the input vector
     *
     * \details
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

        KalmanFilterBase(){}

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


        /// The type of the jacobian df/dx
        Amatrix getAmatrixConstant(double c) const;

        Amatrix getAmatrixRandom() const;

        Amatrix getAmatrixZero() const;

        Amatrix getAmatrixIdentity() const;

        bool checkAmatrix(const Amatrix & ) const;


        /// The type of the jacobian dh/dx
        Cmatrix getCmatrixConstant(double c) const;

        Cmatrix getCmatrixRandom() const;

        Cmatrix getCmatrixZero() const;

        bool checkCmatrix(const Cmatrix &) const;

        /// The type of the covariance matrix of the process noise v
        Qmatrix getQmatrixConstant(double c) const;

        Qmatrix getQmatrixRandom() const;

        Qmatrix getQmatrixZero() const;

        Qmatrix getQmatrixIdentity() const;

        bool checkQmatrix(const Qmatrix &) const;

        /// The type of the covariance matrix of the measurement noise w
        Rmatrix getRmatrixConstant(double c) const;

        Rmatrix getRmatrixRandom() const;

        Rmatrix getRmatrixZero() const;

        Rmatrix getRmatrixIdentity() const;

        bool checkRmatrix(const Rmatrix &) const;

        /// The type of the covariance matrix of the state estimation error.
        Pmatrix getPmatrixConstant(double c) const;

        Pmatrix getPmatrixRandom() const;

        Pmatrix getPmatrixZero() const;

        Pmatrix getPmatrixIdentity() const;

        bool checkPmatrix(const Pmatrix & ) const;


        virtual void setStateSize(unsigned n);

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


#include <state-observer/kalman-filter-base.hxx>

}

#endif //KALMANFILTERBASEHPP
