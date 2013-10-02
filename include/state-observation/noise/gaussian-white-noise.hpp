/**
 * \file      gaussian-white-noise.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief
 *
 *
 *
 */



#ifndef SENSORSIMULATIONGAUSSIANWHITENOISEHPP
#define SENSORSIMULATIONGAUSSIANWHITENOISEHPP

#include <state-observation/noise/noise-base.hpp>




namespace stateObservation
{

    /**
     * \class  GaussienWhiteNoise
     * \brief
     *
     * \details
     *
     */

    class GaussianWhiteNoise : public NoiseBase
    {
    public:
        ~GaussianWhiteNoise(){}

        GaussianWhiteNoise(unsigned dimension);

        GaussianWhiteNoise();

        virtual Vector addNoise(const Vector & state, const Vector & measurement, unsigned timeIndex);

        virtual void setStandardDeviation(const Matrix & std);

        virtual void setCovarianceMatrix(const Matrix & cov);

        virtual void setBias(const Vector & bias);

        virtual void setDimension(unsigned dim_);

        virtual unsigned getDimension() const;

    protected:
        virtual void checkMatrix_(const Matrix & m) const ;

        virtual void checkVector_(const Vector & v) const;

        unsigned dim_;

        Matrix std_;

        Vector bias_;

    };

}

#endif //SENSORSIMULATIONGAUSSIANWHITENOISEHPP
