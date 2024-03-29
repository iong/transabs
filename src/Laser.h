#ifndef transabs_Laser_h
#define transabs_Laser_h

#include <cstdlib>
#include <cmath>

using namespace std;

class LaserPulse
{
protected:
    double  intensity, frequency, fwhm, zero, phase, fmax;
public:
    LaserPulse(void) {}

    LaserPulse ( double i, double w, double T ) :
        intensity ( i ), frequency ( w ), fwhm ( T ), zero ( 0.0 ), phase(0.5)
    {
        fmax = sqrt ( intensity );
    }

    void init( double i, double w, double T) {
        intensity = i;
        frequency = w ;
        fwhm = T;
        fmax = sqrt ( intensity );
    }

    void setZero(double z) {
        zero = z;
    }

    void setPhase(double p) {
        phase = p;
    }

    virtual double envelope(double t) {
        return fabs(t-zero) < 0.5*fwhm ? 1.0 : 0.0;
    }

    virtual double field ( double t ) {
        double e = envelope ( t );
        return fmax * sin ( frequency * ( t - zero ) + phase*M_PI ) * e;
    }
};

class GaussianLaserPulse : public LaserPulse
{
    double sigma;
public:
    GaussianLaserPulse(void) {}
    
    GaussianLaserPulse ( double i, double w, double T ) : LaserPulse(i, w, T) {
        sigma = T / sqrt(2.0 * log(2.0));
    }

    double envelope(double t)
    {
        double s = (t - zero)/sigma;

        return exp(-s*s);
    }
};

class GaussianKick : public GaussianLaserPulse {
public:
   GaussianKick (double i, double T) : GaussianLaserPulse(i, 0.0, T) {}

   double field ( double t ) {
        double e = envelope ( t );
        return fmax * e;
    }

};

#endif
