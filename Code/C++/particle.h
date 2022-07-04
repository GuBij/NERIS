#ifndef particle_H 
#define particle_H

#include <iostream>
#include <memory>
#include <fstream>
#include <random>
#include "meteorology.h"

 class particle
 {
   meteorology* wind;
   meteorology* eddyDiff;
   ifstream windParamFile;
   ifstream eddyDiffParamFile;
   ifstream* windparams_ptr;
   ifstream* eddydiffparams_ptr;
   const double xMax;
   const double yMax;
   const double zMax;
   const unsigned mfreq;
   const unsigned ppr;
   unsigned nz;
   int dt;
   int dt_z;
   unsigned timeIter;
   long releaseTime;
   double Qp;
   double scaleDose;
   bool dead;
   int* warmUpPeriod;
   bool warmingUp;
   double location[6];
   long seed;
   mt19937 generator;  // Mersenne twister PRNG, initialized with seed from previous random device instance
   normal_distribution<double> distribution;
   int measurementIndex;
   int nom;
   long simTime;
   double corr[3];
   double up[6];
   double Up[3];
   double Kp[3];
   double GradKp[3];
   bool eos_;
   double phi;
   bool useLangevinEqn;
   double GradSig2_;

  public:
   particle (int,unsigned long, double*, meteorology*, meteorology*, int&, bool);
   bool inside();
   bool alive();
   bool measurement();
   void eqn(); 
   void LangevinEqn();
   void advance (); 
   double reflection(double, double, double, double, double);
   inline bool eos()
   {
        return eos_;
   }
   inline double* coord()
   {
	return location;
   }
   inline void close()
   {
	windparams_ptr->close(); eddydiffparams_ptr->close();
   }
   inline long time()
   {
	return simTime;
   }
   inline double relTime()
   {
	return releaseTime;
   }
   inline double particleCharge()
   {
        return Qp;
   }
   inline double scaleFactorDose()
   {
        return scaleDose;
   }
 };

#endif
