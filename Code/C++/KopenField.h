#ifndef KopenField_H
#define KopenField_H

#include "openField.h"

 class KopenField
 :
   public openField
 {
     const string name_;
     double sigma_u;
     double sigma_v;
     double sigma_w;
     double L;
     double propConst;
     double timeStep_;
     double vertTimeStep_;
     double tau_u;
     double tau_v;
     double tau_w;

   public:
     KopenField (const dictionary&);
     void K (double*, const double&);
     void GradK (double*, const double&);
     void updateParams (ifstream*);
     inline string name()
     {
        return name_;
     }

     inline int horTimeStep()
     {
        return timeStep_;
     }

     inline int vertTimeStep()
     {
        return vertTimeStep_;
     }
 };

#endif
