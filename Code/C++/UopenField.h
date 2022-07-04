#ifndef UopenField_H
#define UopenField_H

#include "openField.h"

 class UopenField
 :
   public openField
 {
    const string name_;
    const double pi = acos(-1);
    double uStar;
    double L;
    double heff_;
    double windDirection;
    double elevation;

   public:
    UopenField(const dictionary&);
    void U ( double*, const double& );
    void updateParams ( ifstream* );
    inline string name()
    {
	return name_;
    }
    inline double wdir()
    {
	return windDirection;
    }
    inline double heff()
    {
	return heff_;
    }
 };

#endif
