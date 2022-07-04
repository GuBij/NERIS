#ifndef CONCENTRATION_H
#define CONCENTRATION_H

#include <vector>
#include <list>
#include <string>
#include <tuple>
#include "meteorology.h"

 class concentration
 :
 public dictionary
 {
    double* driftPos; 
    double varPos;
    list<double> distrSigma; 
    const unsigned long ppr;
    const double mfreq;
    const double xMax;
    const double yMax;
    const double zMax;
    const double Xmin;
    const double Xmax;
    const double Ymin;
    const double Ymax;
    const double Zmin;
    const double Zmax;
    const double resX;
    const double resY;
    const double resZ;
    const string pfileName;
    const string ofileName;
    double weight(string);
    double meanDrift[3];
    double* position;
    vector<double> bandWidth;
    double currRelTime;
    unsigned noRecords;
    unsigned long nop;
    int fac(int);
    double integralKernel(double,double);

   public:
    concentration(const dictionary&, meteorology*);
    void prepare(double*,double,double,unsigned long);
    void calculate(double,double,double);
    static void reconstruct(int,string,int,int,int,double,double,double,double,double,double);
 };

#endif

