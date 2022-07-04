 #ifndef radiation_H
 #define radiation_H

 #include <string>
 #include <vector>
 #include <utility>
 using namespace std;

 class radiation
 {
    const double E;
    const double A;
    const int noRecords;
    int nor;
    vector<pair<double,double>> location;
    double* dose; 
    double* neighb; 
    double* weightCoeff;
    double* neighbWeightCoeff;
    double L;
    double mu;
    double a1;
    double a2;
    double sigair;
   public:
    radiation(int);
    void update(double*,int,double);
    void updateBayOpt(double*, int, double, int, double);
    void write(string);
    void writeCoeffBayOpt(string);
    static void reconstruct(int,string);
    static void reconstructWeightsBayOpt(int,string);
 };

 #endif
