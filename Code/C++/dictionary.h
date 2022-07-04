#ifndef dictionary_H
#define dictionary_H

#include <vector>
#include "error.h"

 class dictionary
 {
     string dictToStr;
     const vector<string> entryNames;
     vector<string> entryValues;
     void squeezeString(string&);
     void read(string&);
    public:
     const string fileName;
     dictionary (const vector<string>&);
     dictionary (const dictionary&, string, const vector<string>&);
     void findEntry(string&, const string&, string);
     string get(int) const;
     double getfl(int) const;
     void getSubDict(string&, string) const;
 };

#endif
