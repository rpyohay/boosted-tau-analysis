#ifdef __CINT__

#include <vector>

#pragma link off all class;
#pragma link off all function;
#pragma link off all global;
#pragma link off all typedef;

#pragma link C++ class pair<string, string>+;
#pragma link C++ function pair<string, string>::operator=(const pair<string, string>&);
#pragma link C++ class map<string, string>+;
#pragma link C++ class map<string, pair<string, string> >+;
#pragma link C++ class map<pair<string, string>, pair<string, string> >+;
#pragma link C++ function map<pair<string, string>, pair<string, string> >::operator[](const pair<string, string>&);
#pragma link C++ class map<string, vector<string> >+;
#pragma link C++ function std::reverse(vector<string>::iterator, vector<string>::iterator);
#pragma link C++ class pair<short, short>+;
#pragma link C++ class vector<pair<short, short> >+;
#pragma link C++ function vector<pair<short, short> >::push_back(const pair<short, short>&);
#pragma link C++ class vector<vector<string> >+;
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ function vector<vector<float> >::push_back(const vector<float>&);

#endif
