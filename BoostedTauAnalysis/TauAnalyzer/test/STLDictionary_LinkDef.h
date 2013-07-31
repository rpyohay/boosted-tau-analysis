#ifdef __CINT__

#pragma link off all class;
#pragma link off all function;
#pragma link off all global;
#pragma link off all typedef;

#pragma link C++ class pair<string, string>+;
#pragma link C++ function pair<string, string>::operator=(const pair<string, string>&);
#pragma link C++ class map<string, string>+;
#pragma link C++ class map<string, pair<string, string> >+;
#pragma link C++ class map<string, vector<string> >+;

#endif
