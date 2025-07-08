#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

std::map<std::string, std::string> processSlepOpts(std::string filename);
std::vector<std::tuple<double, double>> readLambdaList(std::string filename);
std::vector<std::tuple<std::string, std::string>> readLambdaListAsStrings(std::string filename);
std::string lambdaLabel(const double arr[], int idx);
std::string lambdaLabel_str(std::string lambda);
std::string writeModelToXMLStream(std::string);