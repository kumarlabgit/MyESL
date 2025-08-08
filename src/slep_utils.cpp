#include "slep_utils.hpp"


map<string, string> processSlepOpts(string filename)
{
  map<string, string> slep_opts;
  std::cout << "Processing SLEP options file: " << filename << "..." << std::endl;

  string line;
  ifstream optsFile (filename);

  if (optsFile.is_open())
  {

    int splitpos;
    string opt_key;
    while (getline(optsFile, line))
    {
      splitpos = line.find("\t");
      if (splitpos != std::string::npos)
      {
        opt_key = line.substr(0, line.find("\t"));
        slep_opts[opt_key] = line.substr(line.find("\t")+1, std::string::npos);
      }
    }
    optsFile.close();
  }

  return slep_opts;
}


std::vector<std::tuple<double, double>> readLambdaList(string filename)
{
  std::vector<std::tuple<double, double>> data;
  std::string line;
  std::ifstream file(filename);

  if (!file) {
    std::cerr << "Unable to open the file." << std::endl;
    return data;
  }



  while (std::getline(file, line)) {
    std::istringstream iss(line);
	float val1, val2;

	if (iss >> val1 && iss >> val2) {
	  data.push_back(std::make_tuple(val1, val2));
    }
  }

  return data;
}


std::vector<std::tuple<std::string, std::string>> readLambdaListAsStrings(string filename)
{
  std::vector<std::tuple<std::string, std::string>> data;
  std::string first, second;
  std::ifstream file(filename);

  if (!file) {
    std::cerr << "Unable to open the file." << std::endl;
    return data;
  }
  while (file.good()) {
    if (std::getline(file, first, '\t') && std::getline(file, second)) {
      data.emplace_back(first, second);
    }
  }

  file.close();

  return data;
}


std::vector<std::tuple<std::string, int>> readXValFile(string filename)
{
  std::vector<std::tuple<std::string, int>> data;
  std::string line;
  std::ifstream file(filename);

  if (!file) {
    std::cerr << "Unable to open the file." << std::endl;
    return data;
  }



  while (std::getline(file, line)) {
    std::istringstream iss(line);
	std:string val1;
	int val2;

	if (iss >> val1 && iss >> val2) {
	  data.push_back(std::make_tuple(val1, val2));
    }
  }

  return data;
}


std::string lambdaLabel(const double arr[], int idx) {


    // Convert the double to a string
    std::ostringstream streamObj;
    streamObj << std::noshowpoint << std::fixed << std::setprecision(14) << arr[idx];
    std::string doubleStr = streamObj.str();
    std::cout << doubleStr << std::endl;
    // Find "0." at the beginning of the string and erase it if present
    if(doubleStr.substr(0, 2) == "0.") {
        doubleStr.erase(0, 2);
    }

    size_t lastNonZero = doubleStr.rfind('0', doubleStr.size()-1);
    if (lastNonZero == std::string::npos) {  // no zeroes found
        return "_" + doubleStr;
    }

    size_t pos = lastNonZero;
    while (pos > 0 && doubleStr[pos-1] == '0') {
        --pos;
    }
    std::cout << "_" + doubleStr.substr(0, pos) << std::endl;
    return "_" + doubleStr.substr(0, pos);
}


std::string lambdaLabel_str(std::string lambda) {
    if (lambda.substr(0, 2) == "0.") {
        lambda.erase(0, 2);
    }
    return "_" + lambda;
}