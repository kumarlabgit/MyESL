#include <set>
#include <argparse.hpp>
#include <armadillo>
#include <iomanip>
#include "gl_logisticr.hpp"
#include "slep_utils.hpp"

using namespace std;
using namespace arma;


int main(int argc, char *argv[]) {
  argparse::ArgumentParser program("sg_lasso");

  program.add_argument("-f", "--features")
    .required()
    .help("Specify the input features file.");

  program.add_argument("-n", "--groups")
    .required()
    .help("Specify the group indices file.");

  program.add_argument("-r", "--response")
    .required()
    .help("Specify the response file.");

  program.add_argument("-w", "--output")
    .required()
    .help("specify the output file.");

  program.add_argument("-s", "--slep")
    .default_value(std::string("-"))
    .help("Specify a file of key/value SLEP options.");

  program.add_argument("-x", "--xval")
    .default_value(std::string("-"))
    .help("Specify a file of cross-validation partition indices.");

  program.add_argument("-l", "--lambda_list")
    .default_value(std::string("-"))
    .help("Specify a file of lambda value pairs.");

  program.add_argument("-z", "--lambda1")
    .default_value(0.1)
    .help("Specify individual feature sparsity.")
    .scan<'g', double>();

  program.add_argument("-y", "--lambda2")
    .default_value(0.1)
    .help("Specify group feature sparsity.")
    .scan<'g', double>();

  program.add_argument("--model_format")
    .default_value(std::string("-"))
    .help("Specify an output model format of either xml or flat (defaults to xml)");

  program.add_argument("-c", "--gene_count_threshold")
    .default_value(0)
    .help("Specify gene selection cutoff for grid-based model building.")
    .scan<'i', int>();

  try {
    program.parse_args(argc, argv);
  }
  catch (const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    std::exit(1);
  }

  double lambda[2] = {program.get<double>("lambda1"), program.get<double>("lambda2")};
  bool omit_zeroes;
  string model_ext;
  fmat features;
  mat opts_ind;
  frowvec responses;

  int auto_cancel_threshold = program.get<int>("gene_count_threshold");

  std::ifstream feature_map_filestream;

  if (program.get<std::string>("model_format") == "flat"){
    omit_zeroes = true;
    model_ext = ".txt";
  } else {
    omit_zeroes = false;
    model_ext = ".xml";
  }

  if (omit_zeroes)
  {
    std::string feature_map_file = program.get<std::string>("features");
    std::string from = "feature_";
    std::string to = "feature_mapping_";
    size_t startPos = feature_map_file.find(from);
    if(startPos != std::string::npos)
    {
      feature_map_file.replace(startPos, from.length(), to);
    }
    feature_map_filestream.open(feature_map_file);
    if (!feature_map_filestream)
    {
  	  throw std::invalid_argument("\nFeature mapping file(" + feature_map_file + ") not detected.\n");
    }
  }

  bool features_transposed = true;

  if (features_transposed)
  {
	  fmat features_t;
	  features_t.load(csv_name(program.get<std::string>("features"),csv_opts::trans));
	  //features_t.load(csv_name(program.get<std::string>("features"),csv_opts::trans));
	  features = features_t.t();

	  //features = features_t.t();
	  //features = features_t;
  }
  else
  {
	  //features.load(csv_name(program.get<std::string>("features"),csv_opts::semicolon));
      features.load(csv_name(program.get<std::string>("features"),csv_opts::trans));
  }

  //responses.load(csv_name(program.get<std::string>("response"),csv_opts::semicolon));
  responses.load(csv_name(program.get<std::string>("response"),csv_opts::trans));

  if (responses.n_cols != features.n_rows)
  {
    //Log::Fatal << "The responses must have the same number of columns as the feature set." << endl;
    throw std::invalid_argument("\nThe responses must have the same number of columns as the feature set.\n");
  }

  opts_ind.load(csv_name(program.get<std::string>("groups"),csv_opts::trans));

  if (program.get<std::string>("xval") != "-") {
    std::cout << "Performing cross validation." << std::endl;
    //std::vector<std::tuple<std::string, int>> xval_idxs = readXValFile(program.get<std::string>("xval"));
    rowvec xval_idxs;
    xval_idxs.load(csv_name(program.get<std::string>("xval"),csv_opts::trans));
    std::set<int> xval_ids;

    for(const auto& item : xval_idxs) {
      xval_ids.insert(item);  // Insert the value into the set (automatically handles uniqueness)
    }
//    xval_idxs.print("");


    for(const auto& xval_id : xval_ids) {
	  std::cout << "Performing cross validation " << xval_id << std::endl;
      RRLogisticR* slepModel = new RRLogisticR(features, responses, opts_ind, lambda, processSlepOpts(program.get<std::string>("slep")), xval_idxs, xval_id);
      ofstream fileStream(program.get<std::string>("output") + "_xval_" + std::to_string(xval_id) + model_ext);
      if (fileStream.is_open())
      {
        //fileStream << sgl->writeModelToXMLStream();
        if (omit_zeroes)
        {
          slepModel->writeSparseMappedWeightsToStream(fileStream, feature_map_filestream);
        } else {
          slepModel->writeModelToXMLStream(fileStream);
        }
        fileStream.close();
      } else {
        std::cout << "Could not open output file for writing." << std::endl;
      }
      delete slepModel;
    }

    //return 0;
  }

  if (program.get<std::string>("lambda_list") != "-") {
    std::vector<std::tuple<double, double>> lambda_list = readLambdaList(program.get<std::string>("lambda_list"));
    std::vector<std::tuple<std::string, std::string>> lambda_list_str = readLambdaListAsStrings(program.get<std::string>("lambda_list"));
    double max_glambda2 = 1.0;
    double min_glambda2 = 1.0;

//    for (const auto& item : lambda_list) {
//	  double glambda[2] = {std::get<0>(item), std::get<1>(item)};
	for (const auto& item : lambda_list_str) {
	  std::string glambda_str[2] = {std::get<0>(item), std::get<1>(item)};
	  double glambda[2] = {std::stod(glambda_str[0]), std::stod(glambda_str[1])};
	  if (min_glambda2 > glambda[1]) {
        min_glambda2 = glambda[1];
      }
	  std::cout << glambda[0] << " - " << glambda[1] << std::endl;
	  if (glambda[1]>max_glambda2) {
        std::cout << "Skipping this lambda value due to gene count thresholding..." << std::endl;
        continue;
      }
	  RRLogisticR* slepModel = new RRLogisticR(features, responses, opts_ind, glambda, processSlepOpts(program.get<std::string>("slep")));

      //TODO: make out filename reflect lambda pair
	  ofstream fileStream(program.get<std::string>("output") + lambdaLabel_str(glambda_str[0]) + lambdaLabel_str(glambda_str[1]) + model_ext);
	  if (fileStream.is_open())
	  {
        if (omit_zeroes)
        {
          slepModel->writeSparseMappedWeightsToStream(fileStream, feature_map_filestream);
        } else {
          slepModel->writeModelToXMLStream(fileStream);
        }
        fileStream.close();
      } else {
        std::cout << "Could not open output file for writing." << std::endl;
      }

      std::cout << "Non-zero gene count: " << slepModel->NonZeroGeneCount() << std::endl;
      if (slepModel->NonZeroGeneCount()<=auto_cancel_threshold){
        if (glambda[1] == min_glambda2) {
          std::cout << "Skipping all further lambda pairs due to gene count thresholding..." << std::endl;
          break;
        }
        max_glambda2 = glambda[1];
      }
      delete slepModel;
    }

    return 0;

  }

//  SGLasso* sgl = new SGLasso(features, responses, opts_ind, lambda, processSlepOpts(program.get<std::string>("slep")));
  RRLogisticR* slepModel = new RRLogisticR(features, responses, opts_ind, lambda, processSlepOpts(program.get<std::string>("slep")));

//std::cout << "9 ############### " << std::endl;

  //std::cout << slepModel->writeModelToXMLStream();
  ofstream fileStream(program.get<std::string>("output") + model_ext);
  if (fileStream.is_open())
  {
    //fileStream << slepModel->writeModelToXMLStream();
    if (omit_zeroes)
    {
	  slepModel->writeSparseMappedWeightsToStream(fileStream, feature_map_filestream);
    } else {
      slepModel->writeModelToXMLStream(fileStream);
    }
    fileStream.close();
  } else {
    std::cout << "Could not open output file for writing." << std::endl;
  }
  delete slepModel;

  return 0;

}