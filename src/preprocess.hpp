#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <iterator>

using namespace std;

class alnData
{
	public:
		alnData();
		void initialize(string speciesFile, string alnFileList);
		void setDelimiter(string delimiter);
		void readTraits(string fileName);
		void printTraits();
		void processFastaFileList(string fileName);
		void readAln(string fileName);
		void processAln();
		void generateResponseFile(string fileName);
		void generateFeatureFile(string fileName);
		void generateMappingFile(string fileName);
		void generateGroupIndicesFile(string baseName);
		void generateMissingFile(string baseName);
		void normalizeFeatures(bool normalize);
		void dropSingletons(bool ignoreSingletons);
		void setCountThreshold(int countThreshold);
		void setUpsampleBalance(bool upsampleBalance);
		void setDownsampleBalance(bool downsampleBalance);
		void setDiskCaching(bool useDiskCache);
		void setIndelFuzzing(bool indelFuzzing);
		void setDataType(string dataType);
		void balanceSample();

	private:
		int featureIndex;
		int geneIndex;
		bool normalize = false;
		bool ignoreSingletons = false;
		bool upsampleBalance = false;
		bool downsampleBalance = false;
		bool useDiskCache = false;
		bool indelFuzzing = false;
		bool caseSensitive = true;
		int countThreshold;
		string dataType = "universal";
		string validChars;
		string currentGene;
		string delimiter;
		vector<string> species;
		vector<string> groups;
		vector<string> missingSeqs;
		vector<vector<int>> groupIndices;
		map<string, int> geneGroupIndex;
		map<string, float> traits;
		map<string, string> seqs;
		map<string, string> validCharSets;
		map<int, string> featureMap;
		map<int, vector<float>> features;
		vector<string> featureCacheFiles;

};
