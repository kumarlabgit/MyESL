#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <iterator>
#include <mutex>
#include <condition_variable>
#include <algorithm>
#include <filesystem>
#include <sstream>
#include <random>
#include <thread>
#include <numeric>
#include <chrono>


using namespace std;

class Semaphore {
public:
    Semaphore(int initial) : value(initial) {}

    void wait() {
        std::unique_lock<std::mutex> lock(mutex);
        condition.wait(lock, [this]() { return value > 0; });
        --value;
    }

    void signal() {
        std::unique_lock<std::mutex> lock(mutex);
        ++value;
        condition.notify_one();
    }

private:
    int value;
    std::mutex mutex;
    std::condition_variable condition;
};

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
		void readTable(string fileName);
		void processAln();
		void processTable();
		void dumpSpeciesMulti(int i, Semaphore& semaphore);
		void dumpSpecies(int i);
		void dumpFeatureCache();
		void generateResponseFile(string fileName);
		void generateFeatureFile(string fileName);
		void generateStatsFile(string fileName);
		void generateMappingFile(string fileName);
		void generateGroupIndicesFile(string baseName);
		void generateMissingFile(string baseName);
		void normalizeFeatures(bool normalize);
		void dropSingletons(bool ignoreSingletons);
		void setCountThreshold(int countThreshold);
		void setUpsampleBalance(bool upsampleBalance);
		void setDownsampleBalance(bool downsampleBalance);
		void setDiskCaching(bool useDiskCache);
		void setCacheDir(string dirName);
		void setIndelFuzzing(bool indelFuzzing);
		void setFlatWeights(bool flatWeights);
		void setDataType(string dataType);
		void setThreads(int threadCount);
		void balanceSample();

	private:
		int featureIndex = 1;
		int cacheFeatureIndex = 1;
		int geneIndex;
		bool normalize = false;
		bool ignoreSingletons = false;
		bool upsampleBalance = false;
		bool downsampleBalance = false;
		bool useDiskCache = false;
		bool indelFuzzing = false;
		bool caseSensitive = true;
		bool numericInput = false;
		bool numericHeaders = false;
		bool flatWeights = false;
		char inputDelimiter = '\t';
		int threads = 1;
		int countThreshold;
		int featureSize = 4; // Size of each feature cell in bytes
		unsigned int workingMemLimit = 4294967294; // Maximum size of features in cache before dumping to disk, default 4GB
		//unsigned int workingMemLimit = 4294967; // 4MB for testing
		string dataType = "universal";
		string validChars;
		string currentGene;
		string outputDelimiter = ",";
		string cacheDir = ".";
		vector<string> species;
		vector<string> groups;
		vector<string> missingSeqs;
		vector<vector<int>> groupIndices;
		map<string, int> geneGroupIndex;
		map<string, float> traits;
		map<string, string> seqs;
		map<string, vector<float>> numericSeqs;
		map<string, string> validCharSets;
		map<int, string> featureMap;
		map<int, vector<float>> features;
		vector<string> featureCacheFiles;

};



