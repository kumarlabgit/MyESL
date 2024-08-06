#include "preprocess.hpp"

using namespace std;

void trim(string& s)
{
   size_t p = s.find_first_not_of(" \t\r\n");
   s.erase(0, p);

   p = s.find_last_not_of(" \t\r\n");
   if (string::npos != p)
      s.erase(p+1);
}

std::random_device rd;     // only used once to initialise (seed) engine
std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
int randint(int min, int max)
{
	std::uniform_int_distribution<int> uni(min,max); // guaranteed unbiased
	return uni(rng);
}

alnData::alnData()
{
	this->normalize = false;
	this->useDiskCache = false;
	this->validCharSets["nucleotide"] = "ATCGU";
	this->validCharSets["protein"] = "ACDEFGHIKLMNPQRSTVWY";
	this->validCharSets["molecular"] = "ACDEFGHIKLMNPQRSTVWYU";
}

void alnData::initialize(string speciesFile, string alnFileList)
{
	//Read species traits/response file
	this->readTraits(speciesFile);
	this->balanceSample();
	//this->printTraits();

	cout << "Processing FASTA files list: " << alnFileList << "..." << endl;
	//Read and process FASTA files
	this->processFastaFileList(alnFileList);
}

void alnData::setDelimiter(string delimiter)
{
	this->outputDelimiter = delimiter;
}

void alnData::setCacheDir(string dirName)
{
	this->cacheDir = dirName;
}

void alnData::normalizeFeatures(bool normalize)
{
	this->normalize = normalize;
}

void alnData::dropSingletons(bool ignoreSingletons)
{
	this->ignoreSingletons = ignoreSingletons;
	this->countThreshold = 2;
}

void alnData::setCountThreshold(int countThreshold)
{
	this->countThreshold = countThreshold;
}

void alnData::setUpsampleBalance(bool upsampleBalance)
{
	this->upsampleBalance = upsampleBalance;
}

void alnData::setDownsampleBalance(bool downsampleBalance)
{
	this->downsampleBalance = downsampleBalance;
}

void alnData::setDiskCaching(bool useDiskCache)
{
	this->useDiskCache = useDiskCache;
}

void alnData::setIndelFuzzing(bool indelFuzzing)
{
	this->indelFuzzing = indelFuzzing;
}

void alnData::setDataType(string dataType)
{
	this->dataType = dataType;
	if (dataType == "nucleotide" || dataType == "protein" || dataType == "molecular")
	{
		cout << "Setting dataType to " << dataType << " alignment." << endl;
		this->caseSensitive = false;
		this->validChars = this->validCharSets[dataType];
	}
	else if (dataType == "numeric")
	{
		cout << "Setting dataType to " << dataType << " table/matrix." << endl;
		this->numericInput = true;
	}
	else if (dataType != "universal")
	{
		cout << "Unrecognized dataType " << dataType << ", defaulting to universal dataType setting." << endl;
		this->dataType = "universal";
	}

}

void alnData::setThreads(int threadCount)
{
	this->threads = threadCount;
}

void alnData::readTraits(string speciesFile)
{
	string speciesTrait;
	ifstream speciesList (speciesFile);
	if (speciesList.is_open())
	{
		char delimiter = this->inputDelimiter;
		string species;
		float trait;
		while (getline(speciesList, speciesTrait))
		{
			trim(speciesTrait);
			if (speciesTrait.length() < 1)
			{
				continue;
			}
			species = speciesTrait.substr(0, speciesTrait.find(delimiter));
			trait = stof(speciesTrait.substr(speciesTrait.find(delimiter), string::npos));
			if (trait != 0)
			{
				this->species.push_back(species);
				this->traits[species] = trait;
			}
		}
		speciesList.close();
	}
}

void alnData::balanceSample()
{
	float traitSum = 0;
	float tempval;
	string tempSeqid;
	string newSeqid;
	int poolSize;
	vector<string> seqidsPos;
	vector<string> seqidsNeg;
	for (auto it = this->traits.begin(); it != this->traits.end(); it++)
	{
		tempval = it->second;
		traitSum = traitSum + tempval;
		if (tempval > 0)
		{
			seqidsPos.push_back(it->first);
		}
		if (tempval < 0)
		{
			seqidsNeg.push_back(it->first);
		}
	}
	//Upsample
	if (this->upsampleBalance)
	{
		if (traitSum<=-1)
		{
			//Upsample trait-positive sequences.
			poolSize = seqidsPos.size();
			for (int i = 0; i > traitSum; i--)
			{
				//Duplicate a trait-positive seqid at random
				tempSeqid = seqidsPos[randint(0, poolSize-1)];
				newSeqid = tempSeqid + "_pos_dup" + std::to_string(i);
				this->species.push_back(newSeqid);
				this->traits[newSeqid] = 1;
			}
		}
		if (traitSum>=1)
		{
			//Upsample trait-negative sequences.
			poolSize = seqidsNeg.size();
			for (int i = 0; i < traitSum; i++)
			{
				//Select a trait-negative seqid at random
				tempSeqid = seqidsNeg[randint(0, poolSize-1)];
				newSeqid = tempSeqid + "_neg_dup" + std::to_string(i);
				this->species.push_back(newSeqid);
				this->traits[newSeqid] = -1;
			}
		}
	}
	//Downsample
	else if (this->downsampleBalance)
	{
		string targetSpecies;
		if (traitSum <= -1)
		{
			//Downsample trait-negative sequences
			for (int i = 0; i > traitSum; i--)
			{
				//Delete a trait-negative seqid at random
				targetSpecies = seqidsNeg[randint(0, static_cast<int>(seqidsNeg.size() - 1))];
				this->species.erase(find(this->species.begin(), this->species.end(), targetSpecies));
				this->traits.erase(targetSpecies);
				seqidsNeg.erase(find(seqidsNeg.begin(), seqidsNeg.end(), targetSpecies));
			}
		}
		if (traitSum >= 1)
		{
			//Downsample trait-positive sequences
			for (int i = 0; i < traitSum; i++)
			{
				//Delete a trait-positive seqid at random
				targetSpecies = seqidsPos[randint(0, static_cast<int>(seqidsPos.size() - 1))];
				this->species.erase(find(this->species.begin(), this->species.end(), targetSpecies));
				this->traits.erase(targetSpecies);
				seqidsPos.erase(find(seqidsPos.begin(), seqidsPos.end(), targetSpecies));
			}
		}
	}
}

void alnData::printTraits()
{
	for (int i = 0; i < this->species.size(); i++)
	{
		cout << this->species[i] << "	" << this->traits[this->species[i]] << endl;
	}
}


void alnData::processFastaFileList(string alnFileList)
{
	string fastaFileGroup, fastaFileName;
	ifstream fileList (alnFileList);
	if (fileList.is_open())
	{
		while (getline(fileList,fastaFileGroup))
		{
			trim(fastaFileGroup);
			if (fastaFileGroup.length() < 1)
			{
				continue;
			}
			this->groups.push_back(fastaFileGroup);
			stringstream ss(fastaFileGroup);
			while(getline(ss, fastaFileName, ','))
			{
				if (this->geneGroupIndex.find(fastaFileName) == this->geneGroupIndex.end())
				{
					this->geneGroupIndex[fastaFileName] = this->geneGroupIndex.size();
					if (!this->numericInput)
					{
						this->readAln(fastaFileName);
					} else {
						this->readTable(fastaFileName);
					}
					//fastaFileName.erase(fastaFileName.find(".fa"), string::npos);
					//fastaFileName.erase(0, fastaFileName.find("/")+1);
					std::filesystem::path fastaFilePath(fastaFileName);
					//this->currentGene = fastaFileName;
					this->currentGene = fastaFilePath.stem().string();
					if (!this->numericInput)
					{
						this->processAln();
						this->seqs.clear();
					} else {
						this->processTable();
						this->numericSeqs.clear();
					}
				}
			}
		}
		fileList.close();
		if (this->useDiskCache && this->cacheFeatureIndex > 1)
		{
			this->dumpFeatureCache();
		}
	} else {
		throw std::invalid_argument("Could not open alignment list "+ alnFileList +" for reading.");
	}
}

void alnData::readTable(string tableFileName)
{
	string line;
	int seqlen = 0;
	std::size_t found;
	vector<string> tempSpecies;
	tempSpecies = this->species;
	ifstream tableFile (tableFileName);
	if (tableFile.is_open())
	{
		cout << "Processing numeric file: " << tableFileName << "..." << endl;
		string seqid = "";
		vector<float> numericSeq;
		vector<float> valList;
		if (this->numericHeaders)
		{
			getline(tableFile,line);
		}
		while (getline(tableFile,line))
		{
			trim(line);
			stringstream ss(line);
			string segment;
			getline(ss, seqid, this->inputDelimiter);
			if (this->traits.find(seqid) != this->traits.end())
			{
				while(getline(ss, segment, this->inputDelimiter))
				{
					if (segment == "-" || segment == "?")
					{
						valList.push_back(0.0);
					}
					else
					{
						valList.push_back(stof(segment));
					}
				}
				this->numericSeqs[seqid] = valList;
				tempSpecies.erase(find(tempSpecies.begin(), tempSpecies.end(), seqid));
				valList.clear();
			}
		}
		tableFile.close();
	}
	else
	{
		cout << "Could not open rable file " << tableFileName << ", exiting..." << endl;
		exit(1);
	}
	while (!tempSpecies.empty())
	{
		//If seqid is a duplication, set its sequence to the original
		if ((tempSpecies.at(0).find("_pos_dup") != std::string::npos || tempSpecies.at(0).find("_neg_dup") != std::string::npos) && this->numericSeqs.find(tempSpecies.at(0).substr(0, found)) != this->numericSeqs.end())
		{
			found = tempSpecies.at(0).find("_dup") - 4;
			//bool found = (std::find(my_list.begin(), my_list.end(), my_var) != my_list.end());
			this->numericSeqs[tempSpecies.at(0)] = this->numericSeqs[tempSpecies.at(0).substr(0, found)];
			tempSpecies.erase(tempSpecies.begin());
		}
		//Else set its sequence to all indels and add it to the missing seqs file
		else
		{
			this->missingSeqs.push_back(tableFileName + "\t" + tempSpecies.at(0));
			this->numericSeqs[tempSpecies.at(0)] = vector<float>(seqlen, 0);
			tempSpecies.erase(tempSpecies.begin());
		}
	}
}
void alnData::readAln(string fastaFileName)
{
	auto isInvalidChar = [&](char c){return this->validChars.find(c) == std::string::npos;};
	string line;
	int seqlen = 0;
	std::size_t found;
	vector<string> tempSpecies;
	tempSpecies = this->species;
	ifstream fastaFile (fastaFileName);
	if (fastaFile.is_open())
	{
		cout << "Processing FASTA file: " << fastaFileName << "..." << endl;
		string seqid = "";
		string seq;
		while (getline(fastaFile, line))
		{
			trim(line);
			if (line[0] == '>')
			{
				if (this->traits.find(seqid) != this->traits.end())
				{
					this->seqs[seqid] = seq;
					if (seq.length() < seqlen)
					{
						cout << "Sequence with ID " << seqid << " in file " << fastaFileName << " has improper length, exiting..." << endl;
						exit(1);
					}
					seqlen = seq.length();
					tempSpecies.erase(find(tempSpecies.begin(), tempSpecies.end(), seqid));
				}
				seqid = line.substr(1,string::npos);
				seq = "";
			} else if (seqid.length() > 0) {
				if (!this->caseSensitive) {
					std::transform(line.begin(), line.end(), line.begin(), ::toupper);
					std::replace_if(line.begin(), line.end(), isInvalidChar, '-');
				}
				std::replace(line.begin(), line.end(), '?', '-');
				seq = seq + line;
			}
		}
		if (this->traits.find(seqid) != this->traits.end())
		{
			this->seqs[seqid] = seq;
			if (seq.length() < seqlen)
			{
				cout << "Sequence with ID " << seqid << " in file " << fastaFileName << " has improper length, exiting..." << endl;
				exit(1);
			}
			tempSpecies.erase(find(tempSpecies.begin(), tempSpecies.end(), seqid));
		}
		fastaFile.close();
	}
	else
	{
		cout << "Could not open alignment file " << fastaFileName << ", exiting..." << endl;
		exit(1);
	}
	while (!tempSpecies.empty())
	{
		//If seqid is a duplication, set its sequence to the original
		if ((tempSpecies.at(0).find("_pos_dup") != std::string::npos || tempSpecies.at(0).find("_neg_dup") != std::string::npos) && this->seqs.find(tempSpecies.at(0).substr(0, found)) != this->seqs.end())
		{
			found = tempSpecies.at(0).find("_dup") - 4;
			//bool found = (std::find(my_list.begin(), my_list.end(), my_var) != my_list.end());
			this->seqs[tempSpecies.at(0)] = this->seqs[tempSpecies.at(0).substr(0, found)];
			tempSpecies.erase(tempSpecies.begin());
		}
		//Else set its sequence to all indels and add it to the missing seqs file
		else
		{
			this->missingSeqs.push_back(fastaFileName + "\t" + tempSpecies.at(0));
			this->seqs[tempSpecies.at(0)] = string(seqlen, '-');
			tempSpecies.erase(tempSpecies.begin());
		}
	}
}

void alnData::processTable()
{
	int seqLen = this->numericSeqs[this->species[1]].size();
	int numSpecies = this->species.size();
	int groupStartIndex = this->featureIndex;
	string featureName;
	set<float> values;
	vector<float> feature;

	for (int i = 0; i < seqLen; i++)
	{
		//check if site is constant
		for (int j = 0; j < numSpecies; j++)
		{
			values.insert(this->numericSeqs[this->species[j]][i]);
		}
		if (values.size() < 2)
		{
			values.clear();
			continue;
		}
		featureName = this->currentGene + "_" + to_string(i);
		for (int k = 0; k < numSpecies; k++)
		{
			feature.push_back(this->numericSeqs[this->species[k]][i]);
		}
		this->featureMap[this->featureIndex] = featureName;
		this->features[featureIndex] = feature;
		this->featureIndex++;
		values.clear();
		feature.clear();
	}
	this->groupIndices.push_back({groupStartIndex,this->featureIndex-1});
}

void alnData::processAln()
{
	int seqLen = this->seqs[this->species[1]].size();
	int numSpecies = this->species.size();
	int groupStartIndex = this->featureIndex;
//	int cacheFeatureIndex = 1;

	for (int i = 0; i < seqLen; i++)
	{
		//check if site is constant
		set<char> bases;
		map<char, int> baseCounts;
		bases.insert('-');
		for (int j = 0; j < numSpecies; j++)
		{
			bases.insert(this->seqs[this->species[j]][i]);
		}
		//process it into features if it isn't
		if (bases.size() > 2)
		{
			if (this->ignoreSingletons)
			{
				for (auto elem : bases)
				{
					baseCounts[elem] = 0;
				}
				for (int j = 0; j < numSpecies; j++)
				{
					baseCounts[this->seqs[this->species[j]][i]]++;
				}
				//Sum non-indel baseCounts, then subtract major count, and check if remainder is greater than count threshold
				int baseCountTotal = 0, maxCount = 0;
				for (auto elem : bases)
				{
					if (elem == '-')
						continue;
					if (baseCounts[elem] > maxCount)
						maxCount = baseCounts[elem];
					baseCountTotal += baseCounts[elem];
				}
				if (baseCountTotal - maxCount < this->countThreshold)
					continue;
			}
			set<char>::iterator baseIter = bases.begin();
			for (int j = 0; j < bases.size(); j++)
			{
				if (*baseIter != '-')
				{
					vector<float> oneHot;
					int featureSum = 0;
					string featureName = this->currentGene + "_" + to_string(i) + "_" + *baseIter;
					for (int k = 0; k < numSpecies; k++)
					{
						if (this->seqs[this->species[k]][i] == *baseIter)
						{
							oneHot.push_back(1);
							featureSum += 1;
						} else if (this->indelFuzzing && this->seqs[this->species[k]][i] == '-') {
							oneHot.push_back(0.5);
						} else {
							oneHot.push_back(0);
						}
					}
					if (featureSum < this->countThreshold)
					{
						baseIter++;
						continue;
					}
					this->featureMap[this->featureIndex] = featureName;
					//If disk cache is being used, overwrite features starting from 1 for each gene
					if (this->useDiskCache)
					{
						this->features[this->cacheFeatureIndex] = oneHot;
					}
					else
					{
						this->features[featureIndex] = oneHot;
					}
					this->featureIndex++;
					this->cacheFeatureIndex++;
				}
				baseIter++;
			}
		}
	}
	this->groupIndices.push_back({groupStartIndex,this->featureIndex-1});
	if (this->useDiskCache && this->cacheFeatureIndex * this->featureSize * this->species.size() >= this->workingMemLimit)
	{
		this->dumpFeatureCache();
	}
}


//void alnData::dumpSpecies(int i, Semaphore& semaphore)
void alnData::dumpSpecies(int i)
{
//	semaphore.wait();
	if (this->featureIndex != this->cacheFeatureIndex)
	{
		cout << "Dumping features to cache for " << this->species[i] << " species." << endl;
	}
	string cacheFileName;
	string cacheLineString;
	ofstream cacheFile;
	vector<string> cacheSegment;
	for (int j = 1; j < this->cacheFeatureIndex; j++)
	{
//		cerr << j << ".1 " << i << ".1 " << this->features[j][i] << endl;
//		cout << this->features[j].size() << endl;
//		cout << this->features[j][0] << endl;
//		cout << this->features[j][1] << endl;
//		cout << this->features[j][i] << endl;
		if (this->features[j][i] == 0.0)
		{
			cacheSegment.push_back(to_string(0));
		}
		else if (this->features[j][i] == 1.0)
		{
			cacheSegment.push_back(to_string(1));
		}
		else
		{
			cacheSegment.push_back(to_string(this->features[j][i]));
		}
//		cout << j << ".2 " << i << ".2" << endl;
	}
	stringstream cacheLine;
	copy(cacheSegment.begin(), cacheSegment.end(), ostream_iterator<string>(cacheLine, this->outputDelimiter.c_str()));
	//copy(cacheSegment.begin(), cacheSegment.end(), ostream_iterator<int>(cacheLine, "	"));
	cacheLineString = cacheLine.str();
	cacheFileName = this->cacheDir + "/.cache_" + this->species[i] + ".txt";
	if (this->featureIndex == this->cacheFeatureIndex)
	{
		std::remove(cacheFileName.c_str());
	}
	else
	{
		cacheLineString = this->outputDelimiter + cacheLineString;
	}
	cacheFile.open(cacheFileName, std::ofstream::app);
	if (!cacheFile.is_open())
	{
				cout << "Could not open disk caching file for writing, quitting..." << endl;
		exit;
	}
	cacheLineString.pop_back();
	cacheFile << cacheLineString.c_str();
	cacheFile.close();
	if (this->featureIndex != this->cacheFeatureIndex)
	{
		cout << "Finished dumping features to cache for " << this->species[i] << " species." << endl;
	}
//	semaphore.signal();
}

void alnData::dumpSpeciesMulti(int i, Semaphore& semaphore)
//void alnData::dumpSpecies(int i)
{
//	semaphore.wait();
	if (this->featureIndex != this->cacheFeatureIndex)
	{
		cout << "Dumping features to cache for " << this->species[i] << " species." << endl;
	}
	string cacheFileName;
	string cacheLineString;
	ofstream cacheFile;
	vector<string> cacheSegment;
	for (int j = 1; j < this->cacheFeatureIndex; j++)
	{
//		cerr << j << ".1 " << i << ".1 " << this->features[j][i] << endl;
//		cout << this->features[j].size() << endl;
//		cout << this->features[j][0] << endl;
//		cout << this->features[j][1] << endl;
//		cout << this->features[j][i] << endl;
		if (this->features[j][i] == 0.0)
		{
			cacheSegment.push_back(to_string(0));
		}
		else if (this->features[j][i] == 1.0)
		{
			cacheSegment.push_back(to_string(1));
		}
		else
		{
			cacheSegment.push_back(to_string(this->features[j][i]));
		}
//		cout << j << ".2 " << i << ".2" << endl;
	}
	stringstream cacheLine;
	copy(cacheSegment.begin(), cacheSegment.end(), ostream_iterator<string>(cacheLine, this->outputDelimiter.c_str()));
	//copy(cacheSegment.begin(), cacheSegment.end(), ostream_iterator<int>(cacheLine, "	"));
	cacheLineString = cacheLine.str();
	cacheFileName = this->cacheDir + "/.cache_" + this->species[i] + ".txt";
	if (this->featureIndex == this->cacheFeatureIndex)
	{
		std::remove(cacheFileName.c_str());
	}
	else
	{
		cacheLineString = this->outputDelimiter + cacheLineString;
	}
	cacheFile.open(cacheFileName, std::ofstream::app);
	if (!cacheFile.is_open())
	{
				cout << "Could not open disk caching file for writing, quitting..." << endl;
		exit;
	}
	cacheLineString.pop_back();
	cacheFile << cacheLineString.c_str();
	cacheFile.close();
	if (this->featureIndex != this->cacheFeatureIndex)
	{
		cout << "Finished dumping features to cache for " << this->species[i] << " species." << endl;
	}
	semaphore.signal();
}


void alnData::dumpFeatureCache()
{

	cout << "Dumping " << this->cacheFeatureIndex << " features to disk cache for " << this->species.size() << " species." << endl;

	if (this->threads < 2)
	{
		for (int i = 0; i < this->species.size(); i++)
		{
	//		cerr << "Species IDX starting: " << i << endl;
			this->dumpSpecies(i);
	//		cerr << "Species IDX finished: " << i << endl;
		}
	} else {
//  Static compilation with glibc causes use of condition_variable to generate a seg-fault
		std::vector<std::thread> threads;
		Semaphore semaphore(this->threads);

		for (int i = 0; i < this->species.size(); i++)
		{
			semaphore.wait();
			threads.emplace_back(&alnData::dumpSpeciesMulti, this, i, std::ref(semaphore));
			//semaphore.signal()
		}


		for (auto& thread : threads)
		{
			thread.join();
		}
	}



	this->cacheFeatureIndex = 1;
//	cerr << "2" << endl;
}


void alnData::generateResponseFile(string baseName)
{
	string responseFileName = baseName + "/response_" + baseName + ".txt";
	ofstream responseFile (responseFileName);
	if (responseFile.is_open())
	{
		for (int i = 0; i < this->species.size(); i++)
		{
			responseFile << this->traits[this->species[i]] << endl;
		}
		responseFile.close();
	}
	if (this->downsampleBalance || this->upsampleBalance)
	{
		string resampleFileName = baseName + "/resampled_" + baseName + ".txt";
		ofstream resampleFile (resampleFileName);
		if (resampleFile.is_open())
		{
			for (int i = 0; i < this->species.size(); i++)
			{
				resampleFile << this->species[i] << "\t" << this->traits[this->species[i]] << endl;
			}
			resampleFile.close();
		}
	}
}

void alnData::generateFeatureFile(string baseName)
{
	if (this->useDiskCache)
	{
		string cacheFileName;
		string cacheLine;
		ifstream cacheFile;
		string featuresFileNameNew = baseName + "/feature_" + baseName + ".txt";
		ofstream featuresFileNew (featuresFileNameNew);
		if (!featuresFileNew.is_open())
		{
			cout << "Could not open features output file, quitting..." << endl;
			exit;
		}
		for (int i = 0; i < this->species.size(); i++)
		{
			cacheFileName = this->cacheDir + "/.cache_" + this->species[i] + ".txt";
			cacheFile.open(cacheFileName);
			if (!cacheFile.is_open())
			{
                		cout << "Could not open disk caching file for reading, quitting..." << endl;
				exit;
			}
			getline(cacheFile,cacheLine);
			featuresFileNew << cacheLine << endl;
			cacheFile.close();
			std::remove(cacheFileName.c_str());
		}
		return;
	}
//	string featuresFileName = baseName + "/feature_" + baseName + ".txt";
        //Open as many file handles as possible, then close enough to leave a buffer for writing other files
        vector<ofstream*> outputHandles;
        int bufferHandles = 30;
        int handleCount = 0;
	string tempFileBase = "test_out";
        while (true)
        {
                string ofname = baseName + "/" + tempFileBase + to_string(handleCount) + ".txt";
                outputHandles.push_back(new ofstream(ofname));
//		cout << ofname << endl;
                if (!outputHandles[handleCount]->is_open())
                {
                        break;
                }
		handleCount++;
        }
	int fPerHandle;
	fPerHandle = ceil((float)this->features.size() / (handleCount - bufferHandles));
	bufferHandles = handleCount - ceil((float)this->features.size() / (float)fPerHandle);
//	cout << "0" << endl;
        for (int i = handleCount-1; i >= handleCount - bufferHandles; i--)
        {
                outputHandles[i]->close();
		string handleFName = baseName + "/" + tempFileBase + to_string(i) + ".txt";
		std::remove(handleFName.c_str());
        }
        handleCount = handleCount - bufferHandles;
//	cout << "Features: " << this->features.size() << endl;
//	cout << "Handles: " << handleCount << endl;
//	cout << "Features per handle: " << fPerHandle << endl;
	//transpose features first for efficiency
	//vector<vector<float>> tFeatures;
//	vector<float*> tFeatures;
	vector<float*> featureCache;
//	for (int i = 0; i < this->species.size(); i++)
//	{
//		tFeatures.push_back((float*) malloc(this->features.size() * sizeof(float)));
//	}
        for (int i = 0; i < fPerHandle; i++)
        {
                featureCache.push_back((float*) malloc(this->species.size() * sizeof(float)));
        }
	int handleIdx = 0;
	for (int i = 0; i < this->features.size(); i++)
	{
		vector<float> oneHot = this->features[i+1];
//		if (i==0)
//		{
//			cout << to_string(oneHot[2]) << endl;
//		}
		int cacheIdx = i % fPerHandle;
		handleIdx = i / fPerHandle;
		float sum = 0;
		float val;
//		for (int j = 0; j < this->features[i].size(); j++)
		for (int j = 0; j < oneHot.size(); j++)
		{
//			sum += this->features[i][j];
			sum += oneHot[j];
		}
		val = 1.0;
		if (!this->numericInput)
		{
			if (this->normalize)
			{
				val = 1.0/sum;
			}
	//		if (this->ignoreSingletons && sum == 1.0)
			if (sum < this->countThreshold)
			{
				val = 0.0;
			}
		}
		for (int j = 0; j < oneHot.size(); j++)
		{
			//tFeatures[j].push_back(val * oneHot[j]);
//			tFeatures[j][i] = val * oneHot[j];
			featureCache[cacheIdx][j] = val * oneHot[j];
		}
		if (cacheIdx == fPerHandle - 1 || i == this->features.size() - 1)
		{
			int fPerHandleTemp = fPerHandle;
			if (i == this->features.size() - 1 && this->features.size() % fPerHandle > 0)
			{
				fPerHandleTemp = this->features.size() % fPerHandle;
			}
			for (int j = 0; j < this->species.size(); j++)
			{
				for (int k = 0; k < fPerHandleTemp; k++)
				{
					if (featureCache[k][j] == 0.0)
					{
						*outputHandles[handleIdx] << to_string(0) << this->outputDelimiter;
					}
					else
					{
						*outputHandles[handleIdx] << to_string(featureCache[k][j]) << this->outputDelimiter;
					}
				}
				*outputHandles[handleIdx] << endl;

			}
			outputHandles[handleIdx]->close();
		}
	}
        string featuresFileNameNew = baseName + "/feature_" + baseName + ".txt";
        ofstream featuresFileNew (featuresFileNameNew);
	if (!featuresFileNew.is_open())
        {
                cout << "Could not open features output file, quitting..." << endl;
		exit;
        }
	vector<ifstream*> inputHandles;
	for (int i=0; i < handleCount; i++)
	{
		string ifname = baseName + "/" + tempFileBase + to_string(i) + ".txt";
                inputHandles.push_back(new ifstream(ifname));
	}
	string fragment;
	string featureLineNew;
	while (std::getline(*inputHandles[0], fragment))
	{
//		featuresFileNew << fragment << '\t';
//		cout << "First fragment length: " << fragment.length() << endl;
		featureLineNew = fragment;
		for (int i=1; i < handleCount; i++)
		{
			std::getline(*inputHandles[i], fragment);
			featureLineNew = featureLineNew + fragment;
//			featuresFileNew << fragment << '\t';
		}
//		featuresFileNew << endl;
		featureLineNew.pop_back();
		featuresFileNew << featureLineNew << endl;
	}
	for (int i=0; i < handleCount; i++)
	{
		inputHandles[i]->close();
		string handleFName = baseName + "/" + tempFileBase + to_string(i) + ".txt";
		std::remove(handleFName.c_str());
	}
	featuresFileNew.close();
}

void alnData::generateStatsFile(string baseName)
{
	string statsFileName = baseName + "/feature_stats_" + baseName + ".txt";
	ofstream statsFile (statsFileName);
	if (statsFile.is_open())
	{
		statsFile << "Samples\t" << this->species.size() << endl;
		statsFile << "Features\t" << this->featureIndex << endl;
	}
	statsFile.close();
}

void alnData::generateMappingFile(string baseName)
{
	string mappingFileName = baseName + "/feature_mapping_" + baseName + ".txt";
	ofstream mappingFile (mappingFileName);
	if (mappingFile.is_open())
	{
		for (int i = 0; i < this->featureMap.size(); i++)
		{
			mappingFile << i << "\t" << this->featureMap[i] << endl;
		}
		mappingFile.close();
	}
}

void alnData::generateGroupIndicesFile(string baseName)
{
	string groupIndicesFileName = baseName + "/group_indices_" + baseName + ".txt";
	string fieldFileName = baseName + "/field_" + baseName + ".txt";
	ofstream groupIndicesFile (groupIndicesFileName);
	ofstream fieldFile (fieldFileName);
	string indStarts, indEnds, weights, weightBuffer, gene;
	double weight;
	int geneStart, geneEnd;
	int groupStart = 1;
	int groupEnd = 0;
	if (groupIndicesFile.is_open() and fieldFile.is_open())
	{
		for (int i = 0; i < this->groups.size(); i++)
		{
			stringstream ss(this->groups[i]);
			while(getline(ss, gene, ','))
			{
				geneStart = groupIndices[geneGroupIndex[gene]][0];
				geneEnd = groupIndices[geneGroupIndex[gene]][1];
				for (int j = geneStart; j <= geneEnd; j++)
				{
					fieldFile << to_string(j) + this->outputDelimiter;
				}
				groupEnd += geneEnd-geneStart+1;
			}
			indStarts.append(to_string(groupStart) + this->outputDelimiter);
			indEnds.append(to_string(groupEnd) + this->outputDelimiter);
			weight = sqrt(1+(groupEnd - groupStart));
			weights.append(to_string(weight) + this->outputDelimiter);
			groupStart = groupEnd + 1;
		}
		indStarts.pop_back();
		indEnds.pop_back();
		weights.pop_back();

		groupIndicesFile << indStarts << endl;
		groupIndicesFile << indEnds << endl;
		groupIndicesFile << weights << endl;

		groupIndicesFile.close();
		fieldFile.close();
	}

}

void alnData::generateMissingFile(string baseName)
{
	string missingSeqsFileName = baseName + "/missing_seqs_" + baseName + ".txt";
	ofstream missingSeqsFile (missingSeqsFileName);
	if (missingSeqsFile.is_open())
	{
		for (int i = 0; i < this->missingSeqs.size(); i++)
		{
			missingSeqsFile << this->missingSeqs[i] << endl;
		}
	}
	missingSeqsFile.close();
}



