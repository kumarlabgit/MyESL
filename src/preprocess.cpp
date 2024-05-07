#include "preprocess.hpp"
#include <algorithm>
#include <filesystem>
#include <sstream>
#include <random>

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
	this->featureIndex = 1;
	this->normalize = false;
	this->useDiskCache = false;
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
	this->delimiter = delimiter;
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

void alnData::readTraits(string speciesFile)
{
	string speciesTrait;
	ifstream speciesList (speciesFile);
	if (speciesList.is_open())
	{
		string delimiter = "	";
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
					this->readAln(fastaFileName);
					//fastaFileName.erase(fastaFileName.find(".fa"), string::npos);
					//fastaFileName.erase(0, fastaFileName.find("/")+1);
					std::filesystem::path fastaFilePath(fastaFileName);
					//this->currentGene = fastaFileName;
					this->currentGene = fastaFilePath.stem().string();
					this->processAln();
					this->seqs.clear();
				}
			}
		}
		fileList.close();
	} else {
		throw std::invalid_argument("Could not open alignment list "+ alnFileList +" for reading.");
	}
}

void alnData::readAln(string fastaFileName)
{
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

void alnData::processAln()
{
	int seqLen = this->seqs[this->species[1]].size();
	int numSpecies = this->species.size();
	int groupStartIndex = this->featureIndex;
	int cacheFeatureIndex = 1;

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
						this->features[cacheFeatureIndex] = oneHot;
					}
					else
					{
						this->features[featureIndex] = oneHot;
					}
					this->featureIndex++;
					cacheFeatureIndex++;
				}
				baseIter++;
			}
		}
	}
	this->groupIndices.push_back({groupStartIndex,this->featureIndex-1});
	if (this->useDiskCache)
	{
		string cacheFileName;
		string cacheLineString;
		ofstream cacheFile;
		for (int i = 0; i < this->species.size(); i++)
		{
			vector<int> cacheSegment;
			for (int j = 1; j < cacheFeatureIndex; j++)
			{
				cacheSegment.push_back(this->features[j][i]);
			}
			stringstream cacheLine;
			copy(cacheSegment.begin(), cacheSegment.end(), ostream_iterator<int>(cacheLine, "	"));
			cacheLineString = cacheLine.str();
			cacheFileName = ".cache_" + this->species[i] + ".txt";
			if (this->geneGroupIndex.size() == 1)
			{
				std::remove(cacheFileName.c_str());
			}
			else
			{
				cacheLineString = "	" + cacheLineString;
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
		}
	}
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
			cacheFileName = ".cache_" + this->species[i] + ".txt";
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
		if (this->normalize)
		{
			val = 1.0/sum;
		}
//		if (this->ignoreSingletons && sum == 1.0)
		if (sum < this->countThreshold)
		{
			val = 0.0;
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
						*outputHandles[handleIdx] << to_string(0) << this->delimiter;
					}
					else
					{
						*outputHandles[handleIdx] << to_string(featureCache[k][j]) << this->delimiter;
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


	//Open as many file handles as possible, then close enough to leave a buffer for writing other files
	/*
	vector<ofstream*> outputHandles;
	int bufferHandles = 20;
	int handleCount = 0;
	while (true)
	{
		string ofname = baseName + "/test_out" + to_string(handleCount) + ".txt";
		outputHandles.push_back(new ofstream(ofname));
		handleCount++;
		cout << ofname << endl;
		if (!outputHandles[handleCount-1]->is_open())
		{
			break;
		}
	}
	cout << "0" << endl;
	for (int i = handleCount-1; i > handleCount - bufferHandles; i--)
	{
		outputHandles[i]->close();
	}
	handleCount = handleCount - bufferHandles;
	int fPerHandle;
	cout << "1" << endl;
	//featuresFile.close();
	//write to file
	*/
/*
	ofstream featuresFile (featuresFileName);
        if (!featuresFile.is_open())
        {
                cout << "Could not open features output file, quitting..." << endl;
        }
	cout << "2" << endl;
//	fPerHandle = ceil(this->features.size() / handleCount);
	if (featuresFile.is_open())
	{
		for (int i = 0; i < tFeatures.size(); i++)
		{
			string featureLine;
			int handleIdx = 0;
			for (int j = 0; j < this->features.size(); j++)
			{
//				if (handleIdx != j / fPerHandle)
//				{
//					//and line break to last handle
//					*outputHandles[handleIdx] << endl;
//					handleIdx = j / fPerHandle;
//				}
				if (tFeatures[i][j] == 0.0)
				{
					featureLine.append(to_string(0) + '\t');
//					*outputHandles[handleIdx] << to_string(0) + '\t';
				}
				else
				{
					featureLine.append(to_string(tFeatures[i][j]) + '\t');
//					*outputHandles[handleIdx] << to_string(tFeatures[i][j]) + '\t';
				}
			}
//			*outputHandles[handleIdx] << endl;
			featureLine.pop_back();
			featuresFile << featureLine << endl;
		}
		featuresFile.close();
	}
	cout << "3" << endl;
*/

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
					fieldFile << to_string(j) + this->delimiter;
				}
				groupEnd += geneEnd-geneStart+1;
			}
			indStarts.append(to_string(groupStart) + this->delimiter);
			indEnds.append(to_string(groupEnd) + this->delimiter);
			weight = sqrt(1+(groupEnd - groupStart));
			weights.append(to_string(weight) + this->delimiter);
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



