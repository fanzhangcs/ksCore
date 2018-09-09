/*
Copyright (c) of 2018 by Fan Zhang <fanzhang.cs@gmail.com>
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<vector>
#include<time.h>
#include<iostream>
#include<algorithm>
#include<cstdlib>
#include<ctime>
#include <fstream>
#include "khash.h"

#define testProg 1//0:automation, 1:manual test
#define algorithmChoice 2//0:given k and s, 1:with (k,k-1)-core decomp index, 2: (k,k-1)-core decomp (k-fami), 3: compute k-core, 4: core decomp, 5: compute k-truss
#define computeScore 1//compute active user rate, modularity and cc
#define assesss 0//1:open
#define outputYelp 0
#define outputDBLP 0

using namespace std;

long anchorID;
long maxFollowerNum = -1;
double totalFollowerNum = 0;
long bestAnchorEdgeSize = 0;
char fedges[200], fidnames[150], fstati[200], fmaxCommunityNames[150], fmaxCommunityEdges[150];
long inputK, inputS, datasetID, inikcoreSize, inidegreekNeiSize, visitNumber = 0;
long iniKtrussSize, iniKtrussEdgeSize, iniKm1TrussSize, iniKm1TrussEdgeSize, iniKm2CoreEdgeSize, kcoreVerSize, kcoreEdgeSize, anchoredKtrussSize, anchoredKtrussEdgeSize, ksCoreSize, anchoredKm1TrussEdgeSize;
long triedAncKtrussEdgeSize;
long maxLayerNum, ancInKtNum = 0;
double followerNumber = 0, compTrussTime = 0, constructTimeTag = 0, preprocessingTime = 0, layerBylayerTime = 0, khashConstructTime = 0, edgesBetwNeighTime = 0, edgeBetwNeighRecoTime = 0, layerSearchNeiTime = 0, layerSearchCandTime = 0;
double kCoreTimeTag, kCoreTime = 0, edgeSupportTimeTag, edgeSupportTime = 0, ksCoreTimeTag, km1TrussOnCoreTime = 0, ksCoreTime = 0, kTrussTimeTag, kTrussTime = 0, findBestTimeTag, findBestTime = 0;
double algStartTimeTag, oneKm1TrussTime = 0, oneKtrussOnKm1Time = 0, indexTime = 0, triaCountTime = 0, constructTime = 0;
long triangleNum = 0, triedVerNum = 0, lesskIniNum = 0;
long iniK, insK;

vector<long> degreekVertices, degreekNeighbors, degreekNeiTag, kcoreDeletedTag, degreekVerticesNew, followerTag;//set P and T
vector<long> verSetTag, verOrigIDs, kcoreVerSetIDs, verTag, kcoreTag, verDegree, lesskSet, followersRecord, verSetKcoreIDs;
vector<long> collapserIDs, numFollowersRecord, anchorIDs;
vector<long> neighborInTriangleIndex;
vector<vector<char> > verNames;
vector<vector<long> > verSet0, verSet, kcoreSet, kcoreSetNeiIndex, kcoreSetSupport, km2coreSetOneCount, km2coreEdgeIDSet, km2coreEdgeNeiSet, verNameIDs, trussEdgeToVerID;
vector<vector<long> > km1ShellVerLayers, km1ShellEdges, km1TrussSupport, km1TrussToDelete, candidateSetLayerNum, candidateSetSearchTag;
vector<long> edgeLayerNumForToDelete;
vector<vector<vector<long> > > km1ShellEdgeLayers, candidateTriaSet, candidateTriaSet1, candidatePatTriaSet, candidatePatTriaSet1, candidateHorTriaSet, candidateHorTriaSet1, candidateSupTriaSet, candidateSupTriaSet1;
vector<long> candNeiCompIndex, kcoreDegree, km1TrussVerices, km1ShellVertices, kTrussVerices, candidateDegree, candidateTrussDegree;
vector<vector<long> > anchoredEdgeSet, candidateEdgeAnchorTag, candidateEdgeTag, candidateEdgeTrussTag, candidateSet, candidateTrussSet, candidateTrussEdgeTag, candidateSetNeiIndex, candidateSetSupport, candidateTrussSetSupport, km1TrussNeiEdgeSet, candidateToDelete, candidateTrussToDelete;
vector<long> km1ShellTriaNeighrbors, km1ShellTriaNeigIndex, km1TrussDegree, candidateVer, candidateTag, toDeleteEdgeSetVid, candidateTrussTag, km2coreOrigTag;
vector<long> degSmallVer, degreeVerDegree;
vector<long> degreeVertices, bestFollowers;
vector<vector<long> > kcoreEdgeTag;
vector<long> coreDecomTag, supBlockTag;
vector<vector<long> > supIncOrder, supOrderIndex;
vector<long> kkcoreDecompTag;
vector<long> degreeIncOrder, blockTag, orderIndex;
vector<long> kkdegreeIncOrder, kkblockTag, kkorderIndex;

vector<long> coreTag, coreDegree, coreVertices, coreIndex;
vector<vector<long> > inikcoreSupport, coreSet, coreEdgeSupport, coreEdgeTag, coreNeiIndex;
vector<vector<long> > toDeleteEdgeSet;

KHASH_MAP_INIT_INT(32, long)
vector<khash_t(32)*> kcoreTable, coreTable;

string infile, outfile;

void dataInput()
{
	double time0 = (double)clock() / CLOCKS_PER_SEC;

	//read edges, build verSet //need first row ordered data
	long vertexID, neighborID;
	vector<long> verSetInsertion;
	ifstream fin;
	fin.open(infile.c_str());

	if (!fin) 
	{
		cout << "Open Fail";
		exit(1); 
	}
	else
	{
		long verid = -1, vid, nid;
		char fech = 's';
		while (!fin.eof())
		{
			fin >> vertexID >> neighborID;
			vid = verSetTag[vertexID];
			nid = verSetTag[neighborID];
			if (vid < 0)
			{
				verid++;
				verSetTag[vertexID] = verid; //vertexID -> vid
				verOrigIDs.push_back(vertexID); //vid -> vertexID
				vid = verid;
			}
			if (nid < 0)
			{
				verid++;
				verSetTag[neighborID] = verid;
				verOrigIDs.push_back(neighborID);
				nid = verid;
			}
			if (verSet0.size() == 0)
			{
				verSetInsertion.clear();
				verSetInsertion.push_back(vid);
				verSetInsertion.push_back(nid);
				verSet0.push_back(verSetInsertion);
			}
			else
			{
				if (vid == verSet0[verSet0.size() - 1][0])
				{
					verSet0[verSet0.size() - 1].push_back(nid);
				}
				else
				{
					verSetInsertion.clear();
					verSetInsertion.push_back(vid);
					verSetInsertion.push_back(nid);
					verSet0.push_back(verSetInsertion);
				}
			}
		}
	}
	fin.close();

	//1. remove duplicates£¬ create verSet; 2. make unpaired edge to be paired ; 3. remove self loop
	vector<vector<long> > verNei;
	verNei.resize(verOrigIDs.size());
	for (long i = 0; i < verSet0.size(); i++)
	{
		long id = verSet0[i][0];
		for (long j = 1; j < verSet0[i].size(); j++)
		{
			long nid = verSet0[i][j];
			verNei[id].push_back(nid); 
			verNei[nid].push_back(id); //duplicate
		}
	}
	long eNum = 0;
	vector<long> vTag, nTag;
	vTag.resize(verOrigIDs.size());
	nTag.resize(verOrigIDs.size(), -1);
	verSet.resize(verOrigIDs.size());
	long nT = 0;
	for (long i = 0; i < verNei.size(); i++)
	{
		for (long j = 0; j < verNei[i].size(); j++)
		{
			long nid = verNei[i][j];
			if (nTag[nid] != nT && nid != i) //remove duplicate neighbor and self loop
			{
				verSet[i].push_back(nid);
				nTag[nid] = nT;
				eNum++;
			}
		}
		nT++;
	}
	vector<long>().swap(vTag);
	vector<long>().swap(nTag);
	vector<vector<long> >().swap(verNei);

	cout << "Vertices: " << verSet.size() << " Edges: " << eNum/2 << " Avg. Degree: " << double(eNum)/2/verSet.size();
	double time1 = (double)clock() / CLOCKS_PER_SEC;
	cout << "\nRead time: " << time1 - time0 << "\n";
}

void computeKcore0()
{
	vector<long> verTag1 = verTag, lesskSet1 = lesskSet, verDegree1 = verDegree;

	//kcorePruNum = 0;
	for (long i = 0; i < lesskSet1.size(); i++)
	{
		long id = lesskSet1[i];
		if (verTag1[id] > 0)
		{
			//kcorePruNum++;
			verTag1[id] = 0; //cur vertex deleted
			for (long j = 0; j < verSet[id].size(); j++)
			{
				long nid = verSet[id][j];
				if (verTag1[nid] == 1)
				{
					verDegree1[nid]--; //neighbor degree - 1
					if (verDegree1[nid] < iniK) //new candidate for computing
					{
						lesskSet1.push_back(nid);
						verTag1[nid] = 2;
					}
				}
			}
		}
	}

	vector<long> kcoreInsertion;
	long ki = -1;
	for (long i = 0; i < verTag1.size(); i++)
	{
		if (verTag1[i])
		{
			ki++;
			verSetKcoreIDs[i] = ki;
		}
		else
		{
			verSetKcoreIDs[i] = -1;
		}
	}

	kcoreTag.clear();
	kcoreSet.clear();
	kcoreVerSetIDs.clear();
	kcoreEdgeSize = 0;

	for (long i = 0; i < verTag1.size(); i++)
	{
		if (verTag1[i])
		{
			kcoreInsertion.clear();
			for (long j = 0; j < verSet[i].size(); j++)
			{
				long nid = verSet[i][j];
				if (verTag1[nid])
				{
					kcoreInsertion.push_back(verSetKcoreIDs[nid]);
				}
			}
			kcoreTag.push_back(1); //initialize 
			kcoreSet.push_back(kcoreInsertion); //vid-nids set
			kcoreVerSetIDs.push_back(i); //km2core id -> verSet id
			kcoreEdgeSize += kcoreInsertion.size();
		}
	}
	kcoreEdgeSize /= 2;
}

void computeEdgeSupport0()
{
	int ret, is_missing;
	khash_t(32) *h;
	khiter_t kr;
	vector<long> setInsertion, setInsertion1;
	kcoreSetSupport.clear();
	kcoreSetNeiIndex.clear();
	kcoreTable.clear();
	kcoreEdgeTag.clear();
	double khashStart, khashEnd;

	for (long i = 0; i < kcoreSet.size(); i++)
	{
		setInsertion.clear();
		setInsertion1.clear();
		h = kh_init(32);
		for (long j = 0; j < kcoreSet[i].size(); j++)
		{
			long nid = kcoreSet[i][j];
			setInsertion.push_back(-1);
			setInsertion1.push_back(1);
			khashStart = (double)clock() / CLOCKS_PER_SEC;
			kr = kh_put(32, h, nid, &ret);
			kh_val(h, kr) = j;
			khashEnd = (double)clock() / CLOCKS_PER_SEC;
			khashConstructTime += khashEnd - khashStart;
		}
		kcoreEdgeTag.push_back(setInsertion1);
		kcoreSetSupport.push_back(setInsertion);
		kcoreSetNeiIndex.push_back(setInsertion);
		kcoreTable.push_back(h);
	}

	//compute edge support
	lesskSet.clear();
	kcoreDegree.clear();
	long sDegree;
	long matchIdInIndexVer, indexIdInMatchVer;
	for (long i = 0; i < kcoreSet.size(); i++)
	{
		sDegree = 0;
		for (long j = 0; j < kcoreSet[i].size(); j++)
		{
			if (kcoreSetSupport[i][j] < 0)
			{
				long nid = kcoreSet[i][j];
				long indexId, matchId;
				long support = 0;
				if (kcoreSet[i].size() < kcoreSet[nid].size())
				{
					indexId = i;
					matchId = nid;
				}
				else
				{
					indexId = nid;
					matchId = i;
				}
				//find common neighbors
				h = kcoreTable[matchId];
				kr = kh_get(32, h, indexId);
				indexIdInMatchVer = kh_val(h, kr);
				for (long k = 0; k < kcoreSet[indexId].size(); k++)
				{
					long knid = kcoreSet[indexId][k];
					kr = kh_get(32, h, knid);
					if (kr != kh_end(h))
					{
						support++;
					}
				}
				h = kcoreTable[indexId];
				kr = kh_get(32, h, matchId);
				matchIdInIndexVer = kh_val(h, kr);

				kcoreSetSupport[indexId][matchIdInIndexVer] = support;
				kcoreSetSupport[matchId][indexIdInMatchVer] = support;
				kcoreSetNeiIndex[indexId][matchIdInIndexVer] = indexIdInMatchVer;
				kcoreSetNeiIndex[matchId][indexIdInMatchVer] = matchIdInIndexVer;
				if (support > inputS - 1)
				{
					sDegree++;
				}
				triangleNum += support;
			}
			else if (kcoreSetSupport[i][j] > inputS - 1)
			{
				sDegree++;
			}
		}
		kcoreDegree.push_back(sDegree);
		if (sDegree < iniK)
		{
			lesskSet.push_back(i);
			kcoreTag[i] = 2;
		}
	}
	triangleNum /= 3;
}

void computeKScore0()
{
	khash_t(32) *h;
	khiter_t kr;
	long indexId, matchId;
	vector<long> toDeleteInsertion;
	long vidInMatchVer;

	for (long i = 0; i < lesskSet.size(); i++)
	{
		long id = lesskSet[i];
		kcoreTag[id] = 0;
		kcoreDegree[id] = 0;
		for (long j = 0; j < kcoreSet[id].size(); j++)
		{
			long nid = kcoreSet[id][j];
			kcoreEdgeTag[id][j] = 0;
			long nnid = kcoreSetNeiIndex[id][j];
			kcoreEdgeTag[nid][nnid] = 0;

			if (kcoreTag[nid] == 1)
			{
				if (kcoreSetSupport[id][j] > inputS - 1)
				{
					kcoreDegree[nid]--;
					if (kcoreDegree[nid] < iniK && kcoreTag[nid] == 1)
					{
						lesskSet.push_back(nid);
						kcoreTag[nid] = 2;
					}
				}

				//find common neighbors
				if (kcoreSet[id].size() < kcoreSet[nid].size())
				{
					indexId = id;
					matchId = nid;
				}
				else
				{
					indexId = nid;
					matchId = id;
				}
				h = kcoreTable[matchId];
				for (long k = 0; k < kcoreSet[indexId].size(); k++)
				{
					long knid = kcoreSet[indexId][k];
					if (kcoreTag[knid] == 1)
					{
						kr = kh_get(32, h, knid);
						vidInMatchVer = kh_value(h, kr);
						if (kr != kh_end(h))
						{
							if (kcoreEdgeTag[indexId][k] && kcoreEdgeTag[matchId][vidInMatchVer]) //correctness! knid-indexid and knid-matchid exist
							{
								if (kcoreSetSupport[indexId][k] > inputS - 1)
								{
									//update neighboring edges of indexId-knid
									kcoreSetSupport[indexId][k]--;
									long vidInKnid = kcoreSetNeiIndex[indexId][k];
									kcoreSetSupport[knid][vidInKnid]--;
									if (kcoreSetSupport[indexId][k] == inputS - 1)
									{
										kcoreDegree[knid]--;
										if (kcoreDegree[knid] < iniK && kcoreTag[knid] == 1)
										{
											lesskSet.push_back(knid);
											kcoreTag[knid] = 2;
										}
										kcoreDegree[indexId]--;
										if (kcoreDegree[indexId] < iniK && kcoreTag[indexId] == 1)
										{
											lesskSet.push_back(indexId);
											kcoreTag[indexId] = 2;
										}
									}
								}

								if (kcoreSetSupport[matchId][vidInMatchVer] > inputS - 1)
								{
									//update neighboring edges of matchId-knid
									kcoreSetSupport[matchId][vidInMatchVer]--;
									long vidInKnid = kcoreSetNeiIndex[matchId][vidInMatchVer];
									kcoreSetSupport[knid][vidInKnid]--;
									if (kcoreSetSupport[matchId][vidInMatchVer] == inputS - 1)
									{
										kcoreDegree[knid]--;
										if (kcoreDegree[knid] < iniK && kcoreTag[knid] == 1)
										{
											lesskSet.push_back(knid);
											kcoreTag[knid] = 2;
										}
										kcoreDegree[matchId]--;
										if (kcoreDegree[matchId] < iniK && kcoreTag[matchId] == 1)
										{
											lesskSet.push_back(matchId);
											kcoreTag[matchId] = 2;
										}
									}

								}
							}
						}
					}
				}
			}
		}
	}

	ksCoreSize = kcoreSet.size() - lesskSet.size();
}

void algorithm()
{
	if (inputS + 1 > inputK)
	{
		iniK = inputS + 1;
	}
	else
	{
		iniK = inputK;
	}

	long degree;
	lesskSet.clear();
	for (long i = 0; i < verSet.size(); i++)
	{
		degree = verSet[i].size();
		verDegree.push_back(degree);
		if (degree < iniK)
		{
			verTag.push_back(2);
			lesskSet.push_back(i);
		}
		else
		{
			verTag.push_back(1);
		}
	}
	verSetKcoreIDs.resize(verOrigIDs.size(), -1);

	//compute k-core as the base
	double iterationStartTag = (double)clock() / CLOCKS_PER_SEC;
	computeKcore0();
	kCoreTimeTag = (double)clock() / CLOCKS_PER_SEC;
	kCoreTime = kCoreTimeTag - iterationStartTag;
	printf("k-core vertices:%ld\nk-core time:%.4lf\n", kcoreSet.size(), kCoreTime);
	kcoreVerSize = kcoreSet.size();

	//compute support for each edge in k-core
	computeEdgeSupport0();
	edgeSupportTimeTag = (double)clock() / CLOCKS_PER_SEC;
	edgeSupportTime = edgeSupportTimeTag - kCoreTimeTag;
	printf("edge support time:%.4lf\n", edgeSupportTime);
	lesskIniNum = lesskSet.size();

	//compute k,s-truss based on k-core
	computeKScore0();
	ksCoreTimeTag = (double)clock() / CLOCKS_PER_SEC;
	ksCoreTime = ksCoreTimeTag - iterationStartTag;
	printf("ks-core vertices:%ld\nks-core time:%.4lf\n", ksCoreSize, ksCoreTime);
}

void dataOutput(double &runtime)
{
	char *writeTemp;
	char record[100];
	FILE *fs;
	fs = fopen(outfile.c_str(), "a");
	char fsch = 0;

	if (fs == NULL)
	{
		printf("ERROR!");
		exit(1);
	}
	else
	{
		sprintf(record, "k%ld", inputK);
		fwrite(record, sizeof(*record), strlen(record), fs);
		fsch = putc('\t', fs);

		sprintf(record, "s%ld", inputS);
		fwrite(record, sizeof(*record), strlen(record), fs);
		fsch = putc('\t', fs);

		sprintf(record, "%.2lf", runtime);
		fwrite(record, sizeof(*record), strlen(record), fs);
		fsch = putc('\t', fs);

		sprintf(record, "kc%ld", kcoreVerSize);
		fwrite(record, sizeof(*record), strlen(record), fs);
		fsch = putc('\t', fs);

		sprintf(record, "ks%ld", ksCoreSize);
		fwrite(record, sizeof(*record), strlen(record), fs);
		fsch = putc('\t', fs);
	}
}

int main(int argc, char *argv[])
{
	//configuration
	long maxVerID = 100000000; //max vertex id
	verSetTag.resize(maxVerID, -1);
	infile = "dataset.txt";
	outfile = "result.txt";
	scanf("%ld %ld", &inputK, &inputS);

	//input data 
	dataInput();

	//algorithm start
	algStartTimeTag = (double)clock() / CLOCKS_PER_SEC;

	//compute ks-core
	algorithm();

	double runtime = (double)clock() / CLOCKS_PER_SEC - algStartTimeTag;
	printf("time: %lf\n", runtime);

	//write
	dataOutput(runtime);

	return 0;
}
