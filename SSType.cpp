#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <tr1/unordered_map>
#include "SSType.h"
#include "Utilities.h"

using namespace std;

void parseStride(string stride, vector< vector<int> >& sse, int pdbOffset)
{
	ifstream file;
	string line;
	file.open(stride.c_str());
	string curSSType = ""; // current ss type
	int curStartRes = -1;
	int resnum=0; // current res num
	string temp;
	int idCounter = 0;

	tr1::unordered_map< pair<int,int>, int > strands; // [strand id1 < strand id2, count]

	if (file.is_open())
	{
		while (file.good())
		{
			getline(file,line);
			trim(line);
			if (line.length() == 0)
				continue;
			stringstream tok(line);
			string key;
			tok >> key;
			if (key == "ASG")
			{
				string ssType;
				tok >> temp; // aa type
				tok >> temp; // chain id
				tok >> temp; // pdb res num
				resnum = atoi(temp.c_str())+pdbOffset;
				tok >> temp; // ordinal res num
				tok >> ssType; // ss type
				if (ssType != curSSType)
				{
					// end previous ss type;
					// add if helix (alpha,310,pi) or sheet; and satisfy min length
					if (curStartRes != -1 &&
							(curSSType == "H" ||
							 curSSType == "G" ||  // 3-10 helix
							 curSSType == "I" ||  // pi helix
							 curSSType == "E"))
					{
						int len = resnum-curStartRes;
						if ( (curSSType == "E" && len >  SHEETLOWERBOUNDLEN) ||
								(len > HELIXLOWERBOUNDLEN) )
						{
							// add
							vector<int> ss;
							if (curSSType != "E")
								ss.push_back(HELIX); // helix
							else
								ss.push_back(STRAND); // strand
							ss.push_back(curStartRes);
							ss.push_back(resnum-1);
							ss.push_back(idCounter);
							if (curSSType == "E")
								ss.push_back(DEFAULTID); // no sheet id yet
							sse.push_back(ss);
							idCounter++;
						}
					}
					// start new ss element
					curStartRes = resnum;
					curSSType = ssType;
				}
				// else continue reading this ss element
			}
			else
			{
				// add the last SS element
				if (curStartRes != -1 && key == "REM")
				{
					if (curStartRes != -1 &&
							(curSSType == "H" ||
							 curSSType == "G" ||
							 curSSType == "I" ||
							 curSSType == "E"))
					{
						int len = resnum-curStartRes;
						if ( (curSSType == "E" && len >  SHEETLOWERBOUNDLEN) ||
								(len > HELIXLOWERBOUNDLEN) )
						{
							// add
							vector<int> ss;
							if (curSSType != "E")
								ss.push_back(HELIX); // helix
							else
								ss.push_back(STRAND); // strand
							ss.push_back(curStartRes);
							ss.push_back(resnum-1);
							ss.push_back(idCounter);
							if (curSSType == "E")
								ss.push_back(DEFAULTID); // no sheet id yet
							sse.push_back(ss);
							idCounter++;
						}
					}
					curStartRes = -1; // to skip this if block next time
					continue;
				}
				if (key == "DNR" || key == "ACC")
				{
					tok >> temp; // resname1
					tok >> temp; // chain1
					tok >> temp; // resnum1
					int res1 = atoi(temp.c_str())+pdbOffset;
					tok >> temp; // ordinal1
					tok >> temp; // arrow
					tok >> temp; // resname2
					tok >> temp; // chain2
					tok >> temp; // resnum2
					int res2 = atoi(temp.c_str())+pdbOffset;
					// check if res1, res2 in different strands
					int ssType1 = LOOP;
					int ssType2 = LOOP;

					if (res1 > res2)
						std::swap(res1,res2);

					int ssid1=-1;
					int ssid2=-2;;
					// count strand interactions
					for (vector< vector<int> >::iterator it = sse.begin(); it != sse.end(); ++it)
					{
						vector<int>& ss = *it;
						int r1 = ss[STARTRESINDEX];
						int r2 = ss[ENDRESINDEX];
						if (res1 >= r1 && res1 <= r2)
						{
							ssType1 = ss[TYPEINDEX];
							ssid1 = ss[SSIDINDEX];
							if (ssType1 != STRAND)
								break;
						}
						if (res2 >= r1 && res2 <= r2)
						{
							ssType2 = ss[TYPEINDEX];
							ssid2 = ss[SSIDINDEX];
							if (ssType2 != STRAND)
								break;
						}
					}
					if (ssType1 == STRAND && ssType2 == STRAND && ssid1 != ssid2)
					{
						// increment the count of these pair of strands
						pair<int,int> strandPair(ssid1,ssid2);
						tr1::unordered_map< pair<int,int>, int >::iterator it = strands.find(strandPair);
						if (it != strands.end())
							it->second++;
						else
							strands[strandPair] = 1;
					}
				}
			}
		} // end while file
		// collect pairwise strands into sheets
		for (tr1::unordered_map< pair<int,int>, int >::iterator it = strands.begin(); it != strands.end(); ++it)
		{
			int count = it->second;
			if (count >= MINHBONDS)
			{
				const pair<int,int>& p = it->first;
				int ssid1 = p.first;
				int ssid2 = p.second;
				vector<int> sheet;
				sheet.push_back(SHEET);
				sheet.push_back(2); // 2 strands
				sheet.push_back(DEFAULTID); // sheet index to be set later
				sheet.push_back(ssid1);
				sheet.push_back(ssid2);
				sse.push_back(sheet);
			}
		} // end group strands pairwise
		// combine overlapping sheets; delete the sheet with the higher index when combining
		// assign sheet ids and then update the sheet id field of the strands
      for (vector< vector<int> >::iterator it1 = sse.begin(); it1 != sse.end(); ++it1)
		{
			vector<int>& ss1 = *it1;
			if (ss1[TYPEINDEX] != SHEET)
				continue;
			for (vector< vector<int> >::iterator it2 = sse.begin(); it2 != sse.end(); ++it2)
			{
				vector<int>& ss2 = *it2;
				if (ss2[TYPEINDEX] != SHEET)
					continue;
				if (ss1 == ss2)
					continue;
				bool combine = false; // true if should merge these two sheets
				for (int i = 0; i < ss1[NUMSTRANDSINDEX]; i++)
				{
					int ssid1 = ss1[STRANDSTARTINDEX+i];
					for (int j = 0; j < ss2[NUMSTRANDSINDEX]; j++)
					{
						int ssid2 = ss2[STRANDSTARTINDEX+j];
						if (ssid1 == ssid2)
						{
							combine = true;
							break;
						}
					}
					if (combine)
						break;
				}
				if (combine)
				{
					// add strands in ss2 to ss1 if not already there
					for (int j = 0; j < ss2[NUMSTRANDSINDEX]; j++)
					{
						int ssid2 = ss2[STRANDSTARTINDEX+j];
						bool containsSSID2 = false;
						for (int i = 0; i < ss1[NUMSTRANDSINDEX]; i++)
						{
							int ssid1 = ss1[STRANDSTARTINDEX+i];
							if (ssid1 == ssid2)
							{
								containsSSID2 = true;
								break;
							}
						}
						if (!containsSSID2)
						{
							ss1.push_back(ssid2);
							ss1[NUMSTRANDSINDEX]++;
						}
					}
					ss2[TYPEINDEX] = INVALIDSS; // mark for removal
					break;
				}
			} // end for ss2
		} // end for ss1
		// delete entries marked for removal;
		for (vector< vector<int> >::iterator it = sse.begin(); it != sse.end();)
		{
			vector<int>& ss = *it;
			if (ss[TYPEINDEX] == INVALIDSS)
			{
				// not efficient due to swapping, but number of strands is limited, so not a big issue
				// a swapping with last valid element approach will be better for a large # of strands
				it = sse.erase(it);
			}
			else
				++it;
		}
		// assign sheet ids
		for (vector< vector<int> >::iterator it = sse.begin(); it != sse.end(); ++it)
		{
			vector<int>& ss = *it;
			if (ss[TYPEINDEX] == SHEET)
			{
				int sheetIndex = idCounter;
				idCounter++;
				ss[SSIDINDEX_SHEET] = sheetIndex;
				for (int i = 0; i < ss[NUMSTRANDSINDEX]; i++)
				{
					int ssid = ss[STRANDSTARTINDEX+i];
					sse[ssid][SHEETIDINDEX] = sheetIndex;
				}
			}
		}
		// end collect strands
		// print sse including sheets
		printSSE(sse);
	}
	else
	{
		printf("Unable to open file %s\n",stride.c_str());
		file.close();
		exit(-1);
	}
	file.close();
}


void parseSSTypes(tr1::array<SSTYPE,MAXPROSIZE>& ssTypes, vector< vector<int> >& sse)
{
	SSTYPE curSSType = UNKNOWNSS; // for now, loops are ignored
	int startRes = 0; // starting resnum of the current ss type
	int resnum = 1;
	int ssIndex = 0;
	for(tr1::array<SSTYPE,MAXPROSIZE>::iterator it = ssTypes.begin(); it != ssTypes.end(); ++it)
	{
		SSTYPE type = *it;
		if (type != curSSType)
		{
			// possibly add ss
			if (curSSType != UNKNOWNSS)
			{
				int len = resnum-startRes;
				if ( (curSSType == STRAND && len >  SHEETLOWERBOUNDLEN) ||
						(curSSType == HELIX && len > HELIXLOWERBOUNDLEN) )
				{
					vector<int> ss;
					ss.push_back(curSSType); // helix
					ss.push_back(startRes);
					ss.push_back(resnum-1);
					ss.push_back(ssIndex);
					if (curSSType == STRAND)
						ss.push_back(DEFAULTID); // no sheet id yet
					sse.push_back(ss);
					ssIndex++;
				}
			}
			// start next ss
			if (type != UNKNOWNSS && type != LOOP)
			{
				startRes = resnum;
				curSSType = type;
			}
			else
				curSSType = UNKNOWNSS; // invalide current SS type, so won't add
		}
		resnum++;
	}
	// add last sse element
	if (curSSType != UNKNOWNSS)
	{
		int len = resnum-startRes;
		if ( (curSSType == STRAND && len >  SHEETLOWERBOUNDLEN) ||
				(curSSType == HELIX && len > HELIXLOWERBOUNDLEN) )
		{
			vector<int> ss;
			ss.push_back(curSSType); // helix
			ss.push_back(startRes);
			ss.push_back(resnum-1);
			ss.push_back(ssIndex);
			if (curSSType == STRAND)
				ss.push_back(DEFAULTID); // no sheet id yet
			sse.push_back(ss);
			ssIndex++;
		}
	}
	printSSE(sse);
}

void getContactSSType(vector< vector<int> >& sse, int res1, int res2, int& ssid1, int& ssid2)
{
	bool isStrand1 = false;
	bool isStrand2 = false;
	ssid1 = DEFAULTID;
	ssid2 = DEFAULTID;
	for (vector< vector<int> >::iterator itSS = sse.begin(); itSS != sse.end(); ++itSS)
	{
		vector<int>& ss = *itSS;
		if (ss[TYPEINDEX] == STRAND)
		{
			if (res1 >= ss[STARTRESINDEX] && res1 <= ss[ENDRESINDEX])
			{
				ssid1 = ss[SSIDINDEX];
				isStrand1 = true;
			}
			if (res2 >= ss[STARTRESINDEX] && res2 <= ss[ENDRESINDEX])
			{
				ssid2 = ss[SSIDINDEX];
				isStrand2 = true;
			}
		}
		else if (ss[TYPEINDEX] == HELIX)
		{
			if (res1 >= ss[STARTRESINDEX] && res1 <= ss[ENDRESINDEX])
				ssid1 = ss[SSIDINDEX];
			if (res2 >= ss[STARTRESINDEX] && res2 <= ss[ENDRESINDEX])
				ssid2 = ss[SSIDINDEX];
		}
		// loops are not considered; sheets will be considered below
	}
	if (ssid1 == -1 || ssid2 == -1)
		return;
	// if both strand and in same sheet, compare the strands
	// else treat any strands as sheets (if it has been assigned to a sheet)
	if (!isStrand1 || !isStrand2 || sse[ssid1][SHEETIDINDEX] != sse[ssid2][SHEETIDINDEX])
	{
		// if one is a strand, get it's sheet if != DEFAULTID
		if (isStrand1 && sse[ssid1][SHEETIDINDEX] != DEFAULTID)
			ssid1 = sse[ssid1][SHEETIDINDEX];
		if (isStrand2 && sse[ssid2][SHEETIDINDEX] != DEFAULTID)
			ssid2 = sse[ssid2][SHEETIDINDEX];
	} // else both in same sheet, so compare the strands
}


void getResidueList(vector< vector<int> >& sse, int ssid, list<int>& residues)
{
	vector<int>& ss = sse[ssid];
	if (ss[TYPEINDEX] == STRAND || ss[TYPEINDEX] == HELIX)
	{
		int startRes = ss[STARTRESINDEX];
		int endRes = ss[ENDRESINDEX];
		for (int j = startRes; j <= endRes; j++)
			residues.push_back(j);
	}
	else if (ss[TYPEINDEX] == SHEET)
	{
		for (int i = 0; i < ss[NUMSTRANDSINDEX]; i++)
		{
			int strandID = ss[STRANDSTARTINDEX+i];
			int startRes = sse[strandID][STARTRESINDEX];
			int endRes = sse[strandID][ENDRESINDEX];
			for (int j = startRes; j <= endRes; j++)
				residues.push_back(j);
		}
	}
}

int getSSLength(vector< vector<int> >& sse, int ssid)
{
	int len = 0;
	if (sse[ssid][TYPEINDEX] == SHEET)
	{
		for (int i = 0; i < sse[ssid][NUMSTRANDSINDEX]; i++)
		{
			int id = sse[ssid][STRANDSTARTINDEX+i];
			len += sse[id][ENDRESINDEX]-sse[id][STARTRESINDEX]+1;
		}
	}
	else
	{
		len += sse[ssid][ENDRESINDEX]-sse[ssid][STARTRESINDEX]+1;
	}
	return len;
}

void printSSE(vector< vector<int> >& sse)
{
	for (vector< vector<int> >::iterator it = sse.begin(); it != sse.end(); ++it)
	{
		vector<int>& ss = *it;
		if (ss[TYPEINDEX] == SHEET)
		{
			printf("ID: %d\n\tSHEET\n\tNumstrands: %d\n",ss[SSIDINDEX_SHEET],ss[NUMSTRANDSINDEX]);
			for (int i=0; i < ss[NUMSTRANDSINDEX]; i++)
			{
				printf("\tStrandID: %d\n",ss[STRANDSTARTINDEX+i]);
			}
		}
		else if (ss[TYPEINDEX] == HELIX)
		{
			printf("ID: %d\n\tHELIX\n\t%d-%d\n",ss[SSIDINDEX],ss[STARTRESINDEX],ss[ENDRESINDEX]);
		}
		else if (ss[TYPEINDEX] == STRAND)
		{
			printf("ID: %d\n\tSTRAND\n\t%d-%d\n\tSheetID: %d\n",ss[SSIDINDEX],ss[STARTRESINDEX],ss[ENDRESINDEX],ss[SHEETIDINDEX]);
		}
		else
			printf("Skipping loop\n");
	}
}
