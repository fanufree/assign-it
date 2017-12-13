#include "Main.h"
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <bitset>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <list>
#include <set>
#include <tr1/array>
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include "Score.h"
#include "SSType.h"
#include "NOE.h"
#include "NOECluster.h"
#include "CSProtein.h"
#include "Utilities.h"
#include "Atom.h"
#include "Residue.h"
#include "Rectangle.h"
#include "Assignment.h"
#include "PCAssignment.h"
#include "XplorFormat.h"
#include "ContactMap.h"
#include "ContactMapIterator.h"
#include "IfElseCondition.h"
#include "TestCondition.h"
#include "AndCondition.h"
#include "OrCondition.h"
#include "ScoreCountCondition1.h"
#include "ScoreCountCondition2.h"
#include "ScoreTermCondition1.h"
#include "ScoreTermCondition2.h"
#include "InStructureCondition.h"
#include "InStructureConditionFract.h"
#include "InStructureConditionFractTemp.h"
#include "AssignedCondition.h"
#include "AssignedCondition2.h"
#include "TrueCondition.h"
#include "FalseCondition.h"
#include "SeqSepCondition.h"
#include "NotCondition.h"
#include "InStructureCutCondition.h"
#include "CorrectFilter.h"
#include "InStructureWindow.h"
#include "Interval.h"
#include "NoiseCondition.h"
#include "NumAssignmentsCondition.h"


ILOSTLBEGIN

using namespace std;

// distance calibration constants
// 12 constants
// among the pairs backbone protein, beta proton, methyl proton, other sidechain protein
const int NUM_CALIB_TYPES = 10;
enum CALIBRATION_TYPE {CALIB_BBBB, CALIB_BBBE, CALIB_BEBE,
	CALIB_BBME, CALIB_BBSC,
	CALIB_BEME, CALIB_BESC,
	CALIB_MEME, CALIB_MESC,
	CALIB_SCSC};
std::string CALIBRATION_TYPE_NAMES[] = {"CALIB_BBBB", "CALIB_BBBE", "CALIB_BEBE",
		"CALIB_BBME", "CALIB_BBSC",
		"CALIB_BEME", "CALIB_BESC",
		"CALIB_MEME", "CALIB_MESC",
		"CALIB_SCSC"};

const int NUM_CALIB_HALF = 4;
enum CALIBRATION_TYPE_HALF {
	BB_CALIB,
	BE_CALIB,
	ME_CALIB,
	SC_CALIB
};
std::string CALIBRATION_TYPE_HALF_NAMES[] = {
	"BB_CALIB",
	"BE_CALIB",
	"ME_CALIB",
	"SC_CALIB"
};

// used by expectedAssign table below
class ExpectedAssignKey
{
public:
	int r1; // resnum or AA type depending on the type of map
	string h1; // proton name
	int r2;
	string h2;
	ExpectedAssignKey() : r1(-1), h1(""), r2(-1), h2("")
	{
	}
	// if order is true, will ensure that r1 <= r2 and if r1==r2, h1 <= h2 alphabetically
	ExpectedAssignKey(int res1, string p1, int res2, string p2, bool order) : r1(res1), h1(p1), r2(res2), h2(p2)
	{
		if (order)
		{
			if (res1 > res2)
			{
				std::swap(r1,r2);
				std::swap(h1,h2);
			}
			else if (res1 == res2 && h1 > h2)
			{
				std::swap(h1,h2);
			}
		}
	}
	ExpectedAssignKey(int res1, string p1, int res2, string p2) : r1(res1), h1(p1), r2(res2), h2(p2)
	{
	}
	// use default copy and assignment operator
	bool operator==(const ExpectedAssignKey& k) const
	{
		return r1 == k.r1 && h1 == k.h1 && r2 == k.r2 && h2 == k.h2;
	}
	bool operator!=(const ExpectedAssignKey& k) const
	{
		return !(*this == k);
	}
};
namespace std
{
	namespace tr1
	{
		template <>
		class hash<ExpectedAssignKey>
		{
			public:
			size_t operator()(const ExpectedAssignKey& k) const
			{
				return hash<int>()(k.r1) ^ hash<string>()(k.h1) ^ hash<int>()(k.r2) ^ hash<string>()(k.h2);
			}
		};
	}
}

// used by scscAvgDist;
class SCSCDistKey
{
public:
	int r1; // resnum or AA type depending on the type of map
	string h1;
	int r2;
	string h2;
	int bin;
	SCSCDistKey() : r1(-1), h1(""), r2(-1), h2(""), bin(-1)
	{
	}
	SCSCDistKey(int res1, string p1, int res2, string p2, int distBin) : r1(res1), h1(p1), r2(res2), h2(p2), bin(distBin)
	{
	}
	// use default copy and assignment operator
	bool operator==(const SCSCDistKey& k) const
	{
		return r1 == k.r1 && h1 == k.h1 && r2 == k.r2 && h2 == k.h2 && bin == k.bin;
	}
	bool operator!=(const SCSCDistKey& k) const
	{
		return !(*this == k);
	}
};
namespace std
{
	namespace tr1
	{
		template <>
		class hash<SCSCDistKey>
		{
			public:
			size_t operator()(const SCSCDistKey& k) const
			{
				return hash<int>()(k.r1) ^ hash<string>()(k.h1) ^ hash<int>()(k.r2) ^ hash<string>()(k.h2) ^ hash<int>()(k.bin);
			}
		};
	}
}

/*
 * BEGIN GLOBAL VARIABLES
 * TODO: place the input stuff in its own class, so easier to pass to other classes
 */
const int LRSEQSEP=6; // see Main.h
double LBDIST = 1.8; // lower bound on proton-proton distance for setting restraints; for methyls it is LBDIST-0.5; default value 1.8
double DISTCUTOFF = 6.0; // contact distance cutoff; set in parseConfigFile

bitset<MAXPROSIZE> isDynamic; // [resnum-1] = 1 if resnum is dynamic based on TALOS+
tr1::array<SSTYPE,MAXPROSIZE> ssTypes; // predicted SS type from chem shifts; index = [resnum-1]
// secondary structure of template structures; these are calculated by stride from the 3D pdb structure
vector< vector< vector<int> > > sses; // Each vector element consists of the ss type info for one template structure
                  // see parseStride in SSType.h for description of the data structure (sse) containing the ss type info

// default chem shift tolerances
const double tolx = 0.35;
const double tolhx = 0.035;
const double tolh = 0.035; // direct dim
const double tolx2 = 0.35; // for 4D NOE
// chemical shift error scaling factor for scoring function;
const double factorX1=1.96/(0.5*sqrt(2.0)); // an error of 0.5 corresponds to a 0.05 score
const double factorH1=1.96/(0.05*sqrt(2.0));
const double factorX2=1.96/(0.5*sqrt(2.0));
const double factorH2=1.96/(0.05*sqrt(2.0));
tr1::array< list<NOE*>, NOE::NUMTYPES > peakLists; // indexed by PEAKLISTTYPE
vector<NOE::PEAKLISTTYPE> peakListTypes;
tr1::array< tr1::array<double,4>, NOE::NUMTYPES > tols;  // tolerances for each peak list (to get the assignment possibilities from the chemical shift assignment)
                             // indexed by PEAKLISTTYPE; for 3D order is x,hx,INVALIDSHIFTVAL,h; for 4D order is x1 hx1 x2 hx2;
                             // for 2D it is INVALIDSHIFTVAL,h1,INVALIDSHIFTVAL,h2, where entries with INVALIDSHIFTVAL are ignored

int startRes; // residue range for filtering the assigned contacts (inclusive); default is the length of the sequence
int endRes;
string seq; // amino acid sequence; set by parseConfigFile
Residue::AA3 seqAAType[MAXPROSIZE]; // 0-based, set by parseConfigFile

// chemical shift assignment
CSProtein* bmrbProtein = NULL;

// template 3D protein structure
list<CSProtein*> structures;

// reference PDB if available
list<CSProtein*> refPDB;

// TODO: Water band handling (currently not implemented)
const double WATERLOWER = 4.3;
const double WATERUPPER = 5.1;

// distance bins
const int NUMDISTBINS = 3; // used by getResidueBasedDistanceIntensity()
const double UPPERDISTBINS[] = {3.0, 4.5, 6.0}; // used by getResidueBasedDistanceIntensity()
const double INTENSITY_DISTCUT = 1.5; // cutoff for distance bin violation; used in score.intensity
const double STDEV_DISTCUT = 0.3; // cutoff for standard deviation of distances in the structures; for a given peak
                          // the distances in the assignment possibilities must have stdev less than this for
                         // this peak to be used for calibration
const double INVALIDDISTANCE = 9999999; // see Main.h

// contact maps
ContactMap contactMap; // from input structures, store fraction of
                  // structures with (directed) contact (w/ and w/o intersectSphere); minimum distance among the structures (can be INVALIDDISTANCE if no coords exist)
                 // avg and stdev of distance (INVALIDDISTANCE if coords do not exist in one of the models)
                // the map is symmetric (both directions are included); use r1 < r2 && if r1==r2 h1 < h2 to iterate over all pairs of protons; methyls represented by 1 proton (the first one, e.g. ALA HB1)
ContactMap contactMapExpected; // Similar to contactMap, but takes into account chemical shift completeness and peak list types.
                              // Contacts from template are not added if this contact cannot be generated from the experimental data
                              // This map is NOT symmetric because the peak lists may not be.
int numLRContacts; // >= LRSEQSEPfrom and < 4.0A; initialized in contactMapExpected; each contact treated as undirected, so if have contact and symmetric contact, then they count as one
                   //  (r2-r1) >= LRSEQSEP && minDist < 4.0 && entry[FRACSTRUCINDEX] >= 0.199)
                   // [NOT THAT IMPORTANT; can delete]

ContactMap refContactMap; // contact map of reference structure if known; used for debugging

// number of inter-residue contacts stats
tr1::unordered_map<ExpectedAssignKey,int> expectedAssign; // [resnum1 h1 resnum2 h2] = expected # of assignments between resnum1 and resnum2 given h1, h2 in contact;
                              // depends on peak list types and cs completeness; used in interres score
int expectedAssignAvg[MAXPROSIZE][MAXPROSIZE]; // [resnum1-1][resnum2-1] = avg # of contacts between res1 and res2 based on expectedAssign stats (excluding the entries with zeros)
                                               // basically it is the avg over h1 h2 of [r1 h1 r2 h2] in expectedAssign; used in netCS score

// side-side chain distance stats
tr1::unordered_map<SCSCDistKey,double> scscAvgDist; // [resnum1 h1 resnum2 h2 distbin] avg {SC,n,ca,cb}-{SC,n,ca,b} distance + SCSC_STDEV_FACTOR*max(stdev,SCSC_MIN_STDEV)
                                                    // given h1, h2 in contact within distbin;
                                            // dist bins <3, 4.5, <6 (same as UPPERDISTBINS)
// double scscContactDist[MAXPROSIZE][MAXPROSIZE]; // average sc-sc distance in the templates (excludes CA, CB); for GLY CA returned, and for ALA CB coords returned.
double SCSC_STDEV_FACTOR=2.0; // used to initialize scscAvgDist
double SCSC_MIN_STDEV = 0.5;  // used to initialize scscAvgDist

// network score parameters
const int NETWORKWINDOW = 2; // number of residues flanking residue i to consider for network score; 1 = consider sequential neighor +/- 1 residues

// assignment filtering parameters [TO BE DELETED; replaced by parameters in configuration file]
double FILTERNETPOSS = 0.12; // for filtering assignment possibilities by net score
double FILTERINTERPOSS = 0.15; // for filtering assignment possibilities by inter res score
double LONGDISTCUTFRAC = 0.4; // for filtering assignment possibilities. Uses current fraction of long range assignment possibilities without any structure evidence to decide on cutoff to prune
double FILTERNET = 0.2;  // for filtering assigned contacts by net score
double FILTERINTER = 0.2; // for filtering assigned contacts by inter res score
double REFINE_FRACTSTRCUT = 0.199; // 0.04 = 2/50; fraction of structures that assigned contact must satisfy within DISTCUTOFF.
								    // Used with above criteria to filter assignment output before it can be included in noe_refine.tbl (xplor)
                                   // or for noe.tbl (I-TASSER); this filtering is done after assignment step, but before the output step
double INIT_FRACTSTRCUT = 0;  // minimum fraction of input structures containing each contact, which satisfies DISTCUTOFF. Used to filter out assignment possibilities (OBSOLETE DO NOT USE)
bool REFINE_NETCUT = true;// for assigned contacts not satisfying the above criteria; if true, such contacts can satisfy network criteria for inclusion in new.tbl (xplor) or ambig.tbl (I-TASSER)
double TEMPFRACFILTER = 0.3; // fraction of templates containing long-range assigned contact; used to filter out assignments based on first sorting assignments by score in decreasing order
                          // then selecting the assignments in decreasing score order until this fraction of assigned contacts are not satisified by any template. Selection then stops here.
                          // the smaller the fraction, the more filtering is done
double TEMPFRACEVAL = 0.1;  // similar to TEMPFRACFILTER, but used to evaluate fit of templates to the NOE data; this sum of the scores of these assignments are outputted

// output parameters common to both xplor and I-TASSER
double WEAK_RESTRAINT_WT = 0.333; // scaling factor on the scores for restraints with weak confidence; such restraints will be set as ambiguous restraints
double AMBIG_SCORE_CUT = 1.0; // residue pairs with score sum < than this will be treated as an ambiguous restraint

// parameters for setting the distance in the restraints
double WT_STRDIST = -1; // if < 0, it will be set dynamically based on stdev of dist in templates; the distance is set to WT_STRDIST*<avg dist in templates> + WT_INTDIST*<intensity calibrated distance>
double WT_INTDIST = -1;

// distance calibration from peak volume/intensity parameter
double MAX_CALIB_ERROR = 0.1; // fraction of distances in structure that violate the calibrated intensity-based distance upperbounds

// xplor-only output parameters
double REFINE_FRACTTEMPCUT = 0.399; // fraction of structures containing long-range contact (LRSEQSEP), not already assigned for inclusion in template.tbl
double TEMPLATE_DISTCUT = 5.0; // distance cutoff for restraints from templates for template.tbl output

// parameters for outputting ambiguous assignments
double AMBIG_FILTER = 999999.0;  // assignment possibilities are sorted by score for each NOE, and then we keep assignments until the fraction of the total score exceeds this cutoff
int AMBIG_FILTER_MAX = 1000000; // maximum number of assignment possibilities as set by AMBIG_FILTER (used only if AMBIG_FILTER_BY_COUNT is false)
bool AMBIG_FILTER_BY_COUNT = false; // if true filter by count instead of by score
double POOL_FILTER = -1; // ambiguous assignments with score sum < this are placed in a pool for constraint combination; set to < 0 to disable filter
bool POOL_FILTER_BY_COUNT = false; // if true do constraint combination by count instead of by score

// balance between old and new value of score terms interres, net, and netstr; used in recompute score step
double OLDWT = 0; // weight on the old value; weight on new value is 1.0-OLDWT; default OLDWT is 0
double NEWWT = 1.0;

// other parameters
// Weight training
string weightTrainFile = ""; // if non-null, NOE assignment is skipped and assignment possibilities are outputted to a
                             // user specified file; requires refPDB to be non-empty
int initialQID = 1; // used by writeTrainingData for svm_rank software; the first qid to output

TestCondition* filterAss = NULL; // for filtering assignment output
string filterAssStr = "";
TestCondition* filterAss2 = NULL; // for filtering assignment output
string filterAssStr2 = "";
TestCondition* filterAssPos = NULL; // for filtering assignment possibilities
string filterAssPosStr = "";
TestCondition* filterAmbPos = NULL; // for filtering ambiguous assignment possibilities
string filterAmbPosStr = "";

// assignment stats
int numContactsWithNOEMatchLR = 0; // # number of long range contacts in templates with at least one NOE assignment possibility

// assignment output
int assignMatrix[MAXPROSIZE][MAXPROSIZE] = {{0}}; // # of assignments between resnumA-1 to resnumB-1; all elements initialized to 0
                                                 // matrix is symmetric
// asssignment possibilities
// # of assignment possibilities between resnumA-1 to resnumB-1; all elements initialized to 0
// we count each contact between atoms i and j from resnumA and resnumB at most once even though i j may have
// many NOE peaks that match i j
// initialized by initializeAssignmentFilter, and initializeAmbigAssignmentFilter
// matrix is symmetric
int assignPossibMatrix[MAXPROSIZE][MAXPROSIZE] = {{0}}; // for filtering

// maximum number of iteration cycles in assign()
// if equal to 1, then only a single assignment step is done
int MAXITER = 1;

// Used in parseScoreCount/Term
list<PCAssignment*>* assignPossibPtr = NULL;
vector<PCAssignment*>* assignedPtr = NULL;
//

/*
 *
 * END GLOBAL VARIABLES
 *
 */


class PseudoContactKey // for hash table of PseudoContacts
{
public:
	NOE* noe;
	int  resnum1;
	int  resnum2;
	PseudoContactKey() : noe(NULL), resnum1(0), resnum2(0)
	{
	}
	PseudoContactKey(NOE* n, int r1, int r2) :  noe(n), resnum1(r1), resnum2(r2)
	{
	}
	// use default copy constructor and assignment operator
	bool operator==(const PseudoContactKey& k) const
	{
		if (noe != k.noe || resnum1 != k.resnum1 || resnum2 != k.resnum2)
			return false;
		else
			return true;
	}
	bool operator!=(const PseudoContactKey& k) const
	{
		return !(*this == k);
	}
};
namespace std
{
	namespace tr1
	{
		template <>
		class hash<PseudoContactKey>
		{
			public:
			size_t operator()(const PseudoContactKey& k) const
			{
				return hash<long>()(long(k.noe)) ^ k.resnum1 ^ k.resnum2;
			}
		};
	}
}

// Note: There is no distinction in order, so CALIB_BBSC can mean h1 is SC or h2 is SC
// used by calibrateDistance
CALIBRATION_TYPE getCalibContactType(const Contact& c)
{
	Residue* r1 = c.r1;
	Residue* r2 = c.r2;
	Atom* x1 = c.x1;
	Atom* x2 = c.x2;

	if (x1->isBackbone())
	{
		if (x2->isBackbone())
			return CALIB_BBBB;
		else if (x2->isCB() && r2->type != Residue::ALA)
			return CALIB_BBBE;
		else if (x2->isMethyl())
			return CALIB_BBME;
		else
			return CALIB_BBSC;
	}
	else if (x1->isCB() && r1->type != Residue::ALA)
	{
		if (x2->isBackbone())
			return CALIB_BBBE;
		else if (x2->isCB() && r2->type != Residue::ALA)
			return CALIB_BEBE;
		else if (x2->isMethyl())
			return CALIB_BEME;
		else
			return CALIB_BESC;
	}
	else if (x1->isMethyl())
	{
		if (x2->isBackbone())
			return CALIB_BBME;
		else if (x2->isCB() && r2->type != Residue::ALA)
			return CALIB_BEME;
		else if (x2->isMethyl())
			return CALIB_MEME;
		else
			return CALIB_MESC;
	}
	else // sidechain
	{
		if (x2->isBackbone())
			return CALIB_BBSC;
		else if (x2->isCB() && r2->type != Residue::ALA)
			return CALIB_BESC;
		else if (x2->isMethyl())
			return CALIB_MESC;
		else
			return CALIB_SCSC;
	}
	return CALIB_SCSC;
}

/**
 * Distance calibration step
 * Sets the constants used for setting the upperbound intensity-based distance
 * This is done for each peak list type
 */
void calibrateDistances(list<NOE*>& noes, vector<PCAssignment*>& seed,
		double coeffs[NOE::NUMTYPES][NUM_CALIB_TYPES])
{
	double MINUS_1_6 = -1.0/6.0;
	const unsigned int MIN_CALIB_COUNT = 3; // for each calibration type, need at least this many contacts to use allCoeffs

	// List of all coefficient values for each calibration type. The list is sorted in ascending order
	// [peak list type][CALIBRATION_TYPE][vector of coefficents from distances in template structures]
	tr1::array< tr1::array< vector<double>, NUM_CALIB_TYPES >, NOE::NUMTYPES > allCoeffs;
	for (int t = 0; t < NOE::NUMTYPES; t++)
	{
		for (int i = 0; i < NUM_CALIB_TYPES; i++)
		{
			vector<double> lst;
			allCoeffs[t][i] = lst;
			coeffs[t][i] = 0; // use 0 to denote no calibration constant
		};
	}
	// from seed assignment, use only the local and sequential contacts for calibration
	// since they are expected to have low variance
	for (vector<PCAssignment*>::iterator it = seed.begin(); it != seed.end(); ++it)
	{
		PCAssignment* a = *it;
		NOE* noe = a->noe;
		double intensity = noe->volume;
		if (intensity <= 0)
			continue;
		Contact first = a->pc.getOneContact();
		int seqSep = abs(first.r1->num - first.r2->num);
		if (seqSep < 2)
		{
			tr1::unordered_set<Contact>& contacts = a->pc.getContacts();
			double minDist = 999999;
			Contact minDistContact;
			for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
			{
				const Contact& c = *itC;
				double dist = contactMap[c][ContactMap::MINDISTINDEX];
				double stdev = contactMap[c][ContactMap::STDEVINDEX];
				if (dist <= DISTCUTOFF && stdev <= 0.5)
				{
					if (dist < minDist)
					{
						minDist = dist;
						minDistContact = c;
					}
				}
				else
					continue;
			}
			if (minDist <= DISTCUTOFF)
			{
				CALIBRATION_TYPE contactType = getCalibContactType(minDistContact);
				double k = pow(minDist,-6.0)/intensity; // calibration factor for isolated spin approximation
				allCoeffs[noe->type][contactType].push_back(k);
			}
		} // end if seqSep < 2
	} // end for each seed assignment

	// sort each vector in allCoeffs
	for (int t = 0; t < NOE::NUMTYPES; t++)
	{
		if (allCoeffs[t].size() > 0)
		{
			for (int i = 0; i < NUM_CALIB_TYPES; i++)
			{
				if (allCoeffs[t][i].size() > 0)
				{
					sort(allCoeffs[t][i].begin(), allCoeffs[t][i].end()); // ascending order
				}
			}
		}
	}
	tr1::array< tr1::array< double, NUM_CALIB_HALF >, NOE::NUMTYPES > increment; // the step size for searching for best coefficents; for half coefficients
	tr1::array< tr1::array< double, NUM_CALIB_TYPES >, NOE::NUMTYPES > stepSize; // the step size for searching for best coefficents; for full coefficients
	// find full coefficients kij first, then use half kii, kjj if not enough data on full coefficients

	for (int t = 0; t < NOE::NUMTYPES; t++)
	{
		for (int i = 0; i < NUM_CALIB_HALF; i++)
		{
			increment[t][i] = 0;
		}
		for (int i = 0; i < NUM_CALIB_TYPES; i++)
		{
			stepSize[t][i] = 0;
		}
	}
	double kCalibCoeffs[NOE::NUMTYPES][NUM_CALIB_HALF];
	for (int t = 0; t < NOE::NUMTYPES; t++)
	{
		for (int i = 0; i < NUM_CALIB_HALF; i++)
			kCalibCoeffs[t][i] = 0;
	}
	const int MAXITER = 40; // maximum number of steps for searching for best calibration constant
	double FRACTITER = 1.0/double(MAXITER);
	double FRACTMIN = 1.0/50.0; // used for setting the increment[][] below

	// for determining bb, be, me, sc
	CALIBRATION_TYPE toGet[] = {CALIB_BBBB, CALIB_BEBE, CALIB_MEME, CALIB_SCSC};
	CALIBRATION_TYPE_HALF toGetHalf[] = {BB_CALIB, BE_CALIB, ME_CALIB, SC_CALIB}; // same order as in Main.h
	CALIBRATION_TYPE toGetRelated[][3] = { // same order as toGet
			{CALIB_BBBE, CALIB_BBME, CALIB_BBSC},
			{CALIB_BBBE, CALIB_BEME, CALIB_BESC},
			{CALIB_BBME, CALIB_BEME, CALIB_MESC},
			{CALIB_BBSC, CALIB_BESC, CALIB_MESC}};
	CALIBRATION_TYPE_HALF toGetRelatedHalf[][3] = { // should match with toGetRelated
			{BE_CALIB, ME_CALIB, SC_CALIB},
			{BB_CALIB, ME_CALIB, SC_CALIB},
			{BB_CALIB, BE_CALIB, SC_CALIB},
			{BB_CALIB, BE_CALIB, ME_CALIB}};
	for (int t = 0; t < NOE::NUMTYPES; t++) // try to set kCalibCoeffs for BB_CALIB and BE_CALIB
	{
		if (peakLists[t].size() == 0)
			continue;
		printf("Starting distance calibration of peak list %s\n",NOE::PEAKTYPENAMES[t].c_str());
		for (int i = 0; i < NUM_CALIB_TYPES; i++)
		{
			vector<double>& iCoeffs = allCoeffs[t][i];
			if (iCoeffs.size() >= MIN_CALIB_COUNT)
			{
				int midIndex = floor(iCoeffs.size()/2);
				double midValue = iCoeffs[midIndex];
				double minValue = iCoeffs[0];
				coeffs[t][i] = midValue;
				stepSize[t][i] = FRACTITER*(midValue-FRACTITER*minValue);
			}
		}
		for (int i = 0; i < NUM_CALIB_HALF; i++)
		{
			vector<double>& iCoeffs = allCoeffs[t][toGet[i]];
			if (iCoeffs.size() > 0)
			{
				int midIndex = floor(iCoeffs.size()/2);
				double midValue = iCoeffs[midIndex];
				int minIndex = (midIndex > 1 ? 1 : 0);
				double minValue = iCoeffs[minIndex];
				kCalibCoeffs[t][toGetHalf[i]] = midValue;
				if (midIndex != minIndex)
					increment[t][toGetHalf[i]] = FRACTITER*(midValue-FRACTMIN*minValue);
				else
					increment[t][toGetHalf[i]] = FRACTITER*midValue;
			}
		}
		// if bb, be, me, or sc cannot be set by bbbb, bebe, meme, scsc, try other contact types
		for (int i = 0; i < NUM_CALIB_HALF; i++)
		{
			if (kCalibCoeffs[t][toGetHalf[i]] == 0)
			{
				// pick first one that works
				for (int j = 0; j < NUM_CALIB_HALF-1; j++)
				{
					CALIBRATION_TYPE type = toGetRelated[i][j];
					vector<double>& iCoeffs = allCoeffs[t][type];
					if (iCoeffs.size() > 0)
					{
						CALIBRATION_TYPE_HALF typeHalf = toGetRelatedHalf[i][j];
						int midIndex = floor(iCoeffs.size()/2);
						double midValue = iCoeffs[midIndex];
						int minIndex = (midIndex > 1 ? 1 : 0);
						double minValue = iCoeffs[minIndex];
						double midValue2 = midValue*midValue/kCalibCoeffs[t][typeHalf];
						kCalibCoeffs[t][toGetHalf[i]] = midValue2;
						if (midIndex != minIndex)
						{
							double minValue2 = minValue*minValue/kCalibCoeffs[t][typeHalf];
							increment[t][toGetHalf[i]] = FRACTITER*(midValue2-FRACTMIN*minValue2);
						}
						else
							increment[t][toGetHalf[i]] = FRACTITER*midValue2;
						break;
					}
				}
			}
		}
		// check if missing increment for some entries in toGet[]
		for (int i = 0; i < 4; i++)
		{
			if (increment[t][toGetHalf[i]] == 0)
			{
				increment[t][toGetHalf[i]] = FRACTITER*kCalibCoeffs[t][toGetHalf[i]];
			}
		}
		// from bb,be,me,sc constants, get the joint constants if does not already exist
		for (int i = 0; i < NUM_CALIB_TYPES; i++)
		{
			if (coeffs[t][i] > 0)
				continue; // already set
			switch (i)
			{
				case CALIB_BBBB:
					if (kCalibCoeffs[t][BB_CALIB] > 0) // kCalib is also 0 if no peak volume in peak list
						coeffs[t][CALIB_BBBB] = kCalibCoeffs[t][BB_CALIB];
					break;
				case CALIB_BBBE:
					if (kCalibCoeffs[t][BB_CALIB] > 0 && kCalibCoeffs[t][BE_CALIB] > 0)
						coeffs[t][CALIB_BBBE] = sqrt(kCalibCoeffs[t][BB_CALIB]*kCalibCoeffs[t][BE_CALIB]);
					break;
				case CALIB_BEBE:
					if (kCalibCoeffs[t][BE_CALIB] > 0)
						coeffs[t][CALIB_BEBE] = kCalibCoeffs[t][BE_CALIB];
					break;
				case CALIB_BBME:
					if (kCalibCoeffs[t][BB_CALIB] > 0 && kCalibCoeffs[t][ME_CALIB] > 0)
						coeffs[t][CALIB_BBME] = sqrt(kCalibCoeffs[t][BB_CALIB]*kCalibCoeffs[t][ME_CALIB]);
					break;
				case CALIB_BBSC:
					if (kCalibCoeffs[t][BB_CALIB] > 0 && kCalibCoeffs[t][SC_CALIB] > 0)
						coeffs[t][CALIB_BBSC] = sqrt(kCalibCoeffs[t][BB_CALIB]*kCalibCoeffs[t][SC_CALIB]);
					break;
				case CALIB_BEME:
					if (kCalibCoeffs[t][BE_CALIB] > 0 && kCalibCoeffs[t][ME_CALIB] > 0)
						coeffs[t][CALIB_BEME] = sqrt(kCalibCoeffs[t][BE_CALIB]*kCalibCoeffs[t][ME_CALIB]);
					break;
				case CALIB_BESC:
					if (kCalibCoeffs[t][BE_CALIB] > 0 && kCalibCoeffs[t][SC_CALIB] > 0)
						coeffs[t][CALIB_BESC] = sqrt(kCalibCoeffs[t][BE_CALIB]*kCalibCoeffs[t][SC_CALIB]);
					break;
				case CALIB_MEME:
					if (kCalibCoeffs[t][ME_CALIB] > 0)
						coeffs[t][CALIB_MEME] = kCalibCoeffs[t][ME_CALIB];
					break;
				case CALIB_MESC:
					if (kCalibCoeffs[t][ME_CALIB] > 0 && kCalibCoeffs[t][SC_CALIB] > 0)
						coeffs[t][CALIB_MESC] = sqrt(kCalibCoeffs[t][ME_CALIB]*kCalibCoeffs[t][SC_CALIB]);
					break;
				case CALIB_SCSC:
					if (kCalibCoeffs[t][SC_CALIB] > 0)
						coeffs[t][CALIB_SCSC] = kCalibCoeffs[t][SC_CALIB];
					break;
				default:
					break;
			}
		}

		// fix some potential calibration problems
		if (kCalibCoeffs[t][BB_CALIB] > 0 && kCalibCoeffs[t][ME_CALIB] > 0.33333*kCalibCoeffs[t][BB_CALIB])
			kCalibCoeffs[t][ME_CALIB] = 0.33333*kCalibCoeffs[t][BB_CALIB];
		if (kCalibCoeffs[t][BE_CALIB] > 0 && kCalibCoeffs[t][SC_CALIB] > 0.66667*kCalibCoeffs[t][BE_CALIB])
			kCalibCoeffs[t][SC_CALIB] = 0.66667*kCalibCoeffs[t][BE_CALIB];

		// test the constants for errors; decrease constants by increment if does not satisfy fractError <= MAX_CALIB_ERROR
		double fractError = 1.0;
		int numIter = 0;
		while (numIter < MAXITER)
		{
			int seedCountedType[NUM_CALIB_TYPES];
			int numErrorsType[NUM_CALIB_TYPES];
			int distCounts[3]; // [0] < 3A, 4.5A, 6A
			distCounts[0]=0;
			distCounts[1] = 0;
			distCounts[2] = 0;
			for (int i = 0; i < NUM_CALIB_TYPES; i++)
			{
				seedCountedType[i] = 0;
				numErrorsType[i] = 0;
			}
			double minCalibDist = 99999;
			int numErrorsNative = 0;

			// check fit to templates
			int numErrors = 0;
			int totalSeedCounted = 0;
			for (vector<PCAssignment*>::iterator it = seed.begin(); it != seed.end(); ++it)
			{
				PCAssignment* a = *it;
				NOE* noe = a->noe;
				if (noe->type != t)
					continue; // ignore seed assignments for other peak lists
				double intensity = noe->volume;
				if (intensity <= 0)
					continue; // cannot calibrate if no intensity
				tr1::unordered_set<Contact>& contacts = a->pc.getContacts();
				double calibratedDist = 0; // the intensity-based distance upperbound
				int count = 0;
				double minDist = 999999;

				double minDistNative = 999999;
				double avgStdev = 0;
				double avgStdevCount = 0;

				for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
				{
					const Contact& c = *itC;
					double dist = contactMap[c][ContactMap::MINDISTINDEX];
					if (dist > DISTCUTOFF)
						continue;
					CALIBRATION_TYPE contactType = getCalibContactType(c);
					double ki = coeffs[noe->type][contactType];
					if (ki > 0)
					{
						calibratedDist += ki;
						count++;
					}
					else
						continue;
					if (dist < minDist)
						minDist = dist;

					double stdev = contactMap[c][ContactMap::STDEVINDEX];
					avgStdev += stdev;
					avgStdevCount++;
					if (refPDB.size() > 0)
					{
						double nativeDist = refContactMap[c][ContactMap::MINDISTINDEX];
						if (nativeDist < minDistNative)
							minDistNative = nativeDist;
					}
				}
				if (count == 0)
					continue;
				if (minDist > DISTCUTOFF)
					continue;
				totalSeedCounted++;

				CALIBRATION_TYPE contactType = getCalibContactType(a->pc.getOneContact());
				seedCountedType[contactType]++;

				double wt = 1.0/double(count); // based on equal weight for each contact in pc
				calibratedDist = pow(calibratedDist*intensity*wt, MINUS_1_6);

				if (minDistNative < 999999 && minDistNative > calibratedDist)
					numErrorsNative++;
				avgStdev = avgStdev/avgStdevCount;

				if (avgStdev > 1.0)
					avgStdev = 1.0;
				if (minDist+avgStdev > calibratedDist) // minDist > calibratedDist)
				{
					numErrors++;
					numErrorsType[contactType]++;
				}
				if (calibratedDist < 3.0)
					distCounts[0]++;
				else if (calibratedDist < 4.5)
					distCounts[1]++;
				else if (calibratedDist <= DISTCUTOFF)
					distCounts[2]++;

				if (calibratedDist < minCalibDist) {
					minCalibDist = calibratedDist;
				}
			} // end for each seed assignment
			if (totalSeedCounted == 0)
			{
				numIter++;
				continue;
			}
			fractError = double(numErrors)/(double)(totalSeedCounted);
			if (fractError <= MAX_CALIB_ERROR)
			{
				if (refPDB.size() > 0)
				{
					double fractErrorNative = double(numErrorsNative)/(double)(totalSeedCounted);
					printf("Iter %d   fractNativeError: %5.3f    FractError %5.3f   Min Calibrated Dist: %5.2f\n",
							numIter,fractErrorNative,fractError, minCalibDist);
				}
				else
					printf("Iter %d   FractError %5.3f   Min Calibrated Dist: %5.2f\n",numIter, fractError, minCalibDist);

//				double fractErrorType[NUM_CALIB_TYPES];
//				for (int i = 0; i < NUM_CALIB_TYPES; i++)
//				{
//					fractErrorType[i] = (seedCountedType[i] > 0 ? double(numErrorsType[i])/double(seedCountedType[i]) : 0);
//					printf("%s Error=%f\n",CALIBRATION_TYPE_NAMES[i].c_str(),fractErrorType[i]);
//				}
				for (int d = 0; d < 3; d++)
					printf("DistCounts %d is %d\n",d,distCounts[d]);

				break;
			}
			else
			{
				// adjust coefficients
				for (int i = 0; i < NUM_CALIB_TYPES; i++)
				{
					if (stepSize[t][i] > 0)
					{
						coeffs[t][i] = coeffs[t][i] - stepSize[t][i];
						if (coeffs[t][i] <= 0)
							coeffs[t][i] = coeffs[t][i] + stepSize[t][i]; // go back to previous value
					}
				}

				kCalibCoeffs[t][BB_CALIB] = kCalibCoeffs[t][BB_CALIB] - increment[t][BB_CALIB];
				if (kCalibCoeffs[t][BB_CALIB] <= 0)
					kCalibCoeffs[t][BB_CALIB] = kCalibCoeffs[t][BB_CALIB] + increment[t][BB_CALIB]; // go back to previous value
				kCalibCoeffs[t][BE_CALIB] = kCalibCoeffs[t][BE_CALIB] - increment[t][BE_CALIB];
				if (kCalibCoeffs[t][BE_CALIB] <= 0)
					kCalibCoeffs[t][BE_CALIB] = kCalibCoeffs[t][BE_CALIB] + increment[t][BE_CALIB];
				kCalibCoeffs[t][ME_CALIB] = kCalibCoeffs[t][ME_CALIB] - increment[t][ME_CALIB];
				if (kCalibCoeffs[t][ME_CALIB] <= 0)
					kCalibCoeffs[t][ME_CALIB] = kCalibCoeffs[t][ME_CALIB] + increment[t][ME_CALIB];
				kCalibCoeffs[t][SC_CALIB] = kCalibCoeffs[t][SC_CALIB] - increment[t][SC_CALIB];
				if (kCalibCoeffs[t][SC_CALIB] <= 0)
					kCalibCoeffs[t][SC_CALIB] = kCalibCoeffs[t][SC_CALIB] + increment[t][SC_CALIB];

		        // fix some potential calibration problems
				if (kCalibCoeffs[t][BB_CALIB] > 0 && kCalibCoeffs[t][ME_CALIB] > 0.33333*kCalibCoeffs[t][BB_CALIB])
					kCalibCoeffs[t][ME_CALIB] = 0.33333*kCalibCoeffs[t][BB_CALIB];
				if (kCalibCoeffs[t][BE_CALIB] > 0 && kCalibCoeffs[t][SC_CALIB] > 0.66667*kCalibCoeffs[t][BE_CALIB])
					kCalibCoeffs[t][SC_CALIB] = 0.66667*kCalibCoeffs[t][BE_CALIB];

				for (int i = 0; i < NUM_CALIB_TYPES; i++)
				{
					if (allCoeffs[t][i].size() >= MIN_CALIB_COUNT)
						continue; // already set
					switch (i)
					{
						case CALIB_BBBB:
							if (kCalibCoeffs[t][BB_CALIB] > 0)
								coeffs[t][CALIB_BBBB] = kCalibCoeffs[t][BB_CALIB];
							break;
						case CALIB_BBBE:
							if (kCalibCoeffs[t][BB_CALIB] > 0 && kCalibCoeffs[t][BE_CALIB] > 0)
								coeffs[t][CALIB_BBBE] = sqrt(kCalibCoeffs[t][BB_CALIB]*kCalibCoeffs[t][BE_CALIB]);
							break;
						case CALIB_BEBE:
							if (kCalibCoeffs[t][BE_CALIB] > 0)
								coeffs[t][CALIB_BEBE] = kCalibCoeffs[t][BE_CALIB];
							break;
						case CALIB_BBME:
							if (kCalibCoeffs[t][BB_CALIB] > 0 && kCalibCoeffs[t][ME_CALIB] > 0)
								coeffs[t][CALIB_BBME] = sqrt(kCalibCoeffs[t][BB_CALIB]*kCalibCoeffs[t][ME_CALIB]);
							break;
						case CALIB_BBSC:
							if (kCalibCoeffs[t][BB_CALIB] > 0 && kCalibCoeffs[t][SC_CALIB] > 0)
								coeffs[t][CALIB_BBSC] = sqrt(kCalibCoeffs[t][BB_CALIB]*kCalibCoeffs[t][SC_CALIB]);
							break;
						case CALIB_BEME:
							if (kCalibCoeffs[t][BE_CALIB] > 0 && kCalibCoeffs[t][ME_CALIB] > 0)
								coeffs[t][CALIB_BEME] = sqrt(kCalibCoeffs[t][BE_CALIB]*kCalibCoeffs[t][ME_CALIB]);
							break;
						case CALIB_BESC:
							if (kCalibCoeffs[t][BE_CALIB] > 0 && kCalibCoeffs[t][SC_CALIB] > 0)
								coeffs[t][CALIB_BESC] = sqrt(kCalibCoeffs[t][BE_CALIB]*kCalibCoeffs[t][SC_CALIB]);
							break;
						case CALIB_MEME:
							if (kCalibCoeffs[t][ME_CALIB] > 0)
								coeffs[t][CALIB_MEME] = kCalibCoeffs[t][ME_CALIB];
							break;
						case CALIB_MESC:
							if (kCalibCoeffs[t][ME_CALIB] > 0 && kCalibCoeffs[t][SC_CALIB] > 0)
								coeffs[t][CALIB_MESC] = sqrt(kCalibCoeffs[t][ME_CALIB]*kCalibCoeffs[t][SC_CALIB]);
							break;
						case CALIB_SCSC:
							if (kCalibCoeffs[t][SC_CALIB] > 0)
								coeffs[t][CALIB_SCSC] = kCalibCoeffs[t][SC_CALIB];
							break;
						default:
							break;
					}
				} // end for each calib contact type
			} // end if fract error > cut
			numIter++;
		} // end while numIter
		if (numIter == MAXITER)
		{
			NOE* n1 = peakLists[t].front();
			NOE* n2 = peakLists[t].back();
			if (n1->volume <= 0 && n2->volume <= 0)
				printf("WARNING: Max iter %d reached for distance calibration. FractError %5.3f. Peak list has no intensities?\n",MAXITER,fractError);
			else
				printf("WARNING: Max iter %d reached for distance calibration. FractError %5.3f\n",MAXITER,fractError);
		}
	} // end for each NOE peak list type
	// for contact types missing calibration coefficients; use average value of the other coefficients
	// print the calibration coefficients
	for (int t = 0; t < NOE::NUMTYPES; t++)
	{
		if (peakLists[t].size() == 0)
			continue;
		double minval = 999999;
		for (int i = 0; i < NUM_CALIB_TYPES; i++)
		{
			if (coeffs[t][i] > 0 && coeffs[t][i] < minval)
			{
				minval = coeffs[t][i];
			}
		}
		for (int i = 0; i < NUM_CALIB_TYPES; i++)
		{
			if (coeffs[t][i] == 0)
				coeffs[t][i] = minval;
			// printf("Calibration constant for peak list %s, contact type %s is %e\n",
			//	NOE::PEAKTYPENAMES[t].c_str(), CALIBRATION_TYPE_NAMES[i].c_str(), coeffs[t][i]);
		}
	}
}

/**
 * Calibrates X, HX (adds offsets to NOE peak list chemical shifts)
 * Used by filterCalibratePeakLists()
 * TODO: 4D, 2D; this only supports 3D peak lists
 */
void calibrateNOEPeakList(list<NOE*>& noes, NOE::PEAKLISTTYPE type)
{
	// check if valid chemical shift
	// must have both C,HC
	// count methyls only once
	vector< tr1::array<double,2> > ref; // stores C,HC or N,HN chem shifts depending on peak list type
	vector<double> refH; // stores H chemical shifts; peak list type dependent
	vector< tr1::array<double,2> > input; // list of x,hx from NOEs
	vector<double> inputH; // list of H from NOEs

	char typeChar = 'C';
	if (type == NOE::N15_3D)
		typeChar = 'N'; // else CALI_3D or CARO_3D
	for (int r1 = 1; r1 <= bmrbProtein->size; r1++)
	{
		Residue* res1 = (*bmrbProtein)[r1];
		for (AtomIterator itX1 = res1->begin(typeChar); itX1 != res1->end(typeChar); itX1++)
		{
			Atom* x1Atom = *itX1;
			if (x1Atom->getNumProtons() < 1)
				continue; // no proton
			if (!x1Atom->hasChemShift())
				continue;
			if (type == NOE::CALI_3D)
			{
				if (res1->isCaro(x1Atom->name))
					continue;
			}
			else if (type == NOE::CARO_3D)
			{
				if (!res1->isCaro(x1Atom->name))
					continue;
			} // else not applicable for N15_3D
			bool isMethylX1 = x1Atom->isMethyl(); // NH3 counts as methyl in our definition
			Atom* hx1Atom = NULL;
			for (HIterator itHX1 = x1Atom->beginH(); itHX1 != x1Atom->endH(); itHX1++)
			{
				if (hx1Atom != NULL && isMethylX1)
					break; // count the protons of a methyl group as one pseudoproton
				hx1Atom = *itHX1;
				if (!hx1Atom->hasChemShift())
					continue;
				if (type == NOE::CALI_3D)
				{
					if (hx1Atom->cs > 4.0)
						continue;
				}
				else if (type == NOE::CARO_3D)
				{
					if (hx1Atom->cs < 6.0)
						continue;
				}
				else
				{
					if (hx1Atom->cs < 7.0)
						continue;
				} // TODO 4D
				tr1::array<double,2> ch;
				ch[0] = x1Atom->cs;
				ch[1] = hx1Atom->cs;
				ref.push_back(ch);
				refH.push_back(hx1Atom->cs);
			}
		}
	}
// BEGIN CLUSTER VERSION
	// cluster the NOEs into overlapping 2D groups
	list<NOECluster*> clusters;
	double tolx = tols[type][0];
	double tolhx = tols[type][1];
	for (list<NOE*>::iterator it = noes.begin(); it != noes.end(); ++it)
	{
		NOE* np = *it;

		if (type == NOE::CALI_3D)
		{
			if (np->hx > 4.0)
				continue;
		}
		else if (type == NOE::CARO_3D)
		{
			if (np->hx < 6.0)
				continue;
		}
		else
		{
			if (np->hx < 7.0)
				continue;
		} // TODO 4D

		bool hasSimilarCluster = false; // true if another cluster exists that matches this peak
		for (list<NOECluster*>::iterator itC = clusters.begin(); itC != clusters.end(); ++itC)
		{
			NOECluster* nc = *itC;
			if (nc->testAndAdd(np,tolx,tolhx))
				hasSimilarCluster = true;
		}
		if (!hasSimilarCluster)
			clusters.push_back(new NOECluster(np));
	}
	// delete small clusters
	for (list<NOECluster*>::iterator it = clusters.begin(); it != clusters.end();)
	{
		NOECluster* nc = *it;
		if (nc->size() < 3)
		{
			delete nc;
			it = clusters.erase(it);
		}
		else
			++it;
	}
	// do splitting; overlaps might disappear when do the splitting
	for (list<NOECluster*>::iterator it = clusters.begin(); it != clusters.end();)
	{
		it = NOECluster::split(it,clusters,tolx/2.0,tolhx/2.0); // split adjusts the iterator it;
		                                               // split may delete *it (the NOECluster object)
	}
	// remove overlapped peaks
	for (list<NOECluster*>::iterator it1 = clusters.begin(); it1 != clusters.end(); ++it1)
	{
		NOECluster* nc1 = *it1;
		for (list<NOECluster*>::iterator it2 = clusters.begin(); it2 != clusters.end(); ++it2)
		{
			NOECluster* nc2 = *it2;
			if (nc1 < nc2)
			{
				nc1->remOverlap(nc2);
			}
		}
	}

	// filter again by size; add to input those that are not filtered out
	for (list<NOECluster*>::iterator it = clusters.begin(); it != clusters.end();)
	{
		NOECluster* nc = *it;
		if (nc->size() < 3)
		{
			delete nc;
			it = clusters.erase(it);
		}
		else
		{   // take avg x,hx of cluster
			++it;
			tr1::array<double,2> ch;
			double x = 0;
			double hx = 0;
			nc->getAvg(x,hx);
			ch[0] = x;
			ch[1] = hx;
			input.push_back(ch);

			for (list<NOE*>::iterator it = nc->noes.begin(); it != nc->noes.end(); ++it)
			{
				NOE* np = *it;
				inputH.push_back(np->h);
			}
		}
	}

	// clean up; don't need clusters anymore
	for (list<NOECluster*>::iterator it = clusters.begin(); it != clusters.end();)
	{
		NOECluster* nc = *it;
		delete nc;
		it = clusters.erase(it);
	}
// END CLUSTER VERSION

    double searchN = 0.7; // 1.0
	double searchHN = 0.07; // 0.1
	double searchH = 0.07; // 0.1
	double matchTolN = tolx;
	double matchTolHN = tolhx;
	double matchTolH = tolh;
	double offsetN = 0;
	double offsetHN = 0;
	double offsetH = 0;
	// do not calibrate if not enough clusters

	Rectangle::calibrate(ref, input,matchTolN,matchTolHN,searchN,searchHN, offsetN, offsetHN, false);

	sort(refH.begin(), refH.end());
	sort(inputH.begin(), inputH.end());
	Interval::calibrate(refH, inputH, matchTolH, searchH, offsetH ,false); // refH, inputH must be in ascending order

	printf("Calibrate %s X,HX,H: %4.3f %5.4f %5.4f\n",NOE::PEAKTYPENAMES[type].c_str(), offsetN,offsetHN,offsetH);
	for (list<NOE*>::iterator it = noes.begin(); it != noes.end(); ++it)
	{
		NOE* np = *it;
		np->reference(offsetN,offsetHN);
		np->reference(offsetH);
	}
}

// adds offsets to the chemical shift values of the peak lists to make the lists align with each other
// removes waterband at 4.3-5.1ppm for CALI
// if all peaks have negative intensity, change them to positive
// remove duplicate peaks, which may be due to artifacts from automated peak picking
// remove self (diagonal) peaks (we only care about the cross peaks)
void filterCalibratePeakLists()
{
	vector< tr1::array<double,2> > diagonalXH; // diagonal X, HX chemical shifts
	for (vector<NOE::PEAKLISTTYPE>::iterator it = peakListTypes.begin(); it != peakListTypes.end(); ++it)
	{
		NOE::PEAKLISTTYPE peakListType = *it;
		list<NOE*> diagPeaks;
		tr1::array<double,4>& tol = tols[peakListType]; // 3D: x,hx,INVALIDSHIFTVAL,h, 4D: x1,hx1,x2,hx2, 2D: h1 h2
		double tolH = 0.03;
		if (peakListType == NOE::N15_3D || peakListType == NOE::CALI_3D || peakListType == NOE::CARO_3D)
			tolH = tol[3];
		else
			continue; // TODO: 4D, 2D peak lists not yet supported
		list<NOE*>& noes = peakLists[peakListType];
		// test sign of NOE volumes; if mostly negative sign, then invert sign
		// ignore value -1, which stands for no volume

		int numPositive = 0;
		int numNegative = 0;
		int numIter = 0; // check only the first 20 peaks
		for (list<NOE*>::iterator itNOE = noes.begin(); itNOE != noes.end(); ++itNOE)
		{
			NOE* np = *itNOE;
			if (np->volume != -1)
			{
				if (np->volume < 0)
					numNegative++;
				else if (np->volume > 0)
					numPositive++;
			}
			numIter++;
			if (numIter == 20)
				break;
		}
		if (numNegative > numPositive)
		{
			for (list<NOE*>::iterator itNOE = noes.begin(); itNOE != noes.end(); ++itNOE)
			{
				NOE* np = *itNOE;
				np->volume = -np->volume;
			}
		}

		// remove duplicate NOEs; TODO: do this for 4D and 2D peak lists too
		if (peakListType == NOE::N15_3D || peakListType == NOE::CALI_3D || peakListType == NOE::CARO_3D)
			NOE::removeDuplicates3D(noes, 0.1, 0.01, 0.01);

		// calibrate x, hx, h
		// TODO: 4D, 2D
		calibrateNOEPeakList(noes, peakListType);

		// remove self/diagonal NOEs; if CALI_3D remove waterband between 4.3-5.1
		// also remove aliased peaks with negative intensity (TODO: Do unaliasing)

		for (list<NOE*>::iterator it = noes.begin(); it != noes.end(); )
		{
			NOE* np = *it;
//			if (peakListType == NOE::CALI_3D)
//			{
//				if (np->hx > 4.3 && np->hx < 5.1) // waterband
//				{
//					delete np;
//					it = noes.erase(it);
//					continue;
//				}
//			}
			double dhx = np->hx-np->h; // TODO: make this work for 4D, 2D peak lists

			if (abs(dhx) <= tolH) // potential diagonal peak (mark for deletion, but don't delete)
			{
				tr1::array<double,2> xh = {np->x, np->hx};
				diagonalXH.push_back(xh);
				delete np;
				it = noes.erase(it);
			}
			else
			{
				if (np->volume != -1 && np->volume < 0)
				{
					delete np;
					it = noes.erase(it);
				}
				else
					it++;
			}
		} // end remove diagonal peaks from noes
	} // end for each peak list type
	// use diagonal peaks to filter out peaks if have enough diagonal peaks
	if (diagonalXH.size() > (unsigned)3*bmrbProtein->size)
	{
		int numPeaks = 0;
		int numFiltered = 0;
		for (vector<NOE::PEAKLISTTYPE>::iterator it = peakListTypes.begin(); it != peakListTypes.end(); ++it)
		{
			NOE::PEAKLISTTYPE peakListType = *it;
			tr1::array<double,4>& tol = tols[peakListType]; // 3D: x,hx,INVALIDSHIFTVAL,h, 4D: x1,hx1,x2,hx2, 2D: h1 h2
			double tolH = 0.03;
			if (peakListType == NOE::N15_3D || peakListType == NOE::CALI_3D || peakListType == NOE::CARO_3D)
				tolH = tol[3];
			else
				continue; // TODO: 4D, 2D peak lists not yet supported
			list<NOE*>& noes = peakLists[peakListType];
			for (list<NOE*>::iterator itNOE = noes.begin(); itNOE != noes.end(); ++itNOE)
			{
				NOE* np = *itNOE;
				double h = np->h;
				bool match = false;
				// check diagonal
				for (vector< tr1::array<double,2> >::iterator itd = diagonalXH.begin(); itd != diagonalXH.end(); ++itd)
				{
					tr1::array<double,2>& cs = *itd;
					if (abs(h-cs[1]) <= tolH)
					{
						match = true;
						break;
					}
				}
				if (!match)
				{
					// remove peak
					delete np;
					itNOE = noes.erase(itNOE);
					numFiltered++;
				}
				numPeaks++;
			}
		}
		printf("Filtering with diagonal peaks. NumFiltered %d   NumPeaks %d\n",numFiltered, numPeaks);
	}
	else
	{
		printf("Not enough diagonal peaks %zu. Skipping filtering with diagonal peaks\n",diagonalXH.size());
	}
}

/**
 * True if ca1->h1 sphere intersects ca2->h2 sphere, where radii of sphere is 1/2 distance between ca and h
 * Score is between 0 and 1
 */
bool intersectSphere(Atom* ca1, Atom* h1, Atom* ca2, Atom* h2, double& score)
{
	score = 0;
	if (!ca1->hasCoordinates() || !h1->hasCoordinates() || !ca2->hasCoordinates() || !h2->hasCoordinates())
		return false;
	double r1 = ca1->getDistance(h1); // /2.0; // /4.0;
	double r2 = ca2->getDistance(h2); // /2.0; // /4.0;
	double a = ca1->getDistance(ca2);
	if (r1 > a || r2 > a)
	{
		score = 1.0;
		return true; // one engulfs the other by >= half
	}
	double x = (a*a+r1*r1-r2*r2)/(2*a); // x >= 0
	if (x < r1)
	{
		score = -1+(r1+r2)/a; // max((r1-x)/r1,(x-a+r2)/r2); // for intersection, have [0--r2--x--r1--a]
		// r1+r2 must be >= a otherwise no intersection
		return true;
	}
	else
		return false;
}

// if refPDBFlag is true, the contact map for the reference structures is made instead of for the template structures
// does not take into acccount chemical shift completeness. See makeContactMapExpected instead
void makeContactMap(bool refPDBFlag)
{
	CSProtein& csa = *bmrbProtein; // chemical shift assignment
	int NUMSTRUCTURES = 0;
	if (!refPDBFlag)
	{
		NUMSTRUCTURES = structures.size();
	}
	else
		NUMSTRUCTURES = refPDB.size();

	for (int r1 = 1; r1 <= csa.size; r1++)
	{
		Residue* res1 = csa[r1];
		for (AtomIterator itX1 = res1->begin('X'); itX1 != res1->end('X'); itX1++)
		{
			Atom* x1Atom = *itX1;
			if (x1Atom->getNumProtons() < 1)
				continue; // no proton
			bool isMethylX1 = x1Atom->isMethyl();
			Atom* hx1Atom = NULL; // used to ensure we count methyl groups only once
			for (HIterator itHX1 = x1Atom->beginH(); itHX1 != x1Atom->endH(); itHX1++)
			{
				if (hx1Atom != NULL && isMethylX1)
					break; // count the protons of a methyl group as one pseudoproton
				hx1Atom = *itHX1;
				for (int r2 = r1; r2 <= csa.size; r2++)
				{
					Residue* res2 = csa[r2];
					// NOEs can be in same residue
					// iterator over all heavy atoms, then get their H's
					for (AtomIterator itX2 = res2->begin('X'); itX2 != res2->end('X'); itX2++)
					{
						Atom* x2Atom = *itX2;
						if (x2Atom->getNumProtons() < 1)
							continue;
						bool isMethylX2 = x2Atom->isMethyl();
						Atom* hx2Atom = NULL;
						for (HIterator itHX2 = x2Atom->beginH(); itHX2 != x2Atom->endH(); itHX2++)
						{
							if (hx2Atom != NULL && isMethylX2)
								break;
							hx2Atom = *itHX2;
							if (r1 == r2 && hx1Atom >= hx2Atom)
								continue; // avoid double counting and self contacts

							Contact contact(res1,x1Atom,hx1Atom,res2,x2Atom,hx2Atom);
							int contactCount = 0;
							double contactCountSphere = 0;
							double minDist = INVALIDDISTANCE; // if this dist does not change, then coordinates do not exist
							double* distances = new double[NUMSTRUCTURES];
							int distCount = 0; // number of entries in distances
							if (!refPDBFlag)
							{
								for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
								{
									CSProtein* structure = *itS;
									double dist = 0;
									if (structure->inContact(contact,DISTCUTOFF,dist)) // inContact considers methyl's by taking min dist; -1 if no coords
									{
										contactCount++;
										contactCountSphere += 1.0;
									}
									else if (dist > 0 && dist < 15.0)
									{
										// try intersect sphere
										Residue* r1 = (*structure)[res1->num];
										Residue* r2 = (*structure)[res2->num];
										if (r1 != NULL && r2 != NULL)
										{
											Atom* ca1 = r1->getC("CA");
											Atom* ca2 = r2->getC("CA");
											Atom* x1 = (x1Atom->name[0] != 'N' ? r1->getC(x1Atom->name) : r1->getN(x1Atom->name));
											Atom* hx1 = x1->getH(hx1Atom->name);
											Atom* x2 = (x2Atom->name[0] != 'N' ? r2->getC(x2Atom->name) : r2->getN(x2Atom->name));
											Atom* hx2 = x2->getH(hx2Atom->name);
											double score = 0;
											if (intersectSphere(ca1, hx1, ca2, hx2,score)) // if coords do not exist, false is returned
											{
												contactCountSphere += 0.25*score;
											}
										}
									}
									if (dist > 0 && dist < minDist)
										minDist = dist;
									if (dist >= 0)
									{
										distances[distCount] = dist;
										distCount++;
									}
								} // end for each structure

								double fractStruct = double(contactCount)/(double)(structures.size());
								double fractStructSphere = contactCountSphere/(double)(structures.size());
								double avgDist = INVALIDDISTANCE;
								double stdevDist = INVALIDDISTANCE;
								if (distCount > 0)
									avg_stdev(distances,distCount,avgDist,stdevDist);

								tr1::array<double,ContactMap::CMENTRYSIZE> entry = ContactMap::makeEntry(fractStruct,minDist,fractStructSphere,avgDist,stdevDist);
								contactMap.add(contact, entry);
								Contact contactRev(res2,x2Atom,hx2Atom,res1,x1Atom,hx1Atom);
								contactMap.add(contactRev, entry);
							}
							else // refPDB
							{
								for (list<CSProtein*>::const_iterator itS = refPDB.begin(); itS != refPDB.end(); itS++)
								{
									CSProtein* structure = *itS;
									double dist = 0;
									if (structure->inContact(contact,DISTCUTOFF,dist)) // inContact considers methyl's by taking min dist; -1 if no coords
									{
										contactCount++;
										contactCountSphere += 1.0;
									}
									else if (dist > 0 && dist < 15.0)
									{
										// try intersect sphere
										Residue* r1 = (*structure)[res1->num];
										Residue* r2 = (*structure)[res2->num];
										if (r1 != NULL && r2 != NULL)
										{
											Atom* ca1 = r1->getC("CA");
											Atom* ca2 = r2->getC("CA");
											Atom* x1 = (x1Atom->name[0] != 'N' ? r1->getC(x1Atom->name) : r1->getN(x1Atom->name));
											Atom* hx1 = x1->getH(hx1Atom->name);
											Atom* x2 = (x2Atom->name[0] != 'N' ? r2->getC(x2Atom->name) : r2->getN(x2Atom->name));
											Atom* hx2 = x2->getH(hx2Atom->name);
											double score = 0;
											if (intersectSphere(ca1, hx1, ca2, hx2,score)) // if coords do not exist, false is returned
											{
												contactCountSphere += 0.25*score;
											}
										}
									}
									if (dist > 0 && dist < minDist)
										minDist = dist;
									if (dist >= 0)
									{
										distances[distCount] = dist;
										distCount++;
									}
								} // end for each structure

								double fractStruct = double(contactCount)/(double)(refPDB.size());
								double fractStructSphere = contactCountSphere/(double)(refPDB.size());
								double avgDist = INVALIDDISTANCE;
								double stdevDist = INVALIDDISTANCE;
								if (distCount > 0)
									avg_stdev(distances,distCount,avgDist,stdevDist);
								tr1::array<double,ContactMap::CMENTRYSIZE> entry = ContactMap::makeEntry(fractStruct,minDist,fractStructSphere,avgDist,stdevDist);
								refContactMap.add(contact, entry);
								Contact contactRev(res2,x2Atom,hx2Atom,res1,x1Atom,hx1Atom);
								refContactMap.add(contactRev, entry);
							} // end else refPDB
							delete[] distances;
						} // end for each HX atom in residue j
					} // end for each X atom in residue j
				} // end for each residue j for H
			} // end for each HN
		} // end for each N atom in residue i
	} // end for each residue i for N,HN
}

// similar to makeContactMap, but takes into account chemical shift completeness
// and peak list types; if a contact is not possible due to missing chem shifts or
// NOESY experiment was not performed then the contact will not exist in the map
// TODO: consider 4D peak lists
void makeContactMapExpected()
{
	numLRContacts = 0;
	bool has15N = false; // true if peak list type exists
	bool hasCali = false;
	bool hasCaro = false;
	for (vector<NOE::PEAKLISTTYPE>::iterator it = peakListTypes.begin(); it != peakListTypes.end(); ++it)
	{
		NOE::PEAKLISTTYPE type = *it;
		if (type == NOE::N15_3D)
			has15N = true;
		if (type == NOE::CALI_3D)
			hasCali = true;
		if (type == NOE::CARO_3D)
			hasCaro = true;
	}
	CSProtein& csa = *bmrbProtein; // chemical shift assignment
	for (int r1 = 1; r1 <= csa.size; r1++)
	{
		Residue* res1 = csa[r1];
		for (AtomIterator itX1 = res1->begin('X'); itX1 != res1->end('X'); itX1++)
		{
			Atom* x1Atom = *itX1;
			if (x1Atom->getNumProtons() < 1)
				continue; // no proton
			bool isMethylX1 = x1Atom->isMethyl();
			Atom* hx1Atom = NULL; // used to ensure we count methyl groups only once
			for (HIterator itHX1 = x1Atom->beginH(); itHX1 != x1Atom->endH(); itHX1++)
			{
				if (hx1Atom != NULL && isMethylX1)
					break; // count the protons of a methyl group as one pseudoproton
				hx1Atom = *itHX1;
				for (int r2 = r1; r2 <= csa.size; r2++)
				{
					Residue* res2 = csa[r2];
					// NOEs can be in same residue
					// iterator over all heavy atoms, then get their H's
					for (AtomIterator itX2 = res2->begin('X'); itX2 != res2->end('X'); itX2++)
					{
						Atom* x2Atom = *itX2;
						if (x2Atom->getNumProtons() < 1)
							continue;
						bool isMethylX2 = x2Atom->isMethyl();
						Atom* hx2Atom = NULL;
						for (HIterator itHX2 = x2Atom->beginH(); itHX2 != x2Atom->endH(); itHX2++)
						{
							if (hx2Atom != NULL && isMethylX2)
								break;
							hx2Atom = *itHX2;
							if (r1 == r2 && hx1Atom >= hx2Atom)
								continue; // avoid double counting and self contacts

							if (!hx1Atom->hasChemShift() || !hx2Atom->hasChemShift())
								continue;

							Contact contact(res1,x1Atom,hx1Atom,res2,x2Atom,hx2Atom);
							int contactCount = 0;
							double contactCountSphere = 0;
							double minDist = INVALIDDISTANCE; // if this dist does not change, then coordinates do not exist
							bool skipAvg = false; // true if missing coordinates
							int i = 0;

							double distances[structures.size()];
							for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
							{
								CSProtein* structure = *itS;
								double dist = 0;
								if (structure->inContact(contact,DISTCUTOFF,dist)) // inContact considers methyl's by taking min dist; -1 if no coords
								{
									contactCount++;
									contactCountSphere += 1.0;
								}
								else if (dist > 0 && dist < 15.0)
								{
									// try intersect sphere
									Residue* r1 = (*structure)[res1->num];
									Residue* r2 = (*structure)[res2->num];
									if (r1 != NULL && r2 != NULL)
									{
										Atom* ca1 = r1->getC("CA");
										Atom* ca2 = r2->getC("CA");
										Atom* x1 = (x1Atom->name[0] != 'N' ? r1->getC(x1Atom->name) : r1->getN(x1Atom->name));
										Atom* hx1 = x1->getH(hx1Atom->name);
										Atom* x2 = (x2Atom->name[0] != 'N' ? r2->getC(x2Atom->name) : r2->getN(x2Atom->name));
										Atom* hx2 = x2->getH(hx2Atom->name);
										double score = 0;
										if (intersectSphere(ca1, hx1, ca2, hx2,score)) // if coords do not exist, false is returned
										{
											contactCountSphere += 0.25*score;
										}
									}
								}
								if (dist > 0 && dist < minDist)
									minDist = dist;
								if (dist >= 0)
									distances[i] = dist;
								else
									skipAvg = true;
								i++;
							} // end for each structure

							double fractStruct = double(contactCount)/(double)(structures.size());
							double fractStructSphere = contactCountSphere/(double)(structures.size());
							double avgDist = INVALIDDISTANCE;
							double stdevDist = INVALIDDISTANCE;
							if (!skipAvg)
								avg_stdev(distances,structures.size(),avgDist,stdevDist);
							tr1::array<double,ContactMap::CMENTRYSIZE> entry = ContactMap::makeEntry(fractStruct,minDist,fractStructSphere,avgDist,stdevDist);

							bool updateNumLRContacts = false;
							if ( x1Atom->hasChemShift() &&
									((res1->isCali(x1Atom->name) && hasCali) || (res1->is15N(x1Atom->name) && has15N) || (res1->isCaro(x1Atom->name) && hasCaro)) )
							{
								contactMapExpected.add(contact, entry);
								updateNumLRContacts = true;
							}
							if ( x2Atom->hasChemShift() &&
									((res2->isCali(x2Atom->name) && hasCali) || (res2->is15N(x2Atom->name) && has15N) || (res2->isCaro(x2Atom->name) && hasCaro)) )
							{
								Contact contactRev(res2,x2Atom,hx2Atom,res1,x1Atom,hx1Atom);
								contactMapExpected.add(contactRev, entry);
								updateNumLRContacts = true;
							}
							if (updateNumLRContacts && (r2-r1) >= LRSEQSEP && minDist <= DISTCUTOFF && fractStruct > 0)
							{
								numLRContacts++;
							}
						} // end for each HX atom in residue j
					} // end for each X atom in residue j
				} // end for each residue j for H
			} // end for each HN
		} // end for each N atom in residue i
	} // end for each residue i for N,HN
	printf("MakeContactMapExpected: numLRContacts(undirected)=%d\n",numLRContacts);
}

// this also finds all the assignment possibilities in addition to computing the cs, str, db scores
// assignments are store in asspos
void computeCS_Str_DB_Score(list<Assignment*>& asspos)
{
	CSProtein& csa = *bmrbProtein; // chemical shift assignment
	int numCYS = 0;
	for (int r = 1; r <= csa.size; r++)
	{
		Residue* res = csa[r];
		if (res->type == Residue::CYS)
		{
			Atom* cb = res->getC("CB");
			if (cb != NULL && cb->hasChemShift() && cb->cs > 34.0)
				numCYS++;
		}
	}
	double cysPairProb = 0;
	if (numCYS > 1)
	{
		int numCYSPairs = numCYS*(numCYS-1)/2;
		cysPairProb = 1.0/(double)(numCYSPairs);
	}
	for (vector<NOE::PEAKLISTTYPE>::iterator it = peakListTypes.begin(); it != peakListTypes.end(); ++it)
	{
		NOE::PEAKLISTTYPE peakListType = *it;
		list<NOE*>& noes = peakLists[peakListType];
		tr1::array<double,4>& tol = tols[peakListType];
		for (list<NOE*>::const_iterator it = noes.begin(); it != noes.end(); it++)
		{
			NOE* np = *it;
			for (int r1 = 1; r1 <= csa.size; r1++)
			{
				Residue* res1 = csa[r1];
				bool isDyn1 = isDynamic.test(res1->num-1);
				for (AtomIterator itX1 = res1->begin(NOE::HEAVYATOM1TYPE[peakListType]); itX1 != res1->end(NOE::HEAVYATOM1TYPE[peakListType]); itX1++)
				{
					Atom* x1Atom = *itX1;
					if (x1Atom->getNumProtons() < 1)
						continue; // no proton
					if (!x1Atom->hasChemShift())
						continue;
					if (!np->match(x1Atom->cs,INVALIDSHIFTVAL,INVALIDSHIFTVAL,INVALIDSHIFTVAL,
							tol[0],INVALIDSHIFTVAL,INVALIDSHIFTVAL,INVALIDSHIFTVAL))
						continue;
					bool isMethylX1 = x1Atom->isMethyl();
					Atom* hx1Atom = NULL;
					for (HIterator itHX1 = x1Atom->beginH(); itHX1 != x1Atom->endH(); itHX1++)
					{
						if (hx1Atom != NULL && isMethylX1)
							break; // count the protons of a methyl group as one pseudoproton; e.g. HD21 will represent HD21,HD22,HD23
						hx1Atom = *itHX1;
						if (!hx1Atom->hasChemShift())
							continue;
						if (!np->match(INVALIDSHIFTVAL,hx1Atom->cs,INVALIDSHIFTVAL,INVALIDSHIFTVAL,
								INVALIDSHIFTVAL,tol[1],INVALIDSHIFTVAL,INVALIDSHIFTVAL))
							continue;
						for (int r2 = 1; r2 <= csa.size; r2++)
						{
							Residue* res2 = csa[r2];
							bool isDyn2 = isDynamic.test(res2->num-1);
							int seqSep = abs(res1->num-res2->num);
							if ((isDyn1 && isDyn2) && seqSep > 1)
								continue; // only contact with seq sep <= 1 for dynamic parts

							// NOEs can be in same residue
							// iterator over all heavy atoms, then get their H's
							for (AtomIterator itX2 = res2->begin(NOE::HEAVYATOM2TYPE[peakListType]); itX2 != res2->end(NOE::HEAVYATOM2TYPE[peakListType]); itX2++)
							{
								Atom* x2Atom = *itX2;
								if (np->getNumDimensions() > 3 && !x2Atom->hasChemShift())
									continue; // 4D needs to have chem shift
								if (x2Atom->getNumProtons() < 1)
									continue;
								if (!np->match(INVALIDSHIFTVAL,INVALIDSHIFTVAL,x2Atom->cs,INVALIDSHIFTVAL,
										INVALIDSHIFTVAL,INVALIDSHIFTVAL,tol[2],INVALIDSHIFTVAL)) // always returns true if 3D or 2D NOE
									continue;
								bool isMethylX2 = x2Atom->isMethyl();
								Atom* hx2Atom = NULL;
								for (HIterator itHX2 = x2Atom->beginH(); itHX2 != x2Atom->endH(); itHX2++)
								{
									if (hx2Atom != NULL && isMethylX2)
										break;
									hx2Atom = *itHX2;
									if (!hx2Atom->hasChemShift())
										continue;
									if (hx1Atom == hx2Atom)
										continue; // can't be same atom
									if (!np->match(INVALIDSHIFTVAL,INVALIDSHIFTVAL,INVALIDSHIFTVAL,hx2Atom->cs,
											INVALIDSHIFTVAL,INVALIDSHIFTVAL,INVALIDSHIFTVAL,tol[3]))
										continue;
									// NOE chemical shift match score
									double chemShiftProb = 0;
									double dx1,dhx1,dx2,dhx2;
									np->getShiftDiff(x1Atom->cs,hx1Atom->cs,x2Atom->cs,hx2Atom->cs,dx1,dhx1,dx2,dhx2);
									if (dx1 != INVALIDSHIFTVAL)
									{
										chemShiftProb = erfcc(factorX1*dx1);
									}
									if (dhx1 != INVALIDSHIFTVAL)
									{
										if (chemShiftProb > 0)
											chemShiftProb *= erfcc(factorH1*dhx1);
										else
											chemShiftProb = erfcc(factorH1*dhx1);
									}
									if (dx2 != INVALIDSHIFTVAL)
									{
										if (chemShiftProb > 0)
											chemShiftProb *= erfcc(factorX2*dx2);
										else
											chemShiftProb = erfcc(factorX2*dx2);
									}
									if (dhx2 != INVALIDSHIFTVAL)
									{
										if (chemShiftProb > 0)
											chemShiftProb *= erfcc(factorH2*dhx2);
										else
											chemShiftProb = erfcc(factorH2*dhx2);
									}
									// structure score
									Contact con(res1,x1Atom,hx1Atom,res2,x2Atom,hx2Atom);
									double structureProb = contactMap.get(con,ContactMap::FRACSTRUCSPHINDEX);

									if (structureProb < INIT_FRACTSTRCUT)
									{
										continue;
									}
									// TODO: scale by min dist using dist-6 power
									if (structureProb > 0)
									{
										double minDist = contactMap.get(con,ContactMap::MINDISTINDEX);
										if (minDist > 5.0)
											structureProb *= 0.8;
										else if (minDist > 4.0)
											structureProb *= 0.95;
										else if (minDist > 3.0)
											structureProb *= 0.99;
									}

									// bias score
									double pdbProb = 0;
									SSTYPE ss1 = ssTypes[res1->num-1]; // can be UNKNOWNSS
									SSTYPE ss2 = ssTypes[res2->num-1];
									// TODO: Collect and use stats. For now, use best guess for scoring
									// intra-residue > sequential > alpha-alpha (i,i+3,i+4) >
									// disulfide CS > salt-bridge > hydrophobic-hydrophobic
									if (isDyn1 || isDyn2)
									{
										// only intra and seq
										if (seqSep == 0)
											pdbProb = 1.0;
										else if (seqSep == 1)
											pdbProb = 0.75;
									}
									else
									{
										if (seqSep == 0)
											pdbProb = 1.0;
										else if (seqSep == 1)
											pdbProb = 0.75;
										else if (ss1 == HELIX && ss2 == HELIX && (seqSep == 3 || seqSep == 4))
											pdbProb = 0.5;
										else if ( (res1->isPositive() && res2->isNegative()) ||
												  (res1->isNegative() && res2->isPositive()) )
										{
											pdbProb = 0.1;
										}

										if (res1->isCYS() && res2->isCYS())
										{
											Atom* cb1 = res1->getC("CB");
											Atom* cb2 = res2->getC("CB");
											if (cb1 != NULL && cb1->hasChemShift() && cb2 != NULL && cb2->hasChemShift() &&
													cb1->cs > 34.0 && cb2->cs > 34.0)
											{
												pdbProb = max<double>(cysPairProb, pdbProb);
											}
										}
									}
									Score score(chemShiftProb,structureProb,0,0,0,0,0,0,pdbProb);
									Contact contact(res1,x1Atom,hx1Atom,res2,x2Atom,hx2Atom);
									Assignment* as = new Assignment(np,contact,score);
									asspos.push_back(as);
								} // end for each HX atom in residue j
							} // end for each X atom in residue j
						} // end for each residue j for H
					} // end for each HN
				} // end for each N atom in residue i
			} // end for each residue i for N,HN
		} // end for each NOE
	} // end for each peak list type
}

void computeAmbig_Score(list<Assignment*>& asspos)
{
	tr1::unordered_map<NOE*, list<Assignment*>  >* noe2Assigns = new tr1::unordered_map<NOE*, list<Assignment*>  >;
	// setup noe2Assigns
	for (list<Assignment*>::iterator it = asspos.begin(); it != asspos.end(); ++it)
	{
		Assignment* a = *it;
		NOE* n = a->noe;
		tr1::unordered_map<NOE*, list<Assignment*>  >::iterator itC = noe2Assigns->find(n);
		if (itC != noe2Assigns->end())
		{
			list<Assignment*>& list = (*noe2Assigns)[n];
			list.push_back(a);
		}
		else
		{
			list<Assignment*> list;
			list.push_back(a);
			(*noe2Assigns)[n] = list;
		}
	}
	// compute ambig score
	for (tr1::unordered_map<NOE*, list<Assignment*>  >::iterator it = noe2Assigns->begin(); it != noe2Assigns->end(); ++it)
	{
		list<Assignment*>& as = it->second;
		double ambigScore = 1.0/(double)(as.size());
		for (list<Assignment*>::iterator it2 = as.begin(); it2 != as.end(); ++it2)
		{
			Assignment* a = *it2;
			a->score.ambig = ambigScore;
		}
	}
	delete noe2Assigns;
}

void reComputeAmbig_score(list<PCAssignment*>& asspos)
{
	tr1::unordered_map<NOE*, list<PCAssignment*>  >* noe2Assigns = new tr1::unordered_map<NOE*, list<PCAssignment*>  >;
	// setup noe2Assigns
	for (list<PCAssignment*>::iterator it = asspos.begin(); it != asspos.end(); ++it)
	{
		PCAssignment* a = *it;
		NOE* n = a->noe;
		tr1::unordered_map<NOE*, list<PCAssignment*>  >::iterator itC = noe2Assigns->find(n);
		if (itC != noe2Assigns->end())
		{
			list<PCAssignment*>& list = (*noe2Assigns)[n];
			list.push_back(a);
		}
		else
		{
			list<PCAssignment*> list;
			list.push_back(a);
			(*noe2Assigns)[n] = list;
		}
	}
	// compute ambig score
	for (tr1::unordered_map<NOE*, list<PCAssignment*>  >::iterator it = noe2Assigns->begin(); it != noe2Assigns->end(); ++it)
	{
		list<PCAssignment*>& as = it->second;
		double ambigScore = 1.0/(double)(as.size());
		for (list<PCAssignment*>::iterator it2 = as.begin(); it2 != as.end(); ++it2)
		{
			PCAssignment* a = *it2;
			a->score.ambig = ambigScore;
		}
	}
	delete noe2Assigns;
}

// returns true if this directed contact (from hx1 to hx2) is expected to have
// a symmetric peak (from hx2 to hx1) given
// the peak list types and chemical shift assignments
bool hasExpectedSymmetricPeak(Residue* r1, Atom* x1, Atom* hx1, Residue* r2, Atom* x2, Atom* hx2,
		bool has15N, bool hasCali, bool hasCaro)
{
	if (r2->isCali(x2->name) && hasCali && x2->hasChemShift())
	{
		return true;
	}
	else if (r2->is15N(x2->name) && has15N && x2->hasChemShift())
	{
		return true;
	}
	else if (r2->isCaro(x2->name) && hasCaro && x2->hasChemShift())
	{
		return true;
	}
	return false;
}

/* used by computeSym_Inter_NetScore2 (for computing net and interres scores)
 * returns the scale factor for each pair of residues (not symmetric)
 * r1 = resnum of X,HX and r2 = resnum of H
 * assumes only 3d peak lists
 * */
double getNetScaleFactor(int r1, int r2)
{
	static double scale[MAXPROSIZE][MAXPROSIZE] = {0}; // declared static for lazy initialization
	if (scale[r1-1][r2-1] == 0)
	{
		CSProtein& csa = *bmrbProtein; // chemical shift assignment
		int sum = 0;
		int count = 0;
		for (vector<NOE::PEAKLISTTYPE>::iterator it = peakListTypes.begin(); it != peakListTypes.end(); ++it)
		{
			NOE::PEAKLISTTYPE peakListType = *it;
			tr1::array<double,4>& tol = tols[peakListType]; // [0] = X1, [1] = HX1, [2] = X2, [3] = HX2
			Residue* res1 = csa[r1];
			for (AtomIterator itX1 = res1->begin(NOE::HEAVYATOM1TYPE[peakListType]); itX1 != res1->end(NOE::HEAVYATOM1TYPE[peakListType]); itX1++)
			{
				Atom* x1Atom = *itX1;
				if (x1Atom->getNumProtons() < 1)
					continue; // no proton
				if (!x1Atom->hasChemShift())
					continue;
				bool isMethylX1 = x1Atom->isMethyl();
				Atom* hx1Atom = NULL;
				for (HIterator itHX1 = x1Atom->beginH(); itHX1 != x1Atom->endH(); itHX1++)
				{
					if (hx1Atom != NULL && isMethylX1)
						break; // count the protons of a methyl group as one pseudoproton; e.g. HD21 will represent HD21,HD22,HD23
					hx1Atom = *itHX1;
					if (!hx1Atom->hasChemShift())
						continue;

					int redundancya = 1;
					// get number of other x,hx atoms in r1 with similar chem shifts to hx1Atom
					for (AtomIterator itX1b = res1->begin(NOE::HEAVYATOM1TYPE[peakListType]); itX1b != res1->end(NOE::HEAVYATOM1TYPE[peakListType]); itX1b++)
					{
						Atom* x1Atomb = *itX1b;
						if (x1Atomb->getNumProtons() < 1)
							continue; // no proton
						if (!x1Atomb->hasChemShift())
							continue;
						bool isMethylX1b = x1Atomb->isMethyl();
						Atom* hx1Atomb = NULL;
						for (HIterator itHX1b = x1Atomb->beginH(); itHX1b != x1Atomb->endH(); itHX1b++)
						{
							if (hx1Atomb != NULL && isMethylX1b)
								break; // count the protons of a methyl group as one pseudoproton; e.g. HD21 will represent HD21,HD22,HD23
							hx1Atomb = *itHX1b;
							if (x1Atomb == x1Atom && hx1Atomb == hx1Atom)
								continue; // same atom pair
							if (!hx1Atomb->hasChemShift())
								continue;
							if (abs(x1Atomb->cs-x1Atom->cs) <= tol[0] && abs(hx1Atomb->cs-hx1Atom->cs) <= tol[1])
								redundancya++;
						} // end for each hx1b
					} // end for each x1b
					Residue* res2 = csa[r2];
					// iterator over all other heavy atoms, and get their H's
					for (AtomIterator itX2 = res2->begin(NOE::HEAVYATOM2TYPE[peakListType]); itX2 != res2->end(NOE::HEAVYATOM2TYPE[peakListType]); itX2++)
					{
						Atom* x2Atom = *itX2;
						if (!x2Atom->hasChemShift())
							continue; // 4D needs to have chem shift; always returns true for < 4D
						if (x2Atom->getNumProtons() < 1)
							continue;
						bool isMethylX2 = x2Atom->isMethyl();
						Atom* hx2Atom = NULL;
						for (HIterator itHX2 = x2Atom->beginH(); itHX2 != x2Atom->endH(); itHX2++)
						{
							if (hx2Atom != NULL && isMethylX2)
								break;
							hx2Atom = *itHX2;
							if (!hx2Atom->hasChemShift())
								continue;
							if (hx1Atom == hx2Atom)
								continue; // can't be same atom

							int redundancyb = 1;
							// get number of other h atoms in r2 with similar chem shifts to hx2Atom
							for (AtomIterator itX2b = res2->begin(NOE::HEAVYATOM2TYPE[peakListType]); itX2b != res2->end(NOE::HEAVYATOM2TYPE[peakListType]); itX2b++)
							{
								Atom* x2Atomb = *itX2b;
								if (!x2Atomb->hasChemShift())
									continue; // 4D needs to have chem shift
								if (x2Atomb->getNumProtons() < 1)
									continue;
								bool isMethylX2b = x2Atomb->isMethyl();
								Atom* hx2Atomb = NULL;
								for (HIterator itHX2b = x2Atomb->beginH(); itHX2b != x2Atomb->endH(); itHX2b++)
								{
									if (hx2Atomb != NULL && isMethylX2b)
										break;
									hx2Atomb = *itHX2b;
									if (!hx2Atomb->hasChemShift())
										continue;
									if (hx2Atomb == hx2Atom)
										continue; // can't be same atom
									if (abs(hx2Atomb->cs-hx2Atom->cs) <= tol[3])
										redundancyb++;
								} // end for each hx2b
							} // end for each x2b
							sum += redundancya*redundancyb;
							count++;
						} // end for each HX atom in residue j
					} // end for each X atom in residue j
				} // end for each HN
			} // end for each N atom in residue i
		} // end for each peak list type
		if (count > 0)
			scale[r1-1][r2-1] = double(sum)/double(count);
	} // end if scale initialized
	return scale[r1-1][r2-1];
}

// symmetry and net scores of assignment possibilities based on assigned rather than assignment possibilities
// updates the total score as well
void reComputeSym_Inter_Net_Ambig_Score(list<PCAssignment*>& assignments, vector<PCAssignment*>& assigned)
{
	// TODO: set flags for 4D, 2D
	bool has15N = false; // true if peak list type exists
	bool hasCali = false;
	bool hasCaro = false;
	for (vector<NOE::PEAKLISTTYPE>::iterator it = peakListTypes.begin(); it != peakListTypes.end(); ++it)
	{
		NOE::PEAKLISTTYPE type = *it;
		if (type == NOE::N15_3D)
			has15N = true;
		if (type == NOE::CALI_3D)
			hasCali = true;
		if (type == NOE::CARO_3D)
			hasCaro = true;
	}

	CSProtein& csa = *bmrbProtein; // chemical shift assignment
    // global variable: int expectedAssignMatrix[MAXPROSIZE][MAXPROSIZE] = {{0}}; // expected # of assignments; depends on peak list types
	tr1::unordered_set<Contact>* contactAss = new tr1::unordered_set<Contact>; // based on assigned contacts (directed contact); for finding symmetric peak
	int assignedMatrixDirected[MAXPROSIZE][MAXPROSIZE] = {{0}}; // directed version of assignMatrix for net score
	// setup contactAss; update ambig score, assignedMatrixDirected
	for (vector<PCAssignment*>::iterator it = assigned.begin(); it != assigned.end(); ++it)
	{
		PCAssignment* pca = *it;
		pca->score.ambig = 1.0;
		tr1::unordered_set<Contact>& contacts = pca->pc.getContacts(); // add all contacts to contactAss

		int r1 = 0;
		int r2 = 0;
		pca->pc.getResPair(r1,r2);
		assignedMatrixDirected[r1-1][r2-1]++;

		for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
		{
			const Contact& c = *itC;
			tr1::unordered_set<Contact>::iterator itFind = contactAss->find(c);
			if (itFind == contactAss->end())
				contactAss->insert(c);
		}
	}
//	// scale assignedMatrixDirected
//	for (int r1 = 1; r1 <= csa.size; r1++)
//	{
//		for (int r2 = 1; r2 <= csa.size; r2++)
//		{
//			if (abs(r2-r1) > 1)
//			{
//				double scaleDir = getNetScaleFactor(r1,r2);
//				if (scaleDir < 1.0)
//					scaleDir = 1.0;
//				assignedMatrixDirected[r1-1][r2-1] *= 1.0/scaleDir;
//			}
//		}
//	}

	int structMatrix[MAXPROSIZE][MAXPROSIZE] = {{0}}; // # of contacts between resnumA-1 to resnumB-1 in structures; matrix is NOT symmetric
	                                                  // because matrix depends on contactMapExpected, which is not symmetric
	int assignedStructMatrix[MAXPROSIZE][MAXPROSIZE] = {{0}}; // # of contacts between resnumA-1 to resnumB-1 in structures that have
	                                                        // an assignment; matrix is not NOT symmetric
	// initialize structMatrix
	for (ContactMapIterator it = contactMapExpected.begin(); it != contactMapExpected.end(); ++it)
	{
		const Contact& c = it.first();
		const tr1::array<double,ContactMap::CMENTRYSIZE>& arr = it.second();
		if (arr[ContactMap::MINDISTINDEX] <= DISTCUTOFF)
		{
			int r1 = c.r1->num;
			int r2 = c.r2->num;
			structMatrix[r1-1][r2-1]++;
			tr1::unordered_set<Contact>::const_iterator cFind = contactAss->find(c);
			if (cFind != contactAss->end())
			{
				assignedStructMatrix[r1-1][r2-1]++;
			}
		}
	}

	// residue-level network support
	vector< vector<int> > neighbors; // neighbors[resnum-1] = list of neighboring residues to consider for net score
	for (int r = 1; r <= csa.size; r++)
	{
		SSTYPE ss = ssTypes[r-1];
		vector<int> nList;
		nList.reserve(6);
		nList.push_back(r);
		int numSSMismatch = 0;
		for (int i = r-1; i >= r-NETWORKWINDOW; i--)
		{
			if (i > 0)
			{
				if (ssTypes[i-1] != ss)
				{
					numSSMismatch++;
					if (numSSMismatch > 1)
						break;
				}
				nList.push_back(i);
			}
			else
				break;
		}
		numSSMismatch = 0;
		for (int i = r+1; i <= r+NETWORKWINDOW; i++)
		{
			if (i <= csa.size)
			{
				if (ssTypes[i-1] != ss)
				{
					numSSMismatch++;
					if (numSSMismatch > 1)
						break;
				}
				nList.push_back(i);
			}
			else
				break;
		}
		neighbors.push_back(nList);
	}
	double netCS[MAXPROSIZE][MAXPROSIZE] = {{0}}; // net score
	double netStr[MAXPROSIZE][MAXPROSIZE] = {{0}}; // net str score; initially stores the numerator term
	for (int r1 = 1; r1 <= csa.size; r1++)
	{
		// for overlapped residues ensure r11 <= r22
		vector<int>& n1 = neighbors[r1-1];
		int r1End = -1;
		if (n1.size() > 1)
			r1End = n1[1];
		for (int r2 = r1; r2 <= csa.size; r2++)
		{
			vector<int>& n2 = neighbors[r2-1];
			int r2Begin = 0;
			if (n2.size() > 1)
				r2Begin = n2.back();
			bool overlapFlag = r1End >= r2Begin; // true if n1 overlaps with n2
			int countAss = 0; // sum of num assignments between r11, r22
			int expectedCountAss = 0; // expected sum of num assignments between r11, r22
			int strNumer = 0; // sum of contacts between r11, r22 with assignment
			int strDenom = 0; // sum of contacts between r11, r22

			// parallel orientation
			for (vector<int>::iterator it1 = n1.begin(); it1 != n1.end(); ++it1)
			{
				int r11 = *it1;
				int diff1 = abs(r11-r1);
				for (vector<int>::iterator it2 = n2.begin(); it2 != n2.end(); ++it2)
				{
					int r22 = *it2;
					if (r11 == r1 && r22 == r2)
						continue; // ignore r1-r2
                    // ignore intra-residue to r1 or r2
					if (r11==r22 && (r11 == r1 || r11 == r2))
						continue;
					int diff2 = abs(r22-r2);
					if ( ((r11-r1 >= 0 && r22-r2 >= 0) || (r1-r11 >= 0 && r2-r22 >= 0)) && abs(diff1-diff2) <= 1 )
					{
						if (overlapFlag && r11 >= r2Begin && r11 <= r1End && r22 >= r2Begin && r22 <= r1End && r11 > r22)
							continue; // avoid double counting
						countAss +=  (r11 != r22 ? assignedMatrixDirected[r11-1][r22-1]+assignedMatrixDirected[r22-1][r11-1] : assignedMatrixDirected[r11-1][r22-1]);
						expectedCountAss += 2*expectedAssignAvg[r11-1][r22-1];

						strNumer += (r11 != r22 ? assignedStructMatrix[r11-1][r22-1] + assignedStructMatrix[r22-1][r11-1] : assignedStructMatrix[r11-1][r22-1]);
						strDenom += (r11 != r22 ? structMatrix[r11-1][r22-1] + structMatrix[r22-1][r11-1] : structMatrix[r11-1][r22-1]);
					}
				}
			}
			if (expectedCountAss > 0)
			{
				netCS[r1-1][r2-1] = double(countAss)/double(expectedCountAss);
				if (r1 != r2)
					netCS[r2-1][r1-1] = netCS[r1-1][r2-1];
			}
			else if (countAss > 0)
			{
				netCS[r1-1][r2-1] = 1.0;
				if (r1 != r2)
					netCS[r2-1][r1-1] = 1.0;
			}
			if (strDenom > 0)
			{
				netStr[r1-1][r2-1] = double(strNumer)/double(strDenom);
				if (r1 != r2)
					netStr[r2-1][r1-1] = netStr[r1-1][r2-1];
			}

			// antiparallel
			countAss = 0; // sum of num assignment possibilities between r11, r22
			expectedCountAss = 0; // expected sum of num assignment possibilities between r11, r22
			strNumer = 0;
			strDenom = 0;
			for (vector<int>::iterator it1 = n1.begin(); it1 != n1.end(); ++it1)
			{
				int r11 = *it1;
				int diff1 = abs(r11-r1);
				for (vector<int>::iterator it2 = n2.begin(); it2 != n2.end(); ++it2)
				{
					int r22 = *it2;
					if (r11 == r1 && r22 == r2)
						continue; // ignore r1-r2
					int diff2 = abs(r22-r2);
					if ( ((r11-r1 >= 0 && r2-r22 >= 0) || (r1-r11 >= 0 && r22-r2 >= 0)) && abs(diff1-diff2) <= 1 )
					{
						if (overlapFlag && r11 >= r2Begin && r11 <= r1End && r22 >= r2Begin && r22 <= r1End && r11 > r22)
							continue; // avoid double counting
						countAss += assignedMatrixDirected[r11-1][r22-1]+assignedMatrixDirected[r22-1][r11-1];
						expectedCountAss += 2*expectedAssignAvg[r11-1][r22-1];

						strNumer += (r11 != r22 ? assignedStructMatrix[r11-1][r22-1] + assignedStructMatrix[r22-1][r11-1] : assignedStructMatrix[r11-1][r22-1]);
						strDenom += (r11 != r22 ? structMatrix[r11-1][r22-1] + structMatrix[r22-1][r11-1] : structMatrix[r11-1][r22-1]);
					}
				}
			}
			if (expectedCountAss > 0)
			{
				netCS[r1-1][r2-1] = max(double(countAss)/double(expectedCountAss),netCS[r1-1][r2-1]);
				if (r1 != r2)
					netCS[r2-1][r1-1] = netCS[r1-1][r2-1];
			}
			else if (countAss > 0)
			{
				netCS[r1-1][r2-1] = 1.0;
				if (r1 != r2)
					netCS[r2-1][r1-1] = 1.0;
			}
			if (strDenom > 0)
			{
				netStr[r1-1][r2-1] = max(double(strNumer)/double(strDenom),netStr[r1-1][r2-1]);
				if (r1 != r2)
					netStr[r2-1][r1-1] = netStr[r1-1][r2-1];
			}
		} // end for r2
	} // for r1

	// compute sym, inter, net, netstr score
	for (list<PCAssignment*>::iterator it = assignments.begin(); it != assignments.end(); ++it)
	{
		PCAssignment* pca = *it;
		PseudoContact& pc = pca->pc;
		tr1::unordered_set<Contact>& contacts = pc.getContacts();
		bool isAssigned = false; // true if pca was assigned
		bool symIsAssigned = false; // true if pca has symmetric pca assigned
		int expectedAssignCorrection = 1; // 1 = include this pc and possible +1 if symmetric peak can exist based on peak list types
		double avgExpectedAssign = 0; // average over all contacts in pc
		int avgExpectedAssignCount = 0;
		for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
		{
			const Contact& c = *itC;
			tr1::unordered_set<Contact>::iterator itFind = contactAss->find(c);
			if (itFind != contactAss->end())
				isAssigned = true;

			Residue* res1 = c.r1;
			Residue* res2 = c.r2;
			Atom* x1 = c.x1;
			Atom* x2 = c.x2;
			Atom* hx1 = c.hx1;
			Atom* hx2 = c.hx2;
			ExpectedAssignKey k(res1->num,hx1->name,res2->num,hx2->name);
			tr1::unordered_map<ExpectedAssignKey,int>::iterator itE = expectedAssign.find(k);
			if (itE != expectedAssign.end())
			{
				avgExpectedAssign += 2*itE->second; // *2 to consider both directions
				avgExpectedAssignCount++;
			}

			Contact c2 = c.reverse();
			if (hasExpectedSymmetricPeak(res1,x1,hx1,res2,x2,hx2,has15N,hasCali,hasCaro))
				expectedAssignCorrection = 2;

			itFind = contactAss->find(c2);
			if (itFind != contactAss->end())
				symIsAssigned = true;
		} // end for each contact in pca
		int assignCorrection = 0; // includes pc if assigned and possibly a symmetric contact if also assigned
		if (isAssigned)
			assignCorrection++;
		if (symIsAssigned)
			assignCorrection++;

		if (symIsAssigned)
			pca->score.sym = 1.0;
		else
			pca->score.sym = 0;

		int r1 = 0;
		int r2 = 0;
		pc.getResPair(r1,r2);

		avgExpectedAssign = (avgExpectedAssignCount != 0 ? avgExpectedAssign/avgExpectedAssignCount : 0); // average over all contacts in this pc
		// interres score ignores as (if assigned) and sym if sym assigned
		int assigned12 = (r1 != r2 ? assignedMatrixDirected[r1-1][r2-1] + assignedMatrixDirected[r2-1][r1-1] : assignedMatrixDirected[r1-1][r2-1]);
		if (assigned12-assignCorrection > 0)
		{
			if (avgExpectedAssign-expectedAssignCorrection > 0)
			{
				double newVal = double(assigned12-assignCorrection)/double(avgExpectedAssign-expectedAssignCorrection);
				if (newVal > 1.0)
					newVal = 1.0;
				else if (newVal < 0)
					newVal = 0;
				pca->score.interres = OLDWT*pca->score.interres+NEWWT*newVal;
			}
			else
				pca->score.interres = OLDWT*pca->score.interres+NEWWT; // newVal = 1.0
		}
		else
			pca->score.interres = OLDWT*pca->score.interres; // newVal = 0

		double newNet = netCS[r1-1][r2-1];
		if (newNet > 1.0)
			newNet = 1.0;
		else if (newNet < 0)
		    newNet = 0;
		pca->score.net = OLDWT*pca->score.net + NEWWT*newNet;

		double newNetStr = netStr[r1-1][r2-1];
		if (newNetStr > 1.0)
			newNetStr = 1.0;
		else if (newNetStr < 0)
			newNetStr = 0;
		pca->score.netStr = OLDWT*pca->score.netStr + NEWWT*newNetStr;


		pca->score.setTotal();
	} // end for each possible assignment
	delete contactAss;
}

// computes the intensity compatibility score from the seed assignment
// initializes calibConstants to contain the calibration constants: dist^6 = constant*intensity
// calls pca->score.setTotal() to update the total score
void computeIntensityScore(list<NOE*>& noes, list<PCAssignment*>& assignments, vector<PCAssignment*>& assigned,
		double calibConstants[NOE::NUMTYPES][NUM_CALIB_TYPES])
{
	double MINUS_1_6 = -1.0/6.0;
	calibrateDistances(noes, assigned, calibConstants);
	for (list<PCAssignment*>::iterator it = assignments.begin(); it != assignments.end(); ++it)
	{
		PCAssignment* pca = *it;
		NOE* noe = pca->noe;
		if (noe->volume <= 0)
			continue;
		double expectedDistance = 0; // intensity-based upper bound distance
		PseudoContact& pc = pca->pc;
		tr1::unordered_set<Contact>& contacts = pc.getContacts();
		double minDist = 999999;
		for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
		{
			const Contact& c = *itC;
			double ki = calibConstants[noe->type][getCalibContactType(c)];
			if (ki > 0)
			{
				expectedDistance += ki;
				double dist = contactMap[c][ContactMap::MINDISTINDEX];
				if (dist < minDist)
					minDist = dist;
			} // else no calibration constant exists for this contact
		}
		if (minDist < 999998)
		{
			if (noe->volume > 0)
				expectedDistance = pow(expectedDistance*noe->volume/(double)(contacts.size()), MINUS_1_6);
			else
				expectedDistance = DISTCUTOFF;
			if (minDist <= expectedDistance)
				pca->score.intensity = 1.0;
			else
			{
				double diff = minDist-expectedDistance;
				if (diff < INTENSITY_DISTCUT)
					pca->score.intensity = (INTENSITY_DISTCUT-diff)/INTENSITY_DISTCUT;
				else
					pca->score.intensity = 0;
			}
			pca->score.setTotal();
		}
	} // end for each assignment possibility
}

void computeSym_Inter_Net_Score2(list<Assignment*>& asspos)
{
	tr1::unordered_set<Contact>* contactAss = new tr1::unordered_set<Contact>; // stores contacts with at least one assignment possibility
	                                                                           // [directed contact]; used for symmetry and netstr score
	                                        // the same contact can be assigned to multiple peaks, so use set to count the contact at most once
	double assignPossibMatrixDirected[MAXPROSIZE][MAXPROSIZE] = {{0}}; // directed version of assignPossibMatrix
	          // scaled by getNetScaleFactor
			  // we will also scale assignPossibMatrix
	for (list<Assignment*>::iterator it = asspos.begin(); it != asspos.end(); ++it)
	{
		Assignment* as = *it;
		Contact& c = as->contact;
		tr1::unordered_set<Contact>::iterator itC = contactAss->find(c);
		if (itC == contactAss->end())
			contactAss->insert(c);
	}
	for (tr1::unordered_set<Contact>::iterator it = contactAss->begin(); it != contactAss->end(); ++it)
	{
		const Contact& c = *it;
		int r1 = c.r1->num;
		int r2 = c.r2->num;
		assignPossibMatrixDirected[r1-1][r2-1] += 1.0;
	}
	// scale assignPossibMatrixDirected
	CSProtein& csa = *bmrbProtein; // chemical shift assignment
	for (int r1 = 1; r1 <= csa.size; r1++)
	{
		for (int r2 = 1; r2 <= csa.size; r2++)
		{
			if (abs(r2-r1) > 1)
			{
				double scaleDir = getNetScaleFactor(r1,r2);
				if (scaleDir < 1.0)
					scaleDir = 1.0;
				assignPossibMatrixDirected[r1-1][r2-1] *= 1.0/scaleDir;
			}
		}
	}

	int numNOEs = 0;
	for (vector<NOE::PEAKLISTTYPE>::iterator it = peakListTypes.begin(); it != peakListTypes.end(); ++it)
	{
		NOE::PEAKLISTTYPE peakListType = *it;
		list<NOE*>& noes = peakLists[peakListType];
		numNOEs += noes.size();
	}

	// TODO: set flags for 4D, 2D
	bool has15N = false; // true if peak list type exists
	bool hasCali = false;
	bool hasCaro = false;
	for (vector<NOE::PEAKLISTTYPE>::iterator it = peakListTypes.begin(); it != peakListTypes.end(); ++it)
	{
		NOE::PEAKLISTTYPE type = *it;
		if (type == NOE::N15_3D)
			has15N = true;
		if (type == NOE::CALI_3D)
			hasCali = true;
		if (type == NOE::CARO_3D)
			hasCaro = true;
	}
	// for inter-residue and network support
	// global variable: int expectedAssignMatrix[MAXPROSIZE][MAXPROSIZE] = {{0}}; // expected # of h-h assignments; depends on peak list types
	int structMatrix[MAXPROSIZE][MAXPROSIZE] = {{0}}; // # of contacts between resnumA-1 to resnumB-1 in structures; matrix is NOT symmetric
	                                                  // Is directed
	int assignStructMatrix[MAXPROSIZE][MAXPROSIZE] = {{0}}; // # of contacts between resnumA-1 to resnumB-1 in structures that have at least
	                                                        // one assignment possibility between the residues; matrix is not NOT symmetric
	                                                        // Is directed
	// initialize structMatrix
	for (ContactMapIterator it = contactMapExpected.begin(); it != contactMapExpected.end(); ++it) // contactMapExpected can contain both contact directions
	{
		const Contact& c = it.first();
		const tr1::array<double,ContactMap::CMENTRYSIZE>& arr = it.second();
		if (arr[ContactMap::MINDISTINDEX] <= DISTCUTOFF)
		{
			int r1 = c.r1->num;
			int r2 = c.r2->num;
			structMatrix[r1-1][r2-1]++;
			// check if in contactAss
			tr1::unordered_set<Contact>::const_iterator cFind = contactAss->find(c);
			if (cFind != contactAss->end())
			{
				assignStructMatrix[r1-1][r2-1]++;
			}
		}
	}

	// network score around the residue neighbors of the residues in contact
	// window len depends on NETWORKWINDOW and whether the predicted secondary structures of the other residues
	// are the same (we allow at most one discepancy for the residues around each contact residue)
	vector< vector<int> > neighbors; // neighbors[resnum-1] = list of neighboring residues to consider for net score
	for (int r = 1; r <= csa.size; r++)
	{
		SSTYPE ss = ssTypes[r-1];
		vector<int> nList;
		nList.reserve(6);
		nList.push_back(r);
		int numSSMismatch = 0;
		for (int i = r-1; i >= r-NETWORKWINDOW; i--)
		{
			if (i > 0)
			{
				if (ssTypes[i-1] != ss)
				{
					numSSMismatch++;
					if (numSSMismatch > 1)
						break;
				}
				nList.push_back(i);
			}
			else
				break;
		}
		numSSMismatch = 0;
		for (int i = r+1; i <= r+NETWORKWINDOW; i++)
		{
			if (i <= csa.size)
			{
				if (ssTypes[i-1] != ss)
				{
					numSSMismatch++;
					if (numSSMismatch > 1)
						break;
				}
				nList.push_back(i);
			}
			else
				break;
		}
		neighbors.push_back(nList);
	}
	// compute net scores for all i,j pairs
	double netCS[MAXPROSIZE][MAXPROSIZE] = {{0}}; // net score
	double netStr[MAXPROSIZE][MAXPROSIZE] = {{0}}; // net str score; initially stores the numerator term
	for (int r1 = 1; r1 <= csa.size; r1++)
	{
		// for overlapped residues ensure r11 <= r22
		vector<int>& n1 = neighbors[r1-1];
		int r1End = -1;
		if (n1.size() > 1)
			r1End = n1[1];
		for (int r2 = r1; r2 <= csa.size; r2++)
		{
			vector<int>& n2 = neighbors[r2-1];
			int r2Begin = 0;
			if (n2.size() > 1)
				r2Begin = n2.back();
			bool overlapFlag = r1End >= r2Begin; // true if n1 overlaps with n2
			int countAss = 0; // sum of num assignment possibilities between r11, r22
			int expectedCountAss = 0; // expected sum of num assignment possibilities between r11, r22
			int strNumer = 0; // sum of contacts between r11, r22 with assignment possibility
			int strDenom = 0; // sum of contacts between r11, r22

			// parallel orientation
			for (vector<int>::iterator it1 = n1.begin(); it1 != n1.end(); ++it1)
			{
				int r11 = *it1;
				int diff1 = abs(r11-r1);
				for (vector<int>::iterator it2 = n2.begin(); it2 != n2.end(); ++it2)
				{
					int r22 = *it2;
					if (r11 == r1 && r22 == r2)
						continue; // ignore r1-r2
                    // ignore intra-residue to r1 or r2
					if (r11==r22 && (r11 == r1 || r11 == r2))
						continue;
					int diff2 = abs(r22-r2);
					if ( ((r11-r1 >= 0 && r22-r2 >= 0) || (r1-r11 >= 0 && r2-r22 >= 0)) && abs(diff1-diff2) <= 1 )
					{
						if (overlapFlag && r11 >= r2Begin && r11 <= r1End && r22 >= r2Begin && r22 <= r1End && r11 > r22)
							continue; // avoid double counting

						countAss += (r11 != r22 ? assignPossibMatrixDirected[r11-1][r22-1]+assignPossibMatrixDirected[r22-1][r11-1] : assignPossibMatrixDirected[r11-1][r22-1]);
						expectedCountAss += 2*expectedAssignAvg[r11-1][r22-1]; // 2* to consider both directions (even for r11 == r22)

						strNumer += (r11 != r22 ? assignStructMatrix[r11-1][r22-1] + assignStructMatrix[r22-1][r11-1] : assignStructMatrix[r11-1][r22-1]);
						strDenom += (r11 != r22 ? structMatrix[r11-1][r22-1] + structMatrix[r22-1][r11-1] : structMatrix[r11-1][r22-1]);
					}
				}
			}
			if (expectedCountAss > 0)
			{
				netCS[r1-1][r2-1] = double(countAss)/double(expectedCountAss);
				if (r1 != r2)
					netCS[r2-1][r1-1] = netCS[r1-1][r2-1];
			}
			else if (countAss > 0)
			{
				netCS[r1-1][r2-1] = 1.0;
				if (r1 != r2)
					netCS[r2-1][r1-1] = 1.0;
			}
			if (strDenom > 0)
			{
				netStr[r1-1][r2-1] = double(strNumer)/double(strDenom);
				if (r1 != r2)
					netStr[r2-1][r1-1] = netStr[r1-1][r2-1];
			}

			// antiparallel
			countAss = 0; // sum of num assignment possibilities between r11, r22
			expectedCountAss = 0; // expected sum of num assignment possibilities between r11, r22
			strNumer = 0;
			strDenom = 0;
			for (vector<int>::iterator it1 = n1.begin(); it1 != n1.end(); ++it1)
			{
				int r11 = *it1;
				int diff1 = abs(r11-r1);
				for (vector<int>::iterator it2 = n2.begin(); it2 != n2.end(); ++it2)
				{
					int r22 = *it2;
					if (r11 == r1 && r22 == r2)
						continue; // ignore r1-r2
					int diff2 = abs(r22-r2);
					if ( ((r11-r1 >= 0 && r2-r22 >= 0) || (r1-r11 >= 0 && r22-r2 >= 0)) && abs(diff1-diff2) <= 1 )
					{
						if (overlapFlag && r11 >= r2Begin && r11 <= r1End && r22 >= r2Begin && r22 <= r1End && r11 > r22)
							continue; // avoid double counting
						countAss += (r11 != r22 ? assignPossibMatrixDirected[r11-1][r22-1]+assignPossibMatrixDirected[r22-1][r11-1] : assignPossibMatrixDirected[r11-1][r22-1]);
						expectedCountAss += 2*expectedAssignAvg[r11-1][r22-1]; // 2* to consider both directions (even for r11 == r22)

						strNumer += (r11 != r22 ? assignStructMatrix[r11-1][r22-1] + assignStructMatrix[r22-1][r11-1] : assignStructMatrix[r11-1][r22-1]);
						strDenom += (r11 != r22 ? structMatrix[r11-1][r22-1] + structMatrix[r22-1][r11-1] : structMatrix[r11-1][r22-1]);
					}
				}
			}
			// take max of parallel and antiparallel result
			if (expectedCountAss > 0)
			{
				netCS[r1-1][r2-1] = max(double(countAss)/double(expectedCountAss),netCS[r1-1][r2-1]);
				if (r1 != r2)
					netCS[r2-1][r1-1] = netCS[r1-1][r2-1];
			}
			else if (countAss > 0)
			{
				netCS[r1-1][r2-1] = 1.0;
				if (r1 != r2)
					netCS[r2-1][r1-1] = 1.0;
			}
			if (strDenom > 0)
			{
				netStr[r1-1][r2-1] = max(double(strNumer)/double(strDenom),netStr[r1-1][r2-1]);
				if (r1 != r2)
					netStr[r2-1][r1-1] = netStr[r1-1][r2-1];
			}
		} // end for r2
	} // for r1

	// compute symmetry, inter res, network cs, network str scores
	for (list<Assignment*>::iterator it = asspos.begin(); it != asspos.end(); it++)
	{
		Assignment* as = *it;
		Residue* res1 = as->contact.r1;
		Atom* x1Atom = as->contact.x1;
		Atom* hx1Atom = as->contact.hx1;
		Residue* res2 = as->contact.r2;
		Atom* x2Atom = as->contact.x2;
		Atom* hx2Atom = as->contact.hx2;
		ExpectedAssignKey k(res1->num,hx1Atom->name,res2->num,hx2Atom->name);
		Contact c2(res2,x2Atom,hx2Atom,res1,x1Atom,hx1Atom);
		tr1::unordered_set<Contact>::const_iterator c2Find = contactAss->find(c2);
		int assignCorrection = 1; // includes as->contact and possibly the symmetric contact
		if (c2Find != contactAss->end())
		{
			as->score.sym = 1.0;
			assignCorrection = 2;
		}
		else
			as->score.sym = 0;
		int r1 = res1->num;
		int r2 = res2->num;
		int expectedAssignCorrection = 1; // ignore contact between hx1Atom and hx2Atom
		if (hasExpectedSymmetricPeak(res1,x1Atom,hx1Atom,res2,x2Atom,hx2Atom,has15N,hasCali,hasCaro))
			expectedAssignCorrection = 2;
		// interres score ignores as and sym if sym=1.0

		int expAssign = 2*expectedAssign[k]; // times 2 to count both directions
		int assignPossib12 = (r1 != r2 ? assignPossibMatrixDirected[r1-1][r2-1]+assignPossibMatrixDirected[r2-1][r1-1] : assignPossibMatrixDirected[r1-1][r2-1]);
		if (assignPossib12-assignCorrection > 0)
		{
			if (expAssign-expectedAssignCorrection > 0)
			{
				as->score.interres = double(assignPossib12-assignCorrection)/double(expAssign-expectedAssignCorrection);
			}
			else
				as->score.interres = 1.0;
		}
		else
			as->score.interres = 0;

		as->score.net = netCS[r1-1][r2-1];
		as->score.netStr = netStr[r1-1][r2-1];

		// clamp scores depending on expectedAssign
		if (as->score.interres > 1.0)
			as->score.interres = 1.0;
		else if (as->score.interres < 0)
			as->score.interres = 0;

		if (as->score.net > 1.0)
			as->score.net = 1.0;
		else if (as->score.net < 0)
			as->score.net = 0;
	} // end for each asspos
	delete contactAss;
}


// new version; list version; deletes filtered out PCAssignments
// if deletePCA = true, then filtered PCA in asspos will be deleted, otherwise must delete externally
// Normally filterAssPos sets deltePCA to true and filterAss sets to false since asspos and ass overlap
void filter(TestCondition& filterTest,
		list<PCAssignment*>& asspos, bool deletePCA)
{
	// sort asspos by score
	vector<PCAssignment*> assposV(asspos.begin(),asspos.end());
	sort(assposV.begin(),assposV.end(),PCAssignmentSortByScore()); // ascending order
	reverse(assposV.begin(),assposV.end()); // descending order
	asspos.clear(); // add from assposV to asspos for assignments that pass the filter
	list<PCAssignment*> toRemove;

	int numCorrectPossib = 0;
	int numCorrectPossib6 = 0;
	bool hasRefPDB = (refPDB.size() > 0 ? true : false);

	for (vector<PCAssignment*>::iterator it = assposV.begin(); it != assposV.end(); ++it)
	{
		PCAssignment* pca = *it;
		if (filterTest.test(pca))
			asspos.push_back(pca);
		else
			toRemove.push_back(pca);
		if (hasRefPDB)
		{
			tr1::unordered_set<Contact>& contacts = pca->pc.getContacts();
			bool inContact  = false;
			for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
			{
				const Contact& c = *itC;
				ContactMapIterator itFind = refContactMap.find(c);
				if (itFind != refContactMap.end())
				{
					const tr1::array<double,ContactMap::CMENTRYSIZE>& entry = itFind.second();
					if (entry[ContactMap::MINDISTINDEX] < DISTCUTOFF)
					{
						inContact = true;
						break;
					}
				}
			}
			if (inContact)
			{
				numCorrectPossib++;
				if (pca->pc.getSeqSep() >= LRSEQSEP)
					numCorrectPossib6++;
			}
		}
	}

	int numCorrectlyFiltered = 0;
	int numCorrectlyFiltered6 = 0; // for >= LRSEQSEP
	int numFiltered = 0;
	int numFiltered6 = 0;

	for (list<PCAssignment*>::iterator it = toRemove.begin(); it != toRemove.end(); ++it)
	{
		// check results and delete pca
		PCAssignment* pca = *it;
		if (hasRefPDB)
		{
			tr1::unordered_set<Contact>& contacts = pca->pc.getContacts();
			bool inContact  = false;
			for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
			{
				const Contact& c = *itC;
				ContactMapIterator itFind = refContactMap.find(c);
				if (itFind != refContactMap.end())
				{
					const tr1::array<double,ContactMap::CMENTRYSIZE>& entry = itFind.second();
					if (entry[ContactMap::MINDISTINDEX] < DISTCUTOFF)
					{
						inContact = true;
						break;
					}
				}
			}
			if (!inContact)
			{
				numCorrectlyFiltered++;
				if (pca->pc.getSeqSep() >= LRSEQSEP)
					numCorrectlyFiltered6++;
			}
		}
		numFiltered++;
		if (pca->pc.getSeqSep() >= LRSEQSEP)
			numFiltered6++;
		if (deletePCA)
			delete pca;
	}
	if (hasRefPDB)
	{
		double fr = (numFiltered > 0 ? double(numCorrectlyFiltered)/double(numFiltered) : 1.0);
		printf("numCorrectlyFiltered: %d  numFiltered: %d   fract: %5.3f\n",numCorrectlyFiltered,numFiltered,fr);
		double fr6 = (numFiltered6 > 0 ? double(numCorrectlyFiltered6)/double(numFiltered6) : 1.0);
		printf("numCorrectlyFiltered6: %d  numFiltered6: %d   fract6: %5.3f\n",numCorrectlyFiltered6,numFiltered6,fr6);

		int numInCorrectlyFiltered = numFiltered-numCorrectlyFiltered;
		double frCorrect = (double)numInCorrectlyFiltered/(double)numCorrectPossib;
		printf("numInCorrectlyFiltered: %d  numCorrectPossib: %d   fract: %5.3f\n",numInCorrectlyFiltered,numCorrectPossib,frCorrect);

		int numInCorrectlyFiltered6 = numFiltered6-numCorrectlyFiltered6;
		double frCorrect6 = (double)numInCorrectlyFiltered6/(double)numCorrectPossib6;
		printf("numInCorrectlyFiltered6: %d  numCorrectPossib6: %d   fract6: %5.3f\n",numInCorrectlyFiltered6,numCorrectPossib6,frCorrect6);
	}
	else
		printf("numFiltered: %d\n",numFiltered);
}

// new version; vector version;
// // seed assignments in asspos are filtered out; deletes filtered out PCAssignments if
// deletePCA = true, then filtered PCA in asspos will be deleted, otherwise must delete externally
// Normally filterAssPos sets deltePCA to true and filterAss sets to false since asspos and ass overlap
// if deletePCA=false, deleted will contain the deleted assignments (normally some seed assignments)

void filter(TestCondition& filterTest,
		vector<PCAssignment*>& asspos, bool deletePCA, list<PCAssignment*>& deleted)
{
	// sort asspos by score
	vector<PCAssignment*> assposV(asspos.begin(),asspos.end());
	sort(assposV.begin(),assposV.end(),PCAssignmentSortByScore()); // ascending order
	reverse(assposV.begin(),assposV.end()); // descending order
	asspos.clear(); // add from assposV to asspos for assignments that pass the filter
	list<PCAssignment*> toRemove;

	for (vector<PCAssignment*>::iterator it = assposV.begin(); it != assposV.end(); ++it)
	{
		PCAssignment* pca = *it;
		if (filterTest.test(pca))
		{
			asspos.push_back(pca);
		}
		else
		{
			toRemove.push_back(pca);
		}
	}
	int numCorrectlyFiltered = 0;
	int numFiltered = 0;
	for (list<PCAssignment*>::iterator it = toRemove.begin(); it != toRemove.end(); ++it)
	{
		// check results and delete pca
		PCAssignment* pca = *it;
		if (refPDB.size() > 0)
		{
			tr1::unordered_set<Contact>& contacts = pca->pc.getContacts();
			bool inContact  = false;
			for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
			{
				const Contact& c = *itC;
				ContactMapIterator itFind = refContactMap.find(c);
				if (itFind != refContactMap.end())
				{
					const tr1::array<double,ContactMap::CMENTRYSIZE>& entry = itFind.second();
					if (entry[ContactMap::MINDISTINDEX] < DISTCUTOFF)
					{
						inContact = true;
						break;
					}
				}
			}
			if (!inContact)
				numCorrectlyFiltered++;
		}
		numFiltered++;

		if (deletePCA)
			delete pca;
		else
			deleted.push_back(pca);
	}
	if (refPDB.size() > 0)
	{
		double fr = (numFiltered > 0 ? double(numCorrectlyFiltered)/double(numFiltered) : 1.0);
		printf("numCorrectlyFiltered: %d  numFiltered: %d   fract: %5.3f   NewAssignedCount: %zu\n",numCorrectlyFiltered,
				numFiltered,fr,asspos.size());
	}
	else
		printf("numFiltered: %d\n",numFiltered);
}

// initializes assignMatrix based using assigned
// used by initializeAssignmentFilter
void makeAssignMatrix(vector<PCAssignment*>& assigned)
{
	for (int i = 0; i < MAXPROSIZE; i++)
	{
		for (int j = i; j < MAXPROSIZE; j++)
		{
			assignMatrix[i][j] = 0;
			if (j != i)
				assignMatrix[j][i] = 0;
		}
	}
	for (vector<PCAssignment*>::iterator it = assigned.begin(); it != assigned.end(); ++it)
	{
		PCAssignment* pca = *it;
		PseudoContact& pc = pca->pc;
		int r1 = pc.res1;
		int r2 = pc.res2;
		assignMatrix[r1-1][r2-1]++;
		if (r1 != r2)
			assignMatrix[r2-1][r1-1]++;
	}
}

// initializes assignPossibMatrix based on the assignment possibilities in assignments
// used by initializeAssignmentFilter and initializeAmbigAssignmentFilter
void makeAssignPossibMatrix(list<PCAssignment*>& assignments)
{
	for (int i = 0; i < MAXPROSIZE; i++)
	{
		for (int j = i; j < MAXPROSIZE; j++)
		{
			assignPossibMatrix[i][j] = 0;
			if (j != i)
				assignPossibMatrix[j][i] = 0;
		}
	}
	tr1::unordered_set<PseudoContact>* contactAss = new tr1::unordered_set<PseudoContact>; // stores contacts with at least one assignment possibility
	                                                                           // [directed contact]; used for symmetry and netstr score
	                                        // the same contact can be assigned to multiple peaks, so use set to count the contact at most once
	for (list<PCAssignment*>::iterator it = assignments.begin(); it != assignments.end(); ++it)
	{
		PCAssignment* pca = *it;
		PseudoContact& pc = pca->pc;
		tr1::unordered_set<PseudoContact>::iterator itC = contactAss->find(pc);
		if (itC == contactAss->end())
			contactAss->insert(pc);
	}
	for (tr1::unordered_set<PseudoContact>::iterator it = contactAss->begin(); it != contactAss->end(); ++it)
	{
		const PseudoContact& c = *it;
		int r1 = c.res1;
		int r2 = c.res2;
		assignPossibMatrix[r1-1][r2-1]++;
		if (r1 != r2)
			assignPossibMatrix[r2-1][r1-1]++;
	}

	// scale the matrix by ambiguity
	for (int r1 = 1; r1 <= bmrbProtein->size; r1++)
	{
		for (int r2 = 1; r2 <= bmrbProtein->size; r2++)
		{
			if (abs(r2-r1) > 1 && assignPossibMatrix[r1-1][r2-1] > 0) // don't scale for intra-residue
			{
				double scale = getNetScaleFactor(r1,r2);
				if (scale < 1.0)
					scale = 1.0;
				assignPossibMatrix[r1-1][r2-1] *= 1.0/scale;
				if (assignPossibMatrix[r1-1][r2-1] == 0)
					assignPossibMatrix[r1-1][r2-1] = 1;
			}
		}
	}

	delete contactAss;
}

// helper function for writeAssignmentXplor
// uses dist to decide on ub and lb and clamps dist to 1.8A.
void getDistances(bool isMethyl1, bool isMethyl2, double& dist, double& ub, double& lb)
{
	if (dist < 1.8)
		dist = 1.8;
	if (dist > DISTCUTOFF)
		ub = 0;
	else
		ub = DISTCUTOFF-dist;
	lb = dist-1.8;  // lower bound the same for all bins
	if (lb < 0)
		lb = 0;
	if (isMethyl1)
	{
		ub += 0.5;
		lb += 0.5;
		if (dist-lb < 1.3)
		{
			lb = dist-1.3; // lb distance cannot be less than 1.3
		}
	}
	if (isMethyl2)
	{
		ub += 0.5;
		lb += 0.5;
		if (dist-lb < 1.3)
		{
			lb = dist-1.3;
		}
	}
}

// helper function for writeAssignmentsXplor
string selectionString(tr1::unordered_set<string>& hs, int res, int resOffset)
{
	int size = hs.size();
	string ret = "";
	string prefix = "";
	string suffix = "";
	if (size > 1)
	{
		prefix = "(";
		suffix = ")";
	}

	int i = 0;
	static const int SIZE = 40;
	for (tr1::unordered_set<string>::iterator itH = hs.begin(); itH != hs.end(); ++itH)
	{
		const string& name = *itH;
		char buffer[SIZE];
		if (i < size-1)
			snprintf(buffer,SIZE, "(resid %4d and  name %-4s) OR\n\t",(res+resOffset), name.c_str());
		else
			snprintf(buffer,SIZE, "(resid %4d and  name %-4s)",(res+resOffset), name.c_str());
		i++;
		ret += buffer;
	}

	return prefix+ret+suffix;
}

// from stdev of distance in templates to weight on dist restraint
double getWT_STRDIST(double stdev)
{
	if (stdev <= 0.5)
		return 1.0;
	else if (stdev >= 2.0)
		return 0;
	else
		return stdev*(0.333*stdev-1.5) + 1.6667; // linear from wt=1.0 at stdev=0.5 and wt=0 at stdev=2.0
}


/*
 * does constraint combination
 * output will have atom names in Xplor format
 */
void writeAssignmentsXplor(const string& filenameNOE, const string& filenameNew,
		vector<PCAssignment*>& seed, tr1::unordered_map< NOE*, vector<PCAssignment*> >& ambigAss,
		double calibConstants[NOE::NUMTYPES][NUM_CALIB_TYPES])
{
	XplorFormat::init(); // initializes pdb to psf atom name translation table
	int resOffset = -(startRes-1); // to add to res nums for output
	FILE* out = fopen(filenameNOE.c_str(), "w"); // seed assignments
	FILE* outNew = fopen(filenameNew.c_str(),"w"); // ambiguous assignments
	tr1::unordered_set<Contact> duplicates; // stores all contacts in assigned to avoid duplicates (used for restraints from templates)
	int numRestrNOERefine = 0; // num restraints in noe_refine.tbl
	double scoreNOERefine = 0; // the sum of the scores of the restraints in noe_refine.tbl
	int numRestrNOERefineLR = 0; // num LR (>= LRSEQSEP) restraints in noe_refine.tbl
	double scoreNOERefineLR = 0; // the sum of the scores of the LR restraints in noe_refine.tbl
	           // these assignment possibilities to improve the error tolerance

	// std::srand(7123); // for debugging

	tr1::unordered_map<PCAssignment*,double> pca2MinDist; // min distance in referance map
	if (refPDB.size() > 0)
	{
		for (tr1::unordered_map<NOE*, vector<PCAssignment*>  >::iterator it = ambigAss.begin(); it != ambigAss.end(); ++it)
		{
			vector<PCAssignment*>& ass = it->second;
			for (vector<PCAssignment*>::iterator itV = ass.begin(); itV != ass.end(); ++itV)
			{
				PCAssignment* pca = *itV;
				tr1::unordered_set<Contact>& contacts = pca->pc.getContacts();
				double minDist = 999999;
				for (tr1::unordered_set<Contact>::iterator itC1 = contacts.begin(); itC1 != contacts.end(); ++itC1)
				{
					const Contact& c = *itC1;
					double dd = refContactMap[c][ContactMap::MINDISTINDEX];
					if (dd < minDist)
						minDist = dd;
				}
				pca2MinDist[pca] = minDist;
			}
		}
		for (vector<PCAssignment*>::iterator it = seed.begin(); it != seed.end(); ++it)
		{
			PCAssignment* pca = *it;
			tr1::unordered_set<Contact>& contacts = pca->pc.getContacts();
			double minDist = 999999;
			for (tr1::unordered_set<Contact>::iterator itC1 = contacts.begin(); itC1 != contacts.end(); ++itC1)
			{
				const Contact& c = *itC1;
				double dd = refContactMap[c][ContactMap::MINDISTINDEX];
				if (dd < minDist)
					minDist = dd;
			}
			pca2MinDist[pca] = minDist;
		}
	}

	for (vector<PCAssignment*>::iterator it = seed.begin(); it != seed.end(); ++it)
	{
		PCAssignment* pca = *it;
		PseudoContact& pc = pca->pc;
		tr1::unordered_set<Contact>& contacts = pc.getContacts();
		int res1 = 0;
		int res2 = 0;
		pc.getResPair(res1,res2);
		if (res1 < startRes || res1 > endRes || res2 < startRes || res2 > endRes)
			continue;

		// update stats
		numRestrNOERefine++;
		scoreNOERefine += pca->score.total;
		int seqsep = abs(res2-res1);
		if (seqsep >= LRSEQSEP)
		{
			numRestrNOERefineLR++;
			scoreNOERefineLR += pca->score.total;
		}

		// skip restraints from noes that represent ambiguous intraresidue contacts
		// where atom set 1 overlaps with atom set 2
		// e.g. res 1 hd1 or hd2 to res 1 hd2 or hd1 (which doesn't make sense as a restraint)
		bool containsSelfRef = false;
		for (tr1::unordered_set<Contact>::iterator itC1 = contacts.begin(); itC1 != contacts.end(); ++itC1)
		{
			const Contact& c1 = *itC1;
			if (res1 != res2)
				break; // not intraresidue contact
			if (c1.hx1->name == c1.hx2->name)
			{
				containsSelfRef = true;
				break;
			}
			// check if c1.hx1->name occurs in any hx2->name
			for (tr1::unordered_set<Contact>::iterator itC2 = contacts.begin(); itC2 != contacts.end(); ++itC2)
			{
				if (itC1 == itC2)
					continue;
				const Contact& c2 = *itC2;
				if (c1.hx1->name == c2.hx2->name)
				{
					containsSelfRef = true;
					break;
				}
			}
			if (containsSelfRef)
				break;
		}
		if (containsSelfRef)
			continue;

		// setup seed assignment output (ambiguous if pseudoContact.size > 1)
		tr1::unordered_set<string> h1s; // contacts from proton set 1
		tr1::unordered_set<string> h2s; // to proton set 2
		bool containsMethyl1 = false; // true if exists at least one methyl in contacts
		bool containsMethyl2 = false;
		double bestAvgDist = INVALIDDISTANCE; // set to minimum avg distance among the contacts
		double bestStdev = -1; // set to max standard deviation

		NOE* noe = pca->noe;
		double intensity = noe->volume;
		double intDist = 0; // intensity-based distance

		for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
		{
			const Contact& c = *itC;
			duplicates.insert(c);
			duplicates.insert(c.reverse());
			tr1::array<double,ContactMap::CMENTRYSIZE>& entry = contactMap[c];
			if (entry[ContactMap::AVGINDEX] != INVALIDDISTANCE && entry[ContactMap::AVGINDEX] < bestAvgDist)
				bestAvgDist = entry[ContactMap::AVGINDEX];
			if (entry[ContactMap::STDEVINDEX] != INVALIDDISTANCE && entry[ContactMap::STDEVINDEX] > bestStdev)
				bestStdev = entry[ContactMap::STDEVINDEX];

			bool isMethyl1 = c.x1->isMethyl(); // true if h1 in c is a methyl
			bool isMethyl2 = c.x2->isMethyl();
			if (isMethyl1)
				containsMethyl1 = true;
			if (isMethyl2)
				containsMethyl2 = true;

			string name1 = XplorFormat::translateProton(c.r1,resOffset,isMethyl1,c.hx1);
			string name2 = XplorFormat::translateProton(c.r2,resOffset,isMethyl2,c.hx2);
			h1s.insert(name1);
			h2s.insert(name2);

			if (intensity > 0)
			{
				CALIBRATION_TYPE contactType = getCalibContactType(c);
				double calibConstant = calibConstants[noe->type][contactType];
				intDist += calibConstant;
			}
		} // end for each contact in this pc
		if (intensity > 0)
			intDist = pow(intDist*intensity/(double)(contacts.size()), -1.0/6.0);
		else
			intDist = DISTCUTOFF;
		if (intDist > DISTCUTOFF)
			intDist = DISTCUTOFF;

		if (bestAvgDist > DISTCUTOFF || bestAvgDist == INVALIDDISTANCE)
			bestAvgDist = DISTCUTOFF;

		double dist = DISTCUTOFF;
		double ub = 0;
		double lb = 0;

		if (bestStdev > -1)
		{
			if (WT_STRDIST < 0)  // set wt based on stdev
			{
				double wtStr = getWT_STRDIST(bestStdev);
				double wtInt = 1.0-wtStr;
				dist = wtStr*bestAvgDist + wtInt*intDist;
			}
			else
				dist = WT_STRDIST*bestAvgDist + WT_INTDIST*intDist;
		}
		else
		{
			dist = intDist;
		}
		// set the upper distance bound
		string atoms1 = selectionString(h1s, res1, resOffset);
		string atoms2 = selectionString(h2s, res2, resOffset);
		getDistances(containsMethyl1, containsMethyl2, dist, ub, lb);

		fprintf(out,"assign %s\n\t%s\n\t%3.1f %3.1f %3.1f !! %6.2f %6.3f %6.3f %10.1f %s %4.1f\n",atoms1.c_str(),atoms2.c_str(),
				dist,lb,ub, noe->x, noe->hx, noe->h, noe->volume,pca->score.toString().c_str(),pca2MinDist[pca]);
	} // end for each seed assignment
	// ambiguous restraints
	if (ambigAss.size() > 0)
	{
		tr1::unordered_map< PCAssignment*, pair< string, tr1::array< double, 4 > > > pca2Dist; // pca -> ((atom pairs string), (dist, lb, ub, stdev))
		vector< pair< NOE*, vector< PCAssignment* > > > pool; // low scoring ambiguous assignments are placed here for restraint combination
		pool.reserve(ambigAss.size());
		int poolCount = 0; // number of assignments in entire pool
		double poolSumScore = 0; // sum of scores in entire pool
		for (tr1::unordered_map<NOE*, vector<PCAssignment*>  >::iterator it = ambigAss.begin(); it != ambigAss.end(); ++it)
		{
			NOE* noe = it->first;
			vector<PCAssignment*>& ass = it->second;
			int totalCount = 0; // count of ass
			double totalScore = 0; // sum of scores in ass
			// get the totalScore
			for (vector<PCAssignment*>::iterator itV = ass.begin(); itV != ass.end(); ++itV)
			{
				PCAssignment* pca = *itV;
				PseudoContact& pc = pca->pc;
				int res1 = 0;
				int res2 = 0;
				pc.getResPair(res1,res2);
				if (res1 < startRes || res1 > endRes || res2 < startRes || res2 > endRes)
					continue;
				totalScore += pca->score.total;
				totalCount++;
			}
			poolSumScore += totalScore;
			poolCount += totalCount;

			double intensity = noe->volume;
			double intDist = 0; // intensity derived-distance
			// put the assignments into tbl
			for (vector<PCAssignment*>::iterator itV = ass.begin(); itV != ass.end(); ++itV)
			{
				PCAssignment* pca = *itV;
				PseudoContact& pc = pca->pc;
				tr1::unordered_set<Contact>& contacts = pc.getContacts();
				int res1 = 0;
				int res2 = 0;
				pc.getResPair(res1,res2);
				if (res1 < startRes || res1 > endRes || res2 < startRes || res2 > endRes)
					continue;
				tr1::unordered_set<string> h1s; // contacts from proton set 1
				tr1::unordered_set<string> h2s; // to proton set 2

				bool containsMethyl1 = false; // true if exists at least one methyl in contacts
				bool containsMethyl2 = false;
				double bestAvgDist = INVALIDDISTANCE; // set to minimum distance among the contacts
				double bestStdev = -1; // set to max stdev

				double weight = pca->score.total/(double)(contacts.size()); // from distance calibration

				for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
				{
					const Contact& c = *itC;
					duplicates.insert(c);
					duplicates.insert(c.reverse());
					tr1::array<double,ContactMap::CMENTRYSIZE>& entry = contactMap[c];
					if (entry[ContactMap::AVGINDEX] != INVALIDDISTANCE && entry[ContactMap::AVGINDEX] < bestAvgDist)
						bestAvgDist = entry[ContactMap::AVGINDEX];
					if (entry[ContactMap::STDEVINDEX] != INVALIDDISTANCE && entry[ContactMap::STDEVINDEX] > bestStdev)
						bestStdev = entry[ContactMap::STDEVINDEX];

					bool isMethyl1 = c.x1->isMethyl(); // true if h1 in c is a methyl
					bool isMethyl2 = c.x2->isMethyl();
					if (isMethyl1)
						containsMethyl1 = true;
					if (isMethyl2)
						containsMethyl2 = true;

					string name1 = XplorFormat::translateProton(c.r1,resOffset,isMethyl1,c.hx1);
					string name2 = XplorFormat::translateProton(c.r2,resOffset,isMethyl2,c.hx2);
					h1s.insert(name1);
					h2s.insert(name2);

					if (intensity > 0)
					{
						CALIBRATION_TYPE contactType = getCalibContactType(c);
						double calibConstant = calibConstants[noe->type][contactType];
						intDist += calibConstant*weight;
					}
				} // end for each contact in this pc

				if (bestAvgDist > DISTCUTOFF || bestAvgDist == INVALIDDISTANCE)
					bestAvgDist = DISTCUTOFF;

				double dist = bestAvgDist;
				double ub = 0;
				double lb = 0;
				// set the upper distance bound
				string atoms1 = selectionString(h1s, res1, resOffset);
				string atoms2 = selectionString(h2s, res2, resOffset);
				getDistances(containsMethyl1, containsMethyl2, dist, ub, lb);

				tr1::array<double,4> dists;
				dists[0] = dist;
				dists[1] = lb;
				dists[2] = ub;
				dists[3] = bestStdev;
				pair< string,tr1::array<double,4> > restr(atoms1+"\n\t"+atoms2,dists);
				pca2Dist[pca] = restr;
			} // end for each assignment possibility for this noe
			// setup intensity-based distance
			if (intensity > 0)
			{
				intDist = pow(intDist*intensity/totalScore, -1.0/6.0);
				if (intDist > DISTCUTOFF)
					intDist = DISTCUTOFF;
			}
			else
				intDist = DISTCUTOFF;

			double maxDist = 0;// max dist <= DISTCUTOFF
			double maxStdev = 0;
			double maxUB = 0;
			double minLB = INVALIDDISTANCE;
			// get maxDist, minDist from the structures and combine with intensity-based distance
			for (vector<PCAssignment*>::iterator itV = ass.begin(); itV != ass.end(); ++itV)
			{
				PCAssignment* pca = *itV;
				pair< string,tr1::array<double,4> >& restr = pca2Dist[pca];
				double dist = restr.second[0];
				double stdev = restr.second[3];

				if (stdev > -1)
				{
					if (dist+restr.second[2] > maxUB)
						maxUB = dist+restr.second[2];
					if (dist-restr.second[1] < minLB)
						minLB = dist-restr.second[1];
				}
				if (dist > maxDist)
					maxDist = dist;
				if (stdev > maxStdev)
					maxStdev = stdev;
			}
			double dist = DISTCUTOFF;
			double ub = maxUB-maxDist;
			double lb = maxDist-minLB;
			if (minLB != INVALIDDISTANCE)
			{
				if (WT_STRDIST < 0)  // set wt based on stdev
				{
					double wtStr = getWT_STRDIST(maxStdev);
					double wtInt = 1.0-wtStr;
					dist = wtStr*maxDist + wtInt*intDist;
				}
				else
					dist = WT_STRDIST*maxDist + WT_INTDIST*intDist;
			}
			else
			{
				dist = intDist;
				ub = DISTCUTOFF-dist;
				lb = dist-1.5;
			}
			if (dist-lb < 1.5)
				lb = dist-1.5;
			if ( (!POOL_FILTER_BY_COUNT && totalScore >= POOL_FILTER) || (POOL_FILTER_BY_COUNT && totalCount >= POOL_FILTER) )
			{   // output assignments from tbl if not marked for the pool
				int i = 0;
				for (vector<PCAssignment*>::iterator itV = ass.begin(); itV != ass.end(); ++itV)
				{
					PCAssignment* pca = *itV;
					pair< string,tr1::array<double,4> >& restr = pca2Dist[pca];
					if (i != 0)
					{
						fprintf(outNew,"or\t%s !! %s %4.1f\n",restr.first.c_str(),pca->score.toString().c_str(),pca2MinDist[pca]);
					}
					else
						fprintf(outNew,"assign %s\n\t%3.1f %3.1f %3.1f !! %6.2f %6.3f %6.3f %10.1f %s %4.1f\n",
								restr.first.c_str(), dist, lb, ub, noe->x, noe->hx, noe->h, noe->volume,
								pca->score.toString().c_str(),pca2MinDist[pca]);
					i++;
				}
			}
			else
			{
				// just update dist, lb, ub  for each entry in tbl, then add to pool
				for (vector<PCAssignment*>::iterator itV = ass.begin(); itV != ass.end(); ++itV)
				{
					PCAssignment* pca = *itV;
					pair< string,tr1::array<double,4> >& restr = pca2Dist[pca];
					restr.second[0] = dist;
					restr.second[1] = lb;
					restr.second[2] = ub;
				}
				pair< NOE*, vector< PCAssignment* > > poolEntry(noe,ass);
				pool.push_back(poolEntry);
			}
		} // end for each ambiguous assignment
		// output assignments from pool (constraint combination)
		if (pool.size() > 0 && ((!POOL_FILTER_BY_COUNT && poolSumScore >= POOL_FILTER) || (POOL_FILTER_BY_COUNT && poolCount >= POOL_FILTER))) // skip pool if not enough assignments
		{
			random_shuffle(pool.begin(),pool.end());
			int numPoolRestraints = 0;
			double maxDist = 0;// max dist <= DISTCUTOFF
			double maxUB = 0;
			double minLB = INVALIDDISTANCE;
			vector<string> atomNames;
			vector<PCAssignment*> pcas;
			double currentScore = 0;
			int currentCount = 0;
			int index = 0;
			int poolSize = pool.size();
			bool breakFlag = false; // for adding last restraint
			for (vector< pair< NOE*, vector< PCAssignment* > > >::iterator itV = pool.begin(); itV != pool.end(); ++itV)
			{
				vector< PCAssignment* >& noeAmbigRestr = itV->second;
				bool testDist = true;
				for (vector< PCAssignment* >::iterator itV2 = noeAmbigRestr.begin(); itV2 != noeAmbigRestr.end(); ++itV2)
				{
					PCAssignment* pca = *itV2;
					pair< string,tr1::array<double,4> >& restr = pca2Dist[pca];
					if (testDist)
					{
						double dist = restr.second[0];
						if (dist > maxDist)
							maxDist = dist;
						double ub = restr.second[2]+dist;
						if (ub > maxUB)
							maxUB = ub;
						double lb = dist-restr.second[1];
						if (lb < minLB)
							minLB = lb;
						testDist = false;
					}
					atomNames.push_back(restr.first);
					pcas.push_back(pca);
					currentScore += pca->score.total;
					currentCount++;
				}
				if ( (!POOL_FILTER_BY_COUNT && currentScore >= POOL_FILTER) || (POOL_FILTER_BY_COUNT && currentCount >= POOL_FILTER) ) // add
				{
					double dist = maxDist;
					double ub = maxUB-maxDist;
					double lb = maxDist-minLB;
					if (dist-lb < 1.5)
						lb = dist-1.5;
					int i = 0;
					vector<PCAssignment*>::iterator itP = pcas.begin();
					for (vector<string>::iterator itA = atomNames.begin(); itA != atomNames.end(); ++itA)
					{
						string& s = *itA;
						PCAssignment* p = *itP;
						NOE* n = p->noe;
						if (i != 0)
							fprintf(outNew,"or\t%s !! %6.2f %6.3f %6.3f %10.1f %s %4.1f\n",s.c_str(), n->x, n->hx, n->h,
									n->volume, p->score.toString().c_str(),pca2MinDist[p]);
						else
							fprintf(outNew,"assign %s\n\t%3.1f %3.1f %3.1f !! %6.2f %6.3f %6.3f %10.1f %s %4.1f\n",
									s.c_str(), dist, lb, ub, n->x, n->hx, n->h, n->volume,
									p->score.toString().c_str(),pca2MinDist[p]);
						i++;
						++itP;
					}
					numPoolRestraints++;
					if (breakFlag)
						break;
					maxDist = 0;
					maxUB = 0;
					minLB = INVALIDDISTANCE;
					atomNames.clear();
					pcas.clear();
					currentScore = 0;
					currentCount = 0;
				}
				else if (index == poolSize-1)
				{
					// wrap around
					breakFlag = true;
					itV = pool.begin();
				}
				index++;
			} // end for each entry in pool
			printf("Num Pool Restraints: %d\n",numPoolRestraints);
		} // end if enough restraints in the pool
		else
			printf("Num Pool Restraints: 0\n");
	} // end ambiguous assignments

	fclose(out);
	fclose(outNew);

	//printf("Xplor input stats . . .\n");
	//printf("NumRestrNOERefine: %d\n",numRestrNOERefine);
	//printf("ScoreNOERefine: %f\n",scoreNOERefine);
	//printf("NumRestrNOERefineLR: %d\n",numRestrNOERefineLR);
	//printf("ScoreNOERefineLR: %f\n",scoreNOERefineLR);

//	if (numLRContacts > 0)
//	{
//		double LRRatio = min(double(numRestrNOERefineLR)/double(numLRContacts),1.0);
//		printf("LRRatio: %f\n",LRRatio);
//	}
//	else
//		printf("LRRatio: 0\n");
//	printf("numLRContactsForRatio: %d\n",numLRContacts); // >= LRSEQSEPfrom and < 4.0A; contactMapExpected; each contact treated as undirected, so if have contact and symmetric contact, then they count as one
}

// used by writeITASSER
namespace std
{
	namespace tr1
	{
		template <>
		class hash< pair< pair<int,int>, Contact::ResidueContactType > >
		{
			public:
			size_t operator()(const pair< pair<int,int>, Contact::ResidueContactType>& c) const
			{
				const pair<int,int>& p = c.first;
				return hash<int>()(p.first) ^ hash<int>()(p.second) ^ hash<int>()((int)(c.second));
			}
		};
	}
}

// Used by writeITASSER
// distance is the distance of the atom-atom contact (obtained from the intensity)
// assumes c is not an intraresidue contact
// This distance is converted to a residue-based distance based on a table lookup of scscAvgDist
double getResidueBasedDistanceIntensity(Contact& c,double distance)
{
	string& atom1 = c.hx1->name;
	string& atom2 = c.hx2->name;
	int binIndex = 0;
	for (int i = 0; i < NUMDISTBINS; i++)
	{
		if (distance < UPPERDISTBINS[i]+0.001)
			break;
		binIndex++;
	}
	if (binIndex >= NUMDISTBINS)
		binIndex = NUMDISTBINS-1;
	SCSCDistKey k(c.r1->num,atom1,c.r2->num,atom2,binIndex);
	tr1::unordered_map<SCSCDistKey,double>::iterator itFind = scscAvgDist.find(k);
	if (itFind != scscAvgDist.end())
	{
		return itFind->second;
	}
	return 10.0;
}

// Used by writeITASSER
// if residue r has no coordinates centroids[r-1][0] will return INVALIDCOORD
// centroids[template index in same order as structures forward iterator][resnum-1][0=x, 1=y, 2=z]
double getResidueBasedDistance(Contact& c, vector< vector< tr1::array<double,3> > > centroids)
{
	double avgDist = 0;
	int count = 0;

	// determine contact type from c
	if (c.type == Contact::NN)
	{
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];

			Atom* a1 = r1->getN("N");
			Atom* a2 = r2->getN("N");
			if (a1 != NULL && a2 != NULL)
			{
				avgDist += a1->getDistance(a2);
				count++;
			}
		}
	}
	else if (c.type == Contact::CACA)
	{
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];

			Atom* a1 = r1->getC("CA");
			Atom* a2 = r2->getC("CA");
			if (a1 != NULL && a2 != NULL)
			{
				avgDist += a1->getDistance(a2);
				count++;
			}
		}
	}
	else if (c.type == Contact::CBCB)
	{
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];

			Atom* a1 = r1->getC("CB");
			Atom* a2 = r2->getC("CB");
			if (a1 != NULL && a2 != NULL)
			{
				avgDist += a1->getDistance(a2);
				count++;
			}
		}
	}
	else if (c.type == Contact::SCSC)
	{
		int t = 0;
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];
			tr1::array<double,3>& c1 = centroids[t][r1->num-1];
			tr1::array<double,3>& c2 = centroids[t][r2->num-1];
			if (c1[0] != INVALIDCOORD && c2[0] != INVALIDCOORD)
			{
				double dx = c1[0]-c2[0];
				double dy = c1[1]-c2[1];
				double dz = c1[2]-c2[2];
				dx *= dx;
				dy *= dy;
				dz *= dz;
				avgDist += sqrt(dx+dy+dz);
				count++;
			}
			t++;
		}
	}
	else if (c.type == Contact::NCA)
	{
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];

			Atom* a1 = r1->getN("N");
			Atom* a2 = r2->getC("CA");
			if (a1 != NULL && a2 != NULL)
			{
				avgDist += a1->getDistance(a2);
				count++;
			}
		}
	}
	else if (c.type == Contact::CAN)
	{
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];

			Atom* a1 = r1->getC("CA");
			Atom* a2 = r2->getN("N");
			if (a1 != NULL && a2 != NULL)
			{
				avgDist += a1->getDistance(a2);
				count++;
			}
		}
	}
	else if (c.type == Contact::NCB)
	{
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];

			Atom* a1 = r1->getN("N");
			Atom* a2 = r2->getC("CB");
			if (a1 != NULL && a2 != NULL)
			{
				avgDist += a1->getDistance(a2);
				count++;
			}
		}
	}
	else if (c.type == Contact::CBN)
	{
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];

			Atom* a1 = r1->getC("CB");
			Atom* a2 = r2->getN("N");
			if (a1 != NULL && a2 != NULL)
			{
				avgDist += a1->getDistance(a2);
				count++;
			}
		}
	}
	else if (c.type == Contact::NSC)
	{
		int t = 0;
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];
			Atom* a1 = r1->getN("N");
			tr1::array<double,3>& c2 = centroids[t][r2->num-1];
			if (a1 != NULL && c2[0] != INVALIDCOORD)
			{
				double dx = a1->x-c2[0];
				double dy = a1->y-c2[1];
				double dz = a1->z-c2[2];
				dx *= dx;
				dy *= dy;
				dz *= dz;
				avgDist += sqrt(dx+dy+dz);
				count++;
			}
			t++;
		}
	}
	else if (c.type == Contact::SCN)
	{
		int t = 0;
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];
			tr1::array<double,3>& c1 = centroids[t][r1->num-1];
			Atom* a2 = r2->getN("N");
			if (c1[0] != INVALIDCOORD && a2 != NULL)
			{
				double dx = c1[0]-a2->x;
				double dy = c1[1]-a2->y;
				double dz = c1[2]-a2->z;
				dx *= dx;
				dy *= dy;
				dz *= dz;
				avgDist += sqrt(dx+dy+dz);
				count++;
			}
			t++;
		}
	}
	else if (c.type == Contact::CACB)
	{
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];

			Atom* a1 = r1->getC("CA");
			Atom* a2 = r2->getC("CB");
			if (a1 != NULL && a2 != NULL)
			{
				avgDist += a1->getDistance(a2);
				count++;
			}
		}
	}
	else if (c.type == Contact::CBCA)
	{
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];

			Atom* a1 = r1->getC("CB");
			Atom* a2 = r2->getC("CA");
			if (a1 != NULL && a2 != NULL)
			{
				avgDist += a1->getDistance(a2);
				count++;
			}
		}
	}
	else if (c.type == Contact::CASC)
	{
		int t = 0;
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];
			Atom* a1 = r1->getC("CA");
			tr1::array<double,3>& c2 = centroids[t][r2->num-1];
			if (a1 != NULL && c2[0] != INVALIDCOORD)
			{
				double dx = a1->x-c2[0];
				double dy = a1->y-c2[1];
				double dz = a1->z-c2[2];
				dx *= dx;
				dy *= dy;
				dz *= dz;
				avgDist += sqrt(dx+dy+dz);
				count++;
			}
			t++;
		}
	}
	else if (c.type == Contact::SCCA)
	{
		int t = 0;
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];
			tr1::array<double,3>& c1 = centroids[t][r1->num-1];
			Atom* a2 = r2->getC("CA");
			if (c1[0] != INVALIDCOORD && a2 != NULL)
			{
				double dx = c1[0]-a2->x;
				double dy = c1[1]-a2->y;
				double dz = c1[2]-a2->z;
				dx *= dx;
				dy *= dy;
				dz *= dz;
				avgDist += sqrt(dx+dy+dz);
				count++;
			}
			t++;
		}
	}
	else if (c.type == Contact::CBSC)
	{
		int t = 0;
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];
			Atom* a1 = r1->getC("CB");
			tr1::array<double,3>& c2 = centroids[t][r2->num-1];
			if (a1 != NULL && c2[0] != INVALIDCOORD)
			{
				double dx = a1->x-c2[0];
				double dy = a1->y-c2[1];
				double dz = a1->z-c2[2];
				dx *= dx;
				dy *= dy;
				dz *= dz;
				avgDist += sqrt(dx+dy+dz);
				count++;
			}
			t++;
		}
	}
	else if (c.type == Contact::SCCB)
	{
		int t = 0;
		for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
		{
			CSProtein* structure = *itS;
			Residue* r1 = (*structure)[c.r1->num];
			Residue* r2 = (*structure)[c.r2->num];
			tr1::array<double,3>& c1 = centroids[t][r1->num-1];
			Atom* a2 = r2->getC("CB");
			if (c1[0] != INVALIDCOORD && a2 != NULL)
			{
				double dx = c1[0]-a2->x;
				double dy = c1[1]-a2->y;
				double dz = c1[2]-a2->z;
				dx *= dx;
				dy *= dy;
				dz *= dz;
				avgDist += sqrt(dx+dy+dz);
				count++;
			}
			t++;
		}
	}

	if (count > 0)
	{
		return avgDist/double(count);
	}
	return INVALIDDISTANCE;
}

// used by writeITASSER
// returns the avg dist and stdev of the residue-based contacts in the templates; only non-intraresidue contacts are considered
// the residue contact types are "N","CA","CB","SC"
// The returned maps are symmetric
void makeResidueContactMap(double avgDists[MAXPROSIZE][MAXPROSIZE][Contact::NUMCONTACTTYPES], double stdevDists[MAXPROSIZE][MAXPROSIZE][Contact::NUMCONTACTTYPES])
{
	vector< tr1::array< tr1::array<double,3>,MAXPROSIZE > > scs; // side chain centroids based on the heavy atoms; [template][resnum-1][0=x,1=y,2=z]
	  // ignores N, CA, CB, C, O, S when calculating the centroid; [0]=x, [1]=y, [2]=z
	// get sidechain centroid of templates
	for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
	{
		CSProtein* structure = *itS;
		tr1::array< tr1::array<double,3>, MAXPROSIZE > tempSCs; // the sc of the residues in this template
		for (int r = 1; r <= structure->size; r++)
		{
			Residue* res = (*structure)[r];

			// if residue type GLY or ALA, use CA and CB respectively for the SC
			if (res->type == Residue::GLY)
			{
				Atom* ca = res->getC("CA");
				if (ca == NULL)
					continue;
				tr1::array<double,3> coords;
				coords[0] = ca->x;
				coords[1] = ca->y;
				coords[2] = ca->z;
				tempSCs[r-1] = coords;
				continue;
			}
			else if (res->type == Residue::ALA)
			{
				Atom* cb = res->getC("CB");
				if (cb == NULL)
					continue;
				tr1::array<double,3> coords;
				coords[0] = cb->x;
				coords[1] = cb->y;
				coords[2] = cb->z;
				tempSCs[r-1] = coords;
				continue;
			}
			double x = 0;
			double y = 0;
			double z = 0;
			int count = 0;
			for (AtomIterator itX = res->begin('X'); itX != res->end('X'); itX++)
			{
				Atom* xAtom = *itX;
				if (xAtom->isNMR_SideChainHeavy())
				{
					x += xAtom->x;
					y += xAtom->y;
					z += xAtom->z;
					count++;
				}
			}
			if (count > 0)
			{
				x = x/(double)(count);
				y = y/(double)(count);
				z = z/(double)(count);

				tr1::array<double,3> coords;
				coords[0] = x;
				coords[1] = y;
				coords[2] = z;
				tempSCs[r-1] = coords;
			}
		} // end for each residue
		scs.push_back(tempSCs);
	}

	int counts[MAXPROSIZE][MAXPROSIZE][Contact::NUMCONTACTTYPES];
	for (int i = 0; i < MAXPROSIZE; i++)
	{
		for (int j = 0; j < MAXPROSIZE; j++)
		{
			for (int k = 0; k < Contact::NUMCONTACTTYPES; k++)
			{
				counts[i][j][k] = 0;
				avgDists[i][j][k] = 0;
				stdevDists[i][j][k] = 0;
			}
		}
	}
	int tempIndex = 0;
	string atomTypes[] = {"N","CA","CB","SC"}; // rsidue contact types
	for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
	{
		CSProtein* structure = *itS;
		for (int r1 = 1; r1 <= structure->size; r1++)
		{
			Residue* res1 = (*structure)[r1];
			for (int r2 = r1+1; r2 <= structure->size; r2++)
			{
				Residue* res2 = (*structure)[r2];
				for (int a1 = 0; a1 < 4; a1++)
				{
					string& atom1 = atomTypes[a1];
					double x1,y1,z1;
					if (atom1 == "SC")
					{
						tr1::array<double,3>& coords1 = scs[tempIndex][r1-1];
						x1 = coords1[0];
						y1 = coords1[1];
						z1 = coords1[2];
					}
					else
					{
						Atom* atm1 = res1->getX(atom1);
						if (atm1 == NULL)
							continue;
						x1 = atm1->x;
						y1 = atm1->y;
						z1 = atm1->z;
					}
					for (int a2 = 0; a2 < 4; a2++)
					{
						string& atom2 = atomTypes[a2];
						double x2,y2,z2;
						if (atom2 == "SC")
						{
							tr1::array<double,3>& coords2 = scs[tempIndex][r2-1];
							x2 = coords2[0];
							y2 = coords2[1];
							z2 = coords2[2];
						}
						else
						{
							Atom* atm2 = res2->getX(atom2);
							if (atm2 == NULL)
								continue;
							x2 = atm2->x;
							y2 = atm2->y;
							z2 = atm2->z;
						}
						double dist = distance(x1,y1,z1,x2,y2,z2);
						Contact::ResidueContactType type = Contact::getType(atom1,atom2);
						Contact::ResidueContactType type2 = Contact::getReverseType(type);
						counts[r1-1][r2-1][type]++;
						counts[r2-1][r1-1][type2]++;
						double delta = dist-avgDists[r1-1][r2-1][type];
						avgDists[r1-1][r2-1][type] = avgDists[r1-1][r2-1][type] + delta/(double)(counts[r1-1][r2-1][type]);
						avgDists[r2-1][r1-1][type2] = avgDists[r1-1][r2-1][type];
						stdevDists[r1-1][r2-1][type] = stdevDists[r1-1][r2-1][type] + delta*(dist-avgDists[r1-1][r2-1][type]);
						stdevDists[r2-1][r1-1][type2] = stdevDists[r1-1][r2-1][type];
					} // end for each atom a2
				} // end for each atom a1
			} // end for each residue r2
		} // end for each residue r1
		tempIndex++;
	} // end for each template

	CSProtein* pro = structures.front();
	for (int r1 = 1; r1 <= pro->size; r1++)
	{
		for (int r2 = r1+1; r2 <= pro->size; r2++)
		{
			for (int t = 0; t < Contact::NUMCONTACTTYPES; t++)
			{
				int count = counts[r1-1][r2-1][t];
				if (count > 0)
				{
					stdevDists[r1-1][r2-1][t] = sqrt(stdevDists[r1-1][r2-1][t]/(double)count);
					Contact::ResidueContactType t2 = Contact::getReverseType((Contact::ResidueContactType)t);
					stdevDists[r2-1][r1-1][t2] = stdevDists[r1-1][r2-1][t];
				}
			}
		}
	}
}

// used by writeITASSER
typedef Contact::ResidueContactType CType;
struct ContactInfo {
	int r1;
	int r2;
	CType type;
	double score;
	double dist;
};

// ignores intraresidue, sequential, and 2 res apart contacts
void writeITASSER(const string& filename, const string& filenameAmbig, vector<PCAssignment*>& seed,
		tr1::unordered_map< NOE*, vector<PCAssignment*> >& ambigAss,
		double calibConstants[NOE::NUMTYPES][NUM_CALIB_TYPES])
{
	const double CAERROR = 0.1; // 0.45; // due to using grid for CA
	const double CBERROR = 0.3; // 1.0; // 1.5/2.0; // due to 2 rotamer approximation
	const double NERROR = 0.4; // 1.45; // since CA is used to represent N
	// same order as Residue::AA3; errors due to 2 rotamer approximation of side chain center
	double SCERROR[20] = {
			0,// ALA
			6.0,// ARG
			3.75,// ASN
			3.75,// ASP
			3.0,// CYS
			4.5,// GLN
			4.5,// GLU
			0, // GLY
			4.5,// HIS
			3.375,// ILE
			3.75,// LEU
			5.25,// LYS
			4.5,// MET
			3.75,// PHE
			3.0,// PRO
			3.0,// SER
			3.0,// THR
			5.25,// TRP
			4.5,// TYR
			2.25// VAL
	}; // reduce the above values for better precision
	for (int i = 0; i < 20; i++)
		SCERROR[i] = SCERROR[i]/4.0; // /2.0;

	vector< vector< tr1::array<double,3> > > centroids; // indexed by number of templates
	centroids.reserve(structures.size());
	CSProtein& firstPro = *(structures.front());
	int NUMRES = firstPro.size;

	tr1::array<double,3> invalidCoords;
	invalidCoords[0] = INVALIDCOORD;
	invalidCoords[1] = INVALIDCOORD;
	invalidCoords[2] = INVALIDCOORD;
	// get the sidechain centroids
	for (list<CSProtein*>::iterator it = structures.begin(); it != structures.end(); ++it)
	{
		CSProtein& pro = *(*it);
		vector< tr1::array<double,3> > centroidsPro;
		centroidsPro.reserve(NUMRES);
		for (int r = 1; r <= pro.size; r++)
		{
			Residue* res = pro[r];
			tr1::array<double,3> c;
			c[0] = 0;
			c[1] = 0;
			c[2] = 0;
			if (res->type == Residue::GLY)
			{
				Atom* ca = res->getC("CA");
				if (ca != NULL)
				{
					c[0] = ca->x;
					c[1] = ca->y;
					c[2] = ca->z;
					centroidsPro.push_back(c);
				}
				else
					centroidsPro.push_back(invalidCoords);

				continue;
			}
			else if (res->type == Residue::ALA)
			{
				Atom* cb = res->getC("CB");
				if (cb != NULL)
				{
					c[0] = cb->x;
					c[1] = cb->y;
					c[2] = cb->z;
					centroidsPro.push_back(c);
				}
				else
					centroidsPro.push_back(invalidCoords);

				continue;
			}
			int count = 0;
			for (AtomIterator itX = res->begin('X'); itX != res->end('X'); itX++)
			{
				Atom* xAtom = *itX;
				string& name = xAtom->name;
				if (!starts_with(name,"H") && name != "CA" && name != "CB" && name != "N" && name != "C" && name != "O")
				{
					if (xAtom->x != INVALIDCOORD)
					{
						c[0] += xAtom->x;
						c[1] += xAtom->y;
						c[2] += xAtom->z;
						count++;
					}
				}
			}
			if (count > 0)
			{
				c[0] /= double(count);
				c[1] /= double(count);
				c[2] /= double(count);
				centroidsPro.push_back(c);
			}
			else
				centroidsPro.push_back(invalidCoords);
		} // end for each residue
		centroids.push_back(centroidsPro);
	} // end for each template

	tr1::unordered_map< pair< pair<int,int>, Contact::ResidueContactType >, pair<double,double> > restraints; // seed residue contacts; (res1,res2,type) (score,upperboundDist)
	int resOffset = -(startRes-1);
	FILE* out = fopen(filename.c_str(), "w");
	// get the seed assignments and setup restraints table from seed
	for (vector<PCAssignment*>::iterator it = seed.begin(); it != seed.end(); ++it)
	{
		PCAssignment* pca = *it;
		PseudoContact& pc = pca->pc;
		int res1 = pc.res1;
		int res2 = pc.res2;
		if (res1 < startRes || res1 > endRes || res2 < startRes || res2 > endRes)
			continue;
		if (abs(res1-res2) < 3)
			continue;
		double bestAvgDist = INVALIDDISTANCE;
		double bestFractNumStr = 0;
		tr1::unordered_set<Contact>& contacts = pc.getContacts();
		for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
		{
			const Contact& c = *itC;
			tr1::array<double,ContactMap::CMENTRYSIZE>& entry = contactMap[c];
			if (entry[ContactMap::AVGINDEX] < bestAvgDist)
				bestAvgDist = entry[ContactMap::AVGINDEX];
			if (entry[ContactMap::FRACSTRUCINDEX] > bestFractNumStr)
				bestFractNumStr = entry[ContactMap::FRACSTRUCINDEX];
		}

		double score = pca->score.total/double(contacts.size()); // each contact of pseudocontact contributes a fraction of the restraint score
		NOE* noe = pca->noe;
		double intensity = noe->volume;
		double intDistBase = 0; // intensity-based distance
		if (intensity > 0)
		{
			for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
			{
				Contact c = *itC;
				CALIBRATION_TYPE contactType = getCalibContactType(c);
				double calibConstant = calibConstants[noe->type][contactType];
				intDistBase += calibConstant;
			}
			intDistBase = pow(intDistBase*intensity/(double)(contacts.size()), -1.0/6.0);
		}
		else
		{
			intDistBase = DISTCUTOFF;
		}


		if (intDistBase > DISTCUTOFF)
			intDistBase = DISTCUTOFF;

		for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
		{
			Contact c = *itC;
			tr1::array<double,ContactMap::CMENTRYSIZE>& entry = contactMap[c]; // contactMap is symmetric
			double dist = 0;
			double intDist = getResidueBasedDistanceIntensity(c,intDistBase);
			if (intDist > 15.0)
				intDist = 15.0;
			if (bestFractNumStr > 0)
			{  // take weighted average of intensity-based distance with structure-based distance
				double avgDist = getResidueBasedDistance(c,centroids);
				if (avgDist > 15.0)
					avgDist = 15.0;
				if (WT_STRDIST < 0)
				{  // WT_STRDIST not specified by user, so pick it automatically
					if (avgDist != INVALIDDISTANCE)
					{
						double stdevDist = entry[ContactMap::STDEVINDEX];
						double wtStr = getWT_STRDIST(stdevDist);
						double wtInt = 1.0-wtStr;
						dist = wtStr*avgDist + wtInt*intDist;
					}
					else
					{
						dist = intDist;
					}

				}
				else
				{
					if (avgDist != INVALIDDISTANCE)
						dist = WT_STRDIST*avgDist + WT_INTDIST*intDist;
					else
						dist = intDist;
				}
			}
			else // use only the intensity-based distance
			{
				dist = intDist;
			}

			// update restraints
			pair<int,int> resPair(res1,res2);
			Contact::ResidueContactType type = c.type;
			// For output purpose we output contacts in one direction only (NCA, NCB, NSC, CACB, CASC, CBSC)
			// if the contact is not in any of these directions, then it is reversed so that it will be
			// For symmetric contacts NN, CACA, CBCB, SCSC, set res1 <= res2, so will only output once
			if (type == Contact::CAN || type == Contact::CBN || type == Contact::SCN || type == Contact::CBCA ||
			    type == Contact::SCCA || type == Contact::SCCB)
			{
				resPair = std::make_pair(res2,res1);
				type = c.getReverseType();
			}
			else if (type == Contact::NN || type == Contact::CACA || type == Contact::CBCB || type == Contact::SCSC)
			{
				if (res1 > res2)
					resPair = std::make_pair(res2,res1);
			}
			pair< pair<int,int>, Contact::ResidueContactType > key(resPair,type);
			tr1::unordered_map< pair< pair<int,int>, Contact::ResidueContactType >, pair<double,double> >::iterator itFind = restraints.find(key);
			if (itFind != restraints.end())
			{
				pair<double,double>& scoreDist = itFind->second;
				scoreDist.first += score;
				scoreDist.second = max<double>(scoreDist.second,dist);
			}
			else
			{
				pair<double,double> scoreDist(score,dist);
				restraints[key] = scoreDist;
			}
		} // end for each contact in pc
	} // end for each seed assignment
	// output seed assignment
	int numNonAmbig = 0;
	string outputStr = "";
	outputStr.reserve(6+35*restraints.size()); // upperbound on str length of # restraints and the restraints themselves

	for (tr1::unordered_map< pair< pair<int,int>, Contact::ResidueContactType >, pair<double,double> >::iterator it = restraints.begin(); it != restraints.end(); ++it)
	{
		const pair< pair<int,int>, Contact::ResidueContactType >& key = it->first;
		const pair<int,int>& p = key.first;
		int r1 = p.first;
		int r2 = p.second;
		if (r1 > r2)
			std::swap(r1,r2); // for accessing restraintsResPair only, where r1 <= r2
		pair<int,int> rp(r1,r2);
		const pair<double,double>& scoreDist = it->second;
		numNonAmbig++;
		Contact::ResidueContactType type = key.second;
		double dist = scoreDist.second+2*CAERROR;
		if (type == Contact::CBCB || type == Contact::CBCA || type == Contact::CBN || type == Contact::CBSC)
		{
			dist += CBERROR;
		}
		if (type == Contact::CBCB || type == Contact::CACB || type == Contact::NCB || type == Contact::SCCB)
		{
			dist += CBERROR;
		}
		if (type == Contact::SCCB || type == Contact::SCCA || type == Contact::SCN || type == Contact::SCSC)
		{
			dist += SCERROR[seqAAType[p.first-1]];
		}
		if (type == Contact::CBSC || type == Contact::CASC || type == Contact::NSC || type == Contact::SCSC)
		{
			dist += SCERROR[seqAAType[p.second-1]];
		}
		if (type == Contact::NN)
		{
			dist += 2*NERROR;
		}
		if (type == Contact::CAN || type == Contact::CBN || type == Contact::SCN)
		{
			dist += NERROR;
		}
		if (type == Contact::NCA || type == Contact::NCB || type == Contact::NSC)
		{
			dist += NERROR;
		}
		string typeStr = Contact::CONTACTTYPESTR[type];
		char buffer[35];
		snprintf(buffer,35,"%-4s %5d %5d %9.6f %6.2f\n",typeStr.c_str(),p.first+resOffset,p.second+resOffset,scoreDist.first,dist);
		outputStr += buffer;
	}
	fprintf(out,"%d\n",numNonAmbig);
	fprintf(out,"%s",outputStr.c_str());
	fclose(out);

	// setup ambiguous restraints
	tr1::unordered_map< pair< pair<int,int>, Contact::ResidueContactType >, pair<double,double> > ambigRestraints;   // [(res1,res2),contactType] -> [score,dist]
	double totalPoolScoreSum = 0; // sum of scores in ambigAss
	int totalPoolCount = 0;
	for (tr1::unordered_map<NOE*, vector<PCAssignment*>  >::iterator it = ambigAss.begin(); it != ambigAss.end(); ++it)
	{
		NOE* noe = it->first;
		double intensity = noe->volume;
		double intDistBase = 0;
		vector<PCAssignment*>& ass = it->second;
		double totalScore = 0;
		for (vector<PCAssignment*>::iterator itA = ass.begin(); itA != ass.end(); ++itA)
		{
			PCAssignment* pca = *itA;
			totalScore += pca->score.total;
		}
		totalPoolScoreSum += totalScore;
		totalPoolCount += ass.size();
		if (intensity > 0)
		{
			for (vector<PCAssignment*>::iterator itA = ass.begin(); itA != ass.end(); ++itA)
			{
				PCAssignment* pca = *itA;
				tr1::unordered_set<Contact>& contacts = pca->pc.getContacts();
				double score = pca->score.total/(totalScore*(double)(contacts.size()));
				for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
				{
					Contact c = *itC;
					CALIBRATION_TYPE contactType = getCalibContactType(c);
					double calibConstant = calibConstants[noe->type][contactType];
					intDistBase += score*calibConstant;
				}
			}
			intDistBase = pow(intDistBase*intensity, -1.0/6.0);
		}
		else
		{
			intDistBase = DISTCUTOFF;
		}
		if (intDistBase > DISTCUTOFF)
			intDistBase = DISTCUTOFF;
		// populate ambigRestraints table
		for (vector<PCAssignment*>::iterator itA = ass.begin(); itA != ass.end(); ++itA)
		{
			PCAssignment* pca = *itA;
			PseudoContact& pc = pca->pc;
			int res1 = pc.res1;
			int res2 = pc.res2;
			if (res1 < startRes || res1 > endRes || res2 < startRes || res2 > endRes)
				continue;
			if (abs(res1-res2) < 3)
				continue;
			double bestAvgDist = INVALIDDISTANCE;
			double bestFractNumStr = 0;
			tr1::unordered_set<Contact>& contacts = pc.getContacts();
			for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
			{
				const Contact& c = *itC;
				tr1::array<double,ContactMap::CMENTRYSIZE>& entry = contactMap[c];
				if (entry[ContactMap::AVGINDEX] < bestAvgDist)
					bestAvgDist = entry[ContactMap::AVGINDEX];
				if (entry[ContactMap::FRACSTRUCINDEX] > bestFractNumStr)
					bestFractNumStr = entry[ContactMap::FRACSTRUCINDEX];
			}
			double score = pca->score.total/(double)(contacts.size());
			for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
			{
				Contact c = *itC;
				tr1::array<double,ContactMap::CMENTRYSIZE>& entry = contactMap[c];
				double dist = 0;
				double intDist = getResidueBasedDistanceIntensity(c,intDistBase);
				if (intDist > 15.0)
					intDist = 15.0;
				if (bestFractNumStr > 0)
				{  // take weighted average of intensity-based distance with structure-based distance
					double avgDist = getResidueBasedDistance(c,centroids);
					if (avgDist > 15.0)
						avgDist = 15.0;
					if (avgDist > (intDist+1.0))
						avgDist = intDist+1.0;
					if (WT_STRDIST < 0)
					{  // WT_STRDIST not specified by user, so pick it automatically
						if (avgDist != INVALIDDISTANCE)
						{
							double stdevDist = entry[ContactMap::STDEVINDEX];
							double wtStr = getWT_STRDIST(stdevDist);
							double wtInt = 1.0-wtStr;
							dist = wtStr*avgDist + wtInt*intDist;
						}
						else
						{
							dist = intDist;
						}
					}
					else
					{
						if (avgDist != INVALIDDISTANCE)
							dist = WT_STRDIST*avgDist + WT_INTDIST*intDist;
						else
							dist = intDist;
					}
				}
				else // use only the intensity-based distance
				{
					dist = intDist;
				}
				// update restraints
				pair<int,int> resPair(res1,res2);
				Contact::ResidueContactType type = c.type;
				// For output purpose we output contacts in one direction only (NCA, NCB, NSC, CACB, CASC, CBSC)
				// if the contact is not in any of these directions, then it is reversed so that it will be
				// For symmetric contacts NN, CACA, CBCB, SCSC, set res1 <= res2, so will only output once
				if (type == Contact::CAN || type == Contact::CBN || type == Contact::SCN || type == Contact::CBCA ||
					type == Contact::SCCA || type == Contact::SCCB)
				{
					resPair = std::make_pair(res2,res1);
					type = c.getReverseType();
				}
				else if (type == Contact::NN || type == Contact::CACA || type == Contact::CBCB || type == Contact::SCSC)
				{
					if (res1 > res2)
						resPair = std::make_pair(res2,res1);
				}
				pair< pair<int,int>, Contact::ResidueContactType > key(resPair,type);
				tr1::unordered_map< pair< pair<int,int>, Contact::ResidueContactType >, pair<double,double> >::iterator itFind = ambigRestraints.find(key);
				if (itFind != ambigRestraints.end())
				{
					pair<double,double>& scoreDist = itFind->second;
					scoreDist.first += score;
					scoreDist.second = max<double>(scoreDist.second,dist);
				}
				else
				{
					pair<double,double> scoreDist(score,dist);
					ambigRestraints[key] = scoreDist;
				}
			} // end for each contact in pc
		} // end for each ambiguous assignment for this noe
	} // end for all ambiguous assignments for all noes
	// output ambiguous assignment (in groups of size at least 3, at most 20, and sum of scores at least 1.0)
	FILE* outA = fopen(filenameAmbig.c_str(), "w");
	if (ambigRestraints.size() > 2 && ((!POOL_FILTER_BY_COUNT && totalPoolScoreSum >= POOL_FILTER) || (POOL_FILTER_BY_COUNT && totalPoolCount >= POOL_FILTER)))
	{
		outputStr = "";
		vector<ContactInfo> firstTwo; // contains first 2 entries in ambigRestraints
		int count = 0;
		double scoreSum = 0;
		bool done = false;
		for (tr1::unordered_map< pair< pair<int,int>, Contact::ResidueContactType >, pair<double,double> >::iterator it = ambigRestraints.begin(); it != ambigRestraints.end(); ) // ++it)
		{
			const pair< pair<int,int>, Contact::ResidueContactType >& key = it->first;
			const pair<int,int>& p = key.first;
			int r1 = p.first;
			int r2 = p.second;
			if (r1 > r2)
				std::swap(r1,r2); // for accessing restraintsResPair only, where r1 <= r2
			pair<int,int> rp(r1,r2);
			const pair<double,double>& scoreDist = it->second;
			Contact::ResidueContactType type = key.second;
			double dist = scoreDist.second+2*CAERROR;
			if (type == Contact::CBCB || type == Contact::CBCA || type == Contact::CBN || type == Contact::CBSC)
			{
				dist += CBERROR;
			}
			if (type == Contact::CBCB || type == Contact::CACB || type == Contact::NCB || type == Contact::SCCB)
			{
				dist += CBERROR;
			}
			if (type == Contact::SCCB || type == Contact::SCCA || type == Contact::SCN || type == Contact::SCSC)
			{
				dist += SCERROR[seqAAType[p.first-1]];
			}
			if (type == Contact::CBSC || type == Contact::CASC || type == Contact::NSC || type == Contact::SCSC)
			{
				dist += SCERROR[seqAAType[p.second-1]];
			}
			if (type == Contact::NN)
			{
				dist += 2*NERROR;
			}
			if (type == Contact::CAN || type == Contact::CBN || type == Contact::SCN)
			{
				dist += NERROR;
			}
			if (type == Contact::NCA || type == Contact::NCB || type == Contact::NSC)
			{
				dist += NERROR;
			}

			string typeStr = Contact::CONTACTTYPESTR[type];
			char buffer[35];
			snprintf(buffer,35,"%-4s %5d %5d %9.6f %6.2f\n",typeStr.c_str(),p.first+resOffset,p.second+resOffset,scoreDist.first,dist);
			outputStr += buffer;
			count++;
			scoreSum += scoreDist.first;

			if ( (!POOL_FILTER_BY_COUNT && (count > 2 && (scoreSum > POOL_FILTER || count > 5))) ||
				 (POOL_FILTER_BY_COUNT && count >= POOL_FILTER) )
			{
				char countStr[8];
				snprintf(countStr,8,"%d",count);
				string countS = countStr;
				outputStr = "AMB   "+countS+"\n"+outputStr+"END\n"; // "AMB   3\n"+outputStr+"END\n";
				fprintf(outA,"%s",outputStr.c_str());
				outputStr = "";
				scoreSum = 0;
				count = 0;
				if (done)
					break;
				++it;
			}
			else  // if (count < 3)
			{
				ContactInfo info;
				info.r1 = r1;
				info.r2 = r2;
				info.type = type;
				info.dist = dist;
				info.score = scoreDist.first;
				firstTwo.push_back(info);
				++it;
				if (it == ambigRestraints.end())
				{
					it = ambigRestraints.begin(); // wrap around
					done = true;
				}
			}
		} // end for each ambigRestraints
	}
	else
	{
		printf("Not enough ambiguous restraints. Count=%zu Score=%f. Skipping output (no need to be concerned)\n",
				ambigRestraints.size(),totalPoolScoreSum);
	}
	fclose(outA);
}

/*
 * Merge by PseudoContact the symmetric peaks into assignedNoDup (so assigned contacts are now undirected). The scores are summed.
 * Caller is responsible for deleting the contents of assignedNoDup
 */
void mergeAssignments(vector<PCAssignment*>& assigned, vector<PCAssignment*>& assignedNoDup)
{
	tr1::unordered_set<PCAssignment*> toSkip; // pca's that won't be added to assignedNoDup because it's pca.pc is same as
	                                          // another pca already added
	int i = 0;
	for (vector<PCAssignment*>::iterator it1 = assigned.begin(); it1 != assigned.end(); ++it1)
	{
		PCAssignment* pca1 = *it1;
		tr1::unordered_set<PCAssignment*>::iterator itFind = toSkip.find(pca1);
		if (itFind != toSkip.end())
			continue; // already added to toSkip
		PCAssignment* pca = new PCAssignment(*pca1);
		int j = 0;
		for (vector<PCAssignment*>::iterator it2 = assigned.begin(); it2 != assigned.end(); ++it2)
		{
			PCAssignment* pca2 = *it2;
			if (i < j) // to avoid duplicate comparisons
			{
				tr1::unordered_set<PCAssignment*>::iterator itFind2 = toSkip.find(pca2);
				if (itFind2 != toSkip.end())
					continue; // pca2.pc already counted
				if (pca1->pc == pca2->pc || pca1->pc.areSymmetric(pca2->pc))
				{
					pca->score.setAdd(pca2->score);
					toSkip.insert(pca2); // mark pca2, so won't add to assignedNoDup
				}
			}
			j++;
		}
		assignedNoDup.push_back(pca);
		i++;
	}
	printf("Merging Symmetric Assignments. Before: %zu, After: %zu\n",assigned.size(),assignedNoDup.size());
}

ScoreCountCondition2* parseScoreCount2(stringstream& tok, bool useAssPos, int numAssignableContactsSR, int numAssignableContactsLR)
{
	int countThresholdSR = 0;
	int numTypesSR = 0;
	string temp;
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseScoreCount2\n");
		exit(-1);
	}
	countThresholdSR = atoi(temp.c_str());
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseScoreCount2\n");
		exit(-1);
	}
	numTypesSR = atoi(temp.c_str());
	int* termsSR = new int[numTypesSR];
	double* boundsSR = new double[numTypesSR];
	list<int> negBoundsSR; // index of terms that have negative bounds
	for (int i = 0; i < numTypesSR; i++)
	{
		if (tok.good())
			tok >> temp;
		else
		{
			printf("INPUT ERROR: parseScoreCount2\n");
			exit(-1);
		}
		termsSR[i] = atoi(temp.c_str());
		if (tok.good())
			tok >> temp;
		else
		{
			printf("INPUT ERROR: parseScoreCount2\n");
			exit(-1);
		}
		boundsSR[i] = atof(temp.c_str());
		if (boundsSR[i] < 0)
		{
			negBoundsSR.push_back(i);
		}
	}
	if (negBoundsSR.size() > 0)
	{
		if (useAssPos) // use assignPossibPtr
		{
			vector< vector<double> > scores;
			tr1::unordered_map<unsigned, vector<double> > scoresFlag;
			scores.reserve(10);
			for (int i = 0; i < 10; i++)
			{
				vector<double> temp;
				temp.reserve(assignPossibPtr->size());
				scores.push_back(temp);
			}
			for (list<PCAssignment*>::iterator it = assignPossibPtr->begin(); it != assignPossibPtr->end(); ++it)
			{
				PCAssignment* pca = *it;
				int seqsep = pca->pc.getSeqSep();
				if (seqsep >= LRSEQSEP)
					continue;
				Score& s = pca->score;
				for (list<int>::iterator itB = negBoundsSR.begin(); itB != negBoundsSR.end(); ++itB)
				{
					int index = *itB;
					unsigned int type = termsSR[index];
					if (type > 10000)
					{
						unsigned flag = type - 10000;
						tr1::unordered_map< unsigned, vector<double> >::iterator itFind = scoresFlag.find(flag);
						if (itFind != scoresFlag.end())
							scoresFlag[flag].push_back(s.getScore(flag));
						else
						{
							vector<double> temp;
							temp.push_back(s.getScore(flag));
							scoresFlag[flag] = temp;
						}
						continue;
					}
					switch (type)
					{
					case 0:
						scores[0].push_back(s.cs);
						break;
					case 1:
						scores[1].push_back(s.str);
						break;
					case 2:
						scores[2].push_back(s.intensity);
						break;
					case 3:
						scores[3].push_back(s.sym);
						break;
					case 4:
						scores[4].push_back(s.interres);
						break;
					case 5:
						scores[5].push_back(s.net);
						break;
					case 6:
						scores[6].push_back(s.netStr);
						break;
					case 7:
						scores[7].push_back(s.ambig);
						break;
					case 8:
						scores[8].push_back(s.db);
						break;
					case 9:
						scores[9].push_back(s.total);
						break;
					default:
						printf("Invalid terms\n");
						exit(-1);
					}
				} // end for each negBounds
			} // end for each assignment
			for (int i = 0; i < 10; i++)
			{
				vector<double>& temp = scores[i];
				if (temp.size())
				{
					sort(temp.begin(),temp.end()); // ascending order
					reverse(temp.begin(),temp.end()); // descending order
				}
			}
			for (tr1::unordered_map< unsigned, vector<double> >::iterator it = scoresFlag.begin(); it != scoresFlag.end(); ++it)
			{
				vector<double>& temp = it->second;
				if (temp.size())
				{
					sort(temp.begin(),temp.end());
					reverse(temp.begin(),temp.end());
				}
			}
			// update bounds
			for (list<int>::iterator itB = negBoundsSR.begin(); itB != negBoundsSR.end(); ++itB)
			{
				int index = *itB;
				unsigned int type = termsSR[index];
				vector<double>& temp = (type < 10000 ? scores[type] : scoresFlag[type-10000]);
				double fract = -boundsSR[index];
				unsigned int top = 0;
				if (fract <= 1.0)
					top = round(fract*temp.size());
				else
				{
					top = round(double((fract-1)*numAssignableContactsSR));
					if (top >= temp.size())
						top = temp.size()-1;
				}
				boundsSR[index] = temp[top];
				printf("ScoreTermCount2 SR: Term: %d   NewBounds: %5.3f\n",type,boundsSR[index]);
			}
		}
		else // use assignedPtr
		{
			vector< vector<double> > scores;
			tr1::unordered_map<unsigned, vector<double> > scoresFlag;
			scores.reserve(10);
			for (int i = 0; i < 10; i++)
			{
				vector<double> temp;
				temp.reserve(assignedPtr->size());
				scores.push_back(temp);
			}
			for (vector<PCAssignment*>::iterator it = assignedPtr->begin(); it != assignedPtr->end(); ++it)
			{
				PCAssignment* pca = *it;
				int seqsep = pca->pc.getSeqSep();
				if (seqsep >= LRSEQSEP)
					continue;
				Score& s = pca->score;
				for (list<int>::iterator itB = negBoundsSR.begin(); itB != negBoundsSR.end(); ++itB)
				{
					int index = *itB;
					unsigned int type = termsSR[index];
					if (type > 10000)
					{
						unsigned flag = type - 10000;
						tr1::unordered_map< unsigned, vector<double> >::iterator itFind = scoresFlag.find(flag);
						if (itFind != scoresFlag.end())
							scoresFlag[flag].push_back(s.getScore(flag));
						else
						{
							vector<double> temp;
							temp.push_back(s.getScore(flag));
							scoresFlag[flag] = temp;
						}
						continue;
					}
					switch (type)
					{
					case 0:
						scores[0].push_back(s.cs);
						break;
					case 1:
						scores[1].push_back(s.str);
						break;
					case 2:
						scores[2].push_back(s.intensity);
						break;
					case 3:
						scores[3].push_back(s.sym);
						break;
					case 4:
						scores[4].push_back(s.interres);
						break;
					case 5:
						scores[5].push_back(s.net);
						break;
					case 6:
						scores[6].push_back(s.netStr);
						break;
					case 7:
						scores[7].push_back(s.ambig);
						break;
					case 8:
						scores[8].push_back(s.db);
						break;
					case 9:
						scores[9].push_back(s.total);
						break;
					default:
						printf("Invalid terms\n");
						exit(-1);
					}
				} // end for each negBounds
			} // end for each assignment
			for (int i = 0; i < 10; i++)
			{
				vector<double>& temp = scores[i];
				if (temp.size())
				{
					sort(temp.begin(),temp.end()); // ascending order
					reverse(temp.begin(),temp.end()); // descending order
				}
			}
			for (tr1::unordered_map< unsigned, vector<double> >::iterator it = scoresFlag.begin(); it != scoresFlag.end(); ++it)
			{
				vector<double>& temp = it->second;
				if (temp.size())
				{
					sort(temp.begin(),temp.end());
					reverse(temp.begin(),temp.end());
				}
			}
			// update bounds
			for (list<int>::iterator itB = negBoundsSR.begin(); itB != negBoundsSR.end(); ++itB)
			{
				int index = *itB;
				unsigned int type = termsSR[index];
				vector<double>& temp = (type < 10000 ? scores[type] : scoresFlag[type-10000]);
				double fract = -boundsSR[index];
				unsigned int top = 0;
				if (fract <= 1.0)
					top = round(fract*temp.size());
				else
				{
					top = round(double((fract-1)*numAssignableContactsSR));
					if (top >= temp.size())
						top = temp.size()-1;
				}
				boundsSR[index] = temp[top];
				printf("ScoreTermCount2 SR: Term: %d   NewBounds: %5.3f\n",type,boundsSR[index]);
			}
		}
	}

	// LR
	int countThresholdLR = 0;
	int numTypesLR = 0;
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseScoreCount2\n");
		exit(-1);
	}
	countThresholdLR = atoi(temp.c_str());
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseScoreCount2\n");
		exit(-1);
	}
	numTypesLR = atoi(temp.c_str());
	int* termsLR = new int[numTypesLR];
	double* boundsLR = new double[numTypesLR];
	list<int> negBoundsLR; // index of terms that have negative bounds
	for (int i = 0; i < numTypesLR; i++)
	{
		if (tok.good())
			tok >> temp;
		else
		{
			printf("INPUT ERROR: parseScoreCount2\n");
			exit(-1);
		}
		termsLR[i] = atoi(temp.c_str());
		if (tok.good())
			tok >> temp;
		else
		{
			printf("INPUT ERROR: parseScoreCount2\n");
			exit(-1);
		}
		boundsLR[i] = atof(temp.c_str());
		if (boundsLR[i] < 0)
		{
			negBoundsLR.push_back(i);
		}
	}
	if (negBoundsLR.size() > 0)
	{
		if (useAssPos) // use assignPossibPtr
		{
			vector< vector<double> > scores;
			tr1::unordered_map<unsigned, vector<double> > scoresFlag;
			scores.reserve(10);
			for (int i = 0; i < 10; i++)
			{
				vector<double> temp;
				temp.reserve(assignPossibPtr->size());
				scores.push_back(temp);
			}
			for (list<PCAssignment*>::iterator it = assignPossibPtr->begin(); it != assignPossibPtr->end(); ++it)
			{
				PCAssignment* pca = *it;
				int seqsep = pca->pc.getSeqSep();
				if (seqsep < LRSEQSEP)
					continue;
				Score& s = pca->score;
				for (list<int>::iterator itB = negBoundsLR.begin(); itB != negBoundsLR.end(); ++itB)
				{
					int index = *itB;
					unsigned int type = termsLR[index];
					if (type > 10000)
					{
						unsigned flag = type - 10000;
						tr1::unordered_map< unsigned, vector<double> >::iterator itFind = scoresFlag.find(flag);
						if (itFind != scoresFlag.end())
							scoresFlag[flag].push_back(s.getScore(flag));
						else
						{
							vector<double> temp;
							temp.push_back(s.getScore(flag));
							scoresFlag[flag] = temp;
						}
						continue;
					}
					switch (type)
					{
					case 0:
						scores[0].push_back(s.cs);
						break;
					case 1:
						scores[1].push_back(s.str);
						break;
					case 2:
						scores[2].push_back(s.intensity);
						break;
					case 3:
						scores[3].push_back(s.sym);
						break;
					case 4:
						scores[4].push_back(s.interres);
						break;
					case 5:
						scores[5].push_back(s.net);
						break;
					case 6:
						scores[6].push_back(s.netStr);
						break;
					case 7:
						scores[7].push_back(s.ambig);
						break;
					case 8:
						scores[8].push_back(s.db);
						break;
					case 9:
						scores[9].push_back(s.total);
						break;
					default:
						printf("Invalid terms\n");
						exit(-1);
					}
				} // end for each negBounds
			} // end for each assignment
			for (int i = 0; i < 10; i++)
			{
				vector<double>& temp = scores[i];
				if (temp.size())
				{
					sort(temp.begin(),temp.end()); // ascending order
					reverse(temp.begin(),temp.end()); // descending order
				}
			}
			for (tr1::unordered_map< unsigned, vector<double> >::iterator it = scoresFlag.begin(); it != scoresFlag.end(); ++it)
			{
				vector<double>& temp = it->second;
				if (temp.size())
				{
					sort(temp.begin(),temp.end());
					reverse(temp.begin(),temp.end());
				}
			}
			// update bounds
			for (list<int>::iterator itB = negBoundsLR.begin(); itB != negBoundsLR.end(); ++itB)
			{
				int index = *itB;
				unsigned int type = termsLR[index];
				vector<double>& temp = (type < 10000 ? scores[type] : scoresFlag[type-10000]);
				double fract = -boundsLR[index];
				unsigned int top = 0;
				if (fract <= 1.0)
					top = round(fract*temp.size());
				else
				{
					top = round(double((fract-1)*numAssignableContactsLR));
					if (top >= temp.size())
						top = temp.size()-1;
				}
				boundsLR[index] = temp[top];
				printf("ScoreTermCount2 LR: Term: %d   NewBounds: %5.3f\n",type,boundsLR[index]);
			}
		}
		else // use assignedPtr
		{
			vector< vector<double> > scores;
			tr1::unordered_map<unsigned, vector<double> > scoresFlag;
			scores.reserve(10);
			for (int i = 0; i < 10; i++)
			{
				vector<double> temp;
				temp.reserve(assignedPtr->size());
				scores.push_back(temp);
			}
			for (vector<PCAssignment*>::iterator it = assignedPtr->begin(); it != assignedPtr->end(); ++it)
			{
				PCAssignment* pca = *it;
				int seqsep = pca->pc.getSeqSep();
				if (seqsep < LRSEQSEP)
					continue;
				Score& s = pca->score;
				for (list<int>::iterator itB = negBoundsLR.begin(); itB != negBoundsLR.end(); ++itB)
				{
					int index = *itB;
					unsigned int type = termsLR[index];
					if (type > 10000)
					{
						unsigned flag = type - 10000;
						tr1::unordered_map< unsigned, vector<double> >::iterator itFind = scoresFlag.find(flag);
						if (itFind != scoresFlag.end())
							scoresFlag[flag].push_back(s.getScore(flag));
						else
						{
							vector<double> temp;
							temp.push_back(s.getScore(flag));
							scoresFlag[flag] = temp;
						}
						continue;
					}
					switch (type)
					{
					case 0:
						scores[0].push_back(s.cs);
						break;
					case 1:
						scores[1].push_back(s.str);
						break;
					case 2:
						scores[2].push_back(s.intensity);
						break;
					case 3:
						scores[3].push_back(s.sym);
						break;
					case 4:
						scores[4].push_back(s.interres);
						break;
					case 5:
						scores[5].push_back(s.net);
						break;
					case 6:
						scores[6].push_back(s.netStr);
						break;
					case 7:
						scores[7].push_back(s.ambig);
						break;
					case 8:
						scores[8].push_back(s.db);
						break;
					case 9:
						scores[9].push_back(s.total);
						break;
					default:
						printf("Invalid terms\n");
						exit(-1);
					}
				} // end for each negBounds
			} // end for each assignment
			for (int i = 0; i < 10; i++)
			{
				vector<double>& temp = scores[i];
				if (temp.size())
				{
					sort(temp.begin(),temp.end()); // ascending order
					reverse(temp.begin(),temp.end()); // descending order
				}
			}
			for (tr1::unordered_map< unsigned, vector<double> >::iterator it = scoresFlag.begin(); it != scoresFlag.end(); ++it)
			{
				vector<double>& temp = it->second;
				if (temp.size())
				{
					sort(temp.begin(),temp.end());
					reverse(temp.begin(),temp.end());
				}
			}
			// update bounds
			for (list<int>::iterator itB = negBoundsLR.begin(); itB != negBoundsLR.end(); ++itB)
			{
				int index = *itB;
				unsigned int type = termsLR[index];
				vector<double>& temp = (type < 10000 ? scores[type] : scoresFlag[type-10000]);
				double fract = -boundsLR[index];
				unsigned int top = 0;
				if (fract <= 1.0)
					top = round(fract*temp.size());
				else
				{
					top = round(double((fract-1)*numAssignableContactsLR));
					if (top >= temp.size())
						top = temp.size()-1;
				}
				boundsLR[index] = temp[top];
				printf("ScoreTermCount2: Term: %d   NewBounds: %5.3f\n",type,boundsLR[index]);
			}
		}
	}

	ScoreCountCondition2* scc = new ScoreCountCondition2(termsSR,numTypesSR,boundsSR,countThresholdSR,termsLR,numTypesLR,boundsLR,countThresholdLR); // constructor copies the arrays
	delete [] termsSR;
	delete [] boundsSR;
	delete [] termsLR;
	delete [] boundsLR;
	return scc;
}

ScoreCountCondition1* parseScoreCount1(stringstream& tok, bool useAssPos, int numAssignableContacts)
{
	int countThreshold = 0;
	int numTypes = 0;
	string temp;
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseScoreCount1\n");
		exit(-1);
	}
	countThreshold = atoi(temp.c_str());
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseScoreCount1\n");
		exit(-1);
	}
	numTypes = atoi(temp.c_str());
	int* terms = new int[numTypes];
	double* bounds = new double[numTypes];
	list<int> negBounds; // index of terms that have negative bounds
	for (int i = 0; i < numTypes; i++)
	{
		if (tok.good())
			tok >> temp;
		else
		{
			printf("INPUT ERROR: parseScoreCount1\n");
			exit(-1);
		}
		terms[i] = atoi(temp.c_str());
		if (tok.good())
			tok >> temp;
		else
		{
			printf("INPUT ERROR: parseScoreCount1\n");
			exit(-1);
		}
		bounds[i] = atof(temp.c_str());
		if (bounds[i] < 0)
		{
			negBounds.push_back(i);
		}
	}
	if (negBounds.size() > 0)
	{
		if (useAssPos) // use assignPossibPtr
		{
			vector< vector<double> > scores;
			tr1::unordered_map<unsigned, vector<double> > scoresFlag;
			scores.reserve(10);
			for (int i = 0; i < 10; i++)
			{
				vector<double> temp;
				temp.reserve(assignPossibPtr->size());
				scores.push_back(temp);
			}
			for (list<PCAssignment*>::iterator it = assignPossibPtr->begin(); it != assignPossibPtr->end(); ++it)
			{
				PCAssignment* pca = *it;
				Score& s = pca->score;
				for (list<int>::iterator itB = negBounds.begin(); itB != negBounds.end(); ++itB)
				{
					int index = *itB;
					unsigned int type = terms[index];
					if (type > 10000)
					{
						unsigned flag = type - 10000;
						tr1::unordered_map< unsigned, vector<double> >::iterator itFind = scoresFlag.find(flag);
						if (itFind != scoresFlag.end())
							scoresFlag[flag].push_back(s.getScore(flag));
						else
						{
							vector<double> temp;
							temp.push_back(s.getScore(flag));
							scoresFlag[flag] = temp;
						}
						continue;
					}
					switch (type)
					{
					case 0:
						scores[0].push_back(s.cs);
						break;
					case 1:
						scores[1].push_back(s.str);
						break;
					case 2:
						scores[2].push_back(s.intensity);
						break;
					case 3:
						scores[3].push_back(s.sym);
						break;
					case 4:
						scores[4].push_back(s.interres);
						break;
					case 5:
						scores[5].push_back(s.net);
						break;
					case 6:
						scores[6].push_back(s.netStr);
						break;
					case 7:
						scores[7].push_back(s.ambig);
						break;
					case 8:
						scores[8].push_back(s.db);
						break;
					case 9:
						scores[9].push_back(s.total);
						break;
					default:
						printf("Invalid terms\n");
						exit(-1);
					}
				} // end for each negBounds
			} // end for each assignment
			for (int i = 0; i < 10; i++)
			{
				vector<double>& temp = scores[i];
				if (temp.size())
				{
					sort(temp.begin(),temp.end()); // ascending order
					reverse(temp.begin(),temp.end()); // descending order
				}
			}
			for (tr1::unordered_map< unsigned, vector<double> >::iterator it = scoresFlag.begin(); it != scoresFlag.end(); ++it)
			{
				vector<double>& temp = it->second;
				if (temp.size())
				{
					sort(temp.begin(),temp.end());
					reverse(temp.begin(),temp.end());
				}
			}
			// update bounds
			for (list<int>::iterator itB = negBounds.begin(); itB != negBounds.end(); ++itB)
			{
				int index = *itB;
				unsigned int type = terms[index];
				vector<double>& temp = (type < 10000 ? scores[type] : scoresFlag[type-10000]);
				double fract = -bounds[index];
				unsigned int top = 0;
				if (fract <= 1.0)
					top = round(fract*temp.size());
				else
				{
					top = round(double((fract-1)*numAssignableContacts));
					if (top >= temp.size())
						top = temp.size()-1;
				}
				bounds[index] = temp[top];
				printf("ScoreTermCount1: Term: %d   NewBounds: %5.3f\n",type,bounds[index]);
			}
		}
		else // use assignedPtr
		{
			vector< vector<double> > scores;
			tr1::unordered_map<unsigned, vector<double> > scoresFlag;
			scores.reserve(10);
			for (int i = 0; i < 10; i++)
			{
				vector<double> temp;
				temp.reserve(assignedPtr->size());
				scores.push_back(temp);
			}
			for (vector<PCAssignment*>::iterator it = assignedPtr->begin(); it != assignedPtr->end(); ++it)
			{
				PCAssignment* pca = *it;
				Score& s = pca->score;
				for (list<int>::iterator itB = negBounds.begin(); itB != negBounds.end(); ++itB)
				{
					int index = *itB;
					unsigned int type = terms[index];
					if (type > 10000)
					{
						unsigned flag = type - 10000;
						tr1::unordered_map< unsigned, vector<double> >::iterator itFind = scoresFlag.find(flag);
						if (itFind != scoresFlag.end())
							scoresFlag[flag].push_back(s.getScore(flag));
						else
						{
							vector<double> temp;
							temp.push_back(s.getScore(flag));
							scoresFlag[flag] = temp;
						}
						continue;
					}
					switch (type)
					{
					case 0:
						scores[0].push_back(s.cs);
						break;
					case 1:
						scores[1].push_back(s.str);
						break;
					case 2:
						scores[2].push_back(s.intensity);
						break;
					case 3:
						scores[3].push_back(s.sym);
						break;
					case 4:
						scores[4].push_back(s.interres);
						break;
					case 5:
						scores[5].push_back(s.net);
						break;
					case 6:
						scores[6].push_back(s.netStr);
						break;
					case 7:
						scores[7].push_back(s.ambig);
						break;
					case 8:
						scores[8].push_back(s.db);
						break;
					case 9:
						scores[9].push_back(s.total);
						break;
					default:
						printf("Invalid terms\n");
						exit(-1);
					}
				} // end for each negBounds
			} // end for each assignment
			for (int i = 0; i < 10; i++)
			{
				vector<double>& temp = scores[i];
				if (temp.size())
				{
					sort(temp.begin(),temp.end()); // ascending order
					reverse(temp.begin(),temp.end()); // descending order
				}
			}
			for (tr1::unordered_map< unsigned, vector<double> >::iterator it = scoresFlag.begin(); it != scoresFlag.end(); ++it)
			{
				vector<double>& temp = it->second;
				if (temp.size())
				{
					sort(temp.begin(),temp.end());
					reverse(temp.begin(),temp.end());
				}
			}
			// update bounds
			for (list<int>::iterator itB = negBounds.begin(); itB != negBounds.end(); ++itB)
			{
				int index = *itB;
				unsigned int type = terms[index];
				vector<double>& temp = (type < 10000 ? scores[type] : scoresFlag[type-10000]);
				double fract = -bounds[index];
				unsigned int top = 0;
				if (fract <= 1.0)
					top = round(fract*temp.size());
				else
				{
					top = round(double((fract-1)*numAssignableContacts));
					if (top >= temp.size())
						top = temp.size()-1;
				}
				bounds[index] = temp[top];
				printf("ScoreTermCount1: Term: %d   NewBounds: %5.3f\n",type,bounds[index]);
			}
		}
	}
	ScoreCountCondition1* scc = new ScoreCountCondition1(terms,numTypes,bounds, countThreshold); // constructor copies the arrays
	delete [] terms;
	delete [] bounds;
	return scc;
}

ScoreTermCondition2* parseScoreTerm2(stringstream& tok, bool useAssPos, int numAssignableContactsSR, int numAssignableContactsLR)
{
	int numTypesSR = 0;
	int numTypesLR = 0;
	string temp;
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseScoreTerm2\n");
		exit(-1);
	}
	numTypesSR = atoi(temp.c_str());
	int* termsSR = new int[numTypesSR];
	double* boundsSR = new double[numTypesSR];
	list<int> negBoundsSR; // index of terms that have negative bounds
	for (int i = 0; i < numTypesSR; i++)
	{
		if (tok.good())
			tok >> temp;
		else
		{
			printf("INPUT ERROR: parseScoreTerm2\n");
			exit(-1);
		}
		termsSR[i] = atoi(temp.c_str());
		if (tok.good())
			tok >> temp;
		else
		{
			printf("INPUT ERROR: parseScoreTerm2\n");
			exit(-1);
		}
		boundsSR[i] = atof(temp.c_str());
		if (boundsSR[i] < 0)
		{
			negBoundsSR.push_back(i);
		}
	}
	if (negBoundsSR.size() > 0)
	{
		if (useAssPos) // use assignPossibPtr
		{
			vector< vector<double> > scores;
			tr1::unordered_map<unsigned, vector<double> > scoresFlag; // stores the scores for the flag terms only
			                                                          // index by type-10000
			scores.reserve(10);
			for (int i = 0; i < 10; i++)
			{
				vector<double> temp;
				temp.reserve(assignPossibPtr->size());
				scores.push_back(temp);
			}
			for (list<PCAssignment*>::iterator it = assignPossibPtr->begin(); it != assignPossibPtr->end(); ++it)
			{
				PCAssignment* pca = *it;
				Score& s = pca->score;
				int seqsep = pca->pc.getSeqSep();
				if (seqsep >= LRSEQSEP)
					continue;
				for (list<int>::iterator itB = negBoundsSR.begin(); itB != negBoundsSR.end(); ++itB)
				{
					int index = *itB;
					unsigned int type = termsSR[index];
					if (type > 10000)
					{
						unsigned flag = type - 10000;
						tr1::unordered_map< unsigned, vector<double> >::iterator itFind = scoresFlag.find(flag);
						if (itFind != scoresFlag.end())
							scoresFlag[flag].push_back(s.getScore(flag));
						else
						{
							vector<double> temp;
							temp.push_back(s.getScore(flag));
							scoresFlag[flag] = temp;
						}
						continue;
					}
					switch (type)
					{
					case 0:
						scores[0].push_back(s.cs);
						break;
					case 1:
						scores[1].push_back(s.str);
						break;
					case 2:
						scores[2].push_back(s.intensity);
						break;
					case 3:
						scores[3].push_back(s.sym);
						break;
					case 4:
						scores[4].push_back(s.interres);
						break;
					case 5:
						scores[5].push_back(s.net);
						break;
					case 6:
						scores[6].push_back(s.netStr);
						break;
					case 7:
						scores[7].push_back(s.ambig);
						break;
					case 8:
						scores[8].push_back(s.db);
						break;
					case 9:
						scores[9].push_back(s.total);
						break;
					default:
						printf("Invalid terms\n");
						exit(-1);
					}
				} // end for each negBounds
			} // end for each assignment
			for (int i = 0; i < 10; i++)
			{
				vector<double>& temp = scores[i];
				if (temp.size())
				{
					sort(temp.begin(),temp.end()); // ascending order
					reverse(temp.begin(),temp.end()); // descending order
				}
			}
			for (tr1::unordered_map< unsigned, vector<double> >::iterator it = scoresFlag.begin(); it != scoresFlag.end(); ++it)
			{
				vector<double>& temp = it->second;
				if (temp.size())
				{
					sort(temp.begin(),temp.end());
					reverse(temp.begin(),temp.end());
				}
			}
			// update bounds
			for (list<int>::iterator itB = negBoundsSR.begin(); itB != negBoundsSR.end(); ++itB)
			{
				int index = *itB;
				unsigned int type = termsSR[index];
				vector<double>& temp = (type < 10000 ? scores[type] : scoresFlag[type-10000]);
				double fract = -boundsSR[index];
				unsigned int top = 0;
				if (fract <= 1.0)
					top = round(fract*temp.size());
				else
				{
					top = round(double((fract-1)*numAssignableContactsSR));
					if (top >= temp.size())
						top = temp.size()-1;
				}
				boundsSR[index] = temp[top];
				printf("ScoreTermCount2 SR: Term: %d   NewBounds: %5.3f\n",type,boundsSR[index]);
			}
		}
		else // use assignedPtr
		{
			vector< vector<double> > scores;
			tr1::unordered_map<unsigned, vector<double> > scoresFlag;
			scores.reserve(10);
			for (int i = 0; i < 10; i++)
			{
				vector<double> temp;
				temp.reserve(assignedPtr->size());
				scores.push_back(temp);
			}
			for (vector<PCAssignment*>::iterator it = assignedPtr->begin(); it != assignedPtr->end(); ++it)
			{
				PCAssignment* pca = *it;
				Score& s = pca->score;
				int seqsep = pca->pc.getSeqSep();
				if (seqsep >= LRSEQSEP)
					continue;
				for (list<int>::iterator itB = negBoundsSR.begin(); itB != negBoundsSR.end(); ++itB)
				{
					int index = *itB;
					unsigned int type = termsSR[index];
					if (type > 10000)
					{
						unsigned flag = type - 10000;  ;
						tr1::unordered_map< unsigned, vector<double> >::iterator itFind = scoresFlag.find(flag);
						if (itFind != scoresFlag.end())
							scoresFlag[flag].push_back(s.getScore(flag));
						else
						{
							vector<double> temp;
							temp.push_back(s.getScore(flag));
							scoresFlag[flag] = temp;
						}
						continue;
					}
					switch (type)
					{
					case 0:
						scores[0].push_back(s.cs);
						break;
					case 1:
						scores[1].push_back(s.str);
						break;
					case 2:
						scores[2].push_back(s.intensity);
						break;
					case 3:
						scores[3].push_back(s.sym);
						break;
					case 4:
						scores[4].push_back(s.interres);
						break;
					case 5:
						scores[5].push_back(s.net);
						break;
					case 6:
						scores[6].push_back(s.netStr);
						break;
					case 7:
						scores[7].push_back(s.ambig);
						break;
					case 8:
						scores[8].push_back(s.db);
						break;
					case 9:
						scores[9].push_back(s.total);
						break;
					default:
						printf("Invalid terms\n");
						exit(-1);
					}
				} // end for each negBounds
			} // end for each assignment
			for (int i = 0; i < 10; i++)
			{
				vector<double>& temp = scores[i];
				if (temp.size())
				{
					sort(temp.begin(),temp.end()); // ascending order
					reverse(temp.begin(),temp.end()); // descending order
				}
			}
			for (tr1::unordered_map< unsigned, vector<double> >::iterator it = scoresFlag.begin(); it != scoresFlag.end(); ++it)
			{
				vector<double>& temp = it->second;
				if (temp.size())
				{
					sort(temp.begin(),temp.end());
					reverse(temp.begin(),temp.end());
				}
			}
			// update bounds
			for (list<int>::iterator itB = negBoundsSR.begin(); itB != negBoundsSR.end(); ++itB)
			{
				int index = *itB;
				unsigned int type = termsSR[index];
				vector<double>& temp = (type < 10000 ? scores[type] : scoresFlag[type-10000]);
				double fract = -boundsSR[index];
				unsigned int top = 0;
				if (fract <= 1.0)
					top = round(fract*temp.size());
				else
				{
					top = round(double((fract-1)*numAssignableContactsSR));
					if (top >= temp.size())
						top = temp.size()-1;
				}
				boundsSR[index] = temp[top];
				printf("ScoreTermCount2 SR: Term: %d   NewBounds: %5.3f\n",type,boundsSR[index]);
			}
		}
	}

	// parse LR bounds
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseScoreTerm2\n");
		exit(-1);
	}
	numTypesLR = atoi(temp.c_str());
	int* termsLR = new int[numTypesLR];
	double* boundsLR = new double[numTypesLR];
	list<int> negBoundsLR; // index of terms that have negative bounds

	for (int i = 0; i < numTypesLR; i++)
	{
		if (tok.good())
			tok >> temp;
		else
		{
			printf("INPUT ERROR: parseScoreTerm2\n");
			exit(-1);
		}
		termsLR[i] = atoi(temp.c_str());
		if (tok.good())
			tok >> temp;
		else
		{
			printf("INPUT ERROR: parseScoreTerm2\n");
			exit(-1);
		}
		boundsLR[i] = atof(temp.c_str());
		if (boundsLR[i] < 0)
		{
			negBoundsLR.push_back(i);
		}
	}
	if (negBoundsLR.size() > 0)
	{
		if (useAssPos) // use assignPossibPtr
		{
			vector< vector<double> > scores;
			tr1::unordered_map<unsigned, vector<double> > scoresFlag; // stores the scores for the flag terms only
			                                                          // index by type-10000
			scores.reserve(10);
			for (int i = 0; i < 10; i++)
			{
				vector<double> temp;
				temp.reserve(assignPossibPtr->size());
				scores.push_back(temp);
			}
			for (list<PCAssignment*>::iterator it = assignPossibPtr->begin(); it != assignPossibPtr->end(); ++it)
			{
				PCAssignment* pca = *it;
				Score& s = pca->score;
				int seqsep = pca->pc.getSeqSep();
				if (seqsep < LRSEQSEP)
					continue;
				for (list<int>::iterator itB = negBoundsLR.begin(); itB != negBoundsLR.end(); ++itB)
				{
					int index = *itB;
					unsigned int type = termsLR[index];
					if (type > 10000)
					{
						unsigned flag = type - 10000;
						tr1::unordered_map< unsigned, vector<double> >::iterator itFind = scoresFlag.find(flag);
						if (itFind != scoresFlag.end())
							scoresFlag[flag].push_back(s.getScore(flag));
						else
						{
							vector<double> temp;
							temp.push_back(s.getScore(flag));
							scoresFlag[flag] = temp;
						}
						continue;
					}
					switch (type)
					{
					case 0:
						scores[0].push_back(s.cs);
						break;
					case 1:
						scores[1].push_back(s.str);
						break;
					case 2:
						scores[2].push_back(s.intensity);
						break;
					case 3:
						scores[3].push_back(s.sym);
						break;
					case 4:
						scores[4].push_back(s.interres);
						break;
					case 5:
						scores[5].push_back(s.net);
						break;
					case 6:
						scores[6].push_back(s.netStr);
						break;
					case 7:
						scores[7].push_back(s.ambig);
						break;
					case 8:
						scores[8].push_back(s.db);
						break;
					case 9:
						scores[9].push_back(s.total);
						break;
					default:
						printf("Invalid terms\n");
						exit(-1);
					}
				} // end for each negBounds
			} // end for each assignment
			for (int i = 0; i < 10; i++)
			{
				vector<double>& temp = scores[i];
				if (temp.size())
				{
					sort(temp.begin(),temp.end()); // ascending order
					reverse(temp.begin(),temp.end()); // descending order
				}
			}
			for (tr1::unordered_map< unsigned, vector<double> >::iterator it = scoresFlag.begin(); it != scoresFlag.end(); ++it)
			{
				vector<double>& temp = it->second;
				if (temp.size())
				{
					sort(temp.begin(),temp.end());
					reverse(temp.begin(),temp.end());
				}
			}
			// update bounds
			for (list<int>::iterator itB = negBoundsLR.begin(); itB != negBoundsLR.end(); ++itB)
			{
				int index = *itB;
				unsigned int type = termsLR[index];
				vector<double>& temp = (type < 10000 ? scores[type] : scoresFlag[type-10000]);
				double fract = -boundsLR[index];
				unsigned int top = 0;
				if (fract <= 1.0)
					top = round(fract*temp.size());
				else
				{
					top = round(double((fract-1)*numAssignableContactsLR));
					if (top >= temp.size())
						top = temp.size()-1;
				}
				boundsLR[index] = temp[top];
				printf("ScoreTermCount2 LR: Term: %d   NewBounds: %5.3f\n",type,boundsLR[index]);
			}
		}
		else // use assignedPtr
		{
			vector< vector<double> > scores;
			tr1::unordered_map<unsigned, vector<double> > scoresFlag;
			scores.reserve(10);
			for (int i = 0; i < 10; i++)
			{
				vector<double> temp;
				temp.reserve(assignedPtr->size());
				scores.push_back(temp);
			}
			for (vector<PCAssignment*>::iterator it = assignedPtr->begin(); it != assignedPtr->end(); ++it)
			{
				PCAssignment* pca = *it;
				Score& s = pca->score;
				int seqsep = pca->pc.getSeqSep();
				if (seqsep < LRSEQSEP)
					continue;
				for (list<int>::iterator itB = negBoundsLR.begin(); itB != negBoundsLR.end(); ++itB)
				{
					int index = *itB;
					unsigned int type = termsLR[index];
					if (type > 10000)
					{
						unsigned flag = type - 10000;  ;
						tr1::unordered_map< unsigned, vector<double> >::iterator itFind = scoresFlag.find(flag);
						if (itFind != scoresFlag.end())
							scoresFlag[flag].push_back(s.getScore(flag));
						else
						{
							vector<double> temp;
							temp.push_back(s.getScore(flag));
							scoresFlag[flag] = temp;
						}
						continue;
					}
					switch (type)
					{
					case 0:
						scores[0].push_back(s.cs);
						break;
					case 1:
						scores[1].push_back(s.str);
						break;
					case 2:
						scores[2].push_back(s.intensity);
						break;
					case 3:
						scores[3].push_back(s.sym);
						break;
					case 4:
						scores[4].push_back(s.interres);
						break;
					case 5:
						scores[5].push_back(s.net);
						break;
					case 6:
						scores[6].push_back(s.netStr);
						break;
					case 7:
						scores[7].push_back(s.ambig);
						break;
					case 8:
						scores[8].push_back(s.db);
						break;
					case 9:
						scores[9].push_back(s.total);
						break;
					default:
						printf("Invalid terms\n");
						exit(-1);
					}
				} // end for each negBounds
			} // end for each assignment
			for (int i = 0; i < 10; i++)
			{
				vector<double>& temp = scores[i];
				if (temp.size())
				{
					sort(temp.begin(),temp.end()); // ascending order
					reverse(temp.begin(),temp.end()); // descending order
				}
			}
			for (tr1::unordered_map< unsigned, vector<double> >::iterator it = scoresFlag.begin(); it != scoresFlag.end(); ++it)
			{
				vector<double>& temp = it->second;
				if (temp.size())
				{
					sort(temp.begin(),temp.end());
					reverse(temp.begin(),temp.end());
				}
			}
			// update bounds
			for (list<int>::iterator itB = negBoundsLR.begin(); itB != negBoundsLR.end(); ++itB)
			{
				int index = *itB;
				unsigned int type = termsLR[index];
				vector<double>& temp = (type < 10000 ? scores[type] : scoresFlag[type-10000]);
				double fract = -boundsLR[index];
				unsigned int top = 0;
				if (fract <= 1.0)
					top = round(fract*temp.size());
				else
				{
					top = round(double((fract-1)*numAssignableContactsLR));
					if (top >= temp.size())
						top = temp.size()-1;
				}
				boundsLR[index] = temp[top];
				printf("ScoreTermCount2 LR: Term: %d   NewBounds: %5.3f\n",type,boundsLR[index]);
			}
		}
	}

	ScoreTermCondition2* stc = new ScoreTermCondition2(termsSR,numTypesSR,boundsSR, termsLR, numTypesLR, boundsLR); // constructor copies the arrays
	delete [] termsSR;
	delete [] boundsSR;
	delete [] termsLR;
	delete [] boundsLR;
	return stc;
}

ScoreTermCondition1* parseScoreTerm1(stringstream& tok, bool useAssPos, int numAssignableContacts)
{
	int numTypes = 0;
	string temp;
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseScoreTerm1\n");
		exit(-1);
	}
	numTypes = atoi(temp.c_str());
	int* terms = new int[numTypes];
	double* bounds = new double[numTypes];
	list<int> negBounds; // index of terms that have negative bounds
	for (int i = 0; i < numTypes; i++)
	{
		if (tok.good())
			tok >> temp;
		else
		{
			printf("INPUT ERROR: parseScoreTerm1\n");
			exit(-1);
		}
		terms[i] = atoi(temp.c_str());
		if (tok.good())
			tok >> temp;
		else
		{
			printf("INPUT ERROR: parseScoreTerm1\n");
			exit(-1);
		}
		bounds[i] = atof(temp.c_str());
		if (bounds[i] < 0)
		{
			negBounds.push_back(i);
		}
	}
	if (negBounds.size() > 0)
	{
		if (useAssPos) // use assignPossibPtr
		{
			vector< vector<double> > scores;
			tr1::unordered_map<unsigned, vector<double> > scoresFlag; // stores the scores for the flag terms only
			                                                          // index by type-10000
			scores.reserve(10);
			for (int i = 0; i < 10; i++)
			{
				vector<double> temp;
				temp.reserve(assignPossibPtr->size());
				scores.push_back(temp);
			}
			for (list<PCAssignment*>::iterator it = assignPossibPtr->begin(); it != assignPossibPtr->end(); ++it)
			{
				PCAssignment* pca = *it;
				Score& s = pca->score;
				for (list<int>::iterator itB = negBounds.begin(); itB != negBounds.end(); ++itB)
				{
					int index = *itB;
					unsigned int type = terms[index];
					if (type > 10000)
					{
						unsigned flag = type - 10000;
						tr1::unordered_map< unsigned, vector<double> >::iterator itFind = scoresFlag.find(flag);
						if (itFind != scoresFlag.end())
							scoresFlag[flag].push_back(s.getScore(flag));
						else
						{
							vector<double> temp;
							temp.push_back(s.getScore(flag));
							scoresFlag[flag] = temp;
						}
						continue;
					}
					switch (type)
					{
					case 0:
						scores[0].push_back(s.cs);
						break;
					case 1:
						scores[1].push_back(s.str);
						break;
					case 2:
						scores[2].push_back(s.intensity);
						break;
					case 3:
						scores[3].push_back(s.sym);
						break;
					case 4:
						scores[4].push_back(s.interres);
						break;
					case 5:
						scores[5].push_back(s.net);
						break;
					case 6:
						scores[6].push_back(s.netStr);
						break;
					case 7:
						scores[7].push_back(s.ambig);
						break;
					case 8:
						scores[8].push_back(s.db);
						break;
					case 9:
						scores[9].push_back(s.total);
						break;
					default:
						printf("Invalid terms\n");
						exit(-1);
					}
				} // end for each negBounds
			} // end for each assignment
			for (int i = 0; i < 10; i++)
			{
				vector<double>& temp = scores[i];
				if (temp.size())
				{
					sort(temp.begin(),temp.end()); // ascending order
					reverse(temp.begin(),temp.end()); // descending order
				}
			}
			for (tr1::unordered_map< unsigned, vector<double> >::iterator it = scoresFlag.begin(); it != scoresFlag.end(); ++it)
			{
				vector<double>& temp = it->second;
				if (temp.size())
				{
					sort(temp.begin(),temp.end());
					reverse(temp.begin(),temp.end());
				}
			}
			// update bounds
			for (list<int>::iterator itB = negBounds.begin(); itB != negBounds.end(); ++itB)
			{
				int index = *itB;
				unsigned int type = terms[index];
				vector<double>& temp = (type < 10000 ? scores[type] : scoresFlag[type-10000]);
				double fract = -bounds[index];
				unsigned int top = 0;
				if (fract <= 1.0)
					top = round(fract*temp.size());
				else
				{
					top = round(double((fract-1)*numAssignableContacts));
					if (top >= temp.size())
						top = temp.size()-1;
				}
				bounds[index] = temp[top];
				printf("ScoreTermCount1: Term: %d   NewBounds: %5.3f\n",type,bounds[index]);
			}
		}
		else // use assignedPtr
		{
			vector< vector<double> > scores;
			tr1::unordered_map<unsigned, vector<double> > scoresFlag;
			scores.reserve(10);
			for (int i = 0; i < 10; i++)
			{
				vector<double> temp;
				temp.reserve(assignedPtr->size());
				scores.push_back(temp);
			}
			for (vector<PCAssignment*>::iterator it = assignedPtr->begin(); it != assignedPtr->end(); ++it)
			{
				PCAssignment* pca = *it;
				Score& s = pca->score;
				for (list<int>::iterator itB = negBounds.begin(); itB != negBounds.end(); ++itB)
				{
					int index = *itB;
					unsigned int type = terms[index];
					if (type > 10000)
					{
						unsigned flag = type - 10000;  ;
						tr1::unordered_map< unsigned, vector<double> >::iterator itFind = scoresFlag.find(flag);
						if (itFind != scoresFlag.end())
							scoresFlag[flag].push_back(s.getScore(flag));
						else
						{
							vector<double> temp;
							temp.push_back(s.getScore(flag));
							scoresFlag[flag] = temp;
						}
						continue;
					}
					switch (type)
					{
					case 0:
						scores[0].push_back(s.cs);
						break;
					case 1:
						scores[1].push_back(s.str);
						break;
					case 2:
						scores[2].push_back(s.intensity);
						break;
					case 3:
						scores[3].push_back(s.sym);
						break;
					case 4:
						scores[4].push_back(s.interres);
						break;
					case 5:
						scores[5].push_back(s.net);
						break;
					case 6:
						scores[6].push_back(s.netStr);
						break;
					case 7:
						scores[7].push_back(s.ambig);
						break;
					case 8:
						scores[8].push_back(s.db);
						break;
					case 9:
						scores[9].push_back(s.total);
						break;
					default:
						printf("Invalid terms\n");
						exit(-1);
					}
				} // end for each negBounds
			} // end for each assignment
			for (int i = 0; i < 10; i++)
			{
				vector<double>& temp = scores[i];
				if (temp.size())
				{
					sort(temp.begin(),temp.end()); // ascending order
					reverse(temp.begin(),temp.end()); // descending order
				}
			}
			for (tr1::unordered_map< unsigned, vector<double> >::iterator it = scoresFlag.begin(); it != scoresFlag.end(); ++it)
			{
				vector<double>& temp = it->second;
				if (temp.size())
				{
					sort(temp.begin(),temp.end());
					reverse(temp.begin(),temp.end());
				}
			}
			// update bounds
			for (list<int>::iterator itB = negBounds.begin(); itB != negBounds.end(); ++itB)
			{
				int index = *itB;
				unsigned int type = terms[index];
				vector<double>& temp = (type < 10000 ? scores[type] : scoresFlag[type-10000]);
				double fract = -bounds[index];
				unsigned int top = 0;
				if (fract <= 1.0)
					top = round(fract*temp.size());
				else
				{
					top = round(double((fract-1)*numAssignableContacts));
					if (top >= temp.size())
						top = temp.size()-1;
				}
				bounds[index] = temp[top];
				printf("ScoreTermCount1: Term: %d   NewBounds: %5.3f\n",type,bounds[index]);
			}
		}
	}
	ScoreTermCondition1* stc = new ScoreTermCondition1(terms,numTypes,bounds); // constructor copies the arrays
	delete [] terms;
	delete [] bounds;
	return stc;
}

InStructureCondition* parseInStructure(stringstream& tok)
{
	double distCutoff;
	double fraction;
	int countCutoff;
	string temp;
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseInStructure\n");
		exit(-1);
	}
	distCutoff = atof(temp.c_str());
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseInStructure\n");
		exit(-1);
	}
	fraction = atof(temp.c_str());
	if (numContactsWithNOEMatchLR == 0)
	{
		printf("ERROR parseInStructure: numContactsWithNOEMatch is 0\n");
		exit(-1);
	}
	countCutoff = fraction*numContactsWithNOEMatchLR;
	InStructureCondition* isc = new InStructureCondition(distCutoff,countCutoff,contactMap);
	return isc;
}

InStructureConditionFract* parseInStructureFract(stringstream& tok, int numAssignableContacts)
{
	int minCount;
	double distCutoff;
	double fractCutoff;
	double maxCountF;
	int maxCount;
	int seqSep;
	string temp;
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseInStructureFract\n");
		exit(-1);
	}
	minCount = atoi(temp.c_str());
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseInStructureFract\n");
		exit(-1);
	}
	distCutoff = atof(temp.c_str());
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseInStructureFract\n");
		exit(-1);
	}
	fractCutoff = atof(temp.c_str());
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseInStructureFract\n");
		exit(-1);
	}
	maxCountF = atof(temp.c_str());
	if (maxCountF < 0)
		maxCount = round(-maxCountF*numAssignableContacts);
	else
	{
		maxCount = (int)maxCountF;
		if (maxCount == 0)
			maxCount = 100000000; // no limit value
	}
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseInStructureFract\n");
		exit(-1);
	}
	seqSep = atoi(temp.c_str());
	InStructureConditionFract* iscf = new InStructureConditionFract(minCount,
			distCutoff, fractCutoff, contactMap, maxCount, seqSep);
	return iscf;
}

InStructureConditionFractTemp* parseInStructureFractTemp(stringstream& tok)
{
	double distCutoff;
	double fractTempCutoff;
	string temp;
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseInStructureFractTemp\n");
		exit(-1);
	}
	distCutoff = atof(temp.c_str());
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseInStructureFractTemp\n");
		exit(-1);
	}
	fractTempCutoff = atof(temp.c_str());
	InStructureConditionFractTemp* is = new InStructureConditionFractTemp(distCutoff, fractTempCutoff, structures);
	return is;
}

// if assposFlag = true, then assignPossibMatrix is used; else assignMatrix is used
AssignedCondition* parseAssignedCount(stringstream& tok, bool assposFlag, int numAssignableContacts)
{
	int countCutoff = 0;
	int windowSize = 0;
	string temp;
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseAssignedCount\n");
		exit(-1);
	}
	countCutoff = atoi(temp.c_str());
	double countFract = atof(temp.c_str());
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseAssignedCount\n");
		exit(-1);
	}
	windowSize = atoi(temp.c_str());
	if (assposFlag)
	{
		if (countFract < 0)
		{
			countFract = -countFract;
			vector<int> counts;
			if (countFract <= 1.0)
			{
				counts.reserve(MAXPROSIZE*(MAXPROSIZE-1)/2);
				for (int i = 0; i < MAXPROSIZE; i++)
				{
					for (int j = i; j < MAXPROSIZE; j++)
					{
						if ( (j-windowSize-i-windowSize) < 6 )
							continue; // don't consider close residue pairs to compute countCutoff
						if (assignPossibMatrix[i][j] > 0) // consider only positive entries to get countCutoff
						{
							int c = 0;
							for (int w1 = -windowSize; w1 <= windowSize; w1++)
							{
								if (i+w1 < 0)
									continue;
								if (i+w1 >= MAXPROSIZE)
									continue;
								for (int w2 = -windowSize; w2 <= windowSize; w2++)
								{
									if (j+w2 < 0)
										continue;
									if (j+w2 >= MAXPROSIZE)
										continue;
									c += assignPossibMatrix[i+w1][j+w2];
								}
							}
							counts.push_back(c);
						}
					}
				}
			}
			else
			{
				counts.reserve(assignPossibPtr->size());
				for (int i = 0; i < MAXPROSIZE; i++)
				{
					for (int j = i; j < MAXPROSIZE; j++)
					{
						if ( (j-windowSize-i-windowSize) < 6 )
							continue; // don't consider close residue pairs to compute countCutoff
						if (assignPossibMatrix[i][j] > 0) // consider only positive entries to get countCutoff
						{
							int c = 0;
							for (int w1 = -windowSize; w1 <= windowSize; w1++)
							{
								if (i+w1 < 0)
									continue;
								if (i+w1 >= MAXPROSIZE)
									continue;
								for (int w2 = -windowSize; w2 <= windowSize; w2++)
								{
									if (j+w2 < 0)
										continue;
									if (j+w2 >= MAXPROSIZE)
										continue;
									c += assignPossibMatrix[i+w1][j+w2];
								}
							}
							for (int k = 0; k < assignPossibMatrix[i][j]; k++)
								counts.push_back(c); // each contact can have many possibilities, so count each one
						}
					}
				}
			}
			sort(counts.begin(),counts.end()); // ascending order
			reverse(counts.begin(),counts.end()); // descending order

			unsigned int top = 0;
			if (countFract <= 1.0)
			{
				top = round(countFract*counts.size());
			}
			else
			{
				top = round(double((countFract-1)*numAssignableContacts));
				if (top >= counts.size())
					top = counts.size()-1;
			}
			// if no count, set cutoff to big value so that false is always returned by AssignedCondition
			countCutoff = counts.size() > 0 ? counts[top] : MAXPROSIZE*MAXPROSIZE;
			if (counts.size() > 0)
				printf("AssignedCount: NewCountCutoff: %d   MaxCount: %d\n",countCutoff,counts[0]);
			else
				printf("AssignedCount: NewCountCutoff: %d\n",countCutoff);
		} // end if < countFract
		AssignedCondition* ac = new AssignedCondition(assignPossibMatrix, countCutoff, windowSize);
		return ac;
	}
	else
	{
		if (countFract < 0)
		{
			countFract = -countFract;
			vector<int> counts;
			if (countFract <= 1.0)
			{
				counts.reserve(MAXPROSIZE*(MAXPROSIZE-1)/2);
				for (int i = 0; i < MAXPROSIZE; i++)
				{
					for (int j = i; j < MAXPROSIZE; j++)
					{
						if ( (j-windowSize-i-windowSize) < 6 )
							continue; // don't consider close residue pairs to compute countCutoff
						if (assignMatrix[i][j] > 0) // consider only positive entries to get countCutoff
						{
							int c = 0;
							for (int w1 = -windowSize; w1 <= windowSize; w1++)
							{
								if (i+w1 < 0)
									continue;
								if (i+w1 >= MAXPROSIZE)
									continue;
								for (int w2 = -windowSize; w2 <= windowSize; w2++)
								{
									if (j+w2 < 0)
										continue;
									if (j+w2 >= MAXPROSIZE)
										continue;
									c += assignMatrix[i+w1][j+w2];
								}
							}
							counts.push_back(c);
						}
					}
				}
			}
			else
			{
				int assignMatrixWin[MAXPROSIZE][MAXPROSIZE] = {0};
				for (int i = 0; i < MAXPROSIZE; i++)
				{
					for (int j = i; j < MAXPROSIZE; j++)
					{
						if ( (j-windowSize-i-windowSize) < 6 )
							continue; // don't consider close residue pairs to compute countCutoff
						if (assignMatrix[i][j] > 0) // consider only positive entries to get countCutoff
						{
							int c = 0;
							for (int w1 = -windowSize; w1 <= windowSize; w1++)
							{
								if (i+w1 < 0)
									continue;
								if (i+w1 >= MAXPROSIZE)
									continue;
								for (int w2 = -windowSize; w2 <= windowSize; w2++)
								{
									if (j+w2 < 0)
										continue;
									if (j+w2 >= MAXPROSIZE)
										continue;
									c += assignMatrix[i+w1][j+w2];
								}
							}
							assignMatrixWin[i][j] = c;
							if (i != j)
								assignMatrixWin[j][i] = c;
						}
					}
				}
				counts.reserve(assignedPtr->size());
				for (vector<PCAssignment*>::iterator it = assignedPtr->begin(); it != assignedPtr->end(); ++it)
				{
					PCAssignment* pca = *it;
					int r1;
					int r2;
					pca->pc.getResPair(r1,r2);
					counts.push_back(assignMatrixWin[r1-1][r2-1]);
				}
			}
			sort(counts.begin(),counts.end()); // ascending order
			reverse(counts.begin(),counts.end()); // descending order

			unsigned int top = 0;
			if (countFract <= 1.0)
			{
				top = round(countFract*counts.size());
			}
			else
			{
				top = round(double((countFract-1)*numAssignableContacts));
				if (top >= counts.size())
					top = counts.size()-1;
			}
			// if no count, set cutoff to big value so that false is always returned by AssignedCondition
			countCutoff = counts.size() > 0 ? counts[top] : MAXPROSIZE*MAXPROSIZE;
			if (counts.size() > 0)
				printf("AssignedCount: NewCountCutoff: %d   MaxCount: %d\n",countCutoff,counts[0]);
			else
				printf("AssignedCount: NewCountCutoff: %d\n",countCutoff);
		} // end if < countFract
		AssignedCondition* ac = new AssignedCondition(assignMatrix, countCutoff, windowSize);
		return ac;
	}
}

AssignedCondition2* parseAssignedCount2(stringstream& tok, bool assposFlag,
		int numAssignableContactsSR, int numAssignableContactsLR)
{
	int countCutoffSR = 0;
	int windowSizeSR = 0;
	int countCutoffLR = 0;
	int windowSizeLR = 0;
	string temp;
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseAssignedCount2\n");
		exit(-1);
	}
	countCutoffSR = atoi(temp.c_str());
	double countFractSR = atof(temp.c_str()); // float version of countCutoffSR
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseAssignedCount2\n");
		exit(-1);
	}
	windowSizeSR = atoi(temp.c_str());
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseAssignedCount2\n");
		exit(-1);
	}
	countCutoffLR = atoi(temp.c_str());
	double countFractLR = atof(temp.c_str());
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseAssignedCount2\n");
		exit(-1);
	}
	windowSizeLR = atoi(temp.c_str());

	if (assposFlag)
	{
		// SR
		if (countFractSR < 0)
		{
			countFractSR = -countFractSR;
			vector<int> counts;
			if (countFractSR <= 1.0)
			{
				counts.reserve(MAXPROSIZE*(MAXPROSIZE-1)/2);
				for (int i = 0; i < MAXPROSIZE; i++)
				{
					for (int j = i; j < MAXPROSIZE; j++)
					{
						if (abs(j-i) >= LRSEQSEP)
							continue;
						if (assignPossibMatrix[i][j] > 0)
						{
							int c = 0;
							for (int w1 = -windowSizeSR; w1 <= windowSizeSR; w1++)
							{
								if (i+w1 < 0)
									continue;
								if (i+w1 >= MAXPROSIZE)
									continue;
								for (int w2 = -windowSizeSR; w2 <= windowSizeSR; w2++)
								{
									if (j+w2 < 0)
										continue;
									if (j+w2 >= MAXPROSIZE)
										continue;
									if ( (i+w1) > (j+w2) )
										continue; // prevent double counting
									c += assignPossibMatrix[i+w1][j+w2];
								}
							}
							counts.push_back(c);
						}
					}
				}
			}
			else
			{
				counts.reserve(assignPossibPtr->size());
				for (int i = 0; i < MAXPROSIZE; i++)
				{
					for (int j = i; j < MAXPROSIZE; j++)
					{
						if (abs(j-i) >= LRSEQSEP)
							continue;
						if (assignPossibMatrix[i][j] > 0)
						{
							int c = 0;
							for (int w1 = -windowSizeSR; w1 <= windowSizeSR; w1++)
							{
								if (i+w1 < 0)
									continue;
								if (i+w1 >= MAXPROSIZE)
									continue;
								for (int w2 = -windowSizeSR; w2 <= windowSizeSR; w2++)
								{
									if (j+w2 < 0)
										continue;
									if (j+w2 >= MAXPROSIZE)
										continue;
									if ( (i+w1) > (j+w2) )
										continue; // prevent double counting
									c += assignPossibMatrix[i+w1][j+w2];
								}
							}
							for (int k = 0; k < assignPossibMatrix[i][j]; k++)
								counts.push_back(c);
						}
					}
				}
			}
			sort(counts.begin(),counts.end()); // ascending order
			reverse(counts.begin(),counts.end()); // descending order

			unsigned int top = 0;
			if (countFractSR <= 1.0)
			{
				top = round(countFractSR*counts.size());
			}
			else
			{
				top = round(double((countFractSR-1)*numAssignableContactsSR));
				if (top >= counts.size())
					top = counts.size()-1;
			}
			// if no count, set cutoff to big value so that false is always returned by AssignedCondition
			countCutoffSR = counts.size() > 0 ? counts[top] : MAXPROSIZE*MAXPROSIZE;
			if (counts.size() > 0)
				printf("AssignedCount2: NewCountCutoffSR: %d   MaxCount: %d\n",countCutoffSR,counts[0]);
			else
				printf("AssignedCount2: NewCountCutoffSR: %d\n",countCutoffSR);
		} // end if < countFractSR

		// LR
		if (countFractLR < 0)
		{
			countFractLR = -countFractLR;
			vector<int> counts;
			if (countFractLR <= 1.0)
			{
				counts.reserve(MAXPROSIZE*(MAXPROSIZE-1)/2);
				for (int i = 0; i < MAXPROSIZE; i++)
				{
					for (int j = i; j < MAXPROSIZE; j++)
					{
						if ( (j-i) < LRSEQSEP )
							continue;
						if (assignPossibMatrix[i][j] > 0)
						{
							int c = 0;
							for (int w1 = -windowSizeLR; w1 <= windowSizeLR; w1++)
							{
								if (i+w1 < 0)
									continue;
								if (i+w1 >= MAXPROSIZE)
									continue;
								for (int w2 = -windowSizeLR; w2 <= windowSizeLR; w2++)
								{
									if (j+w2 < 0)
										continue;
									if (j+w2 >= MAXPROSIZE)
										continue;
									if ( (i+w1) > (j+w2) )
										continue; // prevent double counting
									c += assignPossibMatrix[i+w1][j+w2];
								}
							}
							counts.push_back(c);
						}
					}
				}
			}
			else
			{
				counts.reserve(assignPossibPtr->size());
				for (int i = 0; i < MAXPROSIZE; i++)
				{
					for (int j = i; j < MAXPROSIZE; j++)
					{
						if ( (j-i) < LRSEQSEP )
							continue;
						if (assignPossibMatrix[i][j] > 0)
						{
							int c = 0;
							for (int w1 = -windowSizeLR; w1 <= windowSizeLR; w1++)
							{
								if (i+w1 < 0)
									continue;
								if (i+w1 >= MAXPROSIZE)
									continue;
								for (int w2 = -windowSizeLR; w2 <= windowSizeLR; w2++)
								{
									if (j+w2 < 0)
										continue;
									if (j+w2 >= MAXPROSIZE)
										continue;
									if ( (i+w1) > (j+w2) )
										continue; // prevent double counting
									c += assignPossibMatrix[i+w1][j+w2];
								}
							}
							for (int k = 0; k < assignPossibMatrix[i][j]; k++)
								counts.push_back(c);
						}
					}
				}
			}
			sort(counts.begin(),counts.end()); // ascending order
			reverse(counts.begin(),counts.end()); // descending order

			unsigned int top = 0;
			if (countFractLR <= 1.0)
			{
				top = round(countFractLR*counts.size());
			}
			else
			{
				top = round(double((countFractLR-1)*numAssignableContactsLR));
				if (top >= counts.size())
					top = counts.size()-1;
			}
			// if no count, set cutoff to big value so that false is always returned by AssignedCondition
			countCutoffLR = counts.size() > 0 ? counts[top] : MAXPROSIZE*MAXPROSIZE;
			if (counts.size() > 0)
				printf("AssignedCount2: NewCountCutoffLR: %d   MaxCount: %d\n",countCutoffLR,counts[0]);
			else
				printf("AssignedCount2: NewCountCutoffLR: %d\n",countCutoffLR);
		} // end if < countFractLR

		AssignedCondition2* ac = new AssignedCondition2(assignPossibMatrix, countCutoffSR, windowSizeSR, countCutoffLR, windowSizeLR);
		return ac;
	}
	else
	{
		// assigned
		// SR
		if (countFractSR < 0)
		{
			countFractSR = -countFractSR;
			vector<int> counts;
			if (countFractSR <= 1.0)
			{
				counts.reserve(MAXPROSIZE*(MAXPROSIZE-1)/2);
				for (int i = 0; i < MAXPROSIZE; i++)
				{
					for (int j = i; j < MAXPROSIZE; j++)
					{
						if (abs(j-i) >= LRSEQSEP)
							continue;
						if (assignMatrix[i][j] > 0)
						{
							int c = 0;
							for (int w1 = -windowSizeSR; w1 <= windowSizeSR; w1++)
							{
								if (i+w1 < 0)
									continue;
								if (i+w1 >= MAXPROSIZE)
									continue;
								for (int w2 = -windowSizeSR; w2 <= windowSizeSR; w2++)
								{
									if (j+w2 < 0)
										continue;
									if (j+w2 >= MAXPROSIZE)
										continue;
									if ( (i+w1) > (j+w2) )
										continue; // prevent double counting
									c += assignMatrix[i+w1][j+w2];
								}
							}
							counts.push_back(c);
						}
					}
				}
			}
			else
			{
				counts.reserve(assignedPtr->size());
				for (int i = 0; i < MAXPROSIZE; i++)
				{
					for (int j = i; j < MAXPROSIZE; j++)
					{
						if (abs(j-i) >= LRSEQSEP)
							continue;
						if (assignMatrix[i][j] > 0)
						{
							int c = 0;
							for (int w1 = -windowSizeSR; w1 <= windowSizeSR; w1++)
							{
								if (i+w1 < 0)
									continue;
								if (i+w1 >= MAXPROSIZE)
									continue;
								for (int w2 = -windowSizeSR; w2 <= windowSizeSR; w2++)
								{
									if (j+w2 < 0)
										continue;
									if (j+w2 >= MAXPROSIZE)
										continue;
									if ( (i+w1) > (j+w2) )
										continue; // prevent double counting
									c += assignMatrix[i+w1][j+w2];
								}
							}
							for (int k = 0; k < assignMatrix[i][j]; k++)
								counts.push_back(c);
						}
					}
				}
			}
			sort(counts.begin(),counts.end()); // ascending order
			reverse(counts.begin(),counts.end()); // descending order

			unsigned int top = 0;
			if (countFractSR <= 1.0)
			{
				top = round(countFractSR*counts.size());
			}
			else
			{
				top = round(double((countFractSR-1)*numAssignableContactsSR));
				if (top >= counts.size())
					top = counts.size()-1;
			}
			// if no count, set cutoff to big value so that false is always returned by AssignedCondition
			countCutoffSR = counts.size() > 0 ? counts[top] : MAXPROSIZE*MAXPROSIZE;
			if (counts.size() > 0)
				printf("AssignedCount2: NewCountCutoffSR: %d   MaxCount: %d\n",countCutoffSR,counts[0]);
			else
				printf("AssignedCount2: NewCountCutoffSR: %d\n",countCutoffSR);
		} // end if < countFractSR

		// LR
		if (countFractLR < 0)
		{
			countFractLR = -countFractLR;
			vector<int> counts;
			if (countFractLR <= 1.0)
			{
				counts.reserve(MAXPROSIZE*(MAXPROSIZE-1)/2);
				for (int i = 0; i < MAXPROSIZE; i++)
				{
					for (int j = i; j < MAXPROSIZE; j++)
					{
						if ( (j-i) < LRSEQSEP )
							continue;
						if (assignMatrix[i][j] > 0)
						{
							int c = 0;
							for (int w1 = -windowSizeLR; w1 <= windowSizeLR; w1++)
							{
								if (i+w1 < 0)
									continue;
								if (i+w1 >= MAXPROSIZE)
									continue;
								for (int w2 = -windowSizeLR; w2 <= windowSizeLR; w2++)
								{
									if (j+w2 < 0)
										continue;
									if (j+w2 >= MAXPROSIZE)
										continue;
									if ( (i+w1) > (j+w2) )
										continue; // prevent double counting
									c += assignMatrix[i+w1][j+w2];
								}
							}
							counts.push_back(c);
						}
					}
				}
			}
			else
			{
				counts.reserve(assignedPtr->size());
				for (int i = 0; i < MAXPROSIZE; i++)
				{
					for (int j = i; j < MAXPROSIZE; j++)
					{
						if ( (j-i) < LRSEQSEP )
							continue;
						if (assignMatrix[i][j] > 0)
						{
							int c = 0;
							for (int w1 = -windowSizeLR; w1 <= windowSizeLR; w1++)
							{
								if (i+w1 < 0)
									continue;
								if (i+w1 >= MAXPROSIZE)
									continue;
								for (int w2 = -windowSizeLR; w2 <= windowSizeLR; w2++)
								{
									if (j+w2 < 0)
										continue;
									if (j+w2 >= MAXPROSIZE)
										continue;
									if ( (i+w1) > (j+w2) )
										continue; // prevent double counting
									c += assignMatrix[i+w1][j+w2];
								}
							}
							for (int k = 0; k < assignMatrix[i][j]; k++)
								counts.push_back(c);
						}
					}
				}
			}
			sort(counts.begin(),counts.end()); // ascending order
			reverse(counts.begin(),counts.end()); // descending order

			unsigned int top = 0;
			if (countFractLR <= 1.0)
			{
				top = round(countFractLR*counts.size());
			}
			else
			{
				top = round(double((countFractLR-1)*numAssignableContactsLR));
				if (top >= counts.size())
					top = counts.size()-1;
			}
			// if no count, set cutoff to big value so that false is always returned by AssignedCondition
			countCutoffLR = counts.size() > 0 ? counts[top] : MAXPROSIZE*MAXPROSIZE;
			if (counts.size() > 0)
				printf("AssignedCount2: NewCountCutoffLR: %d   MaxCount: %d\n",countCutoffLR,counts[0]);
			else
				printf("AssignedCount2: NewCountCutoffLR: %d\n",countCutoffLR);
		} // end if < countFractLR
		AssignedCondition2* ac = new AssignedCondition2(assignMatrix, countCutoffSR, windowSizeSR, countCutoffLR, windowSizeLR);
		return ac;
	}
}

TestCondition* parseCondition(string& conditionType, stringstream& tok, int numAssignments, bool assposFlag,
		int numAssignableContactsSR, int numAssignableContactsLR, int numAssignedSR_prev, int numAssignedLR_prev); // forward declaration for parseOr and parseAnd

OrCondition* parseOr(stringstream& tok, int numAssignments, bool assposFlag,
		int numAssignableContactsSR, int numAssignableContactsLR, int numAssignedSR_prev, int numAssignedLR_prev)
{
	string conditionType;
	if (tok.good())
		tok >> conditionType;
	else
	{
		printf("INPUT ERROR: parseOr\n");
		exit(-1);
	}
	TestCondition* a = parseCondition(conditionType,tok, numAssignments, assposFlag,
			numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
	if (tok.good())
		tok >> conditionType;
	else
	{
		printf("INPUT ERROR: parseOr\n");
		exit(-1);
	}
	TestCondition* b = parseCondition(conditionType,tok, numAssignments, assposFlag,
			numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
	return new OrCondition(a,b);
}

AndCondition* parseAnd(stringstream& tok, int numAssignments, bool assposFlag,
		int numAssignableContactsSR, int numAssignableContactsLR, int numAssignedSR_prev, int numAssignedLR_prev)
{
	string conditionType;
	if (tok.good())
		tok >> conditionType;
	else
	{
		printf("INPUT ERROR: parseAnd\n");
		exit(-1);
	}
	TestCondition* a = parseCondition(conditionType,tok, numAssignments, assposFlag,
			numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
	if (tok.good())
		tok >> conditionType;
	else
	{
		printf("INPUT ERROR: parseAnd\n");
		exit(-1);
	}
	TestCondition* b = parseCondition(conditionType,tok, numAssignments, assposFlag,
			numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
	return new AndCondition(a,b);
}

NotCondition* parseNot(stringstream& tok, int numAssignments, bool assposFlag,
		int numAssignableContactsSR, int numAssignableContactsLR, int numAssignedSR_prev, int numAssignedLR_prev)
{
	string conditionType;
	if (tok.good())
		tok >> conditionType;
	else
	{
		printf("INPUT ERROR: parseNot\n");
		exit(-1);
	}
	TestCondition* a = parseCondition(conditionType,tok, numAssignments, assposFlag,
			numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
	return new NotCondition(a);
}

IfElseCondition* parseIf(stringstream& tok, int numAssignments, bool assposFlag,
		int numAssignableContactsSR, int numAssignableContactsLR, int numAssignedSR_prev, int numAssignedLR_prev)
{
	string conditionType;
	if (tok.good())
		tok >> conditionType;
	else
	{
		printf("INPUT ERROR: parseIf\n");
		exit(-1);
	}
	TestCondition* tc = parseCondition(conditionType, tok, numAssignments, assposFlag,
			numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
	if (tok.good())
		tok >> conditionType;
	else
	{
		printf("INPUT ERROR: parseIf\n");
		exit(-1);
	}
	TestCondition* ifBody = parseCondition(conditionType,tok, numAssignments, assposFlag,
			numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
	if (tok.good())
		tok >> conditionType; // else
	else
	{
		printf("INPUT ERROR: parseIf\n");
		exit(-1);
	}
	if (conditionType != "else")
	{
		printf("if condition missing else statement\n");
		exit(-1);
	}
	if (tok.good())
		tok >> conditionType;
	else
	{
		printf("INPUT ERROR: parseIf\n");
		exit(-1);
	}
	TestCondition* elseBody = parseCondition(conditionType,tok, numAssignments, assposFlag,
			numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
	IfElseCondition* ifElse = new IfElseCondition(tc, ifBody, elseBody);
	return ifElse;
}

InStructureCutCondition* parseInStructureCutCondition(stringstream& tok, int numAssignments)
{
	double fractCutoff;
	double fractNumAssignedCut;
	string temp;
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseInStructureCutCondition\n");
		exit(-1);
	}
	fractCutoff = atof(temp.c_str());
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseInStructureCutCondition\n");
		exit(-1);
	}
	fractNumAssignedCut = atof(temp.c_str());
	InStructureCutCondition* isc = new InStructureCutCondition(fractCutoff, fractNumAssignedCut, numAssignments, contactMap);
	return isc;
}

CorrectFilter* parseCorrectFilter(stringstream& tok)
{
	double fractCutoff;
	string temp;
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseCorrectFilter\n");
		exit(-1);
	}
	fractCutoff = atof(temp.c_str());
	if (refPDB.size() == 0)
	{
		printf("ERROR: parseCorrectFilter. refPDB does not exist\n");
		exit(-1);
	}
	CorrectFilter* cf = new CorrectFilter(refContactMap, DISTCUTOFF, fractCutoff);
	return cf;
}

InStructureWindow* parseInStructureWindow(stringstream& tok)
{
	int winLen;
	double fractViol;
	string temp;
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseInStructureWindow\n");
		exit(-1);
	}
	winLen = atoi(temp.c_str());
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseInStructureWindow\n");
		exit(-1);
	}
	fractViol = atof(temp.c_str());
	InStructureWindow* isw = new InStructureWindow(winLen,fractViol,contactMap,DISTCUTOFF);
	return isw;
}

NoiseCondition* parseNoiseCondition(stringstream& tok)
{
	string temp;
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseNoiseCondition\n");
		exit(-1);
	}
	double ratio = atof(temp.c_str());
	if (ratio <= 0)
	{
		printf("INPUT ERROR: parseNoiseCondition. Ratio must be positive\n");
		exit(-1);
	}
	int numPeaks = 0;
	for (vector<NOE::PEAKLISTTYPE>::iterator it = peakListTypes.begin(); it != peakListTypes.end(); ++it)
	{
		NOE::PEAKLISTTYPE peakListType = *it;
		list<NOE*>& noes = peakLists[peakListType];
		numPeaks += noes.size();
	}
	if (numPeaks == 0)
	{
		printf("INPUT ERROR: parseNoiseCondition. NumPeaks must be positive\n");
		exit(-1);
	}
	int numRes = bmrbProtein->size;
	if (numRes == 0)
	{
		printf("INPUT ERROR: parseNoiseCondition. NumRes must be positive\n");
		exit(-1);
	}
	NoiseCondition* nc = new NoiseCondition(numPeaks,numRes, ratio);
	return nc;
}

NumAssignmentsCondition* parseNumAssignmentsCondition(stringstream& tok,
		bool isLR, int numAssignableContactsSR, int numAssignableContactsLR,
		int numAssignedContactsSR_prev, int numAssignedContactsLR_prev)
{
	string temp;
	if (tok.good())
		tok >> temp;
	else
	{
		printf("INPUT ERROR: parseNumAssignmentsCondition\n");
		exit(-1);
	}
	double fract = atof(temp.c_str());
	int numAssigned = (assignedPtr != NULL ? assignedPtr->size() : 0);
	if (numAssigned == 0)
	{
		printf("ERROR: parseNumAssignmentsCondition. Nothing was assigned.\n");
		exit(-1);
	}
	if (isLR)
	{
		return new NumAssignmentsCondition(isLR, numAssignedContactsLR_prev, fract, numAssignableContactsLR);
	}
	else
	{
		return new NumAssignmentsCondition(isLR, numAssignedContactsSR_prev, fract, numAssignableContactsSR);
	}
}

// if assposFlag=true, assignPossibMatrix will be used for AssignedCondition; else assignMatrix is used
TestCondition* parseCondition(string& conditionType, stringstream& tok, int numAssignments, bool assposFlag,
		int numAssignableContactsSR, int numAssignableContactsLR, int numAssignedSR_prev, int numAssignedLR_prev)
{
	TestCondition* tc = NULL;
	int numAssignableContacts = numAssignableContactsSR + numAssignableContactsLR;
	if (conditionType == "scoreterm2")
		tc = parseScoreTerm2(tok, assposFlag, numAssignableContactsSR, numAssignableContactsLR);
	else if (conditionType == "scoreterm1")
		tc = parseScoreTerm1(tok, assposFlag, numAssignableContacts);
	else if (conditionType == "scorecount1")
		tc = parseScoreCount1(tok, assposFlag, numAssignableContacts);
	else if (conditionType == "scorecount2")
		tc = parseScoreCount2(tok, assposFlag, numAssignableContactsSR, numAssignableContactsLR);
	else if (conditionType == "instructure")
		tc = parseInStructure(tok);
	else if (conditionType == "instructurefract" || conditionType == "instructurefrac")
		tc = parseInStructureFract(tok, numAssignableContacts);
	else if (conditionType == "instructurefracttemp" || conditionType == "instructurefractemp")
		tc = parseInStructureFractTemp(tok);
	else if (conditionType == "assignedcount")
		tc = parseAssignedCount(tok, assposFlag, numAssignableContacts);
	else if (conditionType == "assignedcount2")
		tc = parseAssignedCount2(tok, assposFlag, numAssignableContactsSR, numAssignableContactsLR);
	else if (conditionType == "and")
		tc = parseAnd(tok, numAssignments, assposFlag, numAssignableContactsSR, numAssignableContactsLR,
				numAssignedSR_prev, numAssignedLR_prev);
	else if (conditionType == "or")
		tc = parseOr(tok, numAssignments, assposFlag, numAssignableContactsSR, numAssignableContactsLR,
				numAssignedSR_prev, numAssignedLR_prev);
	else if (conditionType == "true")
		tc = new TrueCondition();
	else if (conditionType == "false")
		tc = new FalseCondition();
	else if (conditionType == "seqsep")
	{
		string temp;
		if (tok.good())
			tok >> temp;
		else
		{
			printf("INPUT ERROR: parseSeqSep\n");
			exit(-1);
		}
		int seqSepCutoff = atoi(temp.c_str());
		tc = new SeqSepCondition(seqSepCutoff);
	}
	else if (conditionType == "if")
		tc = parseIf(tok, numAssignments, assposFlag, numAssignableContactsSR, numAssignableContactsLR,
				numAssignedSR_prev, numAssignedLR_prev);
	else if (conditionType == "null")
		tc = new TrueCondition();
	else if (conditionType == "not")
		tc = parseNot(tok, numAssignments, assposFlag, numAssignableContactsSR, numAssignableContactsLR,
				numAssignedSR_prev, numAssignedLR_prev);
	else if (conditionType == "instructurecut")
		tc = parseInStructureCutCondition(tok, numAssignments);
	else if (conditionType == "correctfilter")
		tc = parseCorrectFilter(tok);
	else if (conditionType == "instructurewindow")
		tc = parseInStructureWindow(tok);
	else if (conditionType == "noiseratio")
		tc = parseNoiseCondition(tok);
	else if (conditionType == "numlr")
		tc = parseNumAssignmentsCondition(tok, true, numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
	else if (conditionType == "numsr")
		tc = parseNumAssignmentsCondition(tok, false, numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
	else
	{
		printf("parseCondition: Invalid filter conditionType %s\n",conditionType.c_str());
		exit(-1);
	}
	return tc;
}

// used by parseFilters
// if assposFlag=true, assignPossibMatrix will be used for AssignedCondition; else assignMatrix is used
TestCondition* parseFilterTest(stringstream& tok, int numAssignments,bool assposFlag,
		int numAssignableContactsSR, int numAssignableContactsLR,
		int numAssignedSR_prev, int numAssignedLR_prev)
{
	string filterType;
	if (tok.good())
		tok >> filterType;
	else
	{
		printf("INPUT ERROR: parseFilterTest\n");
		exit(-1);
	}
	TestCondition* ft = parseCondition(filterType, tok, numAssignments, assposFlag,
			numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
	return ft;
}

// used by assign()
// if assposFlag=true, filterAssPos is to be initialized; else filterAss is to be initialized
// numAssignments is either the number of assignment possibilities or the number of assigned
void initializeAssignmentPossibFilter(int numAssignments, list<PCAssignment*>& asspos, vector<PCAssignment*>& ass,
		int numAssignableContactsSR, int numAssignableContactsLR,
		int numAssignedSR_prev, int numAssignedLR_prev)
{
	using namespace std;
	if (filterAssPosStr != "")
	{
		assignPossibPtr = &asspos;
		makeAssignPossibMatrix(asspos);
		transform(filterAssPosStr.begin(), filterAssPosStr.end(), filterAssPosStr.begin(), ::tolower);
		printf("FilterAssPos: %s\n",filterAssPosStr.c_str());
		stringstream tok(filterAssPosStr);
		filterAssPos = parseFilterTest(tok,numAssignments,true,
				numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
		filterAssPos->print(0);
	}
}
void initializeAssignedFilter(int numAssignments, list<PCAssignment*>& asspos, vector<PCAssignment*>& ass,
		int numAssignableContactsSR, int numAssignableContactsLR, int numAssignedSR_prev, int numAssignedLR_prev)
{
	using namespace std;
	if (filterAssStr != "")
	{
		assignedPtr = &ass;
		makeAssignMatrix(ass);
		transform(filterAssStr.begin(), filterAssStr.end(), filterAssStr.begin(), ::tolower);
		printf("FilterAss: %s\n",filterAssStr.c_str());
		stringstream tok(filterAssStr);
		filterAss = parseFilterTest(tok,numAssignments,false, numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
		filterAss->print(0);
	}
}
void initializeAssignedFilter2(int numAssignments, list<PCAssignment*>& asspos, vector<PCAssignment*>& ass,
		int numAssignableContactsSR, int numAssignableContactsLR, int numAssignedSR_prev, int numAssignedLR_prev)
{
	using namespace std;
	if (filterAssStr2 != "")
	{
		assignedPtr = &ass;
		makeAssignMatrix(ass);
		transform(filterAssStr2.begin(), filterAssStr2.end(), filterAssStr2.begin(), ::tolower);
		printf("FilterAss2: %s\n",filterAssStr2.c_str());
		stringstream tok(filterAssStr2);
		filterAss2 = parseFilterTest(tok,numAssignments,false, numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
		filterAss2->print(0);
	}
	else if (filterAssStr != "")
	{
		assignedPtr = &ass;
		makeAssignMatrix(ass);
		transform(filterAssStr.begin(), filterAssStr.end(), filterAssStr.begin(), ::tolower);
		printf("FilterAss2: %s\n",filterAssStr.c_str());
		stringstream tok(filterAssStr);
		filterAss = parseFilterTest(tok,numAssignments,false, numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
		filterAss->print(0);
	}
}
void initializeAmbigAssignmentFilter(int numAssignmentPos, list<PCAssignment*>& asspos,
		int numAssignableContactsSR, int numAssignableContactsLR, int numAssignedSR_prev, int numAssignedLR_prev)
{
	if (filterAmbPosStr != "")
	{
		assignPossibPtr = &asspos;
		makeAssignPossibMatrix(asspos);
		transform(filterAmbPosStr.begin(), filterAmbPosStr.end(), filterAmbPosStr.begin(), ::tolower);
		printf("FilterAmbPos: %s\n",filterAmbPosStr.c_str());
		stringstream tok(filterAmbPosStr);
		filterAmbPos = parseFilterTest(tok,numAssignmentPos, true, // use asspos marix    // false, // use ass matrix
				numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
		filterAmbPos->print(0);
	}
}


void accuracyReport(vector<PCAssignment*>& assigned)
{
	if (refContactMap.getNumContacts() > 0)
	{
		int numCorrect = 0;
		int numUnAmbig = 0;
		int numCorrectLR = 0;
		int numAssignedLR = 0;

		for (vector<PCAssignment*>::iterator it = assigned.begin(); it != assigned.end(); ++it)
		{
			PCAssignment* pca = *it;
			tr1::unordered_set<Contact>& contacts = pca->pc.getContacts();
			int r1;
			int r2;
			pca->pc.getResPair(r1,r2);
			bool correct = false;
			for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
			{
				const Contact& c = *itC;
				tr1::array<double,ContactMap::CMENTRYSIZE>& entry = refContactMap[c];
				if (entry[ContactMap::MINDISTINDEX] <= DISTCUTOFF)
				{
					correct = true;
					break;
				}
			}
			if (correct)
			{
				numCorrect++;
				if (abs(r1-r2) >= LRSEQSEP)
				{
					numCorrectLR++;
				}
			}
			else
			{
				// print reason
				if (pca->score.ambig > 0.999)
					numUnAmbig++;
				// pca->print();
			}
			if (abs(r1-r2) >= LRSEQSEP)
				numAssignedLR++;
		}
		int total = assigned.size();
		double fract = (total > 0 ? double(numCorrect)/double(total) : 1.0);
		// int numIncorrect = total-numCorrect;
		double fractLR = (numAssignedLR > 0 ? double(numCorrectLR)/double(numAssignedLR) : 1.0);
		printf("Acc: %5.3f  %d correct out of %d assigned. AccLR: %5.3f  %d correct out of %d assigned.\n",
				fract,numCorrect,total,fractLR,numCorrectLR,numAssignedLR);
		// double fractUnAmbig = double(numUnAmbig)/double(numIncorrect);
		// printf("FractUnAmbigWrong %5.3f  %d num unambig incorrect out of %d incorrect\n",fractUnAmbig,numUnAmbig,numIncorrect);
	}
	else
	{
		printf("WARNING: Cannot check contact accuracy. No reference pdb was provided\n");
	}
}

/**
 * Stores ambiguous assignments in ambigAss
 * Deleted seed assignments in deletedSeed are added to asspos (filterAmbig will test these and other assignments)
 * Put assignments in pool that have score sum < POOL_FILTER and the other assignments in nonPool
 */
void computeAmbiguousAssignments(vector<PCAssignment*>& assigned,
		list<PCAssignment*>& asspos, // assignment possibilities
		tr1::unordered_map< NOE*, vector<PCAssignment*> >& ambigAss, list<PCAssignment*>& deletedSeed,
		int numAssignableContactsSR, int numAssignableContactsLR,
		int numAssignedSR_prev, int numAssignedLR_prev)
{
	tr1::unordered_set<Contact> assignedContacts; // to prevent adding the seed
	tr1::unordered_set<NOE*> assignedNOEs; // to prevent adding the seed
	// setup assignedContacts, assignedNOEs
	for (vector<PCAssignment*>::iterator it = assigned.begin(); it != assigned.end(); ++it)
	{
		PCAssignment* pca = *it;
		assignedNOEs.insert(pca->noe);
		tr1::unordered_set<Contact>& contacts = pca->pc.getContacts(); // add all contacts to contactAss
		for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
		{
			const Contact& c = *itC;
			tr1::unordered_set<Contact>::iterator itFind = assignedContacts.find(c);
			if (itFind == assignedContacts.end())
				assignedContacts.insert(c);
		}
	}
	// add deletedSeed to ambigAss, but do NOT add to assignedContacts, assignedNOEs
	for (list<PCAssignment*>::iterator it = deletedSeed.begin(); it != deletedSeed.end(); ++it)
	{
		PCAssignment* pca = *it;
		vector<PCAssignment*> listv;
		listv.push_back(pca);
		ambigAss[pca->noe] = listv;
	}
	// setup ambigAss
	int NUMLRASSIGNPOSSIB = 0;
	for (list<PCAssignment*>::iterator it = asspos.begin(); it != asspos.end(); ++it)
	{
		PCAssignment* pca = *it;
		NOE* noe = pca->noe;
		// skip if (1) pc of assignment poss is already assigned and (2) noe is already assigned
		if (assignedNOEs.find(noe) != assignedNOEs.end())
			continue;
		// skip if already added to ambigAss (i.e. assignments from deletedSeed)
		bool alreadyAdded = false;
		for (list<PCAssignment*>::iterator it2 = deletedSeed.begin(); it2 != deletedSeed.end(); ++it2)
		{
			PCAssignment* dPCA = *it2;
			if (dPCA == pca)
			{
				alreadyAdded = true;
				break;
			}
		}
		if (alreadyAdded)
			continue;

		int r1 = 0;
		int r2 = 0;
		pca->pc.getResPair(r1,r2);
		if (abs(r1-r2) >= LRSEQSEP)
			NUMLRASSIGNPOSSIB++;
		tr1::unordered_set<Contact>& contacts = pca->pc.getContacts(); // add all contacts to contactAss
		bool found = false;
		for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
		{
			const Contact& c = *itC;
			tr1::unordered_set<Contact>::iterator itFind = assignedContacts.find(c);
			if (itFind != assignedContacts.end())
			{
				found = true;
				break;
			}
		}
		if (found)
			continue;
		tr1::unordered_map<NOE*, vector<PCAssignment*>  >::iterator itC = ambigAss.find(noe);
		if (itC != ambigAss.end())
		{
			vector<PCAssignment*>& listv = ambigAss[noe];
			listv.push_back(pca);
		}
		else
		{
			vector<PCAssignment*> listv;
			listv.push_back(pca);
			ambigAss[noe] = listv;
		}
	}
	// filter ambigAss based on AMBIG_FILTER and filterAmbPos
	if (filterAmbPosStr != "")
	{
		initializeAmbigAssignmentFilter(NUMLRASSIGNPOSSIB, asspos,
				numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev); // initialize the filters;
	}
	list<NOE*> toRemove; // for removing assignments that fail filters
	for (tr1::unordered_map<NOE*, vector<PCAssignment*>  >::iterator it = ambigAss.begin(); it != ambigAss.end(); ++it)
	{
		vector<PCAssignment*>& ass = it->second;
		vector<PCAssignment*> assCopy(ass.begin(),ass.end());
		ass.clear(); // ass will be re-populated
		sort(assCopy.begin(),assCopy.end(),PCAssignmentSortByScore()); // ascending order
		reverse(assCopy.begin(),assCopy.end()); // descending order
		double scoreSum = 0;
		int count = 0;
		for (vector<PCAssignment*>::iterator itA = assCopy.begin(); itA != assCopy.end(); ++itA)
		{
			PCAssignment* pca = *itA;
			// filter ambig pos
			if (filterAmbPos != NULL)
			{
				if (!filterAmbPos->test(pca))
					continue;
			}
			scoreSum += pca->score.total;
			count++;
			ass.push_back(pca);
			if (AMBIG_FILTER_BY_COUNT)
			{
				if (count >= int(AMBIG_FILTER))
					break;
			}
			else if (scoreSum >= AMBIG_FILTER || count >= AMBIG_FILTER_MAX)
				break;
		}
		if (ass.size() == 0) // failed filterAmbPos
			toRemove.push_back(it->first);
	}
	for (list<NOE*>::iterator it = toRemove.begin(); it != toRemove.end(); ++it)
	{
		NOE* n = *it;
		ambigAss.erase(n);
	}
	if (refPDB.size())
	{
		// check accuracy
		int numCorrect = 0;
		int numAmbig = ambigAss.size();
		int numPC = 0;
		for (tr1::unordered_map<NOE*, vector<PCAssignment*>  >::iterator it = ambigAss.begin(); it != ambigAss.end(); ++it)
		{
			vector<PCAssignment*>& ass = it->second;
			numPC += ass.size();
			bool correct = false;
			for (vector<PCAssignment*>::iterator itA = ass.begin(); itA != ass.end(); ++itA)
			{
				PCAssignment* pca = *itA;
				tr1::unordered_set<Contact>& contacts = pca->pc.getContacts();
				for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
				{
					const Contact& c = *itC;
					tr1::array<double,ContactMap::CMENTRYSIZE>& entry = refContactMap[c];
					if (entry[ContactMap::MINDISTINDEX] <= DISTCUTOFF)
					{
						correct = true;
						break;
					}
				}
				if (correct)
					break;
			}
			if (correct)
				numCorrect++;
		}
		double fractCorrect = (double)numCorrect/(double)numAmbig;
		printf("AmbigAcc: %5.3f  %d correct out of %d.\n",fractCorrect,numCorrect,numAmbig);
		double ambigRatio = (double)numPC/(double)numAmbig;
		printf("AmbigRatio: %5.3f\n",ambigRatio);
	}
}

/*
 * Main assignment function
 */
void assign()
{
	// filter and calibrate the peak lists
	int numNOEsBeforeFilter = 0;
	for (vector<NOE::PEAKLISTTYPE>::iterator it = peakListTypes.begin(); it != peakListTypes.end(); ++it)
	{
		NOE::PEAKLISTTYPE peakListType = *it;
		list<NOE*>& noes = peakLists[peakListType];
		numNOEsBeforeFilter += noes.size();
	}
	printf("Num Peaks Before Filter %d\n",numNOEsBeforeFilter);
	filterCalibratePeakLists(); // adds offsets to the chemical shift values of the peak lists to make the lists align with each other
	                            // if all peaks have negative intensity, change them to positive
	                            // remove duplicate peaks, which may be due to artifacts from automated peak picking
	                            // remove self (diagonal) peaks (we only care about the cross peaks)
	int numNOEsAfterFilter = 0;
	for (vector<NOE::PEAKLISTTYPE>::iterator it = peakListTypes.begin(); it != peakListTypes.end(); ++it)
	{
		NOE::PEAKLISTTYPE peakListType = *it;
		list<NOE*>& noes = peakLists[peakListType];
		numNOEsAfterFilter += noes.size();
	}
	printf("Num Peaks After Filter %d\n",numNOEsAfterFilter);
	int numRes = bmrbProtein->size;
	double noesPerRes = (double)numNOEsAfterFilter/(double)numRes;
	printf("NOEsPerRes: %4.2f\n",noesPerRes);

	makeContactMap(false);
	makeContactMapExpected();
	if (refPDB.size() > 0)
		makeContactMap(true);

	list<Assignment*>* assignmentPossib = new list<Assignment*>; // put on heap since many possible assignments

	// compute assignment possibilities
	computeCS_Str_DB_Score(*assignmentPossib); // the other score terms are initialized to 0
	computeAmbig_Score(*assignmentPossib);
	computeSym_Inter_Net_Score2(*assignmentPossib);
	// set total score
	tr1::unordered_set<Contact> contactsWithNOEMatchLR; // assignment possibilities corresponding to long range contacts in str
	int numContactsWithNOEMatch = 0;

	for (list<Assignment*>::iterator it = assignmentPossib->begin(); it != assignmentPossib->end(); ++it)
	{
		Assignment* a = *it;
		a->score.setTotal();
		Contact& c = a->contact;
		double d = contactMapExpected[c][ContactMap::MINDISTINDEX]; // check if in contact in templates
		if (d <= DISTCUTOFF)
		{
			if (abs(c.r1->num-c.r2->num) >= LRSEQSEP)
			{
				tr1::unordered_set<Contact>::iterator itFind = contactsWithNOEMatchLR.find(c);
				if (itFind == contactsWithNOEMatchLR.end())
					contactsWithNOEMatchLR.insert(c);
			}
			numContactsWithNOEMatch++;
		}
	} // end for each assignment possibility

	numContactsWithNOEMatchLR = contactsWithNOEMatchLR.size();
	printf("NumAssPoss %zu   ContactsWithNOEMatchSizeLR %d\n",assignmentPossib->size(),numContactsWithNOEMatchLR);

	// filter assignment possibilities
	// printf("Num Assignment Possibilities Before Filter %zu\n",assignmentPossib->size());

	int numAssignableContacts = 0;
	int numAssignableContactsLR = 0;
	int numAssignableContactsSR = 0;
	for (ContactMapIterator it = contactMapExpected.begin(); it != contactMapExpected.end(); ++it)
	{
		const Contact& c = it.first();
		const tr1::array<double,ContactMap::CMENTRYSIZE>& entry = it.second();
		if (entry[ContactMap::MINDISTINDEX] <= DISTCUTOFF)
		{
			numAssignableContacts++;
			if (abs(c.r1->num-c.r2->num) < LRSEQSEP)
				numAssignableContactsSR++;
			else
				numAssignableContactsLR++;
		}
	}
	double fractContactsNOEMatch = double(numContactsWithNOEMatch)/double(numAssignableContacts);
	printf("NumAssignableContacts: %d   NumAssignableContactsSR: %d   NumAssignableContactsLR: %d\n",
			numAssignableContacts,numAssignableContactsSR,numAssignableContactsLR);
	printf("NumContactsWithNOEMatch: %d  NumExpectedContacts: %d  Fract: %5.3f\n",numContactsWithNOEMatch,numAssignableContacts,fractContactsNOEMatch);
	double ambigLevel = double(assignmentPossib->size())/double(numNOEsAfterFilter);
	double completeness = double(assignmentPossib->size())/double(numAssignableContacts);
	printf("AmbigLevel: %9.4f\nCompleteness: %9.4f\n",ambigLevel,completeness);
	printf("Num Assignment Possibilities Before Grouping %zu\n",assignmentPossib->size());
	// group assignments together with same NOE assigned to contacts between the same pair of residues (can be same res),
	// but different protons
	tr1::unordered_map<PseudoContactKey, PCAssignment*>* pcaTemp = new tr1::unordered_map<PseudoContactKey, PCAssignment*>; // temporary
	                                                                                                            // key = NOE, res1, res2
	list<PCAssignment*>* assignments = new list<PCAssignment*>; // non-temporary version of pcaTemp, store only the assignments

	int NUMLRASSIGNPOSSIB = 0; // number of assignment possibilities to LR contacts
	for (list<Assignment*>::iterator it = assignmentPossib->begin(); it != assignmentPossib->end(); ++it)
	{
		Assignment* as = *it;
		int r1 = as->contact.r1->num;
		int r2 = as->contact.r2->num;
		if (abs(r1-r2) >= LRSEQSEP)
			NUMLRASSIGNPOSSIB++;
		PseudoContactKey key(as->noe,r1,r2);
		tr1::unordered_map<PseudoContactKey,PCAssignment*>::iterator pcIt = pcaTemp->find(key);
		if (pcIt != pcaTemp->end())
		{
			PCAssignment* pca = pcIt->second;
			// add contact and update score if chem shifts are similar; else create new
			pca->pc.add(as->contact);
			pca->score.setMax(as->score);
		}
		else
		{
			PCAssignment* pca = new PCAssignment(as->noe,PseudoContact(as->contact),as->score);
			(*pcaTemp)[key] = pca;
			assignments->push_back(pca);
		}
	}
	printf("Num Assignment Possibilities After Grouping %zu\n",assignments->size());

	vector<PCAssignment*> dummy;
	if (filterAssPosStr != "")
	{
		initializeAssignmentPossibFilter(NUMLRASSIGNPOSSIB, *assignments, dummy,
				numAssignableContactsSR, numAssignableContactsLR, 0, 0); // initialize the filters
		filter(*filterAssPos, *assignments,true);
		reComputeAmbig_score(*assignments);
	}

	printf("Num Assignment Possibilities After Filter %zu\n",assignments->size());

    // list of possible NOEs for each PseudoContact
	tr1::unordered_map< PseudoContact,tr1::unordered_set<NOE*> >* pc2NOEs =
			new tr1::unordered_map< PseudoContact, tr1::unordered_set<NOE*> >; // used for pseudocontact cardinality constraint
	for (list<PCAssignment*>::iterator it = assignments->begin(); it != assignments->end(); ++it)
	{
		PCAssignment* pca = *it;
		tr1::unordered_map< PseudoContact,tr1::unordered_set<NOE*> >::iterator itFind = pc2NOEs->find(pca->pc);
		if (itFind != pc2NOEs->end())
			itFind->second.insert(pca->noe);
		else
		{
			tr1::unordered_set<NOE*> noes;
			noes.insert(pca->noe);
			(*pc2NOEs)[pca->pc] = noes;
		}
	}

	// intermediate clean up step
	delete pcaTemp;
	for (list<Assignment*>::iterator it = assignmentPossib->begin(); it != assignmentPossib->end(); ++it)
		delete *it;
	delete assignmentPossib;

	// setup constraints and ILP program (actually an LP since matrix is totally unimodular)
	IloEnv env;
	try
	{
		IloModel model(env);
		IloNumVarArray vars(env, assignments->size(), 0.0, 1.0);
		IloNumArray coeffs(env, assignments->size());

		// objective fcn
		int i = 0;
		for (list<PCAssignment*>::iterator it = assignments->begin(); it != assignments->end(); ++it)
		{
			PCAssignment* pca = *it;
			coeffs[i++] = pca->score.total;
		}
		IloObjective objFunc = IloMaximize(env,IloScalProd(coeffs,vars));
		model.add(objFunc);

		// constraints
		tr1::unordered_map<NOE*,IloNumExpr> consN; // each NOE assigned to at most one contact or to noise
		tr1::unordered_map<PseudoContact,IloNumExpr> consPC; // each contact assigned to at most num PseudoContact.size NOE
		i = 0;
		for (list<PCAssignment*>::iterator it = assignments->begin(); it != assignments->end(); ++it)
		{
			PCAssignment* as = *it;
			tr1::unordered_map<NOE*,IloNumExpr>::iterator itFind = consN.find(as->noe);
			if (itFind != consN.end())
				itFind->second += vars[i];
			else
			{
				IloNumExpr expr(env);
				expr += vars[i];
				consN[as->noe] = expr;
			}
			const PseudoContact& pc = as->pc;
			tr1::unordered_map<PseudoContact,IloNumExpr>::iterator itFind2 = consPC.find(pc);
			if (itFind2 != consPC.end())
				itFind2->second += vars[i];
			else
			{
				IloNumExpr expr(env);
				expr += vars[i];
				consPC[pc] = expr;
			}
			i++;
		}

		for (tr1::unordered_map<NOE*,IloNumExpr>::const_iterator it = consN.begin(); it != consN.end(); ++it)
		{
			model.add(IloRange(env,0,it->second,1.0));
		}
		for (tr1::unordered_map<PseudoContact,IloNumExpr>::const_iterator it = consPC.begin(); it != consPC.end(); ++it)
		{
			tr1::unordered_map< PseudoContact,tr1::unordered_set<NOE*> >::iterator itFind = pc2NOEs->find(it->first);
			model.add(IloRange(env,0,it->second,min(it->first.contacts.size(),itFind->second.size())));
		}

		vector<PCAssignment*> bestAssigned;
		list<NOE*> noes;
		for (vector<NOE::PEAKLISTTYPE>::iterator it = peakListTypes.begin(); it != peakListTypes.end(); ++it)
		{  // add all noes to noes
			NOE::PEAKLISTTYPE peakListType = *it;
			list<NOE*>& ns = peakLists[peakListType];
			noes.insert(noes.end(),ns.begin(),ns.end());
		}
		double calibConstants[NOE::NUMTYPES][NUM_CALIB_TYPES]; // dist^-6 = calibConstants*intensity
		list<PCAssignment*> deletedAssignments;
		int numAssignedSR_prev = 0; // num sr assignments after filter (before filter for only the first call to initializeAssignedFilter)
		int numAssignedLR_prev = 0;

		// solve the model
		IloCplex cplex(model); // cplex.setParam(IloCplex::NumParam::EpOpt,0.00000001);
		// cplex.setOut(env.getNullStream());
		cplex.solve();
		env.out() << "Solution status: " << cplex.getStatus() << "\n";
		env.out() << "Solution value: " << cplex.getObjValue() << "\n";

		// fetch result
		int numAssigned = 0;
		vector<PCAssignment*> assigned;
		assigned.reserve(numNOEsAfterFilter);
		i = 0;
		int numAssignedLong = 0;
		for (list<PCAssignment*>::iterator it = assignments->begin(); it != assignments->end(); ++it)
		{
			double val = cplex.getValue(vars[i]);
			if (val > 0.99)
			{
				PCAssignment* as = *it;
				numAssigned++;
				assigned.push_back(as);

				int r1 = 0;
				int r2 = 0;
				as->pc.getResPair(r1,r2);
				if (abs(r1-r2) >= LRSEQSEP)
				{
					numAssignedLong++;
					numAssignedLR_prev++;
				}
				else
					numAssignedSR_prev++;
			}
			i++;
		}
		printf("NumAssigned: %d\n",numAssigned);
		printf("NumAssignedLong: %d\n",numAssignedLong);

		// filter seed assignment
		if (filterAssStr != "")
		{
			initializeAssignedFilter(numAssignedLong, *assignments, assigned, numAssignableContactsSR,
					numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
			deletedAssignments.clear();
			filterAss->reset(); // in case this is not the first iteration
			reComputeSym_Inter_Net_Ambig_Score(*assignments,assigned);
			filter(*filterAss,assigned,false, deletedAssignments);

			numAssignedSR_prev = 0;
			numAssignedLR_prev = 0;

			for (vector<PCAssignment*>::iterator it = assigned.begin(); it != assigned.end(); ++it)
			{
				PCAssignment* pca = *it;
				int seqsep = pca->pc.getSeqSep();
				if (seqsep < LRSEQSEP)
					numAssignedSR_prev++;
				else
					numAssignedLR_prev++;
			}
		}

		// calibrate distances from seed assignment
		computeIntensityScore(noes, *assignments, assigned, calibConstants); // this also calls calibrateDistances

		// recompute seed assignment using intensity scores
		i = 0;
		for (list<PCAssignment*>::iterator it = assignments->begin(); it != assignments->end(); ++it)
		{
			PCAssignment* pca = *it;
			coeffs[i] = pca->score.total;
			objFunc.setLinearCoef(vars[i],coeffs[i]);
			i++;
		}

		cplex.solve();

		env.out() << "Solution status: " << cplex.getStatus() << "\n";
		env.out() << "Solution value: " << cplex.getObjValue() << "\n";
		numAssigned = 0;
		numAssignedLong = 0;
		int numAssignedInContact = 0;
		assigned.clear();
		i = 0;
		for (list<PCAssignment*>::iterator it = assignments->begin(); it != assignments->end(); ++it)
		{
			double val = cplex.getValue(vars[i]);
			if (val > 0.99)
			{
				PCAssignment* as = *it;
				numAssigned++;
				assigned.push_back(as);

				int r1 = 0;
				int r2 = 0;
				as->pc.getResPair(r1,r2);
				if (abs(r1-r2) >= LRSEQSEP)
					numAssignedLong++;

				tr1::unordered_set<Contact>& contacts = as->pc.getContacts();
				for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
				{
					const Contact& c = *itC;
					if (contactMapExpected[c][ContactMap::MINDISTINDEX] <= DISTCUTOFF)
						numAssignedInContact++;
				}
			}
			i++;
		}
		printf("NumAssignedInContact: %d\n",numAssignedInContact);
		printf("NumAssigned: %d\n",numAssigned);
		double fractAssignedInContact = double(numAssignedInContact)/double(numAssigned);
		printf("FractAssignedInContact: %5.3f\n",fractAssignedInContact);
		printf("NumAssignedLong: %d\n",numAssignedLong);

		if (filterAssStr != "" || filterAssStr2 != "")
		{
			printf("NumAssignedSR_prev %d  NumAssignedLR_prev %d\n",numAssignedSR_prev, numAssignedLR_prev);

			initializeAssignedFilter2(numAssignedLong, *assignments, assigned,
					numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);
			deletedAssignments.clear();
			reComputeSym_Inter_Net_Ambig_Score(*assignments,assigned);
			if (filterAss2 != NULL)
			{
				filterAss2->reset(); // in case this is not the first iteration
				filter(*filterAss2,assigned,false, deletedAssignments);
			}
			else if (filterAss != NULL)
			{
				filterAss->reset();
				filter(*filterAss,assigned,false, deletedAssignments);
			}

			numAssignedSR_prev = 0;
			numAssignedLR_prev = 0;

			for (vector<PCAssignment*>::iterator it = assigned.begin(); it != assigned.end(); ++it)
			{
				PCAssignment* pca = *it;
				int seqsep = pca->pc.getSeqSep();
				if (seqsep < LRSEQSEP)
					numAssignedSR_prev++;
				else
					numAssignedLR_prev++;
			}
		}

		bestAssigned.assign(assigned.begin(),assigned.end());

		// compute ambiguous assignments
		tr1::unordered_map< NOE*, vector<PCAssignment*> > ambigAss;

		// each vector entry contains the restraints for the ambiguous assignments of a specific noe,
		// each of which can have many assignment possibilities {[atom pair] -> dist, lb, ub, stdev, score of assignment possibility}
		printf("NumAssigned before compute ambiguous assignments: %zu\n",bestAssigned.size());
		computeAmbiguousAssignments(bestAssigned, *assignments, ambigAss, deletedAssignments,
				numAssignableContactsSR, numAssignableContactsLR, numAssignedSR_prev, numAssignedLR_prev);

		int countAmbigSize = 0;
		int maxAmbigSize = 0;

		for (tr1::unordered_map< NOE*, vector<PCAssignment*> >::iterator it = ambigAss.begin(); it != ambigAss.end(); ++it)
		{
			int s = it->second.size();
			countAmbigSize += s;
			if (s > maxAmbigSize)
				maxAmbigSize = s;
		}
		double avgAmbigSize = (double)countAmbigSize/(double)(ambigAss.size());
		printf("AvgAmbigSize: %7.5f  NumAmbig: %zu  MaxAmbigSize: %d  NumSeedsFilteredOut: %zu\n",
				avgAmbigSize,ambigAss.size(),maxAmbigSize, deletedAssignments.size());

		// prepare to output assignments to files
		vector<PCAssignment*> assignedNoDup; // symmetric peak assignments are merged
		assignedNoDup.reserve(bestAssigned.size());
		mergeAssignments(bestAssigned, assignedNoDup); // must clean up contents of assignedNoDup later
		sort(assignedNoDup.begin(),assignedNoDup.end(),PCAssignmentSortByScore()); // ascending order
		writeAssignmentsXplor("noe_refine.tbl","new.tbl", assignedNoDup, ambigAss, calibConstants);
		writeITASSER("noe.tbl", "ambig.tbl", assignedNoDup, ambigAss, calibConstants);
		accuracyReport(assignedNoDup);
		for (vector<PCAssignment*>::iterator it = assignedNoDup.begin(); it != assignedNoDup.end(); ++it)
		{
			delete *it;
		}
	}
	catch (IloException& e)
	{
		cerr << "Concert exception caught: " << e << endl;
	}
	catch (...)
	{
		cerr << "Unknown exception caught" << endl;
	}
	env.end();

	// clean up local vars
	delete pc2NOEs;

	for (list<PCAssignment*>::iterator it = assignments->begin(); it != assignments->end(); ++it)
	{
		delete *it;
	}
	delete assignments;

	// clean up filters
	if (filterAssPos != NULL)
		delete filterAssPos;
	if (filterAss != NULL)
		delete filterAss;
	if (filterAmbPos != NULL)
		delete filterAmbPos;
}

void cleanup()
{
	delete bmrbProtein;
	for (list<CSProtein*>::iterator it = structures.begin(); it != structures.end(); it++)
	{
		CSProtein* p = *it;
		delete p;
	}
	for (list<CSProtein*>::iterator it = refPDB.begin(); it != refPDB.end(); it++)
	{
		CSProtein* p = *it;
		delete p;
	}
	// delete NOEs
	for (vector<NOE::PEAKLISTTYPE>::iterator it = peakListTypes.begin(); it != peakListTypes.end(); ++it)
	{
		NOE::PEAKLISTTYPE type = *it;
		list<NOE*>& nList = peakLists[type];
		for (list<NOE*>::iterator it2 = nList.begin(); it2 != nList.end(); ++it2)
			delete *it2;
	}
}

// true if file is properly formatted
bool parseConfigFile(const string& configFile)
{
	ifstream file;
	string line;
	file.open(configFile.c_str());
	string filenameLong = "/nfs/amino-home/rjang/bin/longrange.txt"; // default path, can change in config file
	string filenameSeq = "/nfs/amino-home/rjang/bin/sequential.txt";
	string filenameIntra = "/nfs/amino-home/rjang/bin/intra.txt";
	string resContactStatsFile = "/nfs/amino-home/rjang/bin/residueContactStats.txt";

	if (file.is_open())
	{
		while (file.good())
		{
			getline(file,line);
			trim(line);
			if (line.length() == 0)
				continue;
			stringstream tok(line);
			string tempKey;
			tok >> tempKey;
			string key;
			bool isComment = false;
			for (string::iterator i=tempKey.begin(); i != tempKey.end(); ++i)
			{
				char c = *i;
				if (isalpha(c))
					c = tolower(c);
				else
				{
					if (c == '#')
					{
						isComment = true;
						break;
					}
				}
				key.push_back(c);
			}
			if (isComment)
				continue;
			else if (key == "fasta" || key == "bmrb" || key == "template" || key == "ref" ||
					 key == "cali" || key == "caro" || key == "n" || key == "h2d" ||
					 key == "chch" || key == "chnh" || key == "nhnh" ||
					 key == "distcutoff" || key == "weights" ||
					 key == "contactstats" ||
					 key == "refine_fractstrcut" || key == "refine_netcut" || key == "refine_fracttempcut" || key == "wt_strdist" ||
					 key == "template_distcut" ||
					 key == "scscstats" ||
					 key == "ambig_score_cut" ||
					 key == "weak_restraint_wt" ||
					 key == "weighttrainfile" ||
					 key == "initialqid" ||
					 key == "tempfracfilter" ||
					 key == "tempfraceval" ||
					 key == "init_fractstrcut" ||
					 key == "filterinter" ||
					 key == "filternet" ||
					 key == "filterinterposs" ||
					 key == "filternetposs" ||
					 key == "longdistcutfrac" ||
					 key == "filterasspos" ||
					 key == "filterass" ||
					 key == "filterass2" ||
					 key == "filterambpos" ||
					 key == "max_calib_error" ||
					 key == "ambig_filter" ||
					 key == "pool_filter" ||
					 key == "maxiter" ||
					 key == "oldwt")
			{
				vector<string> tokens;
				while (tok.good())
				{
					string temp;
					tok >> temp;
					if (temp.at(0) == '\"')
					{
						string token;
						token.append(temp.substr(1));
						// keeping reading until encounter closing "
						while (tok.good())
						{
							tok >> temp;
							int len = temp.length();
							token.append(" ");
							if (temp[len-1] == '\"')
							{
								token.append(temp.substr(0,len-1));
								break;
							}
							else
								token.append(temp);
						}
						tokens.push_back(token);
					}
					else
					{
						if (temp.at(0) == '#') // comment
							break;
						tokens.push_back(temp);
					}
				}
				if (tokens.size() == 0)
				{
					printf("Not enough arguments for %s\n",key.c_str());
					continue;
				}
				if (key == "fasta")
				{
					string fasta = tokens[0]; // first residue considered as residue #1
					readFasta(fasta,seq);
					if (tokens.size() > 1)
					{
						startRes = atoi(tokens[1].c_str());
						endRes = atoi(tokens[2].c_str());
					}
					else
					{
						startRes = 1;
						endRes = seq.length();
					}
					//printf("%s\n",seq.c_str());
					printf("Seq StartRes: %d   EndRes: %d\n",startRes,endRes);
				}
				else if (key == "bmrb")
				{
					string bmrb = tokens[0];
					int bmrbOffset = 0;
					bool bmrbOffsetFlag = false; // true of offset provided
					if (tokens.size() > 1)
					{
						bmrbOffset = atoi(tokens[1].c_str());
						bmrbOffsetFlag = true;
					}
					bmrbProtein = new CSProtein(seq,1);
					bmrbProtein->getShifts(bmrb,bmrbOffset);

					for (int r = 1; r <= bmrbProtein->size; r++)
					{
						Residue* res = (*bmrbProtein)[r];
						seqAAType[r-1] = (Residue::AA3)res->type;
					}

					// bmrbProtein->write15NChemShift("cs_15n.str",bmrbOffset,true,true);

					string talos;
					int talosOffset = 0;
					unsigned int index = 2;
					if (!bmrbOffsetFlag)
						index = 1;
					if (tokens.size() > index)
					{
						talos = tokens[index];
						index++;
						if (tokens.size() > index)
						{
							talosOffset = atoi(tokens[index].c_str());
							index++;
						}
					}
					else
						talos = "pred.tab";
					if (FILE *file = fopen(talos.c_str(), "r"))
					{
						fclose(file);
					}
					else
					{
						printf("TALOS file %s does not exist\n",talos.c_str());
						return false;
					}
					for (int i = 0; i < MAXPROSIZE; i++)
						isDynamic.set(i,true); // initially everything is dynamic
					ifstream file;
					string line;
					file.open(talos.c_str());
					if (file.is_open())
					{
						while (file.good())
						{
							getline(file,line);
							trim(line);
							if (line.empty())
								continue;
							stringstream tok(line);
							string temp;
							tok >> temp;
							if (temp == "FORMAT")
							{
								while (file.good())
								{
									getline(file,line);
									trim(line);
									if (line.empty())
										continue;
									stringstream tok(line);
									string temp;
									tok >> temp;
									int resnum = talosOffset+atoi(temp.c_str());
									if (resnum < 1)
										continue;
									for (int i = 0; i < 9; i++)
										tok >> temp; // skip all columns except last column
									tok >> temp;
									// assume all res dynamic unless specified otherwise
									if (temp != "Dyn" && temp != "None")
									{
										isDynamic.set(resnum-1,false);
									}
								}
								break; // done reading file
							}
						} // end while file good
						file.close();
					} // end if file is open

					string talosSS; // assumed to have the same talosOffset
					if (tokens.size() > index)
					{
						talosSS = tokens[index];
						index++;
					}
					else
					{
						talosSS = "predSS.tab";
					}
					if (FILE *file = fopen(talosSS.c_str(), "r"))
					{
						fclose(file);
					}
					else
					{
						printf("TALOS SS file %s does not exist\n",talosSS.c_str());
						return false;
					}
					file.open(talosSS.c_str());
					for (int i = 0; i < MAXPROSIZE; i++)
						ssTypes[i] = UNKNOWNSS;
					if (file.is_open())
					{
						while (file.good())
						{
							getline(file,line);
							trim(line);
							if (line.empty())
								continue;
							stringstream tok(line);
							string temp;
							tok >> temp;
							if (temp == "FORMAT")
							{
								while (file.good())
								{
									getline(file,line);
									trim(line);
									if (line.empty())
										continue;
									stringstream tok(line);
									string temp;
									tok >> temp;
									int resnum = talosOffset+atoi(temp.c_str());
									if (resnum < 1)
										continue;
									for (int i = 0; i < 6; i++)
										tok >> temp; // skip all remaining columns except last two columns
									tok >> temp;
									double confidence = atof(temp.c_str());
									if (confidence > 0)
									{
										string ssType;
										tok >> ssType;
										if (ssType == "H")
											ssTypes[resnum-1] = HELIX;
										else if (ssType == "E")
											ssTypes[resnum-1] = STRAND;
										else
											ssTypes[resnum-1] = LOOP;
									}
									else
									{
										ssTypes[resnum-1] = UNKNOWNSS;
									}
								}
								break; // done reading file
							}
						} // end while file good
						file.close();
					} // end if file is open
					printf("BMRB: %s   BMRB_Offset: %d   TALOS: %s   TALOS_Offset: %d   TALOS_SS: %s\n",
							bmrb.c_str(),bmrbOffset,talos.c_str(),talosOffset,talosSS.c_str());
//					for (int i = 0; i < MAXPROSIZE; i++)
//					{
//						if (!isDynamic.test(i))
//							printf("%d\n",(i+1));
//						if (ssTypes[i] != UNKNOWNSS)
//						{
//							if (ssTypes[i] == HELIX)
//								printf("%d H\n",(i+1));
//							else if (ssTypes[i] == STRAND)
//								printf("%d E\n",(i+1));
//							else
//								printf("%d L\n",(i+1));
//						}
//					}
				}
				else if (key == "template")
				{
					string pdb = tokens[0]; // may contain multiple models
					int offset = 0;
					bool offsetFlag = false; // true if offset provided
					if (tokens.size() > 1)
					{
						offset = atoi(tokens[1].c_str());
						offsetFlag = true;
					}
					list<CSProtein*> pdbs;
					if (seq.size() == 0)
					{
						printf("ERROR: Please specify sequence first in config file before the templates\n");
						return false;
					}
					CSProtein::getModels(pdb,pdbs,seq,1,offset);  // offset used to adjust for pdb files that do not start at 1
					for (list<CSProtein*>::iterator it = pdbs.begin(); it != pdbs.end(); ++it)
					{
						CSProtein* cs = *it;
						if (cs->getNumProtons() < 1)
						{
							printf("WARNING: %s contains no protons. Ignoring structure\n",pdb.c_str());
						}
						else
							structures.push_back(*it);
					}
					string stride;
					unsigned int index = 2;
					if (!offsetFlag)
						index = 1;
					if (tokens.size() > index)
					{
						stride = tokens[index];
					}
					else
					{
						stride = pdb.substr(0,pdb.size()-3); // chop off pdb file extension
						stride.append("stride");
					}
					if (FILE *file = fopen(stride.c_str(), "r"))
					{
						fclose(file);
					}
					else
					{
						printf("Stride file %s does not exist\n",stride.c_str());
						return false;
					}
					vector< vector<int> > sse;
					parseStride(stride, sse, offset);
					sses.push_back(sse);

					printf("PDB: %s Offset: %d  Stride: %s\n",pdb.c_str(),offset,stride.c_str());
				}
				else if (key == "ref")
				{
					if (seq.size() == 0)
					{
						printf("ERROR: Please specify sequence first in config file before ref pdb file\n");
						return false;
					}
					string pdb = tokens[0];
					int offset = 0;
					if (tokens.size() > 1)
						offset = atoi(tokens[1].c_str());
					CSProtein::getModels(pdb,refPDB,seq,1,offset);

					printf("RefPDB: %s  Offset: %d\n",pdb.c_str(),offset);
				}
				else if ( (key == "cali") || (key == "caro") || (key == "n") )
				{
					string peakListFile = tokens[0];
					if (tokens.size() < 5)
					{
						printf("ERROR: Peak list format. Missing column indices");
						return false;
					}
					int xIndex = atoi(tokens[1].c_str()); // 1-based index
					int hxIndex = atoi(tokens[2].c_str());
					int hIndex = atoi(tokens[3].c_str());
					int volIndex = atoi(tokens[4].c_str());

					NOE::PEAKLISTTYPE type = NOE::CALI_3D;
					if (key == "caro")
						type = NOE::CARO_3D;
					else if (key == "n")
						type = NOE::N15_3D;
					peakListTypes.push_back(type);

					if (tokens.size() > 5)
					{
						tr1::array<double,4> tol;
						tol[0] = atof(tokens[5].c_str()); // x
						tol[1] = atof(tokens[6].c_str()); // hx
						tol[2] = INVALIDSHIFTVAL;
						tol[3] = atof(tokens[7].c_str()); // h
						tols[type] = tol;
					}
					else
					{
						// defaults
						tr1::array<double,4> tol;
						tol[0] = tolx;
						tol[1] = tolhx;
						tol[2] = INVALIDSHIFTVAL;
						tol[3] = tolh;
						tols[type] = tol;
					}
					list<NOE*> noes;
					NOE::readNOE3D(peakListFile,xIndex,hxIndex,hIndex,volIndex,noes,type);
					peakLists[type] = noes;

					printf("PeakList: %s %s %f %f %f\n",
							key.c_str(),peakListFile.c_str(),tols[type][0],tols[type][1],tols[type][3]);
				}
				else if ( (key == "chch") || (key == "chnh") || (key == "nhnh") )
				{
					string peakListFile = tokens[0];
					if (tokens.size() < 6)
					{
						printf("ERROR: Peak list format. Missing column indices");
						return false;
					}
					printf("WARNING: 4D peak lists not yet implemented\n");
					// TODO: 4D peak list parsing
					continue;
				}
				else if (key == "h2d")
				{
					string peakListFile = tokens[0];
					if (tokens.size() < 4)
					{
						printf("ERROR: Peak list format. Missing column indices");
						return false;
					}
					printf("WARNING: 2D peak lists not yet implemented\n");
					// TODO: 2D peak list parsing
					continue;
				}
				else if (key == "distcutoff")
				{
					DISTCUTOFF = atof(tokens[0].c_str());
					printf("DISTCUTOFF: %f\n",DISTCUTOFF);
				}
				else if (key == "weights")
				{
					if (tokens.size() != 9)
					{
						printf("Weights format: cs str intensity sym interres net netstr ambig db\n");
						return false;
					}
					Score::CS_WT =  atof(tokens[0].c_str());
					Score::STR_WT =  atof(tokens[1].c_str());
					Score::INTENSITY_WT =  atof(tokens[2].c_str());
					Score::SYM_WT =  atof(tokens[3].c_str());
					Score::INTERRES_WT =  atof(tokens[4].c_str());
					Score::NET_WT =  atof(tokens[5].c_str());
					Score::NETSTR_WT =  atof(tokens[6].c_str());
					Score::AMBIG_WT =  atof(tokens[7].c_str());
					Score::DB_WT =  atof(tokens[8].c_str());

					printf("CS=%6.3f  STR=%6.3f  INTENSITY=%6.3f  SYM=%6.3f  INTERRES=%6.3f  NET=%6.3f  NETSTR=%6.3f  AMBIG=%6.3f  DB=%6.3f\n",
							Score::CS_WT,Score::STR_WT,Score::INTENSITY_WT,Score::SYM_WT,Score::INTERRES_WT,Score::NET_WT,Score::NETSTR_WT,
							Score::AMBIG_WT,Score::DB_WT);

					double sum = Score::CS_WT+Score::STR_WT+Score::INTENSITY_WT+Score::SYM_WT+Score::INTERRES_WT+Score::NET_WT+Score::NETSTR_WT+Score::AMBIG_WT+Score::DB_WT;
					if (abs(sum-1.0) > 0.01)
					{
						printf("Scaling weights so that the sum is 1\n");
						Score::CS_WT /= sum;
						Score::STR_WT /= sum;
						Score::INTENSITY_WT /= sum;
						Score::SYM_WT /= sum;
						Score::INTERRES_WT /= sum;
						Score::NET_WT /= sum;
						Score::NETSTR_WT /= sum;
						Score::AMBIG_WT /= sum;
						Score::DB_WT /= sum;
						printf("CS=%6.3f  STR=%6.3f  INTENSITY=%6.3f  SYM=%6.3f  INTERRES=%6.3f  NET=%6.3f  NETSTR=%6.3f  AMBIG=%6.3f  DB=%6.3f\n",
								Score::CS_WT,Score::STR_WT,Score::INTENSITY_WT,Score::SYM_WT,Score::INTERRES_WT,Score::NET_WT,Score::NETSTR_WT,
								Score::AMBIG_WT,Score::DB_WT);
					}
				}
				else if (key == "contactstats")
				{
					if (tokens.size() != 3)
					{
						printf("contactstats format: longrange.txt sequential.txt intra.txt\n");
						return false;
					}
					filenameLong = tokens[0];
					filenameSeq = tokens[1];
					filenameIntra = tokens[2];
				} // end if contactstats
				else if (key == "refine_fractstrcut")
				{
					REFINE_FRACTSTRCUT = atof(tokens[0].c_str());
					printf("REFINE_FRACTSTRCUT %f\n",REFINE_FRACTSTRCUT);
				}
				else if (key == "refine_netcut")
				{
					int netcutFlag = atoi(tokens[0].c_str());
					if (netcutFlag > 0)
					{
						REFINE_NETCUT = true;
						printf("REFINE_NETCUT true\n");
					}
					else
					{
						REFINE_NETCUT = false;
						printf("REFINE_NETCUT false\n");
					}
				}
				else if (key == "refine_fracttempcut")
				{
					REFINE_FRACTTEMPCUT = atof(tokens[0].c_str());
					printf("REFINE_FRACTTEMPCUT %f\n",REFINE_FRACTTEMPCUT);
				}
				else if (key == "wt_strdist")
				{
					WT_STRDIST = atof(tokens[0].c_str());
					WT_INTDIST = 1.0 - WT_STRDIST;
					printf("WT_STRDIST %f  WT_INTDIST %f\n",WT_STRDIST,WT_INTDIST);
				}
				else if (key == "template_distcut")
				{
					TEMPLATE_DISTCUT = atof(tokens[0].c_str());
					printf("TEMPLATE_DISTCUT %f\n",TEMPLATE_DISTCUT);
				}
				else if (key == "scscstats")
				{
					resContactStatsFile = tokens[0];
					if (tokens.size() == 3)
					{
						SCSC_STDEV_FACTOR = atof(tokens[1].c_str());
						SCSC_MIN_STDEV = atof(tokens[2].c_str());
						printf("SCSC_STDEV_FACTOR: %f   SCSC_MIN_STDEV: %f\n",SCSC_STDEV_FACTOR,SCSC_MIN_STDEV);
					}
				}
				else if (key == "ambig_score_cut")
				{
					AMBIG_SCORE_CUT = atof(tokens[0].c_str());
					printf("AMBIG_SCORE_CUT: %f\n",AMBIG_SCORE_CUT);
				}
				else if (key == "weak_restraint_wt")
				{
					WEAK_RESTRAINT_WT = atof(tokens[0].c_str());
					printf("WEAK_RESTRAINT_WT: %f\n",WEAK_RESTRAINT_WT);
				}
				else if (key == "weighttrainfile")
				{
					weightTrainFile = tokens[0];
					printf("Weight Training Mode. Output to file: %s\n",weightTrainFile.c_str());
				}
				else if (key == "initialqid")
				{
					initialQID = atoi(tokens[0].c_str());
					printf("Initial QID: %d\n",initialQID);
				}
				else if (key == "tempfracfilter")
				{
					TEMPFRACFILTER = atof(tokens[0].c_str());
					printf("TEMPFRACFILTER: %f\n",TEMPFRACFILTER);
				}
				else if (key == "tempfraceval")
				{
					TEMPFRACEVAL = atof(tokens[0].c_str());
					printf("TEMPFRACEVAL: %f\n",TEMPFRACEVAL);
				}
				else if (key == "init_fractstrcut")
				{
					INIT_FRACTSTRCUT = atof(tokens[0].c_str());
					printf("INIT_FRACTSTRCUT: %f\n",INIT_FRACTSTRCUT);
				}
				else if (key == "filterinter")
				{
					FILTERINTER = atof(tokens[0].c_str());
					printf("FILTERINTER: %f\n",FILTERINTER);
				}
				else if (key == "filternet")
				{
					FILTERNET = atof(tokens[0].c_str());
					printf("FILTERNET: %f\n",FILTERNET);
				}
				else if (key == "filterinterposs")
				{
					FILTERINTERPOSS = atof(tokens[0].c_str());
					printf("FILTERINTERPOSS: %f\n",FILTERINTERPOSS);
				}
				else if (key == "filternetposs")
				{
					FILTERNETPOSS = atof(tokens[0].c_str());
					printf("FILTERNETPOSS: %f\n",FILTERNETPOSS);
				}
				else if ( key == "longdistcutfrac")
				{
					LONGDISTCUTFRAC = atof(tokens[0].c_str());
					printf("LONGDISTCUTFRAC: %f\n",LONGDISTCUTFRAC);
				}
				else if (key == "filterasspos")
				{
					filterAssPosStr = line.substr(13);
					trim(filterAssPosStr);
				}
				else if (key == "filterass")
				{
					filterAssStr = line.substr(10);
					trim(filterAssStr);
				}
				else if (key == "filterass2")
				{
					filterAssStr2 = line.substr(10);
					trim(filterAssStr2);
				}
				else if (key == "filterambpos")
				{
					filterAmbPosStr = line.substr(13);
					trim(filterAmbPosStr);
				}
				else if (key == "max_calib_error")
				{
					MAX_CALIB_ERROR = atof(tokens[0].c_str());
					printf("MAX_CALIB_ERROR: %f\n",MAX_CALIB_ERROR);
				}
				else if (key == "ambig_filter")
				{
					AMBIG_FILTER = atof(tokens[0].c_str());
					if (tokens.size() > 1)
							AMBIG_FILTER_MAX = atoi(tokens[1].c_str());
					if (AMBIG_FILTER < 0)
					{
						AMBIG_FILTER = -AMBIG_FILTER;
						AMBIG_FILTER_BY_COUNT = true;
						printf("AMBIG_FILTER (by count): %f\n",AMBIG_FILTER);
					}
					else
						printf("AMBIG_FILTER (by score): %f %d\n",AMBIG_FILTER,AMBIG_FILTER_MAX);
				}
				else if (key == "pool_filter")
				{
					POOL_FILTER = atof(tokens[0].c_str());
					if (POOL_FILTER < 0)
					{
						POOL_FILTER = -POOL_FILTER;
						POOL_FILTER_BY_COUNT = true;
						printf("POOL_FILTER (by count): %f\n",POOL_FILTER);
					}
					else
						printf("POOL_FILTER (by score): %f\n",POOL_FILTER);
				}
				else if (key == "maxiter")
				{
					MAXITER = atoi(tokens[0].c_str());
					printf("MAXITER: %d\n",MAXITER);
				}
				else if (key == "oldwt")
				{
					OLDWT = atof(tokens[0].c_str());
					if (OLDWT > 1.0 || OLDWT < 0)
					{
						printf("ERROR INVALID VALUE FOR OLDWT: %f\n",OLDWT);
						exit(0);
					}
					NEWWT = 1.0-OLDWT;
					printf("OLDWT %f  NEWWT %f\n",OLDWT,NEWWT);
				}
			}
			else
			{
				// unknown tag
				printf("WARNING: Unknown tag %s\n",key.c_str());
				continue;
			}
		} // end while file
		// check if have minimum input
		if (seq.size() == 0)
		{
			printf("ERROR: No sequence\n");
			return false;
		}
		if (structures.size() == 0)
		{
			printf("ERROR: No structures in input\n");
			return false;
		}
		if (bmrbProtein == NULL)
		{
			printf("ERROR: No chemical shift assignment\n");
			return false;
		}
		if (peakListTypes.size() == 0)
		{
			printf("ERROR: No peak lists\n");
			return false;
		}
	}
	else
	{
		printf("Unable to open file %s\n",configFile.c_str());
		file.close();
		return false;
	}
	file.close();

	// setup expectedAssign and expectedAssignAvg and scscAvgDist
	// for expectedAssignMatrix
	tr1::unordered_map< ExpectedAssignKey, int > expectedNumLong;
	tr1::unordered_map< ExpectedAssignKey, int > expectedNumSeq;
	tr1::unordered_map< ExpectedAssignKey, int > expectedNumIntra;
	tr1::unordered_map< SCSCDistKey, double > scscTemp;
	ifstream fileLong,fileSeq,fileIntra,fileSCSC;
	fileLong.open(filenameLong.c_str());
	fileSeq.open(filenameSeq.c_str());
	fileIntra.open(filenameIntra.c_str());
	fileSCSC.open(resContactStatsFile.c_str());
	if (!fileLong.is_open())
	{
		printf("ERROR: file %s not found\n",filenameLong.c_str());
		return false;
	}
	if (!fileSeq.is_open())
	{
		printf("ERROR: file %s not found\n",filenameSeq.c_str());
		return false;
	}
	if (!fileIntra.is_open())
	{
		printf("ERROR: file %s not found\n",filenameIntra.c_str());
		return false;
	}
	if (!fileSCSC.is_open())
	{
		printf("ERROR: file %s not found\n",resContactStatsFile.c_str());
		return false;
	}
	while (fileLong.good())
	{
		getline(fileLong,line);
	    trim(line);
		if (line.empty())
			continue;
	    stringstream tok(line);
	    string type1,type2,proton1,proton2,temp;
	    int count;
	    tok >> type1;
	    tok >> type2;
	    tok >> proton1;
	    tok >> proton2;
	    tok >> temp;
	    count = (int)(round(atof(temp.c_str())));
	    int t1 = (int)(Residue::index3(type1.c_str()));
	    int t2 = (int)(Residue::index3(type2.c_str()));

	    ExpectedAssignKey k1(t1,proton1,t2,proton2);
	    ExpectedAssignKey k2(t2,proton2,t1,proton1);
	    expectedNumLong[k1] = count;
	    if (k1 != k2)
	    	expectedNumLong[k2] = count;
	}
	while (fileSeq.good())
	{
		getline(fileSeq,line);
	    trim(line);
		if (line.empty())
			continue;
	    stringstream tok(line);
	    string type1,type2,proton1,proton2,temp;
	    int count;
	    tok >> type1;
	    tok >> type2;
	    tok >> proton1;
	    tok >> proton2;
	    tok >> temp;
	    count = (int)(round(atof(temp.c_str())));
	    int t1 = (int)(Residue::index3(type1.c_str()));
	    int t2 = (int)(Residue::index3(type2.c_str()));

	    ExpectedAssignKey k1(t1,proton1,t2,proton2);
	    ExpectedAssignKey k2(t2,proton2,t1,proton1);
	    expectedNumSeq[k1] = count;
	    if (k1 != k2)
	    	expectedNumSeq[k2] = count;
	}
	while (fileIntra.good())
	{
		getline(fileIntra,line);
	    trim(line);
		if (line.empty())
			continue;
	    stringstream tok(line);
	    string type1,type2,proton1,proton2,temp;
	    int count;
	    tok >> type1;
	    tok >> type2;
	    tok >> proton1;
	    tok >> proton2;
	    tok >> temp;
	    count = (int)(round(atof(temp.c_str())));

	    int t1 = (int)(Residue::index3(type1.c_str()));
	    int t2 = (int)(Residue::index3(type2.c_str()));

	    ExpectedAssignKey k1(t1,proton1,t2,proton2);
	    ExpectedAssignKey k2(t2,proton2,t1,proton1);
	    expectedNumIntra[k1] = count;
	    if (k1 != k2)
	    	expectedNumIntra[k2] = count;
	}
	fileLong.close();
	fileSeq.close();
	fileIntra.close();
	while (fileSCSC.good())
	{
		getline(fileSCSC,line);
		trim(line);
		if (line.empty())
			continue;
	    stringstream tok(line);
	    string type1,type2,proton1,proton2,temp;
	    int bin;
	    double dist;
	    double stdev;
	    tok >> type1;
	    tok >> type2;
	    int t1 = (int)(Residue::index3(type1.c_str()));
	    int t2 = (int)(Residue::index3(type2.c_str()));
	    tok >> proton1;
	    tok >> proton2;
	    tok >> temp;
	    bin = atoi(temp.c_str());
	    tok >> temp;
	    dist = atof(temp.c_str());
	    tok >> temp; // count
	    tok >> temp; // stdev
	    stdev = atof(temp.c_str());
	    // set upper distance bound to avg +1.5*stdev
	    if (stdev < SCSC_MIN_STDEV)
	    	dist = dist + SCSC_MIN_STDEV*SCSC_STDEV_FACTOR;
	    else
	    	dist = dist + SCSC_STDEV_FACTOR*stdev;
	    SCSCDistKey k1(t1,proton1,t2,proton2,bin);
	    SCSCDistKey k2(t2,proton2,t1,proton1,bin);
	    scscTemp[k1] = dist;
	    if (k1 != k2)
	    	scscTemp[k2] = dist;
	    // printf("%s %s %s %s %d %f %f\n",type1.c_str(),type2.c_str(),proton1.c_str(),proton2.c_str(),bin,dist,stdev);
	}
	fileSCSC.close();

	// initialize expectedAssign, expectedAssignAvg, scscstats
	int expectedAssignAvgCounts[MAXPROSIZE][MAXPROSIZE]; // [resnum-1][resnum-1]
	for (int i = 0; i < MAXPROSIZE; i++)
	{
		for (int j = 0; j < MAXPROSIZE; j++)
		{
			expectedAssignAvg[i][j] = 0;
			expectedAssignAvgCounts[i][j] = 0;
		}
	}
	bool has15N = false; // true if peak list type exists
	bool hasCali = false;
	bool hasCaro = false;
	for (vector<NOE::PEAKLISTTYPE>::iterator it = peakListTypes.begin(); it != peakListTypes.end(); ++it)
	{
		NOE::PEAKLISTTYPE type = *it;
		if (type == NOE::N15_3D)
			has15N = true;
		if (type == NOE::CALI_3D)
			hasCali = true;
		if (type == NOE::CARO_3D)
			hasCaro = true;
	}
	CSProtein& csa = *bmrbProtein;
	// compute correction value
	for (int r1 = 1; r1 <= csa.size; r1++)
	{
		Residue* res1 = csa[r1];
		for (int r2 = r1; r2 <= csa.size; r2++)
		{
			Residue* res2 = csa[r2];
			int seqSep = r2-r1;
			int correction = 0; // to determine # of contacts expected between r1,r2, need to take into account cs completeness and peak list type
			              // the correction is subtracted from the expected count(r1,r2,hXAtom1,hXAtom2)
			for (AtomIterator itN1 = res1->begin('X'); itN1 != res1->end('X'); itN1++) // 'X' to check both C and N
			{
				Atom* xAtom1 = *itN1;
				if (xAtom1->numProtons < 1)
					continue;
				bool isMethylx1 = xAtom1->isMethyl();
				Atom* hXAtom1 = NULL;
				for (HIterator ithX1 = xAtom1->beginH(); ithX1 != xAtom1->endH(); ithX1++)
				{
					if (hXAtom1 != NULL && isMethylx1)
						break; // count methyl groups only once
					hXAtom1 = *ithX1;
					for (AtomIterator itN2 = res2->begin('X'); itN2 != res2->end('X'); itN2++) // 'X' to check both C and N
					{
						Atom* xAtom2 = *itN2;
						if (xAtom2->numProtons < 1)
							continue;
						bool isMethylx2 = xAtom2->isMethyl();
						Atom* hXAtom2 = NULL;
						for (HIterator ithX2 = xAtom2->beginH(); ithX2 != xAtom2->endH(); ithX2++)
						{
							if (hXAtom2 != NULL && isMethylx2)
								break; // count methyl groups only once
							hXAtom2 = *ithX2;
							if (r1 == r2 && hXAtom1 >= hXAtom2)  // don't double count intraresidue contacts
								continue;

							// setup scscAvgDist
							double defaultDist=0;// take this to be upperbound if no stats available
							for (int b = 0; b < NUMDISTBINS; b++)
							{
								SCSCDistKey k((int)res1->type,hXAtom1->name,(int)res2->type,hXAtom2->name,b);
								tr1::unordered_map<SCSCDistKey,double>::iterator itFind = scscTemp.find(k);
								if (itFind != scscTemp.end())
								{
									double dist = itFind->second;
									if (dist > defaultDist)
										defaultDist = dist;
								}
							}
							if (defaultDist==0 && r1 != r2)
							{
								defaultDist = 12.0;
							}
							for (int b = 0; b < NUMDISTBINS; b++)
							{
								// distance is 0 if r1 == r2
								double dist = defaultDist;
								if (r1 != r2)
								{
									SCSCDistKey k((int)res1->type,hXAtom1->name,(int)res2->type,hXAtom2->name,b);
									tr1::unordered_map<SCSCDistKey,double>::iterator itFind = scscTemp.find(k);
									if (itFind != scscTemp.end())
									{
										dist = itFind->second;
									}
								}
								else
								{
									dist = 0;
								}
								SCSCDistKey k(r1,hXAtom1->name,r2,hXAtom2->name,b);
								SCSCDistKey k2(r2,hXAtom2->name,r1,hXAtom1->name,b);
								scscAvgDist[k] = dist;
								scscAvgDist[k2] = dist;
							}
							// end setup scscAvgDist

							// correction for completeness, peak list type
							// consider both directions of the contact, so for each contact, correction can
							// increment by at most 2
							if ( !hXAtom1->hasChemShift() ||
									((res1->is15N(xAtom1->name) && !has15N) ||
									 (res1->isCali(xAtom1->name) && !hasCali) ||
									 (res1->isCaro(xAtom1->name) && !hasCaro)) )
							{
								correction++;
							}
							if ( !hXAtom2->hasChemShift() ||
									((res2->is15N(xAtom2->name) && !has15N) ||
									 (res2->isCali(xAtom2->name) && !hasCali) ||
									 (res2->isCaro(xAtom2->name) && !hasCaro)) )
							{
								correction++;
							}
						} // end for each hx2
					} // end for each x2
				} // end for each hx1
			} // end for each x1
			// use correction to update the expected counts
			for (AtomIterator itN1 = res1->begin('X'); itN1 != res1->end('X'); itN1++) // 'X' to check both C and N
			{
				Atom* xAtom1 = *itN1;
				if (xAtom1->numProtons < 1)
					continue;
				bool isMethylx1 = xAtom1->isMethyl();
				Atom* hXAtom1 = NULL;
				for (HIterator ithX1 = xAtom1->beginH(); ithX1 != xAtom1->endH(); ithX1++)
				{
					if (hXAtom1 != NULL && isMethylx1)
						break; // count methyl groups only once
					hXAtom1 = *ithX1;
					if (!hXAtom1->hasChemShift())
						continue;
					for (AtomIterator itN2 = res2->begin('X'); itN2 != res2->end('X'); itN2++) // 'X' to check both C and N
					{
						Atom* xAtom2 = *itN2;
						if (xAtom2->numProtons < 1)
							continue;
						bool isMethylx2 = xAtom2->isMethyl();
						Atom* hXAtom2 = NULL;
						for (HIterator ithX2 = xAtom2->beginH(); ithX2 != xAtom2->endH(); ithX2++)
						{
							if (hXAtom2 != NULL && isMethylx2)
								break; // count methyl groups only once
							hXAtom2 = *ithX2;
							if (!hXAtom2->hasChemShift())
								continue;
							if (r1 == r2 && hXAtom1 >= hXAtom2)  // don't double count intraresidue contacts
								continue;
							ExpectedAssignKey k((int)res1->type,hXAtom1->name,(int)res2->type,hXAtom2->name);
							int expected = 0;
							if (seqSep >= LRSEQSEP)
							{
								tr1::unordered_map<ExpectedAssignKey,int>::iterator itFind = expectedNumLong.find(k);
								if (itFind != expectedNumLong.end())
									expected = itFind->second;
								// else
								//	printf("WARNING no expected num long-range contact stats for %d %d %s %s\n",
								//			r1,r2,hXAtom1->name.c_str(),hXAtom2->name.c_str());
							}
							else if (seqSep > 0)
							{
								tr1::unordered_map<ExpectedAssignKey,int>::iterator itFind = expectedNumSeq.find(k);
								if (itFind != expectedNumSeq.end())
									expected = itFind->second;
								//else
								//	printf("WARNING no expected num sequential contact stats for %d %d %s %s\n",
								//			r1,r2,hXAtom1->name.c_str(),hXAtom2->name.c_str());
							}
							else
							{
								tr1::unordered_map<ExpectedAssignKey,int>::iterator itFind = expectedNumIntra.find(k);
								if (itFind != expectedNumIntra.end())
									expected = itFind->second;
								//else
								//	printf("WARNING no expected num intrares contact stats for %d %d %s %s\n",
								//			r1,r2,hXAtom1->name.c_str(),hXAtom2->name.c_str());
							}
							// int maxCount = max<int>(res1->getNumProtonsMethylOnce()*res2->getNumProtonsMethylOnce(),1);
							int count = expected; // don't use correction // -correction; // (int)round(((double)expected)*(1.0-(double)correction/(double)maxCount));
							if (count < 0)
								count = 0;
							//printf("%d %d %s %s %s %s %d cor=%d mC=%d\n",r1,r2,Residue::AA3S[res1->type],Residue::AA3S[res2->type],hXAtom1->name.c_str(),hXAtom2->name.c_str(),
							//		count,correction,maxCount);
							ExpectedAssignKey kr(r1,hXAtom1->name,r2,hXAtom2->name);
							ExpectedAssignKey kr2(r2,hXAtom2->name,r1,hXAtom1->name);
							expectedAssign[kr] = count;
							expectedAssign[kr2] = count;
							// printf("%d %d %s %s %s %s count=%d expected=%d maxCount=%d correction=%d\n",r1,r2,Residue::AA3S[res1->type],Residue::AA3S[res2->type],hXAtom1->name.c_str(),hXAtom2->name.c_str(),count,expected,maxCount,correction);
							if (count > 0)
							{
								expectedAssignAvg[r1-1][r2-1] += count;
								expectedAssignAvgCounts[r1-1][r2-1]++;
							}
						} // end for each hx2
					} // end for each x2
				} // end for each hx1
			} //end for each x1
			if (expectedAssignAvgCounts[r1-1][r2-1] > 0)
			{
				expectedAssignAvg[r1-1][r2-1] = expectedAssignAvg[r1-1][r2-1]/expectedAssignAvgCounts[r1-1][r2-1];
				expectedAssignAvg[r2-1][r1-1] = expectedAssignAvg[r1-1][r2-1];
			}
		} // end for each r2
	} // end for each r1

	return true;
}

int main(int argc, char** argv)
{
	if (argc == 2)
	{
		// std::srand(std::time(0));
		std::srand(13); // for testing purpose, use the same seeed
		string configFile(argv[1]);
		bool ok = parseConfigFile(configFile);
		if (ok)
		{
			assign();
		}
		cleanup();
		return 0;
	}
	else
	{
		printf("USAGE: ASSIGN-IT <config>\n");
		return 0;
	}
}
