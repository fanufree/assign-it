Compilation



Requires CPLEX library

See CPLEX Library File in email attachment. 
CPLEX installation instructions: https://notesblogblog.wordpress.com/2012/08/25/setting-up-cplex-for-academics-with-eclipse-on-linux/

In subdir.mk (in compile subdirectory), change LINE 146 to the correct library path of CPLEX
g++ -I/home/e4k2/Programs/ILOG/CPLEX124/concert/include -I/home/e4k2/Programs/ILOG/CPLEX124/cplex/include 

In makefile (in compile subdirectory), change LINE 46 to the correct library path of CPLEX
g++ -L/home/e4k2/Programs/ILOG/CPLEX124/concert/lib/x86-64_sles10_4.1/static_pic -L/home/e4k2/Programs/ILOG/CPLEX124/cplex/lib/x86-64_sles10_4.1/static_pic -o "ASSIGN-IT"

Type make to compile. If compiler complains about missing .d file, run make clean first.


Running

compile/ASSIGN-IT example/<config_file>

Outputs: 

noe.tbl and ambig.tbl (can be empty if no ambiguous depending on the level of filtering) are for I-TASSER. They contain the 1-1 and ambiguous assignments respectively.
noe_refine.tbl and new.tbl are for Xplor-NIH refinement. Xplor-NIH refinement is not part of NMR-I-TASSER. It is used only if the structure to be refined is known to be correct.


Config File Format

Please see example directory for example configuration file and format of input files
The format is based on a KEY, values format for each line.

Comments are specified with #. Everything after # on the current line are ignored

FASTA <sequence.fasta> <start_residue_number> <end_residue_number>

- FASTA is used to specify the amino acid sequence in FASTA format and the numbers that follow give
the starting and ending residues, which is useful for modeling protein subsequences. If no numbers are
provided, the sequence is assumed to start from 1 to the length of the sequence.


BMRB <chemical_shift_assignment_BMRB_2.1_format> <chemical_shift_residue_offset> <pred.tab> <pred_residue_offset> <predSS.tab>

- BMRB is used to specify the chemical shift assignments in BMRB format
- Make sure residue numbers in sequence matches the numbering in chemical shift assignment file
- <pred.tab> and <predSS.tab> are from the output of running TALOS+ http://spin.niddk.nih.gov/bax/software/TALOS/
- <chemical_shift_residue_offset> gets added to the residue numbers in the chemical shift assignment to get it to match the sequence
- <pred_residue_offset> gets added to the residue numbers in pred.tab and predSS.tab to get them to match the sequence


Peak lists

<PeakList_Type> <file> <hetero_col1> <proton_col1> <proton_col2> <intensity_col> <hetero_match_tolerance> <proton1_match_tolerance> <proton2_match_tolerance>

Indices are 1-based, first column is 1. They indicate the column number (1-based) of the heavy atom, covalently-attached
hydrogen, the other hydrogen, and peak intensity, respectively. The floating point numbers specify the
chemical shift match tolerances for the heavy atom, covalently-attached hydrogen, and the other hydrogen,
respectively. Only 3D peak lists are supported, and only N, CALI, and CARO types are supported

e.g. 

N           n.peaks    2 3 4 7 0.35 0.035 0.035
CALI        cali.peaks    2 3 4 7 0.35 0.035 0.035
CARO        caro.peaks    2 3 4 7 0.35 0.035 0.035

Weights on Score Terms

Default weight is 1.0. Weights are scaled such that they sum to 1.0.
e.g. 
#             cs   str     intensity     sym   interres      net  netstr     ambig      db
WEIGHTS      0.23   0.71       1.3      0.57       0.05      1.74    1.1      4.39    1.52


Input Structures

To specify input structure(s): 

TEMPLATE <pdb1> <offset1> <stride1>
TEMPLATE <pdb2> <offset2> <stride2>
. . .

- <pdb> PDB file MUST have protons added. To use protons, use HAAD http://zhanglab.ccmb.med.umich.edu/HAAD/ or REDUCE http://kinemage.biochem.duke.edu/software/reduce.php
- <offset> is the offset to add to the residue numbers to get them to match the sequence
- <stride> is the output of STRIDE http://webclu.bio.wzw.tum.de/stride/


Filters

Each filter returns true or false for a given assignment/assignment possibility (ass/asspos). True means keep the assignment, false means discard it.
FILTERASSPOS, FILTERASS, and FILTERAMBPOS are used to specify the filters for filtering the
assignment possibilities, the one-to-one assignments, and the ambiguous assignments, respectively.
Various filters are supported, which are described below. A filter returns true to keep the contact and
false to reject the contact for each assignment or assignment possibility. Filters can be combined using
boolean operations.

e.g. 

FILTERASSPOS OR SEQSEP 2 OR SCORECOUNT 1 2 1 0.001 6 0.001 SCORETERM 1 9 -2.5
FILTERASS IF NOISERATIO 100 OR INSTRUCTUREFRACT 10 6.0 0.05 -1.0 2 OR SEQSEP 1 OR AND OR SCORETERM 1 1 0.2 SCORETERM 1 6 0.2 SCORETERM 1 10441 0.3 SCORETERM 1 10441 0.5 ELSE SCORETERM 1 1 0.2
FILTERAMBPOS AND SCORETERM 1 9 -0.6 OR SEQSEP 2 OR SCORECOUNT 1 3 1 0.1 2 0.5 6 -0.025 SCORECOUNT 2 3 3 0.01 4 -0.025 5 -0.025

* Performance is affected by the filtering. Please try different filters to see which one is the best. 
Also see AMBIG_FILTER, which filters ambiguous assignments after filtering by FILTERAMBPOS.

INSTRUCTUREFRACTTEMP <DISTCUT> <FRACT_TEMPLATE>
: If fraction of templates with distance <= <DISTCUT> and fraction >= <FRACT_TEMPLATE> for 
: assignments/assignment possibilities, then true is returned; else false is returned
: Always returns true for ass/asspos with sequence separation <= 1

INSTRUCTUREFRACT <MIN_KEEP> <DISTCUT> <VIO_FRACT> <MAX_KEEP> <SEQSEP>
: First the ass/asspos are sorted in decreasing order by score
: Let VIOL6 = number of ass/asspos (with sequence separation >= <SEQSEP>) where min
: distance in templates is > <DISTCUT>
: Let COUNT6 = number of ass/asspos (with sequence separation >= <SEQSEP>)
: considered so far, starting from the ones with the highest score and includes
: the current ass/asspos
: If VIOL6/COUNT6 <= <VIO_FRACT> or sequence separation < <SEQSEP>, true is returned;
: else false is returned
: When false is returned for the first time, all subsequent times will always return
: false, even for ass/asspos with sequence separation < <SEQSEP>
: Always returns true for the first <MIN_KEEP> ass/asspos (with sequence
: separation >= <SEQSEP>) regardless of the value of VIOL6 unless <VIO_FRACT> is 0,
: in which as false is always returned no matter the value of <MIN_KEEP>
: If num ass/asspos with true return value is > <MAX_KEEP>, then false will always be
: returned.
: MAX_KEEP is the maximum number of ass/asspos to keep (those that return true).
: If <MAX_KEEP> is 0, then the maximum is set to 100 million, which represents no limit.
: If <MAX_KEEP> is -<float>, then <MAX_KEEP> is set to <float>*<num_assignable_contacts>,
: where the latter is the num contacts that potentially may have an NOE peak.
: e.g. INSTRUCTUREFRACT 10 6.0 0.05 -0.4 2
: Keep minimum 10 ass/asspos and maximum 0.4*<num_assignable_contacts> ass/asspos.
: If fraction of LR violations exceeds 0.05, return false. Only count violations
: for ass/asspos with seqsep >= <SEQSEP>. If seq sep < <SEQSEP>, then return true unless
: false was returned previously (in which case false will always be returned)

INSTRUCTUREWINDOW <WIN> <VIO_FRACT>
: First the assignments/asspos are sorted in decreasing order by score
: Consider the last <WIN> number of ass/asspos (with seq sep >= 2) considered so far
: Let VIOL2 = number of ass/asspos (with sequence separation >= 2) in <WIN> where min
: distance in templates is > distance <DISTCUTOFF> (by default set to 6.0A)
: If VIOL2/<WIN> <= <VIO_FRACT> or sequence separation < 2, true is returned; else false
: is returned
: When false is returned for the first time, all subsequent times will always return
: false, even for assignments/asspos with sequence separation < 2
: Always returns true for the first <WIN> assignments/asspos (with sequence sep >= 2)

INSTRUCTURECUT <FRACT_TEMPLATE> <VIO_FRACT>
: First the ass/asspos are sorted in decreasing order by score
: Let NUM = total number of all ass/asspos (auto-detects whether it is ass or asspos)
: Let VIOL = number of assignments/asspos tested so far (including current ass/asspos),
: where distance is > DISTCUTOFF (default is 6A) in fraction number of
: structures < <FRACT_TEMPLATE>
: If VIOL/NUM <= <VIO_FRACT> true is returned, else false is returned
: All sequence separations are considered

ASSIGNEDCOUNT <COUNT_CUT> <WIN>
: Let count=num ass/asspos within window size <WIN> (inclusive) around residues i & j.
: Returns true if count >= <COUNT_CUT>; e.g. WIN=0 only considers assignments between i
: and j (no window); count includes the current ass/asspos when deciding whether or not
: to keep the current ass/asspos so <CUT> should be set > 1
: If <CUT> is < 0 and is float value, then the assignments/asspos are sorted in
: decreasing order count. The top <CUT> % non-zero values return true and the rest
: return false. Only residues |i-j| >= 6 are used to set the cutoff for the top <CUT> %
: e.g. -0.1 returns top 10% |i-j| >= 6 non-zero assigned count values.
: e.g. ASSIGNEDCOUNT 3 1 returns true for at least 3 assignments;
: e.g. ASSIGNEDCOUNT -0.25 1 returns true for top 25% of all ass/asspos
: If the bound is negative, float, and |.| > 1.0, the top
: (count_cut-1)*<num_assignable_contacts> ass/asspos return true and rest return false.
: <num_assignable_contacts> is the number of contacts in the input 3D structures that
: are expected to have NOE peaks
: e.g. ASSIGNEDCOUNT -3.0 1 will return true for at most the top 2*<num_assignable_contacts>
: ass/asspos
: This should not be used to filter ambiguous assignments because it is not very useful. 
: The assignment possibility matrix is used for filtering ambiguous assignments.

ASSIGNEDCOUNT2 <COUNT_CUT_SR> <WIN_SR> <COUNT_CUT_LR> <WIN_LR>
: Similar to ASSIGNEDCOUNT except that residues i, j are divided into short-range and
: long-range assignments/asspos groups with different cutoffs.
: LR is defined by LRSEQSEP and it is by default equal to 6. For negative <CUT>,
: <num_assignable_contacts_SR> and <num_assignable_contacts_LR> are used.
: e.g. ASSIGNEDCOUNT2 -3.0 1 -2.0 1 returns true for the top 2*<num_assignable_contacts_SR>
: and the top 1*<num_assignable_contacts_LR> ass/asspos
: e.g. ASSIGNEDCOUNT2 -0.25 1 -0.15 1 returns true for the top 25% of all SR ass/asspos and
: the top 15% of all LR ass/asspos

SCORETERM1 <NUMTYPES> <TYPE1> <BOUND1> ... <TYPE_NUMTYPES> <BOUND_NUMTYPES>
: Returns true if all specified score terms >= bounds.
: If the bound is negative, float, and <= 1.0, the top bound scoring ass/asspos return true
: and the rest return false.
: e.g. -0.1 will return true for the top 10% scoring and false otherwise
: If the bound is negative, float, and |value| > 1.0, the top
: (bound-1)*<num_assignable_contacts> assignments/asspos return true and rest return false.
: <num_assignable_contacts> is the number of contacts in the input 3D structures that are
: expected to have NOE peaks
: e.g. -2.5 will return true for at most the top scoring 1.5*<num_assignable_contacts>
: assignments/asspos and -1.5 will return true for 0.5*<num_ass...>
: Score term types: 
: 0=cs
: 1=str
: 2=intensity
: 3=sym
: 4=interres
: 5=net
: 6=netStr
: 7=ambig
: 8=bias
: 9=total
: If TYPE is 10000+flag, then flag is used to turn on and off specific score terms
: 1=cs, 2=str, 4=int, 8=sym, 16=interes, 32=net, 64=netstr, 128=ambig, bias=256
: To turn on all except str, intensity use type equal to
: 10505 = 10000+505 (=1+8+16+32+64+128+256)
: To use all except str, int, netstr use 10000+441
: To use only interres, net, netstr use 16+32+64 = 10112
: For LR, can try sym+interres+bias = 10000+8+16+256=10280 (or 10346 for +str+netstr)
: For SR, can try sym+interres+net+bias = 10000+ 8+16+32+64+256 = 10312 (or 10378 for +str+netstr)
: Note: The use of flag usually only makes sense for NUMTYPES=1
: For filtering ambiguous assignments, the assignment possibilities are used instead of the assigned
: to determine the negative bound values

SCORETERM2 <NUMTYPES_SR> <TYPE1_SR> <BOUND1_SR> ... <TYPE_NUMTYPES_SR>
          <BOUND_NUMTYPES_SR> <NUMTYPES_LR> <TYPE1_LR> <BOUND1_LR> ...
          <TYPE_NUMTYPES_LR> <BOUND_NUMTYPES_LR>
: Similar to SCORETERM1 except we separate between short range and long range contacts
: (with seq sep >= 6).
: Use 0 for NUMTYPES if want to ignore SR or LR
: <num_assignable_contacts> is divided into <num_assignable_contacts_sr> and
: <num_assignable_contacts_lr>
: So <BOUNDi_LR> = -2.5 will return the top scoring 1.5*<num_assignable_contacts_lr>
: ass/asspos
: And <Boundi_LR> = -0.1 will return the top 10% scoring LR ass/asspos
: To skip SR or LR, use 0 for NUMTYPES.
: e.g. SCORETERM2 0 1 9 0.5 ignores SR ass and SCORETERM2 1 9 0.5 0 ignores LR.

SCORECOUNT1 <SUM> <NUMTYPES> <TYPE1> <BOUND1> ... <TYPE_NUMTYPES> <BOUND_NUMTYPES>
: Returns true if at least <SUM> (inclusive) specified score terms >= bounds.
: As in SCORETERM, each bound can be negative and float and flags are supported

SCORECOUNT2 <SUM_SR> <NUMTYPES_SR> <TYPE1_SR> <BOUND1_SR> ... <TYPE_NUMTYPES_SR>
            <BOUND_NUMTYPES_SR> <SUM_LR> <NUMTYPES_LR> <TYPE1_LR> <BOUND1_LR> ...
            <TYPE_NUMTYPES_LR> <BOUND_NUMTYPES_LR>
: Analogous to SCORETERM2. Use 0 0 to skip SR or LR
: e.g. SCORECOUNT2 0 0 1 1 9 0.5 skips SR assignments.

SCORENEIGHBORS <RELATIVE_SCORE_DIFF> <SCORETYPE>
: For filtering 1-to-1 only. If the relative difference between the scores between
: the current 1-1 assignment and the remaining assignments of the NOE peak
: (but skipping contacts already in 1-1 assigned)

SEQSEP <SEQ_CUTOFF>
: Returns true if < sequence cutoff. e.g. SEQSEP 2 means keep only intraresidue
: and sequential contacts i, i+1

TRUE
: Always returns true

FALSE
: Always returns false

NOISERATIO <ratio>
: Returns true if num NOE peaks/num residues <= ratio; else false is returned to
: indicate that the data is too noisy
: Recommended cutoff of 100


Combining Filters

Prefix notation is used for AND, OR, and NOT

AND <filter> <filter>
OR <filter> <filter>
NOT <filter>

IF <filter_test> <filter_body> ELSE <filter_else>
: body and else filter can be null (use the word "null" to specify this);
: ELSE keyword must be present (use ELSE true if have no else statement)


Other Filtering Parameters 

AMBIG_FILTER <val> <max>
: For filtering ambiguous assignments. Bigger float <val> means keep more.
: Set to a very big number to keep all ambig ass up to a max number set by
: integer <max>. Default <val> is 999999.0 and <max> 1000000.
: The assignment possibilities for each NOE peak are first sorted in
: decreasing order by score and filtered by FILTERAMBPOS. The top-scoring
: assignments for each NOE are retained until the sum of the scores
: is >= <val> or <max> assignments are reached. The remaining lower scoring
: assignment possibilities for the NOE are discarded. If score sum for an
: NOE peak is < AMBIG_FILTER, all of its assignment possibilities are retained.
: If <val> is negative, then assignments for each noe are filtered by the number
: of ass instead of by score
: e.g. AMBIG_FILTER -10 99 means to keep at most the top 10 scoring ass
: (<max> is ignored here)
: Larger values results in more ambiguity and less filtering

POOL_FILTER <float>
: Set to -1 to turn off (default). If > AMBIG_FILTER, then all amb assignment
: undergo constraint combination
: After filtering by (1) FILTERAMBPOS and (2) AMBIG_FILTER, the assignment
: possibilities with score sum < POOL_FILTER are placed in a pool for
: constraint combination instead of outputting directly.
: Restraints from two ambiguous assignment sets are merged until the score sum
: is >= POOL_FILTER. Wrap around of the pool (repeating assignments) is done to
: make sure score sum is >= POOL_FILTER.
: If the total sum of all pool assignments is < POOL_FILTER, then the pool
: step is skipped.
: If <float> is < -1, pool selection is based on the number of assignment
: possibilities instead of by score
: e.g. POOL_FILTER -7 means to do constraint combination if the number of
: possibilities is less than 7
: Larger values results in more ambiguity (larger ambiguous assignment groups)
:
: For NMR-I-TASSER output, residue-based contacts are outputted from the set of
: ambiguous proton-proton assignments.
: Supported residue-based contacts type={NN, CACA, CBCB, SCSC, NCA, CAN, NCB,
: CBN, NSC, SCN, CACB, CBCA, CASC, SCCA, CBSC, SCCB}
: Ambiguous assignments are grouped into residue-based contact groups
: (res1, res2, type). The scores are summed for each group. Groups are merged
: so there are at least 3 residue-based contacts and (if <float> > 0) at
: most 6 (inclusive).
: It can be less than 6 if the sum of scores is > <float> (if <float> > 0).

MAX_CALIB_ERROR <float>
: Used for distance calibration from intensities. Let R be the fraction
: of contacts used for calibration that violate the calibrated intensity-based
: distance upperbounds in the input structure. We want R < <float>. Therefore,
: the upperbound cutoffs will be increased if R > <float> (to decrease the
: violations). Default is 1.0, which means do not change the upperbounds

