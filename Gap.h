#ifndef GAP_H_INCLUDED 
#define GAP_H_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "ContigSet.h"
#include "ScaffoldSet.h"
#include "Alignment.h"

using namespace std;

typedef struct GapRegion{
    int gapType;
    long int leftStartPosition;
    long int leftContigEndPosition;
    long int rightStartPosition;
    long int rightContigStartPosition;
    bool leftReverse;
    bool rightReverse;
    long int leftReferenceIndex;
    long int rightReferenceIndex;
    char * leftReferenceName;
    char * rightReferenceName;
}GapRegion;

typedef struct ScaffoldGapRegion{
    GapRegion * gapRegionSet;
}ScaffoldGapRegion;


typedef struct ContigAlignmentIndex{
    long int * alignmentIndex;
    long int alignmentIndexCount;
}ContigAlignmentIndex;

typedef struct ScaffoldAlignmentIndex{
    ContigAlignmentIndex * contigAlignmentIndex;
}ScaffoldAlignmentIndex;



ScaffoldAlignmentIndex * GetScaffoldAlignmentIndex(ScaffoldSetHead * scaffoldSetHead, AlignmentSetHead * alignmentSetHead);
int GetGapRegion(ContigSetHead * referenceContigSetHead, ScaffoldSetHead * scaffoldSetHead, ScaffoldAlignmentIndex * scaffoldAlignmentIndex, AlignmentSetHead * alignmentSetHead, ScaffoldGapRegion * scaffoldGapRegion, long int scaffoldIndex, long int leftContigIndex, long int rightContigIndex, long int gapIndex, bool isReference);
ScaffoldGapRegion * GetScaffoldGapRegionInReference(ScaffoldSetHead * scaffoldSetHead, AlignmentSetHead * alignmentSetHead, ContigSetHead * referenceContigSetHead, bool isReference);
long int SearchReferenceIndexFromName(Contig * reference, long int referenceNumber, char * referenceName);
void OutputScaffoldGapRegion(ScaffoldSetHead * scaffoldSetHead, ScaffoldGapRegion * scaffoldGapRegion, ContigSetHead * referenceContigSetHead, char * gapFileName);
void OutputReferenceGapRegion(ScaffoldSetHead * scaffoldSetHead, ScaffoldGapRegion * scaffoldGapRegion, ContigSetHead * referenceContigSetHead, char * gapFileName);
void OutputReferenceGapRegion1(ScaffoldSetHead * scaffoldSetHead, ScaffoldGapRegion * referenceGapRegion, ScaffoldGapRegion * scaffoldGapRegion, ContigSetHead * referenceContigSetHead, char * gapFileName); 
void CheckFillingGapRegion(ScaffoldSetHead * scaffoldSetHead, ScaffoldGapRegion * fillingGapRegion);
bool ReverseComplement(char * temp);

int OptimizeScaffoldAlignmentIndex(ContigSetHead * referenceContigSetHead, ScaffoldSetHead * scaffoldSetHead, ScaffoldAlignmentIndex * scaffoldAlignmentIndex, AlignmentSetHead * alignmentSetHead);

void OutputGapInScaffoldInformation(ScaffoldSetHead * scaffoldSetHead, ScaffoldGapRegion * scaffoldGapRegion, ScaffoldGapRegion * fillingGapRegion, char * fileName);
long int OutputMinContigAsGap(ScaffoldSetHead * scaffoldSetHead, long int i, long int j);

#endif 
