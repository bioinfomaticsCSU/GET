#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ContigSet.h"
#include "ScaffoldSet.h"
#include "Alignment.h"
#include "Gap.h"

using namespace std;
 
int main(int argc, char *argv[])
{  
    
    ScaffoldSetHead * scaffoldSetHead = GetScaffoldSetFromScaffoldFile(argv[1], atoi(argv[9]));
    
    scaffoldSetHead->minSegmentDistanceEndLength = atoi(argv[10]);
    scaffoldSetHead->minDistanceRelocationLength = atoi(argv[11]);
    scaffoldSetHead->minTimesRelocation = atoi(argv[12]);
    

    AlignmentSetHead * referenceAlignmentSetHead = GetAlignmentSet(argv[2]);
    ContigSetHead * referenceContigSetHead = GetContigSet(argv[3]);
    ScaffoldGapRegion * scaffoldGapRegion = GetScaffoldGapRegionInReference(scaffoldSetHead, referenceAlignmentSetHead, referenceContigSetHead, true);

    AlignmentSetHead * fillingAlignmentSetHead = GetAlignmentSet(argv[5]);
    ContigSetHead * fillingContigSetHead = GetContigSet(argv[6]);
    ScaffoldGapRegion * scaffoldFillingGapRegion = GetScaffoldGapRegionInReference(scaffoldSetHead, fillingAlignmentSetHead, fillingContigSetHead, false);

    OutputReferenceGapRegion(scaffoldSetHead, scaffoldGapRegion, scaffoldFillingGapRegion, referenceContigSetHead, argv[4]);

    OutputScaffoldGapRegion(scaffoldSetHead, scaffoldFillingGapRegion, fillingContigSetHead, argv[7]);

    OutputGapInScaffoldInformation(scaffoldSetHead, scaffoldGapRegion, scaffoldFillingGapRegion, argv[8]);
    
}
