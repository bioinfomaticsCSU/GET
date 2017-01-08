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
    
    ScaffoldSetHead * scaffoldSetHead = GetScaffoldSetFromScaffoldFile(argv[1]);
    cout<<"aa"<<endl;
    AlignmentSetHead * referenceAlignmentSetHead = GetAlignmentSet(argv[2]);
    ContigSetHead * referenceContigSetHead = GetContigSet(argv[3]);
    cout<<"bb00"<<endl;
    ScaffoldGapRegion * scaffoldGapRegion = GetScaffoldGapRegionInReference(scaffoldSetHead, referenceAlignmentSetHead, referenceContigSetHead, true);
    cout<<"bb"<<endl;

    AlignmentSetHead * fillingAlignmentSetHead = GetAlignmentSet(argv[5]);
    ContigSetHead * fillingContigSetHead = GetContigSet(argv[6]);
    ScaffoldGapRegion * scaffoldFillingGapRegion = GetScaffoldGapRegionInReference(scaffoldSetHead, fillingAlignmentSetHead, fillingContigSetHead, false);
    CheckFillingGapRegion(scaffoldSetHead, scaffoldFillingGapRegion);
    cout<<"cc"<<endl;
    OutputReferenceGapRegion1(scaffoldSetHead, scaffoldGapRegion, scaffoldFillingGapRegion, referenceContigSetHead, argv[4]);
    //OutputReferenceGapRegion(scaffoldSetHead, scaffoldGapRegion, referenceContigSetHead, argv[4]);
    cout<<"dd"<<endl;
    OutputScaffoldGapRegion(scaffoldSetHead, scaffoldFillingGapRegion, fillingContigSetHead, argv[7]);
    cout<<"ee"<<endl;

    OutputGapInScaffoldInformation(scaffoldSetHead, scaffoldGapRegion, scaffoldFillingGapRegion, argv[8]);
    
}
