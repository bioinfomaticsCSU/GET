#ifndef ALIGNMENT_H_INCLUDED 
#define ALIGNMENT_H_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

using namespace std;

typedef struct Alignment{
    char * contigName;
    char * referenceName;
    long int referenceStart;
    long int referenceEnd;
    long int contigStart;
    long int contigEnd;
    long int referenceOverlapLength;
    long int contigOverlapLength;
    double identityPercent;
    bool isReverse;
    long int contigLength;
    long int scaffoldIndex;
    long int contigIndex;
}Alignment;

typedef struct AlignmentSetHead{
   Alignment * alignmentSet;
   long int alignmentCount;
}AlignmentSetHead;




AlignmentSetHead * GetAlignmentSet(char * alignmentSetFile);
void OutputAlignmentSet(AlignmentSetHead * alignmentSetHead);
#endif 
