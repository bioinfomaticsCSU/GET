#ifndef SCAFFOLDSET_H_INCLUDED 
#define SCAFFOLDSET_H_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

using namespace std;

typedef struct Scaffold{
   char * scaffold;
   char * scaffoldName;
   long int contigCount;
   long int gapCount;
   long int * gapDistance;
   long int * gapStartPosition;
   long int * contigLength;
   long int * contigStartPosition;
   long int scaffoldLength;
}Scaffold;

typedef struct ScaffoldSetHead{
   Scaffold * scaffoldSet;
   long int scaffoldCount;
   long int contigCount;
   long int gapCount;
   long int minGapDistance;
   long int minContigLength;
}ScaffoldSetHead;


ScaffoldSetHead * GetScaffoldSetFromScaffoldFile(char * scaffoldFileName);
void GetGapInScaffoldSet(ScaffoldSetHead * scaffoldSetHead);
void OutputContigSetOfScaffoldSet(ScaffoldSetHead * scaffoldSetHead, char * contigSetFile);
void WriteScaffoldSet(ScaffoldSetHead * scaffoldSetHead, char * scaffoldSetInforFile);

void OptimizeShortContigToGap(ScaffoldSetHead * scaffoldSetHead);


#endif 
