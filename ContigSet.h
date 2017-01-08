#ifndef CONTIGSET_H_INCLUDED 
#define CONTIGSET_H_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

using namespace std;

typedef struct Contig{
   char * contigName;
   char * contig;
   long int contigLength;
}Contig;

typedef struct ContigSetHead{
   Contig * contigSet;
   long int contigCount;
}ContigSetHead;




ContigSetHead * GetContigSet(char * contigSetFile);

#endif 
