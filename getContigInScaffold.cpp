#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ScaffoldSet.h"

using namespace std;
 
int main(int argc, char *argv[])
{
    ScaffoldSetHead * scaffoldSetHead = GetScaffoldSetFromScaffoldFile(argv[1], atoi(argv[3]));

    OutputContigSetOfScaffoldSet(scaffoldSetHead, argv[2]);
}
