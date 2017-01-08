#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "Alignment.h"

using namespace std;


AlignmentSetHead * GetAlignmentSet(char * alignmentSetFile){
    
    long int i = 0;
    
    FILE * fp;
    
    if((fp = fopen(alignmentSetFile, "r")) == NULL){
        printf("%s, does not exist!", alignmentSetFile);
        exit(0);
    }
    
    long int alignmentCount = 0;
    char * line = (char *)malloc(sizeof(char)*10000);
    
    while(!feof(fp)){
        fgets(line, 10000, fp);
        if(feof(fp))break;
        if(i<4){
            i++;
            continue;
        }
        alignmentCount++;
    }
    
    AlignmentSetHead * alignmentSetHead = (AlignmentSetHead *)malloc(sizeof(AlignmentSetHead));
    alignmentSetHead->alignmentCount = alignmentCount;
    alignmentSetHead->alignmentSet = (Alignment*)malloc(sizeof(Alignment)*alignmentCount);
    
    for(i=0;i<alignmentCount;i++){
        alignmentSetHead->alignmentSet[i].referenceName = NULL;
        alignmentSetHead->alignmentSet[i].contigName = NULL;
        alignmentSetHead->alignmentSet[i].referenceStart = 0;
        alignmentSetHead->alignmentSet[i].referenceEnd = 0;
        alignmentSetHead->alignmentSet[i].contigStart = 0;
        alignmentSetHead->alignmentSet[i].contigEnd = 0;
        alignmentSetHead->alignmentSet[i].referenceOverlapLength = 0;
        alignmentSetHead->alignmentSet[i].contigOverlapLength = 0;
        alignmentSetHead->alignmentSet[i].identityPercent = 0;
        alignmentSetHead->alignmentSet[i].isReverse = false;
        alignmentSetHead->alignmentSet[i].contigLength = 0;
        alignmentSetHead->alignmentSet[i].scaffoldIndex = -1;
        alignmentSetHead->alignmentSet[i].contigIndex = -1;
    }
    
    fclose(fp);
    
    if((fp = fopen(alignmentSetFile, "r")) == NULL){
        printf("%s, does not exist!", alignmentSetFile);
        exit(0);
    }
    
    char * p = NULL;
    char * referenceName = NULL;
    
    i = 0;
    long int index = 0;
    
    while(!feof(fp)){
        fgets(line, 10000, fp);
        if(feof(fp))break;
        if(i<4){
            i++;
            continue;
        }
        
        for(long int j = 0; j<13; j++){
            if(j==0){
                p = strtok(line, "\t");
            }else{
                p = strtok(NULL, "\t");
            }
            if(j==0){
                sscanf(p, "%ld", &alignmentSetHead->alignmentSet[index].referenceStart);
                alignmentSetHead->alignmentSet[index].referenceStart--;
            }
            if(j==1){
                sscanf(p, "%ld", &alignmentSetHead->alignmentSet[index].referenceEnd);
                alignmentSetHead->alignmentSet[index].referenceEnd--;
            }
            if(j==2){
                sscanf(p, "%ld", &alignmentSetHead->alignmentSet[index].contigStart);
                alignmentSetHead->alignmentSet[index].contigStart--;
            }
            if(j==3){
                sscanf(p, "%ld", &alignmentSetHead->alignmentSet[index].contigEnd);
                alignmentSetHead->alignmentSet[index].contigEnd--;
            }
            if(j==4){
                sscanf(p, "%ld", &alignmentSetHead->alignmentSet[index].referenceOverlapLength);
            }
            if(j==5){
                sscanf(p, "%ld", &alignmentSetHead->alignmentSet[index].contigOverlapLength);
            }
            if(j==6){
                sscanf(p, "%lf", &alignmentSetHead->alignmentSet[index].identityPercent);
            }
            if(j==8){
                sscanf(p, "%ld", &alignmentSetHead->alignmentSet[index].contigLength);
            }
            if(j==10){
                int t = 0;
                sscanf(p, "%d", &t);
                if(t == -1){
                    alignmentSetHead->alignmentSet[index].isReverse = true;
                }else{
                    alignmentSetHead->alignmentSet[index].isReverse = false;
                }
            }
            if(j==11){
                alignmentSetHead->alignmentSet[index].referenceName = (char *)malloc(sizeof(char)*(strlen(p)+1));
                strcpy(alignmentSetHead->alignmentSet[index].referenceName, p);
                alignmentSetHead->alignmentSet[index].referenceName[strlen(p)] = '\0';
            }
            if(j==12){
                alignmentSetHead->alignmentSet[index].contigName = (char *)malloc(sizeof(char)*(strlen(p)+1));
                strcpy(alignmentSetHead->alignmentSet[index].contigName, p);
                alignmentSetHead->alignmentSet[index].contigName[strlen(p)] = '\0';
                char * tempName = (char *)malloc(sizeof(char)*(strlen(p)+1));
                strcpy(tempName, p);
                tempName[strlen(p)] = '\0';
                char * p1 = NULL;
                p1 = strtok(tempName, "_");
                p1 = strtok(NULL, "_");
                sscanf(p1, "%ld", &alignmentSetHead->alignmentSet[index].scaffoldIndex);
                p1 = strtok(NULL, "_");
                sscanf(p1, "%ld", &alignmentSetHead->alignmentSet[index].contigIndex);
                free(tempName);
            }
            
        }
        index++;
    }
    
    return alignmentSetHead;
    
}

void OutputAlignmentSet(AlignmentSetHead * alignmentSetHead){
    for(long int i = 0; i < alignmentSetHead->alignmentCount; i++){
        cout<<alignmentSetHead->alignmentSet[i].referenceStart<<"\t";
        cout<<alignmentSetHead->alignmentSet[i].referenceEnd<<"\t";
        cout<<alignmentSetHead->alignmentSet[i].contigStart<<"\t";
        cout<<alignmentSetHead->alignmentSet[i].contigEnd<<"\t";
        cout<<alignmentSetHead->alignmentSet[i].referenceOverlapLength<<"\t";
        cout<<alignmentSetHead->alignmentSet[i].contigOverlapLength<<"\t";
        cout<<alignmentSetHead->alignmentSet[i].identityPercent<<"\t";
        cout<<alignmentSetHead->alignmentSet[i].isReverse<<"\t";
        cout<<alignmentSetHead->alignmentSet[i].contigLength<<"\t";
        cout<<alignmentSetHead->alignmentSet[i].referenceName<<"\t";
        cout<<alignmentSetHead->alignmentSet[i].contigName<<"\t";
        cout<<alignmentSetHead->alignmentSet[i].scaffoldIndex<<"\t";
        cout<<alignmentSetHead->alignmentSet[i].contigIndex<<endl;
    }
}
