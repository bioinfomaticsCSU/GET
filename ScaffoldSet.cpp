#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "ScaffoldSet.h"

using namespace std;


ScaffoldSetHead * GetScaffoldSetFromScaffoldFile(char * scaffoldFileName, long int minContigLength){
    
    long int i = 0;
    long int j = 0;
    
    long int maxSize = 10000;
    char * scaffold = NULL;
    if(NULL == (scaffold = (char*)malloc(sizeof(char)*maxSize))){
        perror("malloc error!");
        exit(1);
    }
    
    long int scaffoldCount = 0;
    
    FILE * fp;
    
    if((fp = fopen(scaffoldFileName, "r")) == NULL){
        printf("%s, does not exist!", scaffoldFileName);
        exit(0);
    }
    
    while((fgets(scaffold, maxSize, fp)) != NULL){ 
       if(scaffold[0] == '>'){  
           scaffoldCount++; 
       }  
    }  
    
    fclose(fp);
    
    ScaffoldSetHead * scaffoldSetHead = NULL;
    if(NULL == (scaffoldSetHead = (ScaffoldSetHead*)malloc(sizeof(ScaffoldSetHead)))){
        perror("ScaffoldSetHead malloc error!");
        exit(1);
    }
    scaffoldSetHead->scaffoldCount = scaffoldCount;
    if(NULL == (scaffoldSetHead->scaffoldSet = (Scaffold*)malloc(sizeof(Scaffold)*scaffoldCount))){
        perror("ScaffoldSet malloc error!");
        exit(1);
    }
    
    scaffoldSetHead->contigCount = 0;
    scaffoldSetHead->gapCount = 0;
    scaffoldSetHead->minGapDistance = 1;
    scaffoldSetHead->minContigLength = minContigLength;
    
    
    for(i = 0; i < scaffoldCount; i++){ 
        scaffoldSetHead->scaffoldSet[i].scaffold = NULL;
        scaffoldSetHead->scaffoldSet[i].contigCount = 0;
        scaffoldSetHead->scaffoldSet[i].scaffoldName = NULL;
        scaffoldSetHead->scaffoldSet[i].scaffoldLength = 0;
        scaffoldSetHead->scaffoldSet[i].gapStartPosition = NULL;
        scaffoldSetHead->scaffoldSet[i].gapDistance = NULL;
        scaffoldSetHead->scaffoldSet[i].gapCount = 0;
        scaffoldSetHead->scaffoldSet[i].contigStartPosition = NULL;
        scaffoldSetHead->scaffoldSet[i].contigLength = NULL;
    }
    
    if((fp = fopen(scaffoldFileName, "r")) == NULL){
        printf("%s, does not exist!", scaffoldFileName);
        exit(0);
    }
    
    long int scaffoldIndex = -1;
    long int allocateLength = 0;
    
    while((fgets(scaffold, maxSize, fp)) != NULL){ 
       
       if(scaffold[0] == '>'){  
           
           if(strlen(scaffold) == maxSize-1){              
               while((fgets(scaffold, maxSize, fp)) != NULL){
                   if(strlen(scaffold) != maxSize-1){
                       break;
                   }
               }        
           }
           scaffoldIndex++;
           allocateLength = 0;
           char * p = strtok(scaffold, " ");
           scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldName = (char *)malloc(sizeof(char)*(strlen(p)-1));
           strncpy(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldName, p+1, strlen(p)-2);
           scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldName[strlen(p)-2] = '\0';
           //sprintf(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldName + strlen(p), "_%ld", scaffoldIndex);
           continue;
           
       }
       
       long int extendLength = strlen(scaffold);
       if(scaffold[extendLength-1] == '\n'){
           extendLength--;
       }
       long int scaffoldLength = 0;
       char * tempScaffold = NULL;
       if(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold != NULL){
           if(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength + extendLength >= allocateLength){
               scaffoldLength = scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength;
               scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold = (char *)realloc(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold, allocateLength + maxSize + 1);
               allocateLength = allocateLength + maxSize + 1;
               
               strncpy(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold + scaffoldLength, scaffold, extendLength);
               scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold[scaffoldLength + extendLength] = '\0';    
               scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength = scaffoldLength + extendLength;
                     
           }else{
               strncpy(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold + scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength, scaffold, extendLength);
               scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold[scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength + extendLength] = '\0';
               scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength = scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength + extendLength;
           }   
           
       }else{
           scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold = (char *)malloc(sizeof(char)*(maxSize+1));
           strncpy(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold, scaffold, extendLength);
           scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold[extendLength] = '\0';
           scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength = extendLength;
           allocateLength = maxSize + 1;
       }      
    }  
    fflush(fp);
    fclose(fp);
    GetGapInScaffoldSet(scaffoldSetHead);
    OptimizeShortContigToGap(scaffoldSetHead);
    return scaffoldSetHead;
    
}

void GetGapInScaffoldSet(ScaffoldSetHead * scaffoldSetHead){
    
    for(long int i = 0; i < scaffoldSetHead->scaffoldCount; i++){
        long int gapDistance = 0;
        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].scaffoldLength; j++){
            
            if(scaffoldSetHead->scaffoldSet[i].scaffold[j] == 'N' || scaffoldSetHead->scaffoldSet[i].scaffold[j] == 'n'){
                gapDistance++;
                continue;
            }else{
                if(gapDistance > 0){
                    if(gapDistance >= scaffoldSetHead->minGapDistance){
                        scaffoldSetHead->scaffoldSet[i].gapCount++;
                    }
                    gapDistance = 0;
                }
            }
        }
        if(gapDistance >= scaffoldSetHead->minGapDistance){
            scaffoldSetHead->scaffoldSet[i].gapCount++;
        }
        if(scaffoldSetHead->scaffoldSet[i].gapCount <=0){
            scaffoldSetHead->scaffoldSet[i].contigCount = 1;
            scaffoldSetHead->contigCount++;
            scaffoldSetHead->scaffoldSet[i].contigStartPosition = (long int*)malloc(sizeof(long int)*scaffoldSetHead->scaffoldSet[i].contigCount);
            scaffoldSetHead->scaffoldSet[i].contigLength = (long int*)malloc(sizeof(long int)*scaffoldSetHead->scaffoldSet[i].contigCount);
            scaffoldSetHead->scaffoldSet[i].contigStartPosition[0] = 0;
            scaffoldSetHead->scaffoldSet[i].contigLength[0] = scaffoldSetHead->scaffoldSet[i].scaffoldLength;
            continue;
        }
        gapDistance = 0;
        scaffoldSetHead->scaffoldSet[i].gapDistance = (long int*)malloc(sizeof(long int)*scaffoldSetHead->scaffoldSet[i].gapCount);
        scaffoldSetHead->scaffoldSet[i].gapStartPosition = (long int*)malloc(sizeof(long int)*scaffoldSetHead->scaffoldSet[i].gapCount);
        long int gapIndex = 0;
        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].scaffoldLength; j++){
            
            if(scaffoldSetHead->scaffoldSet[i].scaffold[j] == 'N' || scaffoldSetHead->scaffoldSet[i].scaffold[j] == 'n'){
                gapDistance++;
                continue;
            }else{
                if(gapDistance > 0){
                    if(gapDistance >= scaffoldSetHead->minGapDistance){
                        scaffoldSetHead->scaffoldSet[i].gapDistance[gapIndex] = gapDistance;
                        scaffoldSetHead->scaffoldSet[i].gapStartPosition[gapIndex] = j - gapDistance;
                        gapIndex++;
                    }
                    gapDistance = 0;
                }
            }
        }
        if(gapDistance >= scaffoldSetHead->minGapDistance){
            scaffoldSetHead->scaffoldSet[i].gapDistance[gapIndex] = gapDistance;
            scaffoldSetHead->scaffoldSet[i].gapStartPosition[gapIndex] = scaffoldSetHead->scaffoldSet[i].scaffoldLength - gapDistance;
        }
        
        scaffoldSetHead->scaffoldSet[i].contigCount = scaffoldSetHead->scaffoldSet[i].gapCount + 1;
        int startGap = 0;
        if(scaffoldSetHead->scaffoldSet[i].gapStartPosition[0] ==0){
            scaffoldSetHead->scaffoldSet[i].contigCount--;
            startGap = 1;
        }
        if(scaffoldSetHead->scaffoldSet[i].gapStartPosition[scaffoldSetHead->scaffoldSet[i].gapCount - 1] + 
            scaffoldSetHead->scaffoldSet[i].gapDistance[scaffoldSetHead->scaffoldSet[i].gapCount - 1] == scaffoldSetHead->scaffoldSet[i].scaffoldLength){
            scaffoldSetHead->scaffoldSet[i].contigCount--;
        }
        scaffoldSetHead->scaffoldSet[i].contigStartPosition = (long int*)malloc(sizeof(long int)*scaffoldSetHead->scaffoldSet[i].contigCount);
        scaffoldSetHead->scaffoldSet[i].contigLength = (long int*)malloc(sizeof(long int)*scaffoldSetHead->scaffoldSet[i].contigCount);
        for(long int p = 0; p < scaffoldSetHead->scaffoldSet[i].contigCount; p++){
            
            if(p + startGap - 1 < 0){
                scaffoldSetHead->scaffoldSet[i].contigStartPosition[p] = 0;
            }else{
                scaffoldSetHead->scaffoldSet[i].contigStartPosition[p] = scaffoldSetHead->scaffoldSet[i].gapStartPosition[p + startGap - 1] +  scaffoldSetHead->scaffoldSet[i].gapDistance[p + startGap - 1];
            }
            if(p + startGap >= scaffoldSetHead->scaffoldSet[i].gapCount){
                scaffoldSetHead->scaffoldSet[i].contigLength[p] = scaffoldSetHead->scaffoldSet[i].scaffoldLength - scaffoldSetHead->scaffoldSet[i].contigStartPosition[p];
            }else{
                scaffoldSetHead->scaffoldSet[i].contigLength[p] = scaffoldSetHead->scaffoldSet[i].gapStartPosition[p + startGap] - scaffoldSetHead->scaffoldSet[i].contigStartPosition[p];
            }
            
        }
        scaffoldSetHead->gapCount = scaffoldSetHead->gapCount + scaffoldSetHead->scaffoldSet[i].gapCount;
        scaffoldSetHead->contigCount = scaffoldSetHead->contigCount + scaffoldSetHead->scaffoldSet[i].contigCount;
        
    }   
}

void OptimizeShortContigToGap(ScaffoldSetHead * scaffoldSetHead){
    for(long int i = 0; i < scaffoldSetHead->scaffoldCount; i++){
        long int shortContigCount = 0;
        if(scaffoldSetHead->scaffoldSet[i].gapCount <= 0){
            continue;
        }
        int startGap = 0;
        if(scaffoldSetHead->scaffoldSet[i].gapStartPosition[0] ==0){
            startGap = 1;
        }
        int endGap = 0;
        if(scaffoldSetHead->scaffoldSet[i].gapStartPosition[scaffoldSetHead->scaffoldSet[i].gapCount - 1] + 
            scaffoldSetHead->scaffoldSet[i].gapDistance[scaffoldSetHead->scaffoldSet[i].gapCount - 1] == scaffoldSetHead->scaffoldSet[i].scaffoldLength){
            endGap = 1;
        }
        
        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].contigCount; j++){
            
            if(scaffoldSetHead->scaffoldSet[i].contigLength[j] < scaffoldSetHead->minContigLength){
                long int previousGapIndex = j + startGap - 1;
                long int nextGapIndex = j + startGap;
                if(previousGapIndex >=0){
                    scaffoldSetHead->scaffoldSet[i].gapDistance[previousGapIndex] = scaffoldSetHead->scaffoldSet[i].gapDistance[previousGapIndex]
                       + scaffoldSetHead->scaffoldSet[i].contigLength[j];
                }else{
                    scaffoldSetHead->scaffoldSet[i].gapDistance[nextGapIndex] = scaffoldSetHead->scaffoldSet[i].gapDistance[nextGapIndex]
                       + scaffoldSetHead->scaffoldSet[i].contigLength[j];
                    scaffoldSetHead->scaffoldSet[i].gapStartPosition[nextGapIndex] = scaffoldSetHead->scaffoldSet[i].gapStartPosition[nextGapIndex]
                       - scaffoldSetHead->scaffoldSet[i].contigLength[j];
                }
                scaffoldSetHead->scaffoldSet[i].contigLength[j] = -1;
                scaffoldSetHead->scaffoldSet[i].contigStartPosition[j] = -1;
                shortContigCount++;
            }
        }
        
        if(shortContigCount == 0){
            continue;
        }
        if(scaffoldSetHead->scaffoldSet[i].contigCount - shortContigCount <= 0){
            scaffoldSetHead->scaffoldSet[i].contigCount = 1;
            scaffoldSetHead->scaffoldSet[i].gapCount = 0;
            continue;
        }
        
        long int index = 0;
        long int * tempContigLength = (long int *)malloc(sizeof(long int)*(scaffoldSetHead->scaffoldSet[i].contigCount - shortContigCount));
        long int * tempContigStartPosition = (long int *)malloc(sizeof(long int)*(scaffoldSetHead->scaffoldSet[i].contigCount - shortContigCount));
        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].contigCount; j++){
            if(scaffoldSetHead->scaffoldSet[i].contigLength[j] != -1){
                tempContigLength[index] = scaffoldSetHead->scaffoldSet[i].contigLength[j];
                tempContigStartPosition[index] = scaffoldSetHead->scaffoldSet[i].contigStartPosition[j];
                index++;
            }
        }
        free(scaffoldSetHead->scaffoldSet[i].contigLength);
        free(scaffoldSetHead->scaffoldSet[i].contigStartPosition);
        scaffoldSetHead->scaffoldSet[i].contigLength = tempContigLength;
        scaffoldSetHead->scaffoldSet[i].contigStartPosition = tempContigStartPosition;
        scaffoldSetHead->scaffoldSet[i].contigCount = scaffoldSetHead->scaffoldSet[i].contigCount - shortContigCount;
        
        long int gapRemoveCount = 0;
        long int nextIndex = 0;
        long int j = 0; 
        while(j < scaffoldSetHead->scaffoldSet[i].gapCount - 1){
            nextIndex = j + 1;
            while(scaffoldSetHead->scaffoldSet[i].gapStartPosition[j] + scaffoldSetHead->scaffoldSet[i].gapDistance[j] 
                == scaffoldSetHead->scaffoldSet[i].gapStartPosition[nextIndex]){
                scaffoldSetHead->scaffoldSet[i].gapDistance[j] = scaffoldSetHead->scaffoldSet[i].gapDistance[j] 
                    + scaffoldSetHead->scaffoldSet[i].gapDistance[nextIndex];
                scaffoldSetHead->scaffoldSet[i].gapStartPosition[nextIndex] = -1;
                nextIndex++;
                gapRemoveCount++;
            }
            j = nextIndex;
        }
        
        index = 0;
        long int * tempGapLength = (long int *)malloc(sizeof(long int)*(scaffoldSetHead->scaffoldSet[i].gapCount - gapRemoveCount));
        long int * tempGapStartPosition = (long int *)malloc(sizeof(long int)*(scaffoldSetHead->scaffoldSet[i].gapCount - gapRemoveCount));
        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++){
            if(scaffoldSetHead->scaffoldSet[i].gapStartPosition[j] != -1){
                tempGapLength[index] = scaffoldSetHead->scaffoldSet[i].gapDistance[j];
                tempGapStartPosition[index] = scaffoldSetHead->scaffoldSet[i].gapStartPosition[j];
                index++;
            }
        }
        free(scaffoldSetHead->scaffoldSet[i].gapDistance);
        free(scaffoldSetHead->scaffoldSet[i].gapStartPosition);
        scaffoldSetHead->scaffoldSet[i].gapDistance = tempGapLength;
        scaffoldSetHead->scaffoldSet[i].gapStartPosition = tempGapStartPosition;
        scaffoldSetHead->scaffoldSet[i].gapCount = scaffoldSetHead->scaffoldSet[i].gapCount - gapRemoveCount;
        
        
        
        
        
    }

    
}

void OutputContigSetOfScaffoldSet(ScaffoldSetHead * scaffoldSetHead, char * contigSetFile){
    
    FILE * fp;
    
    if((fp = fopen(contigSetFile, "w")) == NULL){
        printf("%s, can not be build!", contigSetFile);
        exit(0);
    }
    
    for(long int i = 0; i < scaffoldSetHead->scaffoldCount; i++){
        if(scaffoldSetHead->scaffoldSet[i].gapCount <=0){
            continue;
        }
        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].contigCount; j++){
            fprintf(fp, ">scaffold_%ld_%ld_\%ld_%ld\n", i, j, scaffoldSetHead->scaffoldSet[i].contigStartPosition[j], scaffoldSetHead->scaffoldSet[i].contigLength[j]);
            char * tempContig = (char *)malloc(sizeof(char)*(scaffoldSetHead->scaffoldSet[i].contigLength[j]+1));
            strncpy(tempContig, scaffoldSetHead->scaffoldSet[i].scaffold + scaffoldSetHead->scaffoldSet[i].contigStartPosition[j], scaffoldSetHead->scaffoldSet[i].contigLength[j]);
            tempContig[scaffoldSetHead->scaffoldSet[i].contigLength[j]] = '\0';
            fprintf(fp, "%s\n", tempContig);
        }
    }   
    fflush(fp);
    fclose(fp);
}

void WriteScaffoldSet(ScaffoldSetHead * scaffoldSetHead, char * scaffoldSetInforFile){
    
    FILE *fp;
    
    if((fp = fopen(scaffoldSetInforFile, "w")) == NULL){
        printf("%s, can not be build!", scaffoldSetInforFile);
        exit(0);
    }
    
    for(long int i = 0; i < scaffoldSetHead->scaffoldCount; i++){
        
        fprintf(fp, "%s\t%ld\t%ld\t%ld\n", scaffoldSetHead->scaffoldSet[i].scaffoldName, scaffoldSetHead->scaffoldSet[i].scaffoldLength, scaffoldSetHead->scaffoldSet[i].gapCount, scaffoldSetHead->scaffoldSet[i].contigCount);
        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++){
            fprintf(fp, "%ld\t%ld\n", scaffoldSetHead->scaffoldSet[i].gapStartPosition[j], scaffoldSetHead->scaffoldSet[i].gapDistance[j]);
        }
        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].contigCount; j++){
            fprintf(fp, "%ld\t%ld\n", scaffoldSetHead->scaffoldSet[i].contigStartPosition[j], scaffoldSetHead->scaffoldSet[i].contigLength[j]);
        }
        
    }
    fflush(fp);
    fclose(fp);
    
}



