#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <iomanip>

#include "Gap.h"

using namespace std;

ScaffoldAlignmentIndex * GetScaffoldAlignmentIndex(ScaffoldSetHead * scaffoldSetHead, AlignmentSetHead * alignmentSetHead){
    ScaffoldAlignmentIndex * scaffoldAlignmentIndex = (ScaffoldAlignmentIndex *)malloc(sizeof(ScaffoldAlignmentIndex)*scaffoldSetHead->scaffoldCount);
    for(long int i = 0; i < scaffoldSetHead->scaffoldCount; i++){
        scaffoldAlignmentIndex[i].contigAlignmentIndex = (ContigAlignmentIndex *)malloc(sizeof(ContigAlignmentIndex)*scaffoldSetHead->scaffoldSet[i].contigCount);
        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].contigCount; j++){
            scaffoldAlignmentIndex[i].contigAlignmentIndex[j].alignmentIndex = NULL;
            scaffoldAlignmentIndex[i].contigAlignmentIndex[j].alignmentIndexCount = 0;
        }
    }
    for(long int i = 0; i < alignmentSetHead->alignmentCount; i++){
        scaffoldAlignmentIndex[alignmentSetHead->alignmentSet[i].scaffoldIndex].contigAlignmentIndex[alignmentSetHead->alignmentSet[i].contigIndex].alignmentIndexCount++;
    }
    for(long int i = 0; i < scaffoldSetHead->scaffoldCount; i++){
        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].contigCount; j++){
            if(scaffoldAlignmentIndex[i].contigAlignmentIndex[j].alignmentIndexCount > 0){
                scaffoldAlignmentIndex[i].contigAlignmentIndex[j].alignmentIndex = (long int *)malloc(sizeof(long int)*scaffoldAlignmentIndex[i].contigAlignmentIndex[j].alignmentIndexCount);
                scaffoldAlignmentIndex[i].contigAlignmentIndex[j].alignmentIndexCount = 0;
            }
        }
    }
    for(long int i = 0; i < alignmentSetHead->alignmentCount; i++){
        scaffoldAlignmentIndex[alignmentSetHead->alignmentSet[i].scaffoldIndex].contigAlignmentIndex[alignmentSetHead->alignmentSet[i].contigIndex].alignmentIndex[scaffoldAlignmentIndex[alignmentSetHead->alignmentSet[i].scaffoldIndex].contigAlignmentIndex[alignmentSetHead->alignmentSet[i].contigIndex].alignmentIndexCount] = i;
        scaffoldAlignmentIndex[alignmentSetHead->alignmentSet[i].scaffoldIndex].contigAlignmentIndex[alignmentSetHead->alignmentSet[i].contigIndex].alignmentIndexCount++;
    }
    
    return scaffoldAlignmentIndex;
    
}

int OptimizeScaffoldAlignmentIndex(ContigSetHead * referenceContigSetHead, ScaffoldSetHead * scaffoldSetHead, ScaffoldAlignmentIndex * scaffoldAlignmentIndex, AlignmentSetHead * alignmentSetHead){
    
    for(long int i = 0; i < scaffoldSetHead->scaffoldCount; i++){
        if(scaffoldSetHead->scaffoldSet[i].gapCount <=0){
            continue;
        }
        long int tempUniqueLength = 0;
        long int tempUniqueMaxIndex = -1;
        long int tempLength = 0;
        long int tempMaxIndex = -1;
        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].contigCount; j++){
            if(scaffoldAlignmentIndex[i].contigAlignmentIndex[j].alignmentIndexCount == 0){
                continue;
            }
            if(scaffoldAlignmentIndex[i].contigAlignmentIndex[j].alignmentIndexCount == 1){
                if(scaffoldSetHead->scaffoldSet[i].contigLength[j] > tempUniqueLength){
                    tempUniqueLength = scaffoldSetHead->scaffoldSet[i].contigLength[j];
                    tempUniqueMaxIndex = j;
                }
            }
            if(scaffoldSetHead->scaffoldSet[i].contigLength[j] > tempLength){
                tempLength = scaffoldSetHead->scaffoldSet[i].contigLength[j];
                tempMaxIndex = j;
            }
        }
        if(tempUniqueMaxIndex == -1 && tempMaxIndex == -1){
            continue;
        }
        if(tempUniqueMaxIndex == -1 && tempMaxIndex != -1){
            long int contigOverlap = 0;
            long int tempIndex = -1;
            for(long int p = 0; p < scaffoldAlignmentIndex[i].contigAlignmentIndex[tempMaxIndex].alignmentIndexCount; p++){
                long int alignmentIndex = scaffoldAlignmentIndex[i].contigAlignmentIndex[tempMaxIndex].alignmentIndex[p];
                if(alignmentSetHead->alignmentSet[alignmentIndex].contigOverlapLength > contigOverlap){
                    contigOverlap = alignmentSetHead->alignmentSet[alignmentIndex].contigOverlapLength;
                    tempIndex = alignmentIndex;
                }
            }
            free(scaffoldAlignmentIndex[i].contigAlignmentIndex[tempMaxIndex].alignmentIndex);
            scaffoldAlignmentIndex[i].contigAlignmentIndex[tempMaxIndex].alignmentIndex = (long int *)malloc(sizeof(long int));
            scaffoldAlignmentIndex[i].contigAlignmentIndex[tempMaxIndex].alignmentIndexCount = 1;
            scaffoldAlignmentIndex[i].contigAlignmentIndex[tempMaxIndex].alignmentIndex[0] = tempIndex;
            tempUniqueMaxIndex = tempMaxIndex;
        }
        long int noMapContigCount = 0;
        if(tempUniqueMaxIndex != -1){
            for(long int p = tempUniqueMaxIndex + 1; p < scaffoldSetHead->scaffoldSet[i].contigCount; p++){
                if(scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount == 0){
                    continue;
                }
                long int tempP = p - 1;
                noMapContigCount = 0;
                while(tempP>=tempUniqueMaxIndex){
                    if(scaffoldAlignmentIndex[i].contigAlignmentIndex[tempP].alignmentIndexCount>0){
                        break;
                    }
                    noMapContigCount++;
                    tempP--;
                }
                long int tempUniqueMaxIndexAlignmentIndex = scaffoldAlignmentIndex[i].contigAlignmentIndex[p - 1 - noMapContigCount].alignmentIndex[0];
                if(scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount > 1){
                    long int maxContigOverlapLength = 0;
                    long int alignmentIndex = -1;
                    for(long int q = 0; q < scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount; q++){
                        alignmentIndex = scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q];
                        if(alignmentSetHead->alignmentSet[alignmentIndex].contigOverlapLength > maxContigOverlapLength){
                            maxContigOverlapLength = alignmentSetHead->alignmentSet[alignmentIndex].contigOverlapLength;
                        }
                    }
                    maxContigOverlapLength = alignmentSetHead->alignmentSet[alignmentIndex].contigLength*0.96;
                    long int * tempDistance = (long int *)malloc(sizeof(long int)*scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount);
                    long int tempCount = 0;
                    alignmentIndex = 0;

                    for(long int q = 0; q < scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount; q++){
                        alignmentIndex = scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q];
                        if(SearchReferenceIndexFromName(referenceContigSetHead->contigSet, referenceContigSetHead->contigCount, alignmentSetHead->alignmentSet[alignmentIndex].referenceName) 
                           != SearchReferenceIndexFromName(referenceContigSetHead->contigSet, referenceContigSetHead->contigCount, alignmentSetHead->alignmentSet[tempUniqueMaxIndexAlignmentIndex].referenceName)){
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q] = -1;
                            tempDistance[q] = -1;
                            continue;
                        }
                        if(alignmentSetHead->alignmentSet[alignmentIndex].isReverse != alignmentSetHead->alignmentSet[tempUniqueMaxIndexAlignmentIndex].isReverse){
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q] = -1;
                            tempDistance[q] = -1;
                            continue;
                        }
                        if(alignmentSetHead->alignmentSet[alignmentIndex].contigOverlapLength < maxContigOverlapLength){
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q] = -1;
                            tempDistance[q] = -1;
                            continue;
                        }
                        tempDistance[q] = alignmentSetHead->alignmentSet[alignmentIndex].referenceStart - alignmentSetHead->alignmentSet[tempUniqueMaxIndexAlignmentIndex].referenceEnd;
                        tempCount++;
                    }
                    if(tempCount == 0){
                        free(scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex);
                        scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex = NULL;
                        scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount = 0;
                        continue;
                    }
                    alignmentIndex = -1;
                    long int minDistance = -1;
                    long int alignmentIndex1 = -1;
                    long int maxDistance = 0;
                    for(long int q = 0; q < scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount; q++){
                        if(scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q] == -1){
                            continue;
                        }
                        if(tempDistance[q] >=0 && (tempDistance[q] < minDistance || minDistance == -1)){
                            minDistance = tempDistance[q];
                            alignmentIndex = scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q];
                        }
                        if(tempDistance[q] < 0 && (tempDistance[q] > maxDistance || maxDistance == 0)){
                            maxDistance = tempDistance[q];
                            alignmentIndex1 = scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q];
                        }
                    }
                    free(scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex);
                    scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex = (long int *)malloc(sizeof(long int));
                    scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount = 1;
                    if(alignmentSetHead->alignmentSet[tempUniqueMaxIndexAlignmentIndex].isReverse == false){
                        if(alignmentIndex != -1){
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[0] = alignmentIndex;
                        }else if(alignmentIndex1 != -1){
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[0] = alignmentIndex1;
                        }else{
                            free(scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex);
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex = NULL;
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount = 0;
                        }
                    }else{
                        if(alignmentIndex1 != -1){
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[0] = alignmentIndex1;
                        }else if(alignmentIndex != -1){
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[0] = alignmentIndex;
                        }else{
                            free(scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex);
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex = NULL;
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount = 0;
                        }
                    }
                    free(tempDistance);
                    tempDistance = NULL;
                }
            }
            noMapContigCount = 0;
            for(long int p = tempUniqueMaxIndex - 1; p >= 0; p--){
                
                if(scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount == 0){
                    continue;
                }
                long int tempP = p + 1;
                noMapContigCount = 0;
                while(tempP<=tempUniqueMaxIndex){
                    if(scaffoldAlignmentIndex[i].contigAlignmentIndex[tempP].alignmentIndexCount>0){
                        break;
                    }
                    noMapContigCount++;
                    tempP++;
                }
                
                long int tempUniqueMaxIndexAlignmentIndex = scaffoldAlignmentIndex[i].contigAlignmentIndex[p + 1 + noMapContigCount].alignmentIndex[0];
                
                if(scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount > 1){
                    long int maxContigOverlapLength = 0;
                    long int alignmentIndex = -1;
                    for(long int q = 0; q < scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount; q++){
                        alignmentIndex = scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q];
                        if(alignmentSetHead->alignmentSet[alignmentIndex].contigOverlapLength > maxContigOverlapLength){
                            maxContigOverlapLength = alignmentSetHead->alignmentSet[alignmentIndex].contigOverlapLength;
                        }
                    }
                    maxContigOverlapLength = alignmentSetHead->alignmentSet[alignmentIndex].contigLength*0.96;
                    long int * tempDistance = (long int *)malloc(sizeof(long int)*scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount);
                    long int tempCount = 0;
                    alignmentIndex = 0;
                    for(long int q = 0; q < scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount; q++){
                        alignmentIndex = scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q];
                        if(SearchReferenceIndexFromName(referenceContigSetHead->contigSet, referenceContigSetHead->contigCount, alignmentSetHead->alignmentSet[alignmentIndex].referenceName) 
                           != SearchReferenceIndexFromName(referenceContigSetHead->contigSet, referenceContigSetHead->contigCount, alignmentSetHead->alignmentSet[tempUniqueMaxIndexAlignmentIndex].referenceName)){
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q] = -1;
                            tempDistance[q] = -1;
                            continue;
                        }
                        
                        if(alignmentSetHead->alignmentSet[alignmentIndex].isReverse != alignmentSetHead->alignmentSet[tempUniqueMaxIndexAlignmentIndex].isReverse){

                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q] = -1;
                            tempDistance[q] = -1;
                            continue;
                        }

                        if(alignmentSetHead->alignmentSet[alignmentIndex].contigOverlapLength < maxContigOverlapLength){
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q] = -1;
                            tempDistance[q] = -1;
                            continue;
                        }

                        tempDistance[q] = alignmentSetHead->alignmentSet[tempUniqueMaxIndexAlignmentIndex].referenceStart - alignmentSetHead->alignmentSet[alignmentIndex].referenceEnd;
                        tempCount++;
                    }

                    if(tempCount == 0){
                        free(scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex);
                        scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex = NULL;
                        scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount = 0;
                        continue;
                    }
                    alignmentIndex = -1;
                    long int minDistance = -1;
                    long int alignmentIndex1 = -1;
                    long int maxDistance = 0;
                    for(long int q = 0; q < scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount; q++){
                        if(scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q] == -1){
                            continue;
                        }
                        if(tempDistance[q] >=0 && (tempDistance[q] < minDistance || minDistance == -1)){
                            minDistance = tempDistance[q];
                            alignmentIndex = scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q];
                        }
                        if(tempDistance[q] < 0 && (tempDistance[q] > maxDistance || maxDistance == 0)){
                            maxDistance = tempDistance[q];
                            alignmentIndex1 = scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[q];
                        }
                    }

                    free(scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex);
                    scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex = (long int *)malloc(sizeof(long int));
                    scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount = 1;
                    if(alignmentSetHead->alignmentSet[tempUniqueMaxIndexAlignmentIndex].isReverse == false){
                        if(alignmentIndex != -1){
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[0] = alignmentIndex;
                        }else if(alignmentIndex1 != -1){
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[0] = alignmentIndex1;
                        }else{
                            free(scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex);
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex = NULL;
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount = 0;
                        }
                    }else{
                        if(alignmentIndex1 != -1){
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[0] = alignmentIndex1;
                        }else if(alignmentIndex != -1){
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex[0] = alignmentIndex;
                        }else{
                            free(scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex);
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndex = NULL;
                            scaffoldAlignmentIndex[i].contigAlignmentIndex[p].alignmentIndexCount = 0;
                        }
                    }

                    free(tempDistance);
                    tempDistance = NULL;
                }
            }

        }
        
    }
    return 1;
    
}

long int MaxDistance(long int a, long int b){
    if(a >= b){
        return a;
    }else{
        return b;
    }
}

int GetGapRegion(ContigSetHead * referenceContigSetHead, ScaffoldSetHead * scaffoldSetHead, ScaffoldAlignmentIndex * scaffoldAlignmentIndex, AlignmentSetHead * alignmentSetHead, ScaffoldGapRegion * scaffoldGapRegion, long int scaffoldIndex, long int leftContigIndex, long int rightContigIndex, long int gapIndex, bool isReference){
    
    long int leftContigAlignmentIndexCount = 0;   
    long int rightContigAlignmentIndexCount = 0;
    
    long int * leftContigAlignmentIndex = NULL;
    if(leftContigIndex != -1){
        leftContigAlignmentIndex = scaffoldAlignmentIndex[scaffoldIndex].contigAlignmentIndex[leftContigIndex].alignmentIndex;
        leftContigAlignmentIndexCount = scaffoldAlignmentIndex[scaffoldIndex].contigAlignmentIndex[leftContigIndex].alignmentIndexCount; 
    }
    long int * rightContigAlignmentIndex = NULL;
    if(rightContigIndex != -1){
        rightContigAlignmentIndex = scaffoldAlignmentIndex[scaffoldIndex].contigAlignmentIndex[rightContigIndex].alignmentIndex;
        rightContigAlignmentIndexCount = scaffoldAlignmentIndex[scaffoldIndex].contigAlignmentIndex[rightContigIndex].alignmentIndexCount;
    }    
    
    long int leftAlignmentIndex = -1;
    long int rightAlignmentIndex = -1;
    long int tempPosition = -1;
    double identityPercent = 0; 
    long int tempContigOverlapLength = 0;
    
    for(long int i = 0; i < leftContigAlignmentIndexCount; i++){
        if(alignmentSetHead->alignmentSet[leftContigAlignmentIndex[i]].isReverse == 0){
            if(isReference == true && alignmentSetHead->alignmentSet[leftContigAlignmentIndex[i]].contigEnd < tempPosition){
                continue;
            }
            tempPosition = alignmentSetHead->alignmentSet[leftContigAlignmentIndex[i]].contigEnd;
            tempContigOverlapLength = alignmentSetHead->alignmentSet[leftContigAlignmentIndex[i]].contigOverlapLength;
            leftAlignmentIndex = leftContigAlignmentIndex[i];
            identityPercent = alignmentSetHead->alignmentSet[leftContigAlignmentIndex[i]].identityPercent; 
        }
        if(alignmentSetHead->alignmentSet[leftContigAlignmentIndex[i]].isReverse == 1){
            if(isReference == true && alignmentSetHead->alignmentSet[leftContigAlignmentIndex[i]].contigStart < tempPosition){
                continue;
            }
            tempPosition = alignmentSetHead->alignmentSet[leftContigAlignmentIndex[i]].contigStart;
            tempContigOverlapLength = alignmentSetHead->alignmentSet[leftContigAlignmentIndex[i]].contigOverlapLength;
            leftAlignmentIndex = leftContigAlignmentIndex[i];
            identityPercent = alignmentSetHead->alignmentSet[leftContigAlignmentIndex[i]].identityPercent; 
        }
    }
    
    if(leftAlignmentIndex != -1){
        if(alignmentSetHead->alignmentSet[leftAlignmentIndex].isReverse != 1){
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition = 
                alignmentSetHead->alignmentSet[leftAlignmentIndex].referenceEnd + 1;
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftContigEndPosition = 
                alignmentSetHead->alignmentSet[leftAlignmentIndex].contigEnd;
            if(alignmentSetHead->alignmentSet[leftAlignmentIndex].contigLength - alignmentSetHead->alignmentSet[leftAlignmentIndex].contigEnd - 1 > scaffoldSetHead->minSegmentDistanceEndLength){
                leftAlignmentIndex = -1;
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition = -1;
            }else{
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition = scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition
                    + alignmentSetHead->alignmentSet[leftAlignmentIndex].contigLength - alignmentSetHead->alignmentSet[leftAlignmentIndex].contigEnd - 1;
            }
            
        }else{
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition = 
                alignmentSetHead->alignmentSet[leftAlignmentIndex].referenceStart - 1;
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftContigEndPosition = 
                alignmentSetHead->alignmentSet[leftAlignmentIndex].contigStart;
            if(alignmentSetHead->alignmentSet[leftAlignmentIndex].contigLength - alignmentSetHead->alignmentSet[leftAlignmentIndex].contigStart - 1 > scaffoldSetHead->minSegmentDistanceEndLength){
                leftAlignmentIndex = -1; 
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition = -1;
            }else{
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition = scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition
                   - alignmentSetHead->alignmentSet[leftAlignmentIndex].contigLength + alignmentSetHead->alignmentSet[leftAlignmentIndex].contigStart + 1;
            }
        }
        if(leftAlignmentIndex != -1){
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftReverse = 
                alignmentSetHead->alignmentSet[leftAlignmentIndex].isReverse;
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftReferenceIndex = 
                SearchReferenceIndexFromName(referenceContigSetHead->contigSet, referenceContigSetHead->contigCount, alignmentSetHead->alignmentSet[leftAlignmentIndex].referenceName);
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftReferenceName = (char *)malloc(sizeof(char)*(strlen(alignmentSetHead->alignmentSet[leftAlignmentIndex].referenceName) + 1));
                strcpy(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftReferenceName, alignmentSetHead->alignmentSet[leftAlignmentIndex].referenceName);
        }
        
    }
    
    identityPercent = 0; 
    tempContigOverlapLength = 0;
    tempPosition = scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength;
    for(long int i = 0; i < rightContigAlignmentIndexCount; i++){
        
        int token = 0;
        if(alignmentSetHead->alignmentSet[rightContigAlignmentIndex[i]].isReverse == 0){
            if(token == 0 && isReference == true && alignmentSetHead->alignmentSet[rightContigAlignmentIndex[i]].contigStart > tempPosition){
                continue;
            }
            tempPosition = alignmentSetHead->alignmentSet[rightContigAlignmentIndex[i]].contigStart;
            tempContigOverlapLength = alignmentSetHead->alignmentSet[rightContigAlignmentIndex[i]].contigOverlapLength;
            rightAlignmentIndex = rightContigAlignmentIndex[i];
            identityPercent = alignmentSetHead->alignmentSet[rightContigAlignmentIndex[i]].identityPercent;
        }
        if(alignmentSetHead->alignmentSet[rightContigAlignmentIndex[i]].isReverse == 1){
            if(token == 0 && isReference == true && alignmentSetHead->alignmentSet[rightContigAlignmentIndex[i]].contigEnd > tempPosition){
                continue;
            }
            tempPosition = alignmentSetHead->alignmentSet[rightContigAlignmentIndex[i]].contigEnd;
            tempContigOverlapLength = alignmentSetHead->alignmentSet[rightContigAlignmentIndex[i]].contigOverlapLength;
            rightAlignmentIndex = rightContigAlignmentIndex[i];
            identityPercent = alignmentSetHead->alignmentSet[rightContigAlignmentIndex[i]].identityPercent;
        }
    }
    if(rightAlignmentIndex != -1){
        if(alignmentSetHead->alignmentSet[rightAlignmentIndex].isReverse != 1){
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition = 
                alignmentSetHead->alignmentSet[rightAlignmentIndex].referenceStart - 1;
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightContigStartPosition = 
                alignmentSetHead->alignmentSet[rightAlignmentIndex].contigStart;
            if(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightContigStartPosition > scaffoldSetHead->minSegmentDistanceEndLength){
                rightAlignmentIndex = -1;
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition = -1;
            }else{
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition = scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition
                    - alignmentSetHead->alignmentSet[rightAlignmentIndex].contigStart;
            }
        }else{
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition = 
                alignmentSetHead->alignmentSet[rightAlignmentIndex].referenceEnd + 1;
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightContigStartPosition = 
                alignmentSetHead->alignmentSet[rightAlignmentIndex].contigEnd;
            if(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightContigStartPosition > scaffoldSetHead->minSegmentDistanceEndLength){
                rightAlignmentIndex = -1;
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition = -1;
            }else{
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition = scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition
                    + alignmentSetHead->alignmentSet[rightAlignmentIndex].contigEnd;
            }
            
        }
        if(rightAlignmentIndex != -1){
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightReverse = 
                alignmentSetHead->alignmentSet[rightAlignmentIndex].isReverse;
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightReferenceIndex = 
                SearchReferenceIndexFromName(referenceContigSetHead->contigSet, referenceContigSetHead->contigCount, alignmentSetHead->alignmentSet[rightAlignmentIndex].referenceName);
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightReferenceName = (char *)malloc(sizeof(char)*(strlen(alignmentSetHead->alignmentSet[rightAlignmentIndex].referenceName) + 1));
            strcpy(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightReferenceName, alignmentSetHead->alignmentSet[rightAlignmentIndex].referenceName);
        }
        
    }
    
    
    if(leftContigIndex == -1 && rightAlignmentIndex != -1){
        scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 0;
        scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftReverse = scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightReverse;
        scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftReferenceIndex = scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightReferenceIndex;
        scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftContigEndPosition = 0;    
        scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftReferenceName = scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightReferenceName;
        if(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightReverse != 1){
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition = scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition 
                - scaffoldSetHead->scaffoldSet[scaffoldIndex].gapDistance[gapIndex] + 1;
            if(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition < 0){
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition = 0;
            }
        }else{
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition = scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition 
                + scaffoldSetHead->scaffoldSet[scaffoldIndex].gapDistance[gapIndex] - 1;
            if(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition < 0){
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition = 0;
            }
            
        }
        return 0;
    }
    
    if(rightContigIndex == -1 && leftAlignmentIndex != -1){
        scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 0;
        scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightReverse = scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftReverse;
        scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightReferenceIndex = scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftReferenceIndex;
        scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightContigStartPosition = 0;    
        scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightReferenceName = scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftReferenceName;
        if(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightReverse != 1){
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition = scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition 
                + scaffoldSetHead->scaffoldSet[scaffoldIndex].gapDistance[gapIndex] - 1;
            if(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition < 0){
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition = 0;
            }
        }else{
            scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition = scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition 
                - scaffoldSetHead->scaffoldSet[scaffoldIndex].gapDistance[gapIndex] + 1;
            if(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition < 0){
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition = 0;
            }
            
        }
        return 0;
    }
    
    if(leftAlignmentIndex != -1 && rightAlignmentIndex == -1){
        scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 21;//missing gap;
        return 21;
    }
    
    if(leftAlignmentIndex == -1 && rightAlignmentIndex != -1){
        scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 22;//missing gap;
        return 22;
    }
    
    if(leftAlignmentIndex == -1 && rightAlignmentIndex == -1){
        scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 23;//missing gap;
        return 23;
    }
    
    
    if(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightReferenceIndex != scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftReferenceIndex){
        scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 11;//translocation gap;
        return 11;
    }
    
    if(alignmentSetHead->alignmentSet[leftAlignmentIndex].isReverse != alignmentSetHead->alignmentSet[rightAlignmentIndex].isReverse){
        scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 12;//inversion gap;
        return 12;
    }
    
    if(isReference == true){
        if(alignmentSetHead->alignmentSet[leftAlignmentIndex].isReverse == 0){
        
            if(labs(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition - scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition + 1) <= 
                MaxDistance(scaffoldSetHead->minTimesRelocation*scaffoldSetHead->scaffoldSet[scaffoldIndex].gapDistance[gapIndex], scaffoldSetHead->minDistanceRelocationLength) &&
            labs(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition - scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition + 1) > 0){
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 0;
                return 0;
            }else{
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 1;
                return 1;
            }
            
        }
        if(alignmentSetHead->alignmentSet[leftAlignmentIndex].isReverse == 1){
            
            if(labs(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition - scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition + 1) <= 
                MaxDistance(scaffoldSetHead->minTimesRelocation*scaffoldSetHead->scaffoldSet[scaffoldIndex].gapDistance[gapIndex], scaffoldSetHead->minDistanceRelocationLength) && 
            labs(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition - scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition + 1) > 0){
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 0;
                return 0;
            }else{
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 3;
                return 3;
            }
            
        }
    }else{
        if(alignmentSetHead->alignmentSet[leftAlignmentIndex].isReverse == 0){
        
            if(labs(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition - scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition + 1) > 0){
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 0;
                return 0;
            }else{
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 1;
                return 1;
            }
            
        }
        if(alignmentSetHead->alignmentSet[leftAlignmentIndex].isReverse == 1){
            
            if(labs(scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].leftStartPosition - scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].rightStartPosition + 1) > 0){
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 0;
                return 0;
            }else{
                scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 3;
                return 3;
            }
            
        }
    }
    scaffoldGapRegion[scaffoldIndex].gapRegionSet[gapIndex].gapType = 0;//correct gap;
    return 0;
}

ScaffoldGapRegion * GetScaffoldGapRegionInReference(ScaffoldSetHead * scaffoldSetHead, AlignmentSetHead * alignmentSetHead, ContigSetHead * referenceContigSetHead, bool isReference){

    ScaffoldGapRegion * scaffoldGapRegion = (ScaffoldGapRegion *)malloc(sizeof(ScaffoldGapRegion)*scaffoldSetHead->scaffoldCount);
    for(long int i = 0; i < scaffoldSetHead->scaffoldCount; i++){
        
        if(scaffoldSetHead->scaffoldSet[i].gapCount <=0){
            scaffoldGapRegion[i].gapRegionSet = NULL;
            continue;
        }
        
        scaffoldGapRegion[i].gapRegionSet = (GapRegion *)malloc(sizeof(GapRegion)*scaffoldSetHead->scaffoldSet[i].gapCount);
        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++){
            scaffoldGapRegion[i].gapRegionSet[j].gapType = -1;
            scaffoldGapRegion[i].gapRegionSet[j].leftStartPosition = -1;
            scaffoldGapRegion[i].gapRegionSet[j].leftContigEndPosition = -1;
            scaffoldGapRegion[i].gapRegionSet[j].rightStartPosition = -1;
            scaffoldGapRegion[i].gapRegionSet[j].rightContigStartPosition = -1;
            scaffoldGapRegion[i].gapRegionSet[j].leftReverse = 0;
            scaffoldGapRegion[i].gapRegionSet[j].rightReverse = 0;
            scaffoldGapRegion[i].gapRegionSet[j].leftReferenceIndex = -1;
            scaffoldGapRegion[i].gapRegionSet[j].rightReferenceIndex = -1;
            scaffoldGapRegion[i].gapRegionSet[j].leftReferenceName = NULL;
            scaffoldGapRegion[i].gapRegionSet[j].rightReferenceName = NULL;
        }
    }
    
    ScaffoldAlignmentIndex * scaffoldAlignmentIndex = GetScaffoldAlignmentIndex(scaffoldSetHead, alignmentSetHead);
    
    if(isReference != true){
        OptimizeScaffoldAlignmentIndex(referenceContigSetHead, scaffoldSetHead, scaffoldAlignmentIndex, alignmentSetHead);
    }
    

    for(long int i = 0; i < scaffoldSetHead->scaffoldCount; i++){
        if(scaffoldSetHead->scaffoldSet[i].gapCount <=0){
            continue;
        }
        long int startIndex = 0;
        long int endIndex = 0;
        if(scaffoldSetHead->scaffoldSet[i].gapStartPosition[0] == 0){
            startIndex = -1;
        }

        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount - 1; j++){
            endIndex = startIndex + 1;
            GetGapRegion(referenceContigSetHead, scaffoldSetHead, scaffoldAlignmentIndex, alignmentSetHead, scaffoldGapRegion, i, startIndex, endIndex, j, isReference);
            startIndex++;
        }
        endIndex = startIndex + 1;
        if(endIndex > scaffoldSetHead->scaffoldSet[i].contigCount - 1){
            GetGapRegion(referenceContigSetHead, scaffoldSetHead, scaffoldAlignmentIndex, alignmentSetHead, scaffoldGapRegion, i, startIndex, -1, scaffoldSetHead->scaffoldSet[i].gapCount - 1, isReference);
        }else{
            GetGapRegion(referenceContigSetHead, scaffoldSetHead, scaffoldAlignmentIndex, alignmentSetHead, scaffoldGapRegion, i, startIndex, endIndex, scaffoldSetHead->scaffoldSet[i].gapCount - 1, isReference);
        }
    }
    return scaffoldGapRegion;
}

long int SearchReferenceIndexFromName(Contig * reference, long int referenceNumber, char * referenceName){ 
    
    for(long int i = 0; i < referenceNumber; i++){
        if(strlen(reference[i].contigName) < strlen(referenceName)){
            continue;
        }else if(strlen(reference[i].contigName) == strlen(referenceName)){
            if(strcmp(reference[i].contigName, referenceName) == 0){
                return i;
            }
        }else if(strlen(reference[i].contigName) > strlen(referenceName)){
            long int length = strlen(referenceName);
            long int j = 0;
            for(j = 0; j < length; j++){
                if(reference[i].contigName[j] != referenceName[j]){
                    break;
                }
            }
            if(j != length){
                continue;
            }else{
                return i;
            }
        }
    }
    return -1;
}

void OutputScaffoldGapRegion(ScaffoldSetHead * scaffoldSetHead, ScaffoldGapRegion * scaffoldGapRegion, ContigSetHead * referenceContigSetHead, char * gapFileName){

    FILE * fp;
    if((fp = fopen(gapFileName, "w")) == NULL){
        printf("%s, can not open!", gapFileName);
        exit(0);
    }
    long int gapIndex = 0;
    for(long int i = 0; i < scaffoldSetHead->scaffoldCount; i++){
        if(scaffoldSetHead->scaffoldSet[i].gapCount <=0){
            continue;
        }
        
        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++){
            fprintf(fp, ">scaffold_%ld_%ld_%ld_%d_0\n", i, j, gapIndex, scaffoldGapRegion[i].gapRegionSet[j].gapType);
            if(scaffoldGapRegion[i].gapRegionSet[j].gapType == 0){
                long int referenceIndex = scaffoldGapRegion[i].gapRegionSet[j].leftReferenceIndex;
                long int gapLength = labs(scaffoldGapRegion[i].gapRegionSet[j].rightStartPosition - scaffoldGapRegion[i].gapRegionSet[j].leftStartPosition) + 1;
                
                char * tempGap = (char *)malloc(sizeof(char)*(gapLength + 1));
                if(scaffoldGapRegion[i].gapRegionSet[j].leftReverse != true){
                    strncpy(tempGap, referenceContigSetHead->contigSet[referenceIndex].contig + scaffoldGapRegion[i].gapRegionSet[j].leftStartPosition, gapLength);
                    tempGap[gapLength] = '\0';
                }else{
                    strncpy(tempGap, referenceContigSetHead->contigSet[referenceIndex].contig + scaffoldGapRegion[i].gapRegionSet[j].rightStartPosition, gapLength);
                    tempGap[gapLength] = '\0';
                    ReverseComplement(tempGap);
                }
                fprintf(fp, "%s\n", tempGap);
                free(tempGap);
            }
            gapIndex++;
        }
    }
    fflush(fp);
    fclose(fp);
}


void OutputReferenceGapRegion(ScaffoldSetHead * scaffoldSetHead, ScaffoldGapRegion * referenceGapRegion, ScaffoldGapRegion * scaffoldGapRegion, ContigSetHead * referenceContigSetHead, char * gapFileName){

    FILE * fp;
    if((fp = fopen(gapFileName, "w")) == NULL){
        printf("%s, can not open!", gapFileName);
        exit(0);
    }
    long int gapIndex = 0;
    for(long int i = 0; i < scaffoldSetHead->scaffoldCount; i++){
        if(scaffoldSetHead->scaffoldSet[i].gapCount <=0){
            continue;
        }
        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++){
            long int count = OutputMinContigAsGap(scaffoldSetHead, i, j);
            fprintf(fp, ">scaffold_%ld_%ld_%ld_%d_%ld\n", i, j, gapIndex, referenceGapRegion[i].gapRegionSet[j].gapType, count);
            fflush(stdout);
            if(referenceGapRegion[i].gapRegionSet[j].gapType == 0){
                long int referenceIndex = referenceGapRegion[i].gapRegionSet[j].leftReferenceIndex;
                long int gapLength = labs(referenceGapRegion[i].gapRegionSet[j].rightStartPosition - referenceGapRegion[i].gapRegionSet[j].leftStartPosition) + 1;
                char * tempGap = (char *)malloc(sizeof(char)*(gapLength + 1));
                if(referenceGapRegion[i].gapRegionSet[j].leftReverse != true){
                    strncpy(tempGap, referenceContigSetHead->contigSet[referenceIndex].contig + referenceGapRegion[i].gapRegionSet[j].leftStartPosition, gapLength);
                    tempGap[gapLength] = '\0';
                }else{
                    strncpy(tempGap, referenceContigSetHead->contigSet[referenceIndex].contig + referenceGapRegion[i].gapRegionSet[j].rightStartPosition, gapLength);
                    tempGap[gapLength] = '\0';
                    ReverseComplement(tempGap);
                }
                fprintf(fp, "%s\n", tempGap);
                free(tempGap);
            }else{
                long int fillingGapLength = labs(scaffoldGapRegion[i].gapRegionSet[j].rightStartPosition - scaffoldGapRegion[i].gapRegionSet[j].leftStartPosition) + 1;
                if(referenceGapRegion[i].gapRegionSet[j].leftReferenceIndex != -1){
                    long int referenceIndex = referenceGapRegion[i].gapRegionSet[j].leftReferenceIndex;
                    long int gapLength = fillingGapLength;
                    if(referenceGapRegion[i].gapRegionSet[j].gapType == 2 || referenceGapRegion[i].gapRegionSet[j].gapType == 4){
                        gapLength = labs(referenceGapRegion[i].gapRegionSet[j].rightStartPosition - referenceGapRegion[i].gapRegionSet[j].leftStartPosition) + 1;
                    }else if(scaffoldGapRegion[i].gapRegionSet[j].gapType > 0){
                        gapLength = scaffoldSetHead->scaffoldSet[i].gapDistance[j];
                    }
                
                    char * tempGap = (char *)malloc(sizeof(char)*(gapLength + 1));
                    if(referenceGapRegion[i].gapRegionSet[j].leftReverse != true){
                        strncpy(tempGap, referenceContigSetHead->contigSet[referenceIndex].contig + referenceGapRegion[i].gapRegionSet[j].leftStartPosition, gapLength);
                        tempGap[gapLength] = '\0';
                    }else{
                        strncpy(tempGap, referenceContigSetHead->contigSet[referenceIndex].contig + referenceGapRegion[i].gapRegionSet[j].leftStartPosition - gapLength + 1, gapLength);
                        tempGap[gapLength] = '\0';
                        ReverseComplement(tempGap);
                    }
                    fprintf(fp, "%s\n", tempGap);
                    free(tempGap);
                }
                if(referenceGapRegion[i].gapRegionSet[j].rightReferenceIndex != -1){
                    long int referenceIndex = referenceGapRegion[i].gapRegionSet[j].rightReferenceIndex;
                    long int gapLength = fillingGapLength;
                    if(referenceGapRegion[i].gapRegionSet[j].gapType == 2 || referenceGapRegion[i].gapRegionSet[j].gapType == 4){
                        gapLength = labs(referenceGapRegion[i].gapRegionSet[j].rightStartPosition - referenceGapRegion[i].gapRegionSet[j].leftStartPosition + 1);
                    }else if(scaffoldGapRegion[i].gapRegionSet[j].gapType > 0){
                        gapLength = scaffoldSetHead->scaffoldSet[i].gapDistance[j];
                    }
                    
                    char * tempGap = (char *)malloc(sizeof(char)*(gapLength + 1));
                    if(referenceGapRegion[i].gapRegionSet[j].rightReverse != true){
                        strncpy(tempGap, referenceContigSetHead->contigSet[referenceIndex].contig + referenceGapRegion[i].gapRegionSet[j].rightStartPosition - gapLength + 1, gapLength);
                        tempGap[gapLength] = '\0';
                    }else{
                        strncpy(tempGap, referenceContigSetHead->contigSet[referenceIndex].contig + referenceGapRegion[i].gapRegionSet[j].rightStartPosition, gapLength);
                        tempGap[gapLength] = '\0';
                        ReverseComplement(tempGap);
                    }
                    fprintf(fp, "%s\n", tempGap);
                    free(tempGap);
                }
            }
            fflush(fp);
            gapIndex++;
        }
    }
    fflush(fp);
    fclose(fp);
}

long int OutputMinContigAsGap(ScaffoldSetHead * scaffoldSetHead, long int i, long int j){
    long int count = 0;
    long int start = scaffoldSetHead->scaffoldSet[i].gapStartPosition[j];
    long int len = scaffoldSetHead->scaffoldSet[i].gapDistance[j];
    for(long int p = 0; p < len; p++){
        if(scaffoldSetHead->scaffoldSet[i].scaffold[start + p] != 'N' && scaffoldSetHead->scaffoldSet[i].scaffold[start + p] != 'n'){
            count++;
        }
    }
    return count;
}


bool ReverseComplement(char * temp){
    long int len = strlen(temp);
    char a;
    char b;
    long int m = len/2;
    for(long int i = 0; i < len/2; i++){
        a = temp[i];
        b = temp[len - 1 - i];
        if(a == 'A'){
            temp[len - 1 - i] = 'T';
        }else if(a == 'T'){
            temp[len - 1 - i] = 'A';
        }else if(a == 'G'){
            temp[len - 1 - i] = 'C';
        }else if(a == 'C'){
            temp[len - 1 - i] = 'G';
        }else if(a == 'N'){
            temp[len - 1 - i] = 'N';
        }
        
        if(b == 'A'){
            temp[i] = 'T';
        }else if(b == 'T'){
            temp[i] = 'A';
        }else if(b == 'G'){
            temp[i] = 'C';
        }else if(b == 'C'){
            temp[i] = 'G';
        }else if(b == 'N'){
            temp[i] = 'N';
        }
        
    }
    if(len%2 != 0){
        if(temp[m] == 'A'){
            temp[m] = 'T';
        }else if(temp[m] == 'T'){
            temp[m] = 'A';
        }else if(temp[m] == 'G'){
            temp[m] = 'C';
        }else if(temp[m] == 'C'){
            temp[m] = 'G';
        }else if(temp[m] == 'N'){
            temp[m] = 'N';
        }
        
        if(temp[m] == 'a'){
            temp[m] = 'T';
        }else if(temp[m] == 't'){
            temp[m] = 'A';
        }else if(temp[m] == 'g'){
            temp[m] = 'C';
        }else if(temp[m] == 'c'){
            temp[m] = 'G';
        }else if(temp[m] == 'n'){
            temp[m] = 'N';
        }
    }
    temp[len]='\0';
    return true;
}


void OutputGapInScaffoldInformation(ScaffoldSetHead * scaffoldSetHead, ScaffoldGapRegion * scaffoldGapRegion, ScaffoldGapRegion * fillingGapRegion, char * fileName){
    long int gapCount = 0;
    long int len = 0;
    FILE * fp;
    fp = fopen(fileName, "w");
    fprintf(fp, "1-th column:\tscaffold name;\n");
    fprintf(fp, "2-th column:\tgap start position in the scaffold;\n");
    fprintf(fp, "3-th column:\tgap end position in the scaffold;\n");
    fprintf(fp, "4-th column:\tgap start position in the reference;\n");
    fprintf(fp, "5-th column:\tgap end position in the reference;\n");
    fprintf(fp, "6-th column:\tgap start position in the gap-closed scaffold;\n");
    fprintf(fp, "7-th column:\tgap end position in the gap-closed scaffold;\n");
    
    fprintf(fp, "8-th column:\tname of reference which the left flanking region is aligned to;\n");
    fprintf(fp, "9-th column:\tname of reference which the right flanking region is aligned to;\n");
    fprintf(fp, "10-th column:\taligning strand of the left flanking region for the reference;\n");
    fprintf(fp, "11-th column:\taligning strand of the right flanking region for the reference;\n");
    fprintf(fp, "12-th column:\tgap type about reference;\n");
    
    fprintf(fp, "13-th column:\tname of gap-closed scaffold which the left flanking region is aligned to;\n");
    fprintf(fp, "14-th column:\tname of gap-closed scaffold which the right flanking region is aligned to;\n");
    fprintf(fp, "15-th column:\taligning strand of the left flanking region for the gap-closed scaffold;\n");
    fprintf(fp, "16-th column:\taligning strand of the right flanking region for the gap-closed scaffold;\n");
    //fprintf(fp, "17-th column:\tgap type about gap-closed scaffold;\n");
    fprintf(fp, "\n");
    for(long int i = 0; i < scaffoldSetHead->scaffoldCount; i++){
        if(scaffoldSetHead->scaffoldSet[i].gapCount <=0){
            continue;
        }
        for(long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++){
            
            fprintf(fp, "%-s\t", scaffoldSetHead->scaffoldSet[i].scaffoldName);
            
            fprintf(fp, "%-ld\t", scaffoldSetHead->scaffoldSet[i].gapStartPosition[j]);
            fprintf(fp, "%-ld\t", scaffoldSetHead->scaffoldSet[i].gapStartPosition[j] + scaffoldSetHead->scaffoldSet[i].gapDistance[j] - 1);
            
            fprintf(fp, "%-ld\t", scaffoldGapRegion[i].gapRegionSet[j].leftStartPosition);
            fprintf(fp, "%-ld\t", scaffoldGapRegion[i].gapRegionSet[j].rightStartPosition);
            
            fprintf(fp, "%-ld\t", fillingGapRegion[i].gapRegionSet[j].leftStartPosition);
            fprintf(fp, "%-ld\t", fillingGapRegion[i].gapRegionSet[j].rightStartPosition);
            
            fprintf(fp, "%-s\t", scaffoldGapRegion[i].gapRegionSet[j].leftReferenceName);
            fprintf(fp, "%-s\t", scaffoldGapRegion[i].gapRegionSet[j].rightReferenceName);
            fprintf(fp, "%-d\t", scaffoldGapRegion[i].gapRegionSet[j].leftReverse);
            fprintf(fp, "%-d\t", scaffoldGapRegion[i].gapRegionSet[j].rightReverse);
            if(scaffoldGapRegion[i].gapRegionSet[j].gapType == 0){
                fprintf(fp, "Normal\t");
            }else if(scaffoldGapRegion[i].gapRegionSet[j].gapType == 1 || scaffoldGapRegion[i].gapRegionSet[j].gapType == 3){
                fprintf(fp, "Relocation\t");
            }else if(scaffoldGapRegion[i].gapRegionSet[j].gapType == 21 || scaffoldGapRegion[i].gapRegionSet[j].gapType == 22 || scaffoldGapRegion[i].gapRegionSet[j].gapType == 23){
                fprintf(fp, "Missing\t");
            }else if(scaffoldGapRegion[i].gapRegionSet[j].gapType == 11){
                fprintf(fp, "Translocation\t");
            }else if(scaffoldGapRegion[i].gapRegionSet[j].gapType == 12){
                fprintf(fp, "Inversion\t");
            }else{
                fprintf(fp, "undefined\t");
            }
            
            fprintf(fp, "%-s\t", fillingGapRegion[i].gapRegionSet[j].leftReferenceName);
            fprintf(fp, "%-s\t", fillingGapRegion[i].gapRegionSet[j].rightReferenceName);
            fprintf(fp, "%-d\t", fillingGapRegion[i].gapRegionSet[j].leftReverse);
            fprintf(fp, "%-d\t", fillingGapRegion[i].gapRegionSet[j].rightReverse);
            //fprintf(fp, "%-d\t", fillingGapRegion[i].gapRegionSet[j].gapType);
            
            fprintf(fp, "\n");
            gapCount++;
        }
        
    }

    fflush(fp);
    fclose(fp);
}

