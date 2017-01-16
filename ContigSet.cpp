#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "ContigSet.h"

using namespace std;


ContigSetHead * GetContigSet(char * contigSetFile){
    
    ContigSetHead * contigSetHead = (ContigSetHead *)malloc(sizeof(ContigSetHead));
    contigSetHead->contigSet = NULL;
    contigSetHead->contigCount = 0;
    
    long int maxSize = 10000;
    char * contig = NULL;
    if(NULL == (contig = (char*)malloc(sizeof(char)*maxSize))){
        perror("malloc error!");
        exit(1);
    }
    
    FILE * fp; 
    if((fp = fopen(contigSetFile, "r")) == NULL){
        printf("%s, does not exist!", contigSetFile);
        exit(0);
    }
    while((fgets(contig, maxSize, fp)) != NULL){ 
       if(contig[0] == '>'){  
           contigSetHead->contigCount++; 
       }  
    }  
    fclose(fp);
    
    contigSetHead->contigSet = (Contig *)malloc(sizeof(Contig)*contigSetHead->contigCount);
    for(long int i = 0; i < contigSetHead->contigCount; i++){
        contigSetHead->contigSet[i].contigName = NULL;
        contigSetHead->contigSet[i].contig = NULL;
        contigSetHead->contigSet[i].contigLength = 0;
    }
    
    if((fp = fopen(contigSetFile, "r")) == NULL){
        printf("%s, does not exist!", contigSetFile);
        exit(0);
    }
    
    long int allocateLength = 0;
    long int contigIndex = -1;
    while((fgets(contig, maxSize, fp)) != NULL){ 
       
       if(contig[0] == '>'){  
           if(strlen(contig) == maxSize-1){              
               while((fgets(contig, maxSize, fp)) != NULL){
                   if(strlen(contig) != maxSize-1){
                       break;
                   }
               }        
           }
           contigIndex++;
           
           contigSetHead->contigSet[contigIndex].contigName = (char*)malloc(sizeof(char)*(strlen(contig)-1));
           strncpy(contigSetHead->contigSet[contigIndex].contigName, contig+1, strlen(contig)-2); 
           contigSetHead->contigSet[contigIndex].contigName[strlen(contig)-2] = '\0';
           continue;
           
       }
       
       
       long int extendLength = strlen(contig);
       if(contig[extendLength-1] == '\n'){
           extendLength--;
       }
       long int contigLength = 0;
       char * tempContig = NULL;
       if(contigSetHead->contigSet[contigIndex].contig != NULL){
           if(contigSetHead->contigSet[contigIndex].contigLength + extendLength >= allocateLength){
               contigLength = contigSetHead->contigSet[contigIndex].contigLength;    
               contigSetHead->contigSet[contigIndex].contig = (char *)realloc(contigSetHead->contigSet[contigIndex].contig, allocateLength + maxSize + 1);

               allocateLength = allocateLength + maxSize + 1;
               
               strncpy(contigSetHead->contigSet[contigIndex].contig + contigLength, contig, extendLength);
               contigSetHead->contigSet[contigIndex].contig[contigLength + extendLength] = '\0';    
               contigSetHead->contigSet[contigIndex].contigLength = contigLength + extendLength;
                       
           }else{
               strncpy(contigSetHead->contigSet[contigIndex].contig + contigSetHead->contigSet[contigIndex].contigLength, contig, extendLength);
               contigSetHead->contigSet[contigIndex].contig[contigSetHead->contigSet[contigIndex].contigLength + extendLength] = '\0';
               contigSetHead->contigSet[contigIndex].contigLength = contigSetHead->contigSet[contigIndex].contigLength + extendLength;
           }   
           
       }else{
           contigSetHead->contigSet[contigIndex].contig = (char *)malloc(sizeof(char)*(maxSize+1));
           strncpy(contigSetHead->contigSet[contigIndex].contig, contig, extendLength);
           contigSetHead->contigSet[contigIndex].contig[extendLength] = '\0';
           contigSetHead->contigSet[contigIndex].contigLength = extendLength;
           allocateLength = maxSize + 1;
       }  
    }  
    fflush(fp);
    fclose(fp);
    
    return contigSetHead;
}



