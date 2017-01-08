CC=g++

CPPFLAGS = -g -Wall -O3

getContigInScaffold:	getContigInScaffold.o ScaffoldSet.o
	$(CC) -o $@ $^
	
getGapRegion:	getGapRegion.o ContigSet.o ScaffoldSet.o Alignment.o Gap.o
	$(CC) -o $@ $^
	
getContigInScaffold.o: getContigInScaffold.cpp ScaffoldSet.h
	$(CC) -c getContigInScaffold.cpp
	
getGapRegion.o:	getGapRegion.cpp ContigSet.h ScaffoldSet.h Alignment.h Gap.h
	$(CC) -c getGapRegion.cpp
	
ContigSet.o:	ContigSet.cpp ContigSet.h
	$(CC) -c ContigSet.cpp
	
ScaffoldSet.o:	ScaffoldSet.cpp ScaffoldSet.h
	$(CC) -c ScaffoldSet.cpp
	
Alignment.o:	Alignment.cpp Alignment.h
	$(CC) -c Alignment.cpp
	
Gap.o:	Gap.cpp Gap.h ContigSet.h ScaffoldSet.h Alignment.h
	$(CC) -c Gap.cpp
	
all: getContigInScaffold getGapRegion
	make all -C ./needleman_wunsch_code
	cp ./needleman_wunsch_code/needleman_wunsch ./
	rm -f *.o ./needleman_wunsch_code/needleman_wunsch
clean:
	rm -f *.o
	rm getContigInScaffold getGapRegion needleman_wunsch 

