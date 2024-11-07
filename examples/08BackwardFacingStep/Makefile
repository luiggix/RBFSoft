#////////////////////////////////////////////////////////////////////////////
#//
#// Author      : Luis M. de la Cruz Salas
#// Date        : Mon Apr  7 13:33:35 BST 2008
#// Version     : 0.2.0
#//
#////////////////////////////////////////////////////////////////////////////

include	../rules.in

all : backfacing1 backfacing2

backfacing1: backfacing1.o 
	$(CXX) $< -o $@ $(LIBS)

backfacing1.o: backfacing1.cpp
	$(CXX) $(CXX_PREFLAGS) -D_2D_ $(INCLUDES) -c $<

backfacing2: backfacing2.o 
	$(CXX) $< -o $@ $(LIBS)

backfacing2.o: backfacing2.cpp
	$(CXX) $(CXX_PREFLAGS) -D_2D_ $(INCLUDES) -c $<

clean:
	/bin/rm -f backfacing1 backfacing1.o backfacing2 backfacing2.o
