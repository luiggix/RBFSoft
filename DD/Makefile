#////////////////////////////////////////////////////////////////////////////
#//
#// Author      : Luis M. de la Cruz Salas
#// Date        : Tue Feb 19 12:38:16 GMT 2008
#// Version     : 0.1.0
#//
#////////////////////////////////////////////////////////////////////////////

include	../rulesMPI.in

all : Heat2DParallel

Heat2DParallel: Heat2DParallel.o 
	$(CXX) $< -o $@ $(LIBS)

Heat2DParallel.o: Heat2DParallel.cpp
	$(CXX) $(CXXFLAGS) -D_2D_ $(INCLUDES) -c $<

clean:
	/bin/rm -f Heat2DParallel Heat2DParallel.o
