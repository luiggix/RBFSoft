#////////////////////////////////////////////////////////////////////////////
#//
#// Author      : Luis M. de la Cruz Salas
#// Date        : Wed Dec 19 11:22:11 GMT 2007
#// Version     : 0.1.0
#//
#////////////////////////////////////////////////////////////////////////////

include	../../rules.in

all : heat2DTest

heat2DTest: heat2DTest.o 
	$(CXX) $< -o $@ $(LIBS)

heat2DTest.o: heat2DTest.cpp
	$(CXX) $(CXXFLAGS) -D_2D_ $(INCLUDES) -c $<

clean:
	/bin/rm -f heat2DTest heat2DTest.o
