#////////////////////////////////////////////////////////////////////////////
#//
#// Author      : Luis M. de la Cruz Salas
#// Date        : Sun May  3 18:51:48 CDT 2015
#// Version     : 0.1.0
#//
#////////////////////////////////////////////////////////////////////////////

include	../rules.in

all : invSurf 

invSurf: invSurf.o 
	$(CXX) $< -o $@ $(LIBS)

invSurf.o: invSurf.cpp
	$(CXX) $(CXXFLAGS) -D_2D_ $(INCLUDES) -c $<

clean:
	/bin/rm -f invSurf invSurf.o *.~
