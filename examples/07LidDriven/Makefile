#////////////////////////////////////////////////////////////////////////////
#//
#// Author      : Luis M. de la Cruz Salas
#// Date        : Thu Feb 14 11:16:19 GMT 2008
#// Version     : 0.1.0
#//
#////////////////////////////////////////////////////////////////////////////

include	../rules.in

all : lid_driven1 lid_driven2

lid_driven1: lid_driven1.o 
	$(CXX) $< -o $@ $(LIBS)

lid_driven1.o: lid_driven1.cpp
	$(CXX) $(CXX_PREFLAGS) -D_2D_ $(INCLUDES) -c $<

lid_driven2: lid_driven2.o 
	$(CXX) $< -o $@ $(LIBS)

lid_driven2.o: lid_driven2.cpp
	$(CXX) $(CXX_PREFLAGS) -D_2D_ $(INCLUDES) -c $<

clean:
	/bin/rm -f lid_driven1 lid_driven1.o lid_driven2 lid_driven2.o 