#////////////////////////////////////////////////////////////////////////////
#//
#// Author      : Luis M. de la Cruz Salas
#// Date        : Wed Dec 19 11:22:11 GMT 2007
#// Version     : 0.1.0
#//
#////////////////////////////////////////////////////////////////////////////

include	../rules.in

all : heatsemi

heatsemi: heatsemi.o 
	$(CXX) $< -o $@ $(LIBS)

heatsemi.o: heatsemi.cpp
	$(CXX) $(CXXFLAGS) -D_2D_ $(INCLUDES) -c $<

clean:
	/bin/rm -f heatsemi heatsemi.o 
