#////////////////////////////////////////////////////////////////////////////
#//
#// Author      : Luis M. de la Cruz Salas
#// Date        : Wed Mar 19 10:22:04 GMT 2008
#// Version     : 0.1.0
#//
#////////////////////////////////////////////////////////////////////////////

include	../../rulesMPI.in

all : lid_driven_par 

lid_driven_par: lid_driven_par.o 
	$(CXX) $< -o $@ $(LIBS)

lid_driven_par.o: lid_driven_par.cpp
	$(CXX) $(CXX_PREFLAGS) -D_2D_ $(INCLUDES) -c $<

clean:
	/bin/rm -f lid_driven_par lid_driven_par.o 
