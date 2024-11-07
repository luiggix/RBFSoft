#////////////////////////////////////////////////////////////////////////////
#//
#// Author      : Luis M. de la Cruz Salas
#// Date        : Wed Dec 19 11:22:11 GMT 2007
#// Version     : 0.1.0
#//
#////////////////////////////////////////////////////////////////////////////

include	../rules.in

all : convdiff2D convdiff2D_blocks

convdiff2D: convdiff2D.o 
	$(CXX) $< -o $@ $(LIBS)

convdiff2D.o: convdiff2D.cpp
	$(CXX) $(CXXFLAGS) -D_2D_ $(INCLUDES) -c $<

convdiff2D_blocks: convdiff2D_blocks.o 
	$(CXX) $< -o $@ $(LIBS)

convdiff2D_blocks.o: convdiff2D_blocks.cpp
	$(CXX) $(CXXFLAGS) -D_2D_ $(INCLUDES) -c $<

clean:
	/bin/rm -f convdiff2D convdiff2D.o convdiff2D_blocks convdiff2D_blocks.o 