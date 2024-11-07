#////////////////////////////////////////////////////////////////////////////
#//
#// Author      : Luis M. de la Cruz Salas
#// Date        : Wed Dec 19 11:22:11 GMT 2007
#// Version     : 0.1.0
#//
#////////////////////////////////////////////////////////////////////////////

include	../rules.in

all : heat2D heat2D_blocks

heat2D: heat2D.o 
	$(CXX) $< -o $@ $(LIBS)

heat2D.o: heat2D.cpp
	$(CXX) $(CXXFLAGS) -D_2D_ $(INCLUDES) -c $<

heat2D_blocks: heat2D_blocks.o 
	$(CXX) $< -o $@ $(LIBS)

heat2D_blocks.o: heat2D_blocks.cpp
	$(CXX) $(CXXFLAGS) -D_2D_ $(INCLUDES) -c $<

clean:
	/bin/rm -f heat2D heat2D.o heat2D_blocks heat2D_blocks.o 
