#////////////////////////////////////////////////////////////////////////////
#//
#// Author      : Luis M. de la Cruz Salas
#// Date        : Mon Apr  7 13:33:35 BST 2008
#// Version     : 0.2.0
#//
#////////////////////////////////////////////////////////////////////////////

include	../rules.in

all : poiseuille

poiseuille: poiseuille.o 
	$(CXX) $< -o $@ $(LIBS)

poiseuille.o: poiseuille.cpp
	$(CXX) $(CXX_PREFLAGS) -D_2D_ $(INCLUDES) -c $<

clean:
	/bin/rm -f poiseuille poiseuille.o 
