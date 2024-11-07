#////////////////////////////////////////////////////////////////////////////
#//
#// Author      : Luis M. de la Cruz Salas
#// Date        : Mon Jan 31 18:04:40 CST 2005
#// Version     : 1.0
#// Last Rev.   : Mon Jun  4 13:39:54 CDT 200
#//
#////////////////////////////////////////////////////////////////////////////

include	../rules.in

all : convdiff01

convdiff01: convdiff01.o 
	$(CXX) $< -o $@ $(LIBS)

convdiff01.o: convdiff01.cpp
	$(CXX) $(CXX_PREFLAGS) -O2 -D_1D_ $(INCLUDES) -c $<

clean:
	/bin/rm -f convdiff01 convdiff01.o 
