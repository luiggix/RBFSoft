#////////////////////////////////////////////////////////////////////////////
#//
#// Author      : Luis M. de la Cruz Salas
#// Date        : Mon Jan 31 18:04:40 CST 2005
#// Version     : 1.0
#// Last Rev.   : Mon Jun  4 13:39:54 CDT 200
#//
#////////////////////////////////////////////////////////////////////////////

include	../rules.in
# Path where Boost libraries are installed
BOOST = /home/lcz/Research/MathUtils/boost_1_34_1
INCLUDES = -I$(BOOST)/ 

all : mattest

mattest: mattest.o 
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $@ $(LIBS)

mattest.o: mattest.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<

clean:
	/bin/rm -f mattest mattest.o 
