#////////////////////////////////////////////////////////////////////////////
#//
#// Author      : Luis M. de la Cruz Salas
#// Date        : Mon Dec 17 11:18:29 GMT 2007
#// Version     : 0.1.0
#//
#////////////////////////////////////////////////////////////////////////////

include	../rules.in

all : rectangle semicircle line rectangle_blocks box

rectangle_blocks: rectangle_blocks.o 
	$(CXX) $< -o $@ $(LIBS)

rectangle_blocks.o: rectangle_blocks.cpp
	$(CXX) $(CXXFLAGS) -D_2D_ $(INCLUDES) -c $<

rectangle: rectangle.o 
	$(CXX) $< -o $@ $(LIBS)

rectangle.o: rectangle.cpp
	$(CXX) $(CXXFLAGS) -D_2D_ $(INCLUDES) -c $<

semicircle: semicircle.o 
	$(CXX) $< -o $@ $(LIBS)

semicircle.o: semicircle.cpp
	$(CXX) $(CXXFLAGS) -D_2D_ $(INCLUDES) -c $<


line: line.o 
	$(CXX) $< -o $@ $(LIBS)

line.o: line.cpp
	$(CXX) $(CXXFLAGS) -D_1D_ $(INCLUDES) -c $<

box: box.o 
	$(CXX) $< -o $@ $(LIBS)

box.o: box.cpp
	$(CXX) $(CXXFLAGS) -D_3D_ $(INCLUDES) -c $<

clean:
	/bin/rm -f rectangle rectangle.o semicircle semicircle.o line line.o rectangle_blocks rectangle_blocks.o box box.o 
