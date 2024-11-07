#////////////////////////////////////////////////////////////////////////////
#//
#// Author      : Luis M. de la Cruz Salas
#// Date        : Wed Dec 19 11:22:11 GMT 2007
#// Version     : 0.1.0
#//
#////////////////////////////////////////////////////////////////////////////

include	../rules.in

all : vtkformat

vtkformat: vtkformat.o 
	$(CXX) $< -o $@ $(LIBS)

vtkformat.o: vtkformat.cpp
	$(CXX) $(CXX_PREFLAGS) -D_3D_ $(INCLUDES) -c $<

clean:
	/bin/rm -f vtkformat vtkformat.o 
