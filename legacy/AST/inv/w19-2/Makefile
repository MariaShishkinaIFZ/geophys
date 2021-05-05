CILK = g++
INCADD = 
LFLAGS = -lm -leinspline

GCCOPT = -g -O2 -fno-rtti -fno-exceptions 

s-invert19-2-par : 3d-seis-w19-2.o haar.o submem.o subio.o nrutil.o lsqr.o
	$(CILK) $(INCADD) $(GCCOPT) $(LFLAGS) -o s-invert19-2-par -lm 3d-seis-w19-2.o haar.o submem.o subio.o nrutil.o lsqr.o

3d-seis-w19-2.o : 3d-seis-w19-2.cpp subio.h submem.h nrutil.h lsqr.h haar.h gmtio.h
	$(CILK) $(INCADD) $(GCCOPT) -c -Wall 3d-seis-w19-2.cpp
haar.o : haar.cpp haar.h
	$(CILK) $(INCADD) $(GCCOPT) -c -Wall haar.cpp
submem.o : submem.cpp submem.h
	$(CILK) $(INCADD) $(GCCOPT) -c -Wall submem.cpp
subio.o : subio.cpp subio.h
	$(CILK) $(INCADD) $(GCCOPT) -c -Wall subio.cpp
nrutil.o : nrutil.cpp nrutil.h
	$(CILK) $(INCADD) $(GCCOPT) -c -Wall nrutil.cpp
lsqr.o : lsqr.cpp lsqr.h
	$(CILK) $(INCADD) $(GCCOPT) -c -Wall lsqr.cpp

clean :
	rm s-invert19-2-par *.o
	
install :
	#cp s-invert19-2-par /usr/local/tomo/bin
	cp s-invert19-2-par /common/tomo/bin
	
