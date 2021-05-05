CILK = /usr/local/cilk
INCADD = -I$(CILK)/include -I$(CILK)/examples/include

GCCOPT = -g -O2 -fno-rtti -fno-exceptions 

s-invert18-par : 3d-seis-w19.o haar.o submem.o subio.o nrutil.o lsqr.o
	$(CILK)/bin/cilk++ -DCILKPARALLEL $(INCADD) $(GCCOPT) -o s-invert19-par -lm 3d-seis-w19.o haar.o submem.o subio.o nrutil.o lsqr.o

3d-seis-w19.o : 3d-seis-w19.cpp subio.h submem.h nrutil.h lsqr.h haar.h gmtio.h
	$(CILK)/bin/cilk++ -DCILKPARALLEL $(INCADD) $(GCCOPT) -c -Wall 3d-seis-w19.cpp
haar.o : haar.cpp haar.h
	$(CILK)/bin/cilk++ -DCILKPARALLEL $(INCADD) $(GCCOPT) -c -Wall haar.cpp
submem.o : submem.cpp submem.h
	$(CILK)/bin/cilk++ -DCILKPARALLEL $(INCADD) $(GCCOPT) -c -Wall submem.cpp
subio.o : subio.cpp subio.h
	$(CILK)/bin/cilk++ -DCILKPARALLEL $(INCADD) $(GCCOPT) -c -Wall subio.cpp
nrutil.o : nrutil.cpp nrutil.h
	$(CILK)/bin/cilk++ -DCILKPARALLEL $(INCADD) $(GCCOPT) -c -Wall nrutil.cpp
lsqr.o : lsqr.cpp lsqr.h
	$(CILK)/bin/cilk++ -DCILKPARALLEL $(INCADD) $(GCCOPT) -c -Wall lsqr.cpp

clean :
	rm s-invert19-par *.o
	
install :
	cp s-invert19-par /usr/local/tomo/bin
