specfem2sgy : specfem2sgy.o submem.o subio.o segyio.o
	g++ -o specfem2sgy -lm specfem2sgy.o submem.o subio.o segyio.o nrutil.o

specfem2sgy.o : specfem2sgy.cpp subio.h submem.h segyio.h
	g++ -c -Wall specfem2sgy.cpp
submem.o : submem.cpp submem.h
	g++ -c -Wall submem.cpp
subio.o : subio.cpp subio.h
	g++ -c -Wall subio.cpp
segyio.o : segyio.cpp segyio.h
	g++ -c -Wall segyio.cpp
nrutil.o : nrutil.cpp nrutil.h
	g++ -c -Wall nrutil.cpp

clean :
	rm specfem2sgy *.o
	
install:
	cp specfem2sgy $(ISDG_HOME)

