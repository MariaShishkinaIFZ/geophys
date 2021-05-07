3D_mask : mask.o submem.o subio.o nrutil.o
	g++ -g -o 3D_mask -lm mask.o submem.o subio.o nrutil.o

mask.o : mask.cpp subio.h submem.h nrutil.h gmtio.h
	g++ -c -g -Wall mask.cpp
submem.o : submem.cpp submem.h
	g++ -c -g -Wall submem.cpp
subio.o : subio.cpp subio.h
	g++ -c -g -Wall subio.cpp
nrutil.o : nrutil.cpp nrutil.h
	g++ -c -g -Wall nrutil.cpp

clean :
	rm 3D_mask *.o

install :
	cp 3D_mask /usr/local/tomo/bin/
