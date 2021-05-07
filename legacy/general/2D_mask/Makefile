2D_mask : 2D_mask.o submem.o subio.o nrutil.o
	g++ -g -o 2D_mask -lm 2D_mask.o submem.o subio.o nrutil.o

2D_mask.o : 2D_mask.cpp subio.h submem.h nrutil.h gmtio.h
	g++ -c -g -Wall 2D_mask.cpp
submem.o : submem.cpp submem.h
	g++ -c -g -Wall submem.cpp
subio.o : subio.cpp subio.h
	g++ -c -g -Wall subio.cpp
nrutil.o : nrutil.cpp nrutil.h
	g++ -c -g -Wall nrutil.cpp

clean :
	rm 3D_mask *.o
