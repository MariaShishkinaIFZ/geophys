fine2coarse : fine2coarse.o submem.o subio.o nrutil.o
	g++ -g -o fine2coarse -lm fine2coarse.o submem.o subio.o nrutil.o

fine2coarse.o : fine2coarse.cpp subio.h submem.h nrutil.h gmtio.h
	g++ -c -g -Wall fine2coarse.cpp
submem.o : submem.cpp submem.h
	g++ -c -g -Wall submem.cpp
subio.o : subio.cpp subio.h
	g++ -c -g -Wall subio.cpp
nrutil.o : nrutil.cpp nrutil.h
	g++ -c -g -Wall nrutil.cpp

clean :
	rm fine2coarse *.o
