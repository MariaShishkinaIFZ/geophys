cattraces : cattraces.o submem.o subio.o
	cc -o cattraces -lm cattraces.o submem.o subio.o

cattraces.o : cattraces.c subio.h submem.h
	cc -c -Wall cattraces.c
submem.o : submem.c submem.h
	cc -c -Wall submem.c
subio.o : subio.c subio.h

clean :
	rm cattraces *.o
