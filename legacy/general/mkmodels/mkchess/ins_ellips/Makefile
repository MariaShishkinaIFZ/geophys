ins_ellips : ins_ellips.o submem.o subio.o iogras.o
	cc -o ins_ellips -lm ins_ellips.o submem.o subio.o iogras.o

ins_ellips.o : ins_ellips.c subio.h submem.h iogras.h
	cc -c -Wall ins_ellips.c
submem.o : submem.c submem.h
	cc -c -Wall submem.c
subio.o : subio.c subio.h
	cc -c -Wall subio.c
iogras.o : iogras.c iogras.h
	cc -c -Wall iogras.c

clean :
	rm ins_ellips *.o
	
install:
	cp ins_ellips /usr/local/tomo/bin

