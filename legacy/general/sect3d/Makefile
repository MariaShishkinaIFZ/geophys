sect3d : sect3d.o subio.o submem.o
	cc -o sect3d -lm sect3d.o subio.o submem.o
	
sect3d.o : sect3d.c subio.h submem.h gmt_nan.h
	cc -c sect3d.c

submem.o : submem.c submem.h
	gcc -c -g -Wall submem.c
	
subio.o : subio.c subio.h
	gcc -c -g -Wall subio.c
	


clean :
	rm *.o sect3d

install :
	cp sect3d /usr/local/tomo/bin/
	
