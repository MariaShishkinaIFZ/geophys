stacklays_Y : stacklays_Y.o submem.o subio.o iogras.o
	gcc -o stacklays_Y -lm stacklays_Y.o submem.o subio.o iogras.o
	
stacklays_Y.o : stacklays_Y.c submem.h subio.h  iogras.h
	gcc -c stacklays_Y.c -Wall

submem.o : submem.c submem.h
	gcc -c -g -Wall submem.c
subio.o : subio.c subio.h
	gcc -c -g -Wall subio.c
iogras.o : iogras.c iogras.h
	gcc -c -g -Wall iogras.c

clean :
	rm *.o stacklays_Y

install :
	cp stacklays_Y /common/tomo/bin/
	
