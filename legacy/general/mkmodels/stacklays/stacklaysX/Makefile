stacklays_X : stacklays_X.o submem.o subio.o iogras.o
	gcc -o stacklays_X -lm stacklays_X.o submem.o subio.o iogras.o
	
stacklays_X.o : stacklays_X.c submem.h subio.h  iogras.h
	gcc -c stacklays_X.c -Wall

submem.o : submem.c submem.h
	gcc -c -g -Wall submem.c
subio.o : subio.c subio.h
	gcc -c -g -Wall subio.c
iogras.o : iogras.c iogras.h
	gcc -c -g -Wall iogras.c

clean :
	rm *.o stacklays_X

install :
	cp stacklays_X /common/tomo/bin/
	
