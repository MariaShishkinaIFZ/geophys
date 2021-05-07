stacklays_f : stacklays_f.o submem.o subio.o iogras.o
	gcc -o stacklays_f -lm stacklays_f.o submem.o subio.o iogras.o
	
stacklays_f.o : stacklays_f.c submem.h subio.h  iogras.h
	gcc -c stacklays_f.c -Wall

submem.o : submem.c submem.h
	gcc -c -g -Wall submem.c
subio.o : subio.c subio.h
	gcc -c -g -Wall subio.c
iogras.o : iogras.c iogras.h
	gcc -c -g -Wall iogras.c

clean :
	rm *.o stacklays_f

install :
	cp stacklays_f /usr/local/tomo/bin/
	
