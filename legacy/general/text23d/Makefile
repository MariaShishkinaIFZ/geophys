all 	: text23d.o submem.o subio.o iogras.o
	cc -o text23d -lm text23d.o submem.o subio.o iogras.o


text23d.o : text23d.c subio.h submem.h iogras.h
	cc -c -Wall text23d.c

submem.o : submem.c submem.h
	cc -c -Wall submem.c
subio.o : subio.c subio.h
	cc -c -Wall subio.c
iogras.o : iogras.c iogras.h
	cc -c -Wall iogras.c

clean :
	rm text23d *.o
	
install: 
	cp text23d $(ISDG_HOME)/
