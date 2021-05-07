all 	: 3d2text.o submem.o subio.o iogras.o
	cc -o 3d2text -lm 3d2text.o submem.o subio.o iogras.o


3d2text.o : 3d2text.c subio.h submem.h iogras.h
	cc -c -Wall 3d2text.c

submem.o : submem.c submem.h
	cc -c -Wall submem.c
subio.o : subio.c subio.h
	cc -c -Wall subio.c
iogras.o : iogras.c iogras.h
	cc -c -Wall iogras.c

clean :
	rm 3d2text *.o
	
install: 
	cp 3d2text $(ISDG_HOME)
