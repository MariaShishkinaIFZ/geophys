mkbnd_f : mkbnd_f.o subio.o submem.o iogras.o
	cc -o mkbnd_f -lm mkbnd_f.o subio.o submem.o iogras.o
	
mkbnd_f.o : mkbnd_f.c subio.h submem.h iogras.h
	cc -c mkbnd_f.c
	
mkbnd_f1 : mkbnd_f1.o subio.o submem.o iogras.o
	cc -o mkbnd_f1 -lm mkbnd_f1.o subio.o submem.o iogras.o
	
mkbnd_f1.o : mkbnd_f1.c subio.h submem.h iogras.h
	cc -c mkbnd_f1.c

submem.o : submem.c submem.h
	gcc -c -g -Wall submem.c
subio.o : subio.c subio.h
	gcc -c -g -Wall subio.c
iogras.o : iogras.c iogras.h
	gcc -c -g -Wall iogras.c

clean :
	rm *.o mkbnd_f mkbnd_f1

install :
	cp mkbnd_f /usr/local/tomo/bin/
	cp mkbnd_f1 /usr/local/tomo/bin/
	
