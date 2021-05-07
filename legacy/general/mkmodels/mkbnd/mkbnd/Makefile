mkbnd : mkbnd.o subio.o submem.o iogras.o
	cc -o mkbnd -lm mkbnd.o subio.o submem.o iogras.o
	
mkbnd.o : mkbnd.c subio.h submem.h iogras.h
	cc -c mkbnd.c -lm

submem.o : submem.c submem.h
	gcc -c -g -Wall submem.c -lm
subio.o : subio.c subio.h
	gcc -c -g -Wall subio.c -lm
iogras.o : iogras.c iogras.h
	gcc -c -g -Wall iogras.c -lm 

clean :
	rm *.o mkbnd

install :
	cp mkbnd $(ISDG_HOME)
	
