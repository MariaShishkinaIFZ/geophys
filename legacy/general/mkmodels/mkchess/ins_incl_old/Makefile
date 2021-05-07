ins_incl : ins_incl.o submem.o subio.o iogras.o
	cc -o ins_incl -lm ins_incl.o submem.o subio.o iogras.o

ins_incl.o : ins_incl.c subio.h submem.h iogras.h
	cc -c -Wall ins_incl.c
submem.o : submem.c submem.h
	cc -c -Wall submem.c
subio.o : subio.c subio.h
	cc -c -Wall subio.c
iogras.o : iogras.c iogras.h
	cc -c -Wall iogras.c

clean :
	rm ins_incl *.o
	
install:
	cp ins_incl /usr/local/tomo/bin

