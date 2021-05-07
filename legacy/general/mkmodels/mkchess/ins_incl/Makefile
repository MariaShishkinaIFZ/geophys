ins_incl1 : ins_incl1.o submem.o subio.o iogras.o
	cc -o ins_incl1 -lm ins_incl1.o submem.o subio.o iogras.o

ins_incl1.o : ins_incl1.c subio.h submem.h iogras.h
	cc -c -Wall ins_incl1.c
submem.o : submem.c submem.h
	cc -c -Wall submem.c
subio.o : subio.c subio.h
	cc -c -Wall subio.c
iogras.o : iogras.c iogras.h
	cc -c -Wall iogras.c

clean :
	rm ins_incl1 *.o
	
install:
	cp ins_incl1 $(ISDG_HOME)

