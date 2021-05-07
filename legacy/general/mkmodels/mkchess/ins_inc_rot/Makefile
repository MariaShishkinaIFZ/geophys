ins_incl_rot : ins_incl_rot.o submem.o subio.o iogras.o
	cc -o ins_incl_rot -lm ins_incl_rot.o submem.o subio.o iogras.o

ins_incl_rot.o : ins_incl_rot.c subio.h submem.h iogras.h
	cc -c -Wall ins_incl_rot.c
submem.o : submem.c submem.h
	cc -c -Wall submem.c
subio.o : subio.c subio.h
	cc -c -Wall subio.c
iogras.o : iogras.c iogras.h
	cc -c -Wall iogras.c

clean :
	rm ins_incl_rot *.o
	
install:
	cp ins_incl_rot $(ISDG_HOME)

