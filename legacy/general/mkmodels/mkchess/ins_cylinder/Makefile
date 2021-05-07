ins_cylinder_new : ins_cylinder_new.o submem.o subio.o iogras.o
	cc -o ins_cylinder_new -lm ins_cylinder_new.o submem.o subio.o iogras.o

ins_cylinder_new.o : ins_cylinder_new.c subio.h submem.h iogras.h
	cc -c -Wall ins_cylinder_new.c
submem.o : submem.c submem.h
	cc -c -Wall submem.c
subio.o : subio.c subio.h
	cc -c -Wall subio.c
iogras.o : iogras.c iogras.h
	cc -c -Wall iogras.c

clean :
	rm ins_cylinder_new *.o
	
install:
	cp ins_incl_rot $(ISDG_HOME)

