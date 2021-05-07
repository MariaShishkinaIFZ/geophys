rand_tpk : rand_tpk.o submem.o subio.o gasdev.o ran1.o
	cc -o rand_tpk -lm rand_tpk.o submem.o subio.o gasdev.o ran1.o

rand_tpk.o : rand_tpk.c subio.h submem.h
	cc -c -Wall rand_tpk.c
	
submem.o : submem.c submem.h
	cc -c -Wall submem.c
	
subio.o : subio.c subio.h
	cc -c -Wall subio.c
	
gasdev.o : gasdev.c
	cc -c -Wall gasdev.c

ran1.o : ran1.c
	cc -c -Wall ran1.c

clean :
	rm rand_tpk *.o
