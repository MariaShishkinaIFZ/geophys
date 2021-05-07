mkcrvel_prof_X : mkcrvel_prof_X.o subio.o submem.o
	cc -o mkcrvel_prof_X -lm mkcrvel_prof_X.o subio.o submem.o
	
mkcrvel_prof_X.o : mkcrvel_prof_X.c subio.h submem.h 
	cc -c mkcrvel_prof_X.c

submem.o : submem.c submem.h
	gcc -c -g -Wall submem.c
subio.o : subio.c subio.h
	gcc -c -g -Wall subio.c

clean :
	rm *.o mkcrvel_prof_X

install :
	cp mkcrvel_prof_X /common/tomo/bin/
	
