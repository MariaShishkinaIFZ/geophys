mkcrvel_prof_Y : mkcrvel_prof_Y.o subio.o submem.o
	cc -o mkcrvel_prof_Y -lm mkcrvel_prof_Y.o subio.o submem.o
	
mkcrvel_prof_Y.o : mkcrvel_prof_Y.c subio.h submem.h 
	cc -c mkcrvel_prof_Y.c

submem.o : submem.c submem.h
	gcc -c -g -Wall submem.c
subio.o : subio.c subio.h
	gcc -c -g -Wall subio.c

clean :
	rm *.o mkcrvel_prof_Y

install :
	cp mkcrvel_prof_Y /usr/local/tomo/bin/
	
