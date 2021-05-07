mkcrvel_prof : mkcrvel_prof.o subio.o submem.o
	cc -o mkcrvel_prof -lm mkcrvel_prof.o subio.o submem.o
	
mkcrvel_prof.o : mkcrvel_prof.c subio.h submem.h 
	cc -c mkcrvel_prof.c

submem.o : submem.c submem.h
	gcc -c -g -Wall submem.c
subio.o : subio.c subio.h
	gcc -c -g -Wall subio.c

clean :
	rm *.o mkcrvel_prof

install :
	cp mkcrvel_prof /usr/local/tomo/bin/
	
