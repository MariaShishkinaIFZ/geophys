mkcrvel : mkcrvel.o subio.o submem.o
	cc -o mkcrvel -lm mkcrvel.o subio.o submem.o
	
mkcrvel.o : mkcrvel.c subio.h submem.h 
	cc -c mkcrvel.c

submem.o : submem.c submem.h
	gcc -c -g -Wall submem.c
subio.o : subio.c subio.h
	gcc -c -g -Wall subio.c

clean :
	rm *.o mkcrvel

install :
	cp mkcrvel /usr/local/tomo/bin/
	
