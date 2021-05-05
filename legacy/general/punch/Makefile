punch : punch.o subio.o
	cc -o punch -lm punch.o subio.o
punch.o : punch.c
	cc -c punch.c
subio.o : subio.c subio.h
	cc -c -Wall subio.c

clean :
	rm punch punch.o subio.o
	
install :
	cp punch /usr/local/tomo/bin/
