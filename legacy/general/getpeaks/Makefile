all :  getpeaks-1 

getpeaks-1 : getpeaks-1.o subio.o submem.o
	gcc -o getpeaks-1 -lm getpeaks-1.o subio.o submem.o

getpeaks-1.o : getpeaks-1.c subio.h submem.h
	gcc -c getpeaks-1.c

clean :
	rm *.o  getpeaks-1

install :
	cp getpeaks-1 /usr/local/tomo/bin