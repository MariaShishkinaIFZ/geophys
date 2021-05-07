cnvpeaks-1 : cnvpeaks-1.o subio.o submem.o
	cc -o cnvpeaks-1 -lm cnvpeaks-1.o subio.o submem.o
	
cnvpeaks-1.o : cnvpeaks-1.c subio.h submem.h
	cc -c cnvpeaks-1.c

clean :
	rm *.o cnvpeaks-1
install :
	cp cnvpeaks-1 /usr/local/tomo/bin/
	
