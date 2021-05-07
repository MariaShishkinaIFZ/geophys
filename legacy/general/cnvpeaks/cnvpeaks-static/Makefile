cnvpeaks-static : cnvpeaks-static.o subio.o submem.o
	cc -o cnvpeaks-static -lm cnvpeaks-static.o subio.o submem.o
	
cnvpeaks-static.o : cnvpeaks-static.c subio.h submem.h
	cc -c cnvpeaks-static.c

clean :
	rm *.o cnvpeaks-static

install :
	cp cnvpeaks-static /usr/local/tomo/bin/
	
