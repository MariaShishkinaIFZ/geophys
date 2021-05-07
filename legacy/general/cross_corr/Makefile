cross_corr : cross_corr.o submem.o subio.o
	gcc -o cross_corr -lm cross_corr.o submem.o subio.o

cross_corr.o : cross_corr.c subio.h submem.h
	gcc -c -Wall cross_corr.c
submem.o : submem.c submem.h
	gcc -c -Wall submem.c
subio.o : subio.c subio.h
	gcc -c -Wall subio.c

clean :
	rm cross_corr *.o
