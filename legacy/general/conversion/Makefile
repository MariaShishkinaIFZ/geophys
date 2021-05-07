2dcgmt : 2dcgmt.o subio.o submem.o
	gcc -o 2dcgmt -lm 2dcgmt.o subio.o submem.o
	
2dcgmt.o : 2dcgmt.c subio.h submem.h
	gcc -c 2dcgmt.c

submem.o : submem.c submem.h
	gcc -c -g -Wall submem.c
subio.o : subio.c subio.h
	gcc -c -g -Wall subio.c

clean :
	rm *.o 2dcgmt
	
install :
	cp 2dcgmt /usr/local/tomo/bin/
    