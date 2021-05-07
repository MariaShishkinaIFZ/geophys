all	: 2dXYgmt 2dXZgmt 2dYZgmt

2dXYgmt : 2dXYgmt.o subio.o submem.o
	gcc -o 2dXYgmt -lm 2dXYgmt.o subio.o submem.o
	
2dXYgmt.o : 2dXYgmt.c subio.h submem.h
	gcc -c 2dXYgmt.c

2dXZgmt : 2dXZgmt.o subio.o submem.o
	gcc -o 2dXZgmt -lm 2dXZgmt.o subio.o submem.o
	
2dXZgmt.o : 2dXZgmt.c subio.h submem.h
	gcc -c 2dXZgmt.c


2dYZgmt : 2dYZgmt.o subio.o submem.o
	gcc -o 2dYZgmt -lm 2dYZgmt.o subio.o submem.o
	
2dYZgmt.o : 2dYZgmt.c subio.h submem.h
	gcc -c 2dYZgmt.c

submem.o : submem.c submem.h
	gcc -c -g -Wall submem.c
subio.o : subio.c subio.h
	gcc -c -g -Wall subio.c

clean :
	rm *.o 2dXYgmt 2dXZgmt 2dYZgmt
	
install :
	cp 2dXYgmt /usr/local/tomo/bin/
	cp 2dXZgmt /usr/local/tomo/bin/    
	cp 2dYZgmt /usr/local/tomo/bin/
	