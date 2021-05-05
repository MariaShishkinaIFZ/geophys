trhwgmt : trhwgmt.o  subio.o 
	gcc -g -o trhwgmt -lm trhwgmt.o  subio.o 
	
trhwgmt.o : trhwgmt.c subio.h 
	gcc -c -g -Wall trhwgmt.c
	
subio.o : subio.c subio.h
	gcc -c -g -Wall subio.c
	
clean :
	rm trhwgmt *.o

install :
	cp trhwgmt /usr/local/tomo/bin/
