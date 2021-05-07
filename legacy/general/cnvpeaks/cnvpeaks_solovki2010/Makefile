all :  cnvpeaks_solovki2010 

cnvpeaks_solovki2010 : cnvpeaks_solovki2010.o subio.o submem.o
	cc -o cnvpeaks_solovki2010 -lm cnvpeaks_solovki2010.o subio.o submem.o

cnvpeaks_solovki2010.o : cnvpeaks_solovki2010.c subio.h submem.h
	cc -c cnvpeaks_solovki2010.c

clean :
	rm *.o  cnvpeaks_solovki2010
