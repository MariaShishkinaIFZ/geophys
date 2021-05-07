all :  cnvpeaks_solovki 

cnvpeaks_solovki : cnvpeaks_solovki.o subio.o submem.o
	cc -o cnvpeaks_solovki -lm cnvpeaks_solovki.o subio.o submem.o

cnvpeaks_solovki.o : cnvpeaks_solovki.c subio.h submem.h
	cc -c cnvpeaks_solovki.c

clean :
	rm *.o  cnvpeaks_solovki
