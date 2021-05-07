all :  cnvpeaks_all 

cnvpeaks_all : cnvpeaks_all.o subio.o submem.o
	cc -o cnvpeaks_all -lm cnvpeaks_all.o subio.o submem.o

cnvpeaks_all.o : cnvpeaks_all.c subio.h submem.h
	cc -c cnvpeaks_all.c

clean :
	rm *.o  cnvpeaks_all
