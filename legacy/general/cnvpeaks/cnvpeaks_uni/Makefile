all :  cnvpeaks_uni 

cnvpeaks_uni : cnvpeaks_uni.o subio.o submem.o
	cc -o cnvpeaks_uni -lm cnvpeaks_uni.o subio.o submem.o

cnvpeaks_uni.o : cnvpeaks_uni.c subio.h submem.h
	cc -c cnvpeaks_uni.c

clean :
	rm *.o  cnvpeaks_uni
