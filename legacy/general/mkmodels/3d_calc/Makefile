all 	: 3d_add-1.o 3d_sub-1.o 3d_div-1.o 3d_mul-1.o submem.o subio.o iogras.o
	cc -o 3d_add-1 -lm 3d_add-1.o submem.o subio.o iogras.o
	cc -o 3d_sub-1 -lm 3d_sub-1.o submem.o subio.o iogras.o
	cc -o 3d_div-1 -lm 3d_div-1.o submem.o subio.o iogras.o
	cc -o 3d_mul-1 -lm 3d_mul-1.o submem.o subio.o iogras.o

3d_add-1.o : 3d_add-1.c subio.h submem.h iogras.h
	cc -c -Wall 3d_add-1.c
3d_sub-1.o : 3d_sub-1.c subio.h submem.h iogras.h
	cc -c -Wall 3d_sub-1.c
3d_div-1.o : 3d_div-1.c subio.h submem.h iogras.h
	cc -c -Wall 3d_div-1.c
3d_mul-1.o : 3d_mul-1.c subio.h submem.h iogras.h
	cc -c -Wall 3d_mul-1.c

submem.o : submem.c submem.h
	cc -c -Wall submem.c
subio.o : subio.c subio.h
	cc -c -Wall subio.c
iogras.o : iogras.c iogras.h
	cc -c -Wall iogras.c

clean :
	rm 3d_add-1 *.o
