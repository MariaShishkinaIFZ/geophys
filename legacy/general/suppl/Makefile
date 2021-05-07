all : getpeaks cnvpeaks sl3d sl2d 2dXYgmt 2dXZgmt 2dYZgmt
getpeaks : getpeaks.o subio.o submem.o
	cc -o getpeaks -lm getpeaks.o subio.o submem.o
cnvpeaks : cnvpeaks.o subio.o submem.o
	cc -o cnvpeaks -lm cnvpeaks.o subio.o submem.o
sl3d: sl3d.o subio.o submem.o
	cc -o sl3d -lm sl3d.o subio.o submem.o
sl2d: sl2d.o subio.o submem.o
	cc -o sl2d -lm sl2d.o subio.o submem.o
2dXYgmt: 2dXYgmt.o subio.o submem.o
	cc -o 2dXYgmt -lm 2dXYgmt.o subio.o submem.o
2dXZgmt: 2dXZgmt.o subio.o submem.o
	cc -o 2dXZgmt -lm 2dXZgmt.o subio.o submem.o
2dYZgmt: 2dYZgmt.o subio.o submem.o
	cc -o 2dYZgmt -lm 2dYZgmt.o subio.o submem.o

getpeaks.o : getpeaks.c subio.h submem.h
	cc -c getpeaks.c
cnvpeaks.o : cnvpeaks.c subio.h submem.h
	cc -c cnvpeaks.c
sl3d.o : sl3d.c subio.h submem.h
	cc -c sl3d.c
sl2d.o : sl2d.c subio.h submem.h
	cc -c sl2d.c
2dXYgmt.o : 2dXYgmt.c subio.h submem.h gmt_nan.h
	cc -c 2dXYgmt.c
2dXZgmt.o : 2dXZgmt.c subio.h submem.h gmt_nan.h
	cc -c 2dXZgmt.c
2dYZgmt.o : 2dYZgmt.c subio.h submem.h gmt_nan.h
	cc -c 2dYZgmt.c

clean :
	rm *.o getpeaks cnvpeaks sl3d 2dXYgmt 2dXZgmt 2dYZgmt

install :
	cp  {getpeaks,cnvpeaks,sl3d,2dXYgmt,2dXZgmt,2dYZgmt} $(ISDG_HOME)/
