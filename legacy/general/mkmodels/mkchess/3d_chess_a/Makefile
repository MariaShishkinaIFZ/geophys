3d_chess_rot_z : 3d_chess-3.o submem.o subio.o iogras.o
	cc -o 3d_chess_rot_z -lm 3d_chess-3.o submem.o subio.o iogras.o

3d_chess-3.o : 3d_chess-3.c subio.h submem.h iogras.h
	cc -c -Wall 3d_chess-3.c
submem.o : submem.c submem.h
	cc -c -Wall submem.c
subio.o : subio.c subio.h
	cc -c -Wall subio.c
iogras.o : iogras.c iogras.h
	cc -c -Wall iogras.c

clean :
	rm 3d_chess-3 *.o
	
install:
	cp 3d_chess_rot_z $(ISDG_HOME)

