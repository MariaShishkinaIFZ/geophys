# Inversion Software Development Group
# ISDG makefile
#  
# written by Sergey A. Tikhotski sometime
# modified by Ilya V. Fokin at
# 
# 25.01.2011
#
# Russia, Institute of Physics of the Earth (IPE)
#

MODIFICATION = shahta

all :  cnvpeaks_${MODIFICATION}

cnvpeaks_${MODIFICATION} : cnvpeaks_${MODIFICATION}.o subio.o submem.o
	cc -o cnvpeaks_${MODIFICATION} -lm cnvpeaks_${MODIFICATION}.o subio.o submem.o

cnvpeaks_${MODIFICATION}.o : cnvpeaks_${MODIFICATION}.c subio.h submem.h
	cc -c cnvpeaks_${MODIFICATION}.c

clean :
	rm *.o  cnvpeaks_${MODIFICATION}
