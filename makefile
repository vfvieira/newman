include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

all: runNewman

runNewman: runNewman.o newman.o network.o
	${CLINKER} -o runNewman -O3 runNewman.o newman.o network.o ${PETSC_SYS_LIB} ${PETSC_MAT_LIB} $(LIBS) -lm

clear:
	rm -f *.o runNewman
