CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio1262/cplex
INCLUDE=$(CPLEX_HOME)/include
STATIC_LIB=$(CPLEX_HOME)/lib/x86-64_linux/static_pic
DYNAMIC_LIB=$(CPLEX_HOME)/bin/x86-64_linux
CPLEX_LIB_STATIC=cplex
CPLEX_LIB_DYNAMIC=cplex1262
LFLAGS=-lrt -lm -lpthread
CFLAGS=-DIL_STD -I$(INCLUDE) -O0 -g3 -c -fmessage-length=0 -std=c++0x
CC=g++

BOXES:	main_boxes.o
	$(CC) -L$(STATIC_LIB) -o BOXES ./main_boxes.o -l$(CPLEX_LIB_STATIC) $(LFLAGS)

main_boxes.o:	main_boxes.cpp boxes2.h
	$(CC) $(CFLAGS) main_boxes.cpp

clean:
	rm main_boxes.o BOXES

