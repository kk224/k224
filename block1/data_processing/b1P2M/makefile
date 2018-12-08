CC=g++

SRC_FILES=$(wildcard *.cpp)
OBJ_FILES=$(patsubst %.cpp,%.o,$(SRC_FILES))
M4RI_PATH=/home/kk224/program/m4ri-20130416
CUDD_PATH=/home/kk224/program/cudd-3.0.0
OBJ_LINK_FILES = $(M4RI_PATH)/brilliantrussian.o $(M4RI_PATH)/echelonform.o $(M4RI_PATH)/graycode.o \
	 $(M4RI_PATH)/misc.o $(M4RI_PATH)/mmc.o \
	$(M4RI_PATH)/mzd.o $(M4RI_PATH)/mzp.o $(M4RI_PATH)/ple.o \
	$(M4RI_PATH)/ple_russian.o $(M4RI_PATH)/solve.o $(M4RI_PATH)/strassen.o \
	$(M4RI_PATH)/triangular.o $(M4RI_PATH)/triangular_russian.o

CFLAG=-c -O3 -I/$(M4RI_PATH) -I $(CUDD_PATH)/cudd -I $(CUDD_PATH)/util -I $(CUDD_PATH)/dddmp -I $(CUDD_PATH)/mtr -I $(CUDD_PATH)/epd -I $(CUDD_PATH)/st 
LDFLAG= -lpthread -L $(CUDD_PATH)/cudd/.libs/ -lcudd -lm
TARGET=b1P2M

all:
	$(CC) $(CFLAG) $(SRC_FILES) 
	$(CC) -o $(TARGET)  $(OBJ_FILES) $(OBJ_LINK_FILES) $(LDFLAG) 
.PHONY: clean
clean:
	rm $(OBJ_FILES) $(TARGET) -rf
