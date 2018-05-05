DIR_INC = ./inc
DIR_SRC = ./src
DIR_OBJ = ./obj
BINDIR=/usr/local/bin

CPP_SRC = $(wildcard ${DIR_SRC}/*.cpp)
CPP_OBJ = $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${CPP_SRC}))
D_SRC = $(wildcard ${DIR_SRC}/*.d)
D_OBJ = $(patsubst %.d,${DIR_OBJ}/%.o,$(notdir ${D_SRC}))

SRC = ${CPP_SRC} ${D_SRC}
OBJ = ${CPP_OBJ} ${D_OBJ}

TARGET = fastp

BIN_TARGET = ${TARGET}

CC = g++
CFLAGS = -std=c++11 -g -I${DIR_INC}
DC = dmd
DFLAGS = -g -I${DIR_SRC}

${BIN_TARGET}: ${OBJ}
	$(CC) $(OBJ) -lphobos2 -lz -lpthread -o $@

${DIR_OBJ}/%.o: ${DIR_SRC}/%.cpp make_obj_dir
	$(CC) $(CFLAGS) -O3 -c  $< -o $@

${DIR_OBJ}/%.o: ${DIR_SRC}/%.d make_obj_dir
	$(DC) $(DFLAGS) -O -c  $< -of$@

.PHONY:clean
clean:
	rm obj/*.o
	rm $(TARGET)

make_obj_dir:
	@if test ! -d $(DIR_OBJ) ; \
	then \
		mkdir $(DIR_OBJ) ; \
	fi

install:
	install $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."
