OBJS := main.o memory.o Ax.o utilities.o cg.o

EXECUTABLE := pollution

# choose compiler:
CC := mpiicc
# CC := gcc

# choose flags:
# flags for Intel compiler icc on maya:
CFLAGS := -O3 -std=c99 -Wall -mkl # -qopenmp
# flags for Portland Group compiler pgcc on maya:
# CFLAGS := -O3 -c99 -Minform=warn -fastsse
# flags for GNU compiler gcc anywhere:
# CFLAGS := -O3 -std=c99 -Wall -Wno-unused-variable

DEFS := -DPARALLEL # -DBLAS
INCLUDES :=
LDFLAGS := -lm

%.o: %.c %.h
	$(CC) $(CFLAGS) $(DEFS) $(INCLUDES) -c $< -o $@

$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) $(DEFS) $(INCLUDES) $(OBJS) -o $@ $(LDFLAGS)

clean:
	-rm -f *.o $(EXECUTABLE)

