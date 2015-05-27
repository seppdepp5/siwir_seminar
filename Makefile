CC      = /usr/bin/g++
CFLAGS  = -O3 -Wall -Winline -Wshadow -std=c++11
LDFLAGS = #-lm -pthread
OBJDIR = ./objects
SRCDIR = ./src/

OBJ = $(addprefix $(OBJDIR)/, MGSolver.o Smoother.o Array.o Debug.o Stencil.o main.o )

make: $(OBJ)
	$(CC) $(CFLAGS) -o mgsolve $(OBJ)

$(OBJDIR)/%.o: $(SRCDIR)%.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm $(OBJDIR)/*.o
