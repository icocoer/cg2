CC = clang
CXX = clang++
CFLAGS = -Wall
LDFLAGS = -lglut -lGL -lGLU -lm -lSOIL

main: main.c
	$(CC) $< -o $@ $(LDFLAGS)
demo: demo.c solver.c
	$(CC) $^ -o $@ $(LDFLAGS)