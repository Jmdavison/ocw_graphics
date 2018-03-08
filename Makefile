INCFLAGS = -I /usr/include/GL
INCFLAGS += -I /usr/include/vecmath

LINKFLAGS = -lglut -lGL -lGLU
LINKFLAGS += -lvecmath

CFLAGS    = -O2 -Wall -ansi -DSOLN -ggdb
CC        = g++
SRCS      = main.cpp parse.cpp curve.cpp surf.cpp camera.cpp
OBJS      = $(SRCS:.cpp=.o)
PROG      = a1

all: $(SRCS) $(PROG)

$(PROG): $(OBJS)
	$(CC) $(FRAMEWORKS) $(CFLAGS) $(OBJS) -o $@ $(INCFLAGS) $(LINKFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -c -o $@ $(INCFLAGS)

depend:
	makedepend $(INCFLAGS) -Y $(SRCS)

clean:
	rm $(OBJS) $(PROG)
