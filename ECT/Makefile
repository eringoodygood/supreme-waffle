CC=g++
CFLAGS=$(shell root-config --cflags)
LDFLAGS=$(shell root-config --glibs)
SOURCES=main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=read

all: clean $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(SOURCES) -o $@ $(LDFLAGS)

$(OBJECTS):
	$(CC) $(CFLAGS) $(SOURCES) -o $@ $(LDFLAGS)

clean:
	rm -f *.o
	rm -f *.*~
