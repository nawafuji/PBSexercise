
EXE=exercise2_linux
SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))

LDFLAGS    = -lglut -lGL
CFLAGS     = -c -O2 -DLINUXVERSION=1


default: $(EXE)

$(EXE): $(OBJECTS)
	$(CXX) -o $(EXE) $(OBJECTS) $(LDFLAGS)

%.o : %.cpp
	$(CXX) $(CFLAGS) -o $@ $<

clean:
	 rm -f *.o $(EXE)

