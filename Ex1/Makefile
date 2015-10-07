EXE= ./build/exercise1_mac
OBJDIR= ./build

SOURCES = $(wildcard *.cpp)
OBJECTS = $(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))

# -L... needed for OpenGL with OSX 10.5
LDFLAGS    = -framework CoreServices -framework glut -framework OpenGL -L/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries
CFLAGS     = -c -O2 -DMACVERSION=1

default: $(EXE)

$(EXE): $(OBJECTS)
	$(CXX) -o $(EXE) $(OBJECTS) $(LDFLAGS) 

$(OBJDIR)/%.o : %.cpp
	$(CXX) $(CFLAGS) -o $@ $<

clean:
	 rm -f $(OBJECTS) $(EXE)

