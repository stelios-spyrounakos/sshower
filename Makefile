# compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2
LDFLAGS = -lgsl -lgslcblas -lz

# executable name
EXEC = sshower

# source and object files
SRCS = alphaS.cpp splittings.cpp shower_helpers.cpp showering.cpp \
	kinematic_reco.cpp lhe_io.cpp progress_bar.cpp main.cpp
OBJS = $(SRCS:.cpp=.o)

# headers (for dependency tracking)
HEADERS = alphaS.hpp splittings.hpp shower_helpers.hpp showering.hpp \
	kinematic_reco.hpp lhe_io.hpp progress_bar.hpp constants.hpp

# default target
all: $(EXEC)

# link the object files to create the executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJS) $(LDFLAGS)

# compile each source file into an object file
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# clean up object files and the executable
clean:
	rm -f $(OBJS) $(EXEC)

# target to force recompilation
rebuild: clean all
