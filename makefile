SRCDIR := src
INCDIR := include
OBJDIR := obj
CXX = g++ -std=c++11
FFLAG = -I$(INCDIR)

SRC = $(wildcard $(SRCDIR)/*.cpp $(SRCDIR)/*/*.cpp)
OBJ = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRC))


all: DecayWidth.x


$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(FFLAG) -c $< -o $@

%.x:%.cpp $(OBJ)
	$(CXX) $(FFLAG) -o $@ $< $(OBJ)

.PHONY: clean

clean:
	rm -f *.xenophobic
	rm -f $(OBJDIR)/*.o

