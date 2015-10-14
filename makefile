CXX=g++
CPPFLAGS=-pthread -std=c++11 -O3

all: public

debug: CPPFLAGS += -g
debug: public

public: main.o
	$(CXX) $(CPPFLAGS) -o public $^

clean:
	rm *.o public
