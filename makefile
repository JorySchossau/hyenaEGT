CXX=g++
CPPFLAGS=-pthread -std=c++11 -lcurses

release: CXX += -O3
release: public

debug: CXX += -gdwarf-3 -g3 -pg
debug: public

public: main.o
	$(CXX) $(CPPFLAGS) -o public $^

clean:
	rm *.o public
