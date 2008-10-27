CXX?=g++

TIPSYLIBINC=$(HOME)/zTools/include/tipsylib
TIPSYLIB=$(HOME)/zTools/lib64/

INCLUDE=-I$(TIPSYLIBINC)

CPPFLAGS=-g -O3 -Wall $(WARN) $(LISTING) $(INCLUDE) $(LIBS)
WARN=-Wfloat-equal -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wmissing-noreturn -Winline

tip2gad: 
	$(CXX) $(CPPFLAGS) -o $@ tip2gad.cpp $(TIPSYLIB)/libtipsy.a
	#$(CXX) $(CPPFLAGS) -ftree-vectorize -ftree-vectorizer-verbose=5 -o $@ $< -L. -ltipsy

gad2tip: 
	$(CXX) $(CPPFLAGS) -o $@ gad2tip.cpp $(TIPSYLIB)/libtipsy.a

verify_gadget:
	$(CXX) $(CPPFLAGS) -o $@ verify_gadget.cpp

verify_tipsy:
	$(CXX) $(CPPFLAGS) -o $@ verify_tipsy.cpp $(TIPSYLIB)/libtipsy.a

all: tip2gad gad2tip verify_gadget verify_tipsy

clean:
	rm -f tip2gad gad2tip verify_gadget verify_tipsy *.o

# DO NOT DELETE

#tip2gad: tip2gad.hpp 

