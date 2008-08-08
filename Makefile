CXX?=g++
WARN=-Wfloat-equal -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wmissing-noreturn -Winline
INCLUDE=-I$(HOME)/zTools/include/tipsylib
#LIBS=-L$(HOME)/zTools/lib64
#LIBS=-L$(HOME)/zTools/lib
CPPFLAGS=-g -O3 -Wall $(WARN) $(LISTING) $(INCLUDE) $(LIBS)
#CPPFLAGS=-O3 -Wall $(WARN) $(LISTING) -g -pg #-ftree-vectorize -ftree-vectorizer-verbose=5


tip2gad: 
	$(CXX) $(CPPFLAGS) -o $@ tip2gad.cpp $(HOME)/zTools/lib64/libtipsy.a
	#$(CXX) $(CPPFLAGS) -ftree-vectorize -ftree-vectorizer-verbose=5 -o $@ $< -L. -ltipsy

gad2tip: 
	$(CXX) $(CPPFLAGS) -o $@ gad2tip.cpp $(HOME)/zTools/lib64/libtipsy.a

verify_gadget:
	$(CXX) $(CPPFLAGS) -o $@ verify_gadget.cpp

verify_tipsy:
	$(CXX) $(CPPFLAGS) -o $@ verify_tipsy.cpp $(HOME)/zTools/lib64/libtipsy.a

all: tip2gad gad2tip verify_gadget verify_tipsy

# DO NOT DELETE

#tip2gad: tip2gad.hpp 

