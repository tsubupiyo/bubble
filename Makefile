CXX = clang++
CXXFLAGS = -W -Wall -pedantic -std=c++1z -stdlib=libc++
TARGET = test

$(TARGET): test.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@
