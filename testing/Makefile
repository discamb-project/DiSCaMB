.PHONY: test

test: test.cpp
	bear -- clang++ -I../include -std=c++17 -O0 -g $< -o $@ ./build/build/lib/libdiscambd.a -lsycl
