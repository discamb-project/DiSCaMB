.PHONY: test

test: test.cpp
	bear -- ${CXX} -I../include -std=c++17 -fsycl -fsycl-targets=nvptx64-nvidia-cuda -O0 -g $< -o $@ ./build/build/lib/libdiscambd.a -lsycl
