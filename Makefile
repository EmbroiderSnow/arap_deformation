remove:
	rm -rf build/

build:
	mkdir -p build
	cd build && cmake .. && make VERBOSE=1

rebuild:
	rm -rf build/
	mkdir -p build
	cd build && cmake .. && make VERBOSE=1

