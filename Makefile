remove:
	rm -rf build/

build:
	cd build && make VERBOSE=1

rebuild:
	rm -rf build/
	mkdir -p build
	cd build && cmake .. && make VERBOSE=1

test:
	cd build && ./arap-deformation ../data/deformation.obj 