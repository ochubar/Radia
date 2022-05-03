.PHONY: python clean

MODE ?= 0 

python: radia.so

radia.so:
# MODE=$(MODE) python setup.py build_ext --build-lib='../gcc' 
	MODE=$(MODE) python-sirius setup.py install
	cp ../gcc/radia*.so ../../env/radia_python/radia.so
	rm -rf build

clean:
	rm -rf build
	rm -f ../gcc/radia*.so
	rm -f ../../env/radia_python/radia*.so
