.PHONY: python clean

python: radia.so

radia.so:
	python setup.py build_ext --build-lib='../gcc'
	cp ../gcc/radia*.so ../../env/radia_python/
	rm -rf build

clean:
	rm -rf build
	rm -f ../gcc/radia*.so
	rm -f ../../env/radia_python/radia*.so