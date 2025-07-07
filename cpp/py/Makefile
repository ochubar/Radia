
PACKAGE = Radia_Python_Interface
PREFIX ?=
PIP ?= pip
PYTHON ?= python3
ifeq ($(CONDA_PREFIX),)
	PREFIX = sudo -H
	PIP = pip-sirius
	PYTHON = python-sirius
endif
MODE ?= 0

python: radia.so

radia.so:
	MODE=$(MODE) $(PYTHON) setup.py build_ext --build-lib='../gcc'
	cp ../gcc/radia*.so ../../env/radia_python/radia.so
	rm -rf build

clean:
	rm -rf build
	rm -f ../gcc/radia*.so
	rm -f ../../env/radia_python/radia*.so

install: uninstall python
	$(PREFIX) $(PIP) install --no-deps ./
	$(PREFIX) git clean -fdX

uninstall:
	$(PREFIX) $(PIP) uninstall -y $(PACKAGE)

.PHONY: python clean install uninstall
