.PHONY: default install clean

all: build

build:
	python setup.py build
install:
	python setup.py install
clean:
	rm -r build dist kali.egg-info
