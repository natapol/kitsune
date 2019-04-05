.PHONY: default install clean

all: build

build:
	python setup.py build
install:
	python setup.py install
sdist:
	python setup.py sdist bdist_wheel
	twine upload dist/*
clean:
	rm -r build dist kitsune.egg-info
