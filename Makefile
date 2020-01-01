.PHONY: default clean sdist upload

default: clean sdist upload

sinstall:
	pip install dist/kitsune-*.tar.gz

stest:
	bash test.sh

build:
	python setup.py build

install:
	python setup.py install

sdist:
	python setup.py sdist bdist_wheel

upload:
	twine upload dist/*

clean:
	rm -r build/* dist/* kitsune.egg-info
