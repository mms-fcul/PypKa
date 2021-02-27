.PHONY: clean clean-test clean-pyc clean-build docs help
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

compile mc:
	cd pypka/mc && python3 mc_setup.py build_ext --inplace

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .pytest_cache

lint: ## check style with flake8	
	mypy
	flake8

format:
	isort *.py clean/*.py
	black pypka

test: ## run tests quickly with the default Python
	cd tests && pytest test.py

coverage: ## check code coverage quickly with the default Python
	cd tests && pytest --cov=../pypka test.py
	cd tests && coverage report -m
	cd tests && coverage html
	$(BROWSER) htmlcov/index.html

report-coverage:
	cd tests && bash <(curl -Ls https://coverage.codacy.com/get.sh) report -r coverage.xml

docs: ## generate Sphinx HTML documentation, including API docs
	rm -f docs/pypka.rst
	rm -f docs/modules.rst
	sphinx-apidoc -o docs/ pypka
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html

docs-coverage:
	cd docs/ && sphinx-build -v -b coverage source/ build/	
	cat docs/build/python.txt

release: dist ## package and upload a release
	twine upload dist/*

dist: clean ## builds source and wheel package
	python setup.py sdist
	python setup.py bdist_wheel
	ls -l dist

install: clean ## install the package to the active Python's site-packages
	python setup.py install
