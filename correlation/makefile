
all:
	clear
	swig -c++ -python correlation.i
	python setup_correlation.py build_ext --inplace
	python test_correlation.py
