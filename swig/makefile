python:
	swig -python -c++ -o titanlib_wrap_python.cpp titanlib.i
	g++ -c titanlib.cpp titanlib_wrap_python.cpp -I/usr/include/python2.7 -fPIC -rdynamic
	g++ -shared titanlib.o titanlib_wrap_python.o -o _titanlib.so -lc -lstdc++ -L/usr/lib/gcc/x86_64-linux-gnu/7/

r:
	swig -c++ -r -o titanlib_wrap_r.cpp titanlib.i
	R CMD SHLIB -o titanlib.so titanlib_wrap_r.cpp titanlib.cpp

clean:
	rm -f titanlib*.o titanlib*.so titanlib.py titanlib.R titanlib_wrap* *~ *.pyc .Rhistory
