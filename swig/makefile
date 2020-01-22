python:
	swig -python -c++ example.i
	g++ -c example.cpp example_wrap.cxx -I/usr/include/python2.7 -fPIC -rdynamic
	g++ -shared example.o example_wrap.o -o _example.so -lc -lstdc++ -L/usr/lib/gcc/x86_64-linux-gnu/7/

r:
	swig -c++ -r -o example_wrap.cpp example.i
	R CMD SHLIB -o example.so example_wrap.cpp example.cpp

	#PKG_LIBS="example.cxx" R CMD SHLIB example_wrap.cpp
	#g++  -I/usr/share/R/include -DNDEBUG      -fPIC  -g -O2 -fdebug-prefix-map=/build/r-base-AitvI6/r-base-3.4.4=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c example_wrap.cpp -o example_wrap.o
	#g++ -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions -Wl,-z,relro -o example_wrap.so example_wrap.o example.cxx -L/usr/lib/R/lib -lR -fPIC


	#g++ -c example.cpp example_wrap.cxx -I/home/thomasn/disk1/metnosw/stow/R-2.15.2-amd64/lib/R/include/ -fPIC -rdynamic
	#g++ -shared example.o example_wrap.o -o _example.so -lc -lstdc++ -L/usr/lib/gcc/x86_64-linux-gnu/7/

clean:
	rm example.o *.so example.py *_wrap.* *~ *.pyc
