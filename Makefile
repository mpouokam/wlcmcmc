all:
	swig -c++ -python wlcmclib.i
	python wlcmclib_swigsetup.py build_ext --inplace

install:
	cp _domontecarlosteps.so ..
	cp domontecarlosteps.py ..

clean:
	\rm -rf build *_wrap.cpp *_wrap.cxx


clidebug:
	g++ wlcmclib.cpp random.cpp mtrand.cpp -o wlcmclib
	./wlcmclib