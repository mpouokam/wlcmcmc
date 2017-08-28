install: swig
	pip install numpy && pip install .[all]

devinstall: swig
	pip install numpy && pip install -e .[all]

swig:
	swig -c++ -python wlcmcmc/wlcmcmclib.i

clean:
	\rm -rf build wlcmcmc/*_wrap.cpp wlcmcmc/*_wrap.cxx wlcmcmc/_*.so _wlcmcmclib.so wlcmcmclib wlcmcmc/wlcmcmclib.py


clidebug:
	g++ wlcmcmc/wlcmcmclib.cpp wlcmcmc/mtrand/mtrand.cpp -o wlcmcmclib
	./wlcmcmclib

venv: swig
	# remove virtual env
	(conda env list | grep 'devenv-wlcmcmc' -c && conda env remove --yes -q -n devenv-wlcmcmc)

	# create new virtual env
	conda create --yes -q -n devenv-wlcmcmc python=2.7

	# set the right environment and install packages
	source activate devenv-wlcmcmc && pip install numpy && pip install -e .[all]
