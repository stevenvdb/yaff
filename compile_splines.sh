#export LIBRARY_PATH=/home/steven/build/Splines/lib/
#export LD_LIBRARY_PATH=/home/steven/lib/python/horton/grid/
#export C_INCLUDE_PATH=/home/steven/build/Splines/srcs/
#export CPLUS_INCLUDE_PATH=/home/steven/build/Splines/srcs/:/home/steven/build/Splines/srcs_utils/:/home/steven/build/horton/horton/grid/
#export CC=gcc; python setup.py build_ext -i; export CC=mpic++; python setup.py build_ext -i
#export CC=gcc; python setup.py install --home=~
export CC=g++; python setup.py install --home=~
#export CC=g++; python setup.py build_ext -i
