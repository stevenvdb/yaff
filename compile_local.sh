export LIBRARY_PATH=:/home/steven/build/lammps-ro/src/
export CPLUS_INCLUDE_PATH=/home/steven/build/lammps-ro/src/
export CC=gcc; python setup.py build_ext -i; export CC=mpic++; python setup.py build_ext -i
export CC=gcc; python setup.py install --home=~; export CC=mpic++; python setup.py install --home=~
