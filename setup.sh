export HOOMD_ROOT=/home/steven/bin/hoomd-install
export LD_LIBRARY_PATH=/home/steven/bin/hoomd-install/lib/hoomd/python-module:/home/steven/build/yaff/yaff/pes:$LD_LIBRARY_PATH
export LDFLAGS=`~/bin/hoomd-install/bin/hoomd-config.sh --libs`
export CFLAGS=`~/bin/hoomd-install/bin/hoomd-config.sh --cflags`
python setup.py build_ext -i
nosetests yaff/pes/test/test_hoomd.py:test_hoomd_forces
