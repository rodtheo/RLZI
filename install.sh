CUR_DIR=`pwd`
SDSL_INSTALL_PREFIX="${CUR_DIR}/builds/sdsl-lite"
RLZI_INSTALL_PREFIX="${CUR_DIR}/bin"

mkdir -p "${SDSL_INSTALL_PREFIX}" 2> /dev/null

echo "SDSL-lite library will be installed in '${SDSL_INSTALL_PREFIX}'"



cd "${CUR_DIR}/src/sdsl-lite/"
cd build
if [ $? != 0 ]; then
	exit 1
fi
./clean.sh # clean-up build directory
if [ $? != 0 ]; then
	exit 1
fi

cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX="${SDSL_INSTALL_PREFIX}" .. # run cmake
if [ $? != 0 ]; then
	echo "ERROR: CMake build failed."
	exit 1
fi
make sdsl # run make
if [ $? != 0 ]; then
	echo "ERROR: Build failed."
	exit 1
fi
make install # install library
if [ $? != 0 ]; then
	echo "ERROR: Installation failed."
	exit 1
fi

cd "${CUR_DIR}"

echo "SUCCESS: sdsl was installed successfully!"
echo "The sdsl include files are located in '${SDSL_INSTALL_PREFIX}/include'."
echo "The library files are located in '${SDSL_INSTALL_PREFIX}/lib'."
echo " "
echo "RLZI will be installed in '${RLZI_INSTALL_PREFIX}'"
cd src
make
if [ $? != 0 ]; then
	echo "ERROR: Installation of RLZI failed."
	exit 1
fi
cd "${CUR_DIR}"

echo " "
echo "Tests in the test-directory"
echo "A cheat sheet in the extras/cheatsheet-directory."
echo "Have fun!"
