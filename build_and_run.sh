rm -rf build
mkdir build
cd build/
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
cmake --build .
./run