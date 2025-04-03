rm -rf build
mkdir build
cd build/
cmake ..
cmake --build .
./elliptical_orbit_tests
./circular_orbit_tests
./attitude_tests
./misc_tests
gcovr -r .. --filter ../src/ --filter ../include/ --html-details output.html