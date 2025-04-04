rm -rf build
mkdir build
cd build/
cmake .. -DCMAKE_BUILD_TYPE=DEBUG
cmake --build .
./Satellite_tests
./utils_tests
gcovr -r .. --filter ../src/ --filter ../include/ --html-details output.html
cd ..
rm test_plot.png