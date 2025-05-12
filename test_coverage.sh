rm -rf build
mkdir build
cd build/
cmake .. -DCMAKE_BUILD_TYPE=DEBUG
cmake --build .
./Satellite_tests
./utils_tests
./gs_tests
gcovr -r .. --filter ../src/ --filter ../include/ --json-summary ../tests/test_coverage_summary.json --html-details ../tests/test_coverage_detailed.html
cd ..
rm test_plot.png