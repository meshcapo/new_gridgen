## Commands
- `cmake -B build -S .` - Configure the project (run once, or after CMakeLists.txt changes)
- `cmake --build build` - Build the project (gridexec)
- `cmake --build build --target tests_runner` - Build test suite
- `rm -rf build/` - Clean all build artifacts

## Architecture
- Header files in include/
- Source files in source/
- Tests in tests/
- Files used for regression tests in tests/files

## Conventions
- Use snake\_case for functions and variables
- Check all malloc() return values
- Free all allocated memory before function return
- Use Criterion as the testing framework for writing tests

## Compiler, Flags and Dependencies
- Use clang
- Debug build (default): -g -O0; Release: -O3 -DNDEBUG (set via -DCMAKE_BUILD_TYPE=Release at configure time)
- Include path: include/
- Dependencies: GSL (GNU Scientific Library), Blas/Lapack, math library, criterion, pthread

## Tests
- Directory tests/files/smoothbump contains files for regression testing file read and write functions in file\_operations.c
- These files come from the smoothbump example from SU2
