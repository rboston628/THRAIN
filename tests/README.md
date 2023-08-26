# Unit  Tests

These are some unit tests to ensure correct operation of the program.  They require `CxxTest` to run.  If you do not already have `CxxTest`, macos users can download with 
```
brew install CxxTest
```
Linux users can try the same with their own package manager.  Windows users should consider installing a real operating system.

From the top directory, the tests can be compiled and run using 
```
make tests
./tests/tests.out
```

A lot of output is generated, most of it not meaningful to the tests.  If at the end you see `OK!`, then the tests all ran and passed.