# Unit  Tests

Tests can be run from pixi.  To run all tests,
```
pixi run doctest
```
Note this will also run the "system tests", which can be rather slow.

To run just the unit tests (not the full calculation tests), run
```
pixi run doctest -ts="*[unit]*"
```
