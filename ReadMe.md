# Hierarchical HSP

#### Python Bindings
- python version 3.7.4
- virtual environment
```bash
    python3 -m venv env
    python3 -m pip install numpy
    python3 -m pip install matplotlib
```
- compilation: 
```bash
    mkdir build
    cd build
    cmake -S .. -B .
    make -j
```
- to use:
```python
    from HHSP import HHSP

    dimension = 3
    datasetSize = 10
    alg = HHSP(dimension, datasetSize)
    res = alg.HSP(0)
```