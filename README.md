### sshower: A simple parton shower program

#### Requires for compilation:
- C++11 and newer
- GSL (GNU Scientific Library)
- Eigen
- zlib

(The executable is compiled using gcc 14.2.1 for Manjaro linux)

#### Usage:
The program accepts a .lhe file (gz compressed or not) as input, showers it and outputs a .lhe file. It can be run with two arguments (input and output files) as:
```
./sshower /input/file/path /output/file/path
```

