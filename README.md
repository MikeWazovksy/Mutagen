The idea of ​​the program is to mutate one private key from the previous puzzle to find the solution to the next one.

We take, for example, the key to the 19th puzzle 0000000000000000000000000000000000000000000000000000000000005749f and the address of the next one in hash160 : b907c3a2a3b27789dfb509b730dd47703c272868

and mutate its private key bits, this allows us to reduce the number of options for enumerating options

build : make

# Solve puzzle 38 with 8 threads and override flip count to 21

./puzzle_solver -p 38 -t 8 -f 21

# Solve puzzle 20 with default settings

./puzzle_solver

# Show help

./puzzle_solver -h

# Key Features

Command-line options:

-p/--puzzle: Puzzle number (20-68)

-t/--threads: Number of CPU cores to use

-f/--flips: Override default flip count

-h/--help: Show usage information

# Optimizations:

AVX2 vectorization for cryptographic operations

Multi-threading with OpenMP

Batched processing for better cache utilization

The solver is designed to be both efficient and easy to use, with sensible defaults that can be overridden as needed.

Idea Denevron !

dotation address :

bc1qa3c5xdc6a3n2l3w0sq3vysustczpmlvhdwr8vc

Thanks for the implementation nomachine!

donation address :

bc1qdwnxr7s08xwelpjy3cc52rrxg63xsmagv50fa8
