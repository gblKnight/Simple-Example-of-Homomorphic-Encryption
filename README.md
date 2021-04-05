# Simple Example of Homomorphic Encryption
## Environment Requirement

The required packages are as follows:

* NTL == 11.4.3
* HElib == 1.0.2

## Example to Run the Codes

```
g++ -pthread -g -O2 -std=c++14 -march=native -DFHE_THREADS -DFHE_BOOT_THREADS -o main main.cpp -lhelib -lntl -lgmp -lm
./main
```
