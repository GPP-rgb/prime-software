# Mersenne Prime Lucas-Lehmer Test with FFT

## Description

This project implements the Lucas-Lehmer test for determining the primality of Mersenne numbers using Fast Fourier Transform (FFT) and the GNU Multiple Precision Arithmetic Library (GMP). The code efficiently multiplies large integers and checks if a given Mersenne number M_p = 2^p - 1 is prime.

## Features

- Efficient multiplication of large integers using FFT.
- Implementation of the Lucas-Lehmer test for Mersenne primes.
- Utilizes parallel processing with OpenMP for improved performance.

## Prerequisites

To compile and run this project, ensure you have the following installed:

- [GMP](https://gmplib.org/) - GNU Multiple Precision Arithmetic Library
- [FFTW](http://www.fftw.org/) - Fastest Fourier Transform in the West
- A C++ compiler (e.g., g++)
- [OpenMP](https://www.openmp.org/) - for parallel processing support

## Installation

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd <repository-directory>

2. Install the required libraries. For example, on Ubuntu, you can use:

sudo apt-get install libgmp-dev libfftw3-dev


3. Compile the code using g++:

g++ -o lucas_lehmer_test main.cpp -lgmp -lfftw3 -fopenmp



## Usage

To run the program and test the primality of a Mersenne number, execute the following command:

./lucas_lehmer_test

You can modify the p value in the code to test different Mersenne primes.

## Example Output

Testing primality of 2^398147777 - 1 using Lucas-Lehmer test with FFT...
2^398147777 - 1 is a prime number (Mersenne prime).

## Contributing

Contributions are welcome! If you have suggestions for improvements or found bugs, please create an issue or submit a pull request.
