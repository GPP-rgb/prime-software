#include <gmp.h>
#include <fftw3.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>


// Convert a large GMP number to a coefficient array for FFT
void gmp_to_coefficients(const mpz_t num, fftw_complex* coeffs, size_t size) {
    // Loop through the number bits and assign them to FFT coefficients
    for (size_t i = 0; i < size; ++i) {
        coeffs[i][0] = mpz_tstbit(num, i) ? 1.0 : 0.0;  // real part: binary coefficient
        coeffs[i][1] = 0.0;  // imaginary part is zero
    }
}

// Rebuild a large GMP number from FFT coefficients
void coefficients_to_gmp(mpz_t result, const fftw_complex* coeffs, size_t size) {
    mpz_set_ui(result, 0);

    for (size_t i = 0; i < size; ++i) {
        // Reconstruction with rounded real parts
        unsigned long coeff = (std::round(coeffs[i][0]) > 0.5) ? 1 : 0;  // Binary approximation
        mpz_mul_2exp(result, result, 1);  // Left shift (multiply by 2)
        mpz_add_ui(result, result, coeff);  // Add the coefficient
    }
}

// Ensure that the size is the next power of 2
size_t next_power_of_2(size_t x) {
    return 1 << static_cast<int>(std::ceil(std::log2(x)));
}

// Multiplication function using FFT with parallelization
void fft_multiply(mpz_t result, const mpz_t a, const mpz_t b) {
    // Size of the large numbers
    size_t n = mpz_sizeinbase(a, 2);
    size_t m = mpz_sizeinbase(b, 2);
    size_t size = next_power_of_2(n + m);  // Ensure size is a power of 2

    // Prepare arrays for FFT
    fftw_complex* fa = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
    fftw_complex* fb = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
    fftw_complex* fc = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);

    // Convert the large numbers to coefficients for FFT
    gmp_to_coefficients(a, fa, size);
    gmp_to_coefficients(b, fb, size);

    // Enable parallelism in FFTW
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());

    // Initialize FFT plans with a more patient strategy for better performance
    fftw_plan plan_a = fftw_plan_dft_1d(size, fa, fa, FFTW_FORWARD, FFTW_PATIENT);
    fftw_plan plan_b = fftw_plan_dft_1d(size, fb, fb, FFTW_FORWARD, FFTW_PATIENT);
    fftw_plan plan_c = fftw_plan_dft_1d(size, fc, fc, FFTW_BACKWARD, FFTW_PATIENT);

    // Perform FFT on the coefficients
    fftw_execute(plan_a);
    fftw_execute(plan_b);

    // Multiply the coefficients in the frequency domain
    #pragma omp parallel for
    for (size_t i = 0; i < size; ++i) {
        fc[i][0] = fa[i][0] * fb[i][0] - fa[i][1] * fb[i][1];  // real part
        fc[i][1] = fa[i][0] * fb[i][1] + fa[i][1] * fb[i][0];  // imaginary part
    }

    // Inverse FFT to retrieve the product
    fftw_execute(plan_c);

    // Rebuild the large GMP integer from the coefficients
    coefficients_to_gmp(result, fc, size);

    // Free memory and cleanup FFTW threads
    fftw_destroy_plan(plan_a);
    fftw_destroy_plan(plan_b);
    fftw_destroy_plan(plan_c);
    fftw_free(fa);
    fftw_free(fb);
    fftw_free(fc);
    fftw_cleanup_threads();  // Cleanup threads after execution
}

// Lucas-Lehmer function with FFT and GMP
bool lucas_lehmer_test(unsigned long long p) {
    // Initialize GMP variables
    mpz_t s, Mp, tmp;
    mpz_init_set_ui(s, 4);  // s = 4
    mpz_init(tmp);

    // Compute Mp = 2^p - 1
    mpz_init(Mp);
    mpz_ui_pow_ui(Mp, 2, p);
    mpz_sub_ui(Mp, Mp, 1);

    // Lucas-Lehmer loop
    for (unsigned long long i = 0; i < p - 2; ++i) {
        fft_multiply(tmp, s, s);   // s = s^2 using FFT
        mpz_sub_ui(s, tmp, 2);     // s = s - 2
        mpz_mod(s, s, Mp);         // s = s % Mp
    }

    bool is_prime = (mpz_cmp_ui(s, 0) == 0);  // If s == 0, Mp is prime

    // Free GMP variables
    mpz_clear(s);
    mpz_clear(Mp);
    mpz_clear(tmp);

    return is_prime;
}

int main() {
    unsigned long long p = 398147777 ;

    std::cout << "Testing primality of 2^" << p << " - 1 using Lucas-Lehmer test with FFT..." << std::endl;

    if (lucas_lehmer_test(p)) {
        std::cout << "2^" << p << " - 1 is a prime number (Mersenne prime)." << std::endl;
    } else {
        std::cout << "2^" << p << " - 1 is not a prime number." << std::endl;
    }

    return 0;
}
