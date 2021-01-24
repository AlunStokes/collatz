#include <time.h>
#include <stdio.h>
#include <inttypes.h>
#include <stdint.h>
/* This is a comment */

int64_t power(int64_t base, unsigned int exp) {
    int i = 1;
    int64_t result = 1;
    for (i = 0; i < exp; i++)
        result *= base;
    return result;
 }

int is_C_self_cont(int64_t n) {
  int64_t d = n;
  while (n > 1) {
    if (n % 2 == 0) {
      n = n / 2;
    }
    else {
      n = 3 * n + 1;
    }
    if (n % d == 0) {
      return 1;
    }
  }
  return 0;
}

void *get_C_self_cont(int64_t L, int64_t U, int offset, int skip) {
  int64_t n = 0;
  if (L % 2 == 0) {
    n = L + 2 * offset + 1;
  }
  else {
    n = L + 2 * offset;
  }
  while (n <= U) {
    if (is_C_self_cont(n) == 1) {
      printf("%" PRIu64 "\n", n);
    }
    n = n + 2 * skip;
  }
}

int main()
{
    clock_t begin = clock();

    get_C_self_cont(1, power(10, 8), 0, 1);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Took %fs\n", time_spent);

    return 0;
}
