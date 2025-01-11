#include <stdio.h>
#include <stdlib.h>

int main() {
    int *arr;
    int n = 5;

    // Allocating memory for an array of 5 integers
    arr = (int *)malloc(n * sizeof(int));
    if (arr == NULL) {
        printf("Memory allocation failed!\n");
        return 1;  // Exit if allocation fails
    }

    // Assign values to the array
    for (int i = 0; i < n; i++) {
        arr[i] = i + 1;
    }

    // Print the values of the array
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }

    // Free the allocated memory
    free(arr);
    return 0;
}

