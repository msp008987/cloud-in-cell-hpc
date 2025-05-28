#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define FILENAME "input.bin"  // Output binary filename

void generateBinaryInputFile(int NX, int NY, int NUM_Points, int Maxiter) {
    FILE *file = fopen(FILENAME, "wb"); // Open file in binary write mode
    if (file == NULL) {
        printf("Error: Unable to create file %s\n", FILENAME);
        exit(1);
    }

    // Write grid dimensions (NX, NY)
    fwrite(&NX, sizeof(int), 1, file);
    fwrite(&NY, sizeof(int), 1, file);

    // Write number of Points, OpenMP threads, and max iterations
    fwrite(&NUM_Points, sizeof(int), 1, file);
    fwrite(&Maxiter, sizeof(int), 1, file);

    // Seed random number generator
    srand(time(NULL));

    // Generate and write point positions for `Maxiter` iterations
    for (int iter = 0; iter < Maxiter; iter++) {
        for (int i = 0; i < NUM_Points; i++) {
            double x = (double)rand() / RAND_MAX;  // Random x in [0,1]
            double y = (double)rand() / RAND_MAX;  // Random y in [0,1]
            fwrite(&x, sizeof(double), 1, file);
            fwrite(&y, sizeof(double), 1, file);
        }
    }

    fclose(file);
    printf("Binary file '%s' generated successfully!\n", FILENAME);
}

int main() {
    int NX, NY, NUM_Points, Maxiter;

    // Get user input
    printf("Enter grid dimensions (NX NY): ");
    scanf("%d %d", &NX, &NY);

    printf("Enter number of points: ");
    scanf("%d", &NUM_Points);

    printf("Enter max iterations: ");
    scanf("%d", &Maxiter);

    // Generate the binary input file
    generateBinaryInputFile(NX, NY, NUM_Points, Maxiter);

    return 0;
}
