#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>

typedef struct {
    double x, y;
} Particle;

int GRID_W, GRID_H, X_CELLS, Y_CELLS;
int TOTAL_PARTICLES, NUM_THREADS, NUM_ITERATIONS;
double deltaX, deltaY;
int m, n;

void seedRandomPoints(Particle *particles) {
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (m = 0; m < TOTAL_PARTICLES; m++) {
        #pragma omp critical
        {
            particles[m].x = (double) rand() / RAND_MAX;
            particles[m].y = (double) rand() / RAND_MAX;
        }
    }
}

void loadParticleData(FILE *file, Particle *particles) {
    double *buffer = (double *) malloc(2 * TOTAL_PARTICLES * sizeof(double));
    if (!buffer) {
        printf("Buffer allocation failed\n");
        exit(1);
    }

    fread(buffer, sizeof(double), 2 * TOTAL_PARTICLES, file);

    #pragma omp parallel for num_threads(NUM_THREADS)
    for (m = 0; m < TOTAL_PARTICLES; m++) {
        particles[m].x = buffer[2 * m];
        particles[m].y = buffer[2 * m + 1];
    }

    free(buffer);
}

void dumpGridToFile(double *gridData) {
    FILE *outFile = fopen("Mesh1.out", "w");
    if (!outFile) {
        printf("Error creating output file\n");
        exit(1);
    }

    char *buffer = (char *) malloc(GRID_W * 20 * sizeof(char));
    char line[GRID_W * 20];

    for (m = 0; m < GRID_H; m++) {
        line[0] = '\0';
        for (n = 0; n < GRID_W - 1; n++) {
            sprintf(buffer, "%lf ", gridData[m * GRID_W + n]);
            strcat(line, buffer);
        }
        sprintf(buffer, "%lf\n", gridData[m * GRID_W + (GRID_W - 1)]);
        strcat(line, buffer);
        fputs(line, outFile);
    }

    free(buffer);
    fclose(outFile);
}

void distributeCharge(double *gridData, Particle *particles) {
    #pragma omp parallel for num_threads(NUM_THREADS) schedule(static)
    for (m = 0; m < GRID_W * GRID_H; m++) {
        gridData[m] = 0.0;
    }

    double **localGrids = (double **) malloc(NUM_THREADS * sizeof(double *));
    for (m = 0; m < NUM_THREADS; m++) {
        localGrids[m] = (double *) calloc(GRID_W * GRID_H, sizeof(double));
        if (!localGrids[m]) {
            printf("Allocation failed for thread-local grid\n");
            exit(1);
        }
    }

    const int CHUNK = 256;
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        int tid = omp_get_thread_num();
        double *local = localGrids[tid];

        #pragma omp for schedule(guided, CHUNK)
        for (m = 0; m < TOTAL_PARTICLES; m++) {
            double px = particles[m].x;
            double py = particles[m].y;
            double w = 1.0;

            int gx = (int)(px / deltaX);
            int gy = (int)(py / deltaY);

            double fx = px - gx * deltaX;
            double fy = py - gy * deltaY;

            int idx1 = gy * GRID_W + gx;
            int idx2 = gy * GRID_W + (gx + 1);
            int idx3 = (gy + 1) * GRID_W + gx;
            int idx4 = idx3 + 1;

            double w1 = (deltaX - fx) * (deltaY - fy);
            double w2 = fx * (deltaY - fy);
            double w3 = (deltaX - fx) * fy;
            double w4 = fx * fy;

            local[idx1] += w1 * w;
            local[idx2] += w2 * w;
            local[idx3] += w3 * w;
            local[idx4] += w4 * w;
        }
    }

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(static)
    for (m = 0; m < GRID_W * GRID_H; m++) {
        int j;
        for (j = 0; j < NUM_THREADS; j++) {
            gridData[m] += localGrids[j][m];
        }
    }

    for (m = 0; m < NUM_THREADS; m++)
        free(localGrids[m]);
    free(localGrids);
}

int main(int argc, char **argv) {
    char inputFile[50];
    double t_start, t_end, total_time = 0.0;

    if (argc != 3) {
        printf("Usage: %s <input_file> <num_threads>\n", argv[0]);
        return 1;
    }

    strcpy(inputFile, argv[1]);
    NUM_THREADS = atoi(argv[2]);

    FILE *fp = fopen(inputFile, "rb");
    if (!fp) {
        printf("Unable to open file: %s\n", inputFile);
        return 1;
    }

    fread(&X_CELLS, sizeof(int), 1, fp);
    fread(&Y_CELLS, sizeof(int), 1, fp);
    fread(&TOTAL_PARTICLES, sizeof(int), 1, fp);
    fread(&NUM_ITERATIONS, sizeof(int), 1, fp);
    printf("Read from binary file:\n");
    printf("NX: %d, NY: %d, NUM_Points: %d, Maxiter: %d\n", X_CELLS, Y_CELLS, TOTAL_PARTICLES, NUM_ITERATIONS);


    GRID_W = X_CELLS + 1;
    GRID_H = Y_CELLS + 1;
    deltaX = 1.0 / X_CELLS;
    deltaY = 1.0 / Y_CELLS;

    double *gridData = (double *) calloc(GRID_W * GRID_H, sizeof(double));
    Particle *particles = (Particle *) malloc(TOTAL_PARTICLES * sizeof(Particle));

    omp_set_dynamic(0);
    omp_set_nested(1);
    int step;
    for (step = 0; step < NUM_ITERATIONS; step++) {
        loadParticleData(fp, particles);
        t_start = omp_get_wtime();
        distributeCharge(gridData, particles);
        t_end = omp_get_wtime();
        total_time += t_end - t_start;
    }

    dumpGridToFile(gridData);
    printf("Interpolation execution time = %.6lf seconds\n", total_time);

    free(gridData);
    free(particles);
    fclose(fp);
    return 0;
}
