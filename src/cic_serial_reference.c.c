#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h> 


// Define structure for points
typedef struct {
    double x, y;
}Points;

// Global variables for grid and points
int GRID_X, GRID_Y, NX, NY;
int NUM_Points,NTHR_C,Maxiter;
double dx, dy;  // Grid spacing

// Particle initialization function
// Use if random location of partilces is required
void initializepoints(Points *points) {
    for (int i = 0; i < NUM_Points; i++) {
        points[i].x = (double) rand() / (double)RAND_MAX;
        points[i].y = (double) rand() / (double)RAND_MAX;
    }
}

void readPoints(FILE *file, Points *points) {
    for (int i = 0; i < NUM_Points; i++) {
        fread(&points[i].x, sizeof(double), 1, file);
        fread(&points[i].y, sizeof(double), 1, file);
    }
}


void printmesh(double *meshValue){
	int i,j;
	FILE *fd=fopen("Mesh.out","w");
	if(fd == NULL){
		printf("File not Created: PIC/out/<DIR_NAME>/rho.out\n");
		exit(1);
	}

	for(i=0;i<GRID_Y;i++){
		for(j=0;j<GRID_X - 1;j++){
			fprintf(fd,"%lf ",meshValue[i*GRID_X+j]);
		}
		fprintf(fd,"%lf\n",meshValue[i*GRID_X+j]);
	}
	fclose(fd);
}



// Charge deposition using Cloud-in-a-Cell (CIC)
void interpolation(double *meshValue, Points *points) {
    
    memset(meshValue, 0, GRID_X * GRID_Y * sizeof(double));
    double *privateMeshValue = (double *) calloc(NTHR_C * GRID_X * GRID_Y, sizeof(double));

    #pragma omp parallel for num_threads(NTHR_C)
    for (int i = 0; i < NUM_Points; i++) {
        
        int id = omp_get_thread_num();
        double px = points[i].x;             //x location
        double py = points[i].y;             //y location
        double weighing = 1.0;

        int gx = (int)(px / dx);                //index value of grid point in x direction
        int gy = (int)(py / dy);                //index value of grid point in y direction

        double lx = px - gx * dx;               //fraction of x length inside the cell
        double ly = py - gy * dy;               //fraction of y length inside the cell

        int p1 = gy * GRID_X + gx;              //Grid point in the lower left boundary
        int p2 = gy * GRID_X + (gx + 1);        //Grid point in the lower right boundary
        int p3 = (gy + 1) * GRID_X + gx;        //Grid point in the upper left boundary
        int p4 = (gy + 1) * GRID_X + (gx + 1);  //Grid point in the upper right boundary

        double a1 = (dx - lx) * (dy - ly);      //Area fraction
        double a2 = lx * (dy - ly);             //Area fraction
        double a3 = (dx - lx) * ly;             //Area fraction
        double a4 = lx * ly;                    //Area fraction

        int offset = id * GRID_X * GRID_Y;      //Offset for parallel execution

        privateMeshValue[offset + p1] += a1 * weighing;     //Stores partial sum from each thread
        privateMeshValue[offset + p2] += a2 * weighing;     //Stores partial sum from each thread
        privateMeshValue[offset + p3] += a3 * weighing;     //Stores partial sum from each thread
        privateMeshValue[offset + p4] += a4 * weighing;     //Stores partial sum from each thread

        
    }


    for(int id=0;id<NTHR_C;id++){
        for(int i=0;i<GRID_X*GRID_Y;i++){
            meshValue[i]+=privateMeshValue[id*GRID_X*GRID_Y+i];   //Reduction to get final answer
        }
    }

free(privateMeshValue);
}

// Main function
int main(int argc, char **argv) {

    char filename[50];
    double start_time, end_time, elapsed_time;

    // Check if correct number of arguments is provided
    if (argc != 3) {
        printf("Usage: %s <input_filename> <num_threads>\n", argv[0]);
        return 1;
    }

    strcpy(filename, argv[1]); // copy the filename
    NTHR_C = atoi(argv[2]);    // Convert thread count to integer

    FILE *file = fopen(filename, "rb"); // Open binary file for reading
    if (file == NULL) {
        printf("Error: Unable to open file %s\n", filename);
        exit(1);
    }

    // Read grid dimensions
    fread(&NX, sizeof(int), 1, file);
    fread(&NY, sizeof(int), 1, file);

    // Read number of Points and max iterations
    fread(&NUM_Points, sizeof(int), 1, file);
    fread(&Maxiter, sizeof(int), 1, file);

    printf("Read from binary file:\n");
    printf("NX: %d, NY: %d, NUM_Points: %d, Maxiter: %d\n", NX, NY, NUM_Points, Maxiter);


    // Since Number of points will be 1 more than number of cells
    //****************************************************** */
    GRID_X = NX + 1; 
    GRID_Y = NY + 1;
    dx = 1.0 / NX;
    dy = 1.0 / NY;
    //****************************************************** */  
          

    // Allocate memory for grid and Points
    double *meshValue = (double *) calloc(GRID_X * GRID_Y, sizeof(double));
    Points *points = (Points *) calloc(NUM_Points, sizeof(Points));


    for(int iteration=0; iteration<Maxiter; iteration++){
        // Initialize Points randomly
        //initializepoints(Points);
        //Use it to randomly initalize points inside the domain
        //you can also use this to generate values instead of using input_filemaker.c
        

        //Read scattered points from file
        readPoints(file, points);

        // Start the timer
        start_time = omp_get_wtime();

        // Perform interpolation
        interpolation(meshValue, points);

        // Stop the timer
        end_time = omp_get_wtime();
        elapsed_time = elapsed_time + end_time - start_time;
    }

    printmesh(meshValue);
    printf("Interpolation execution time = %lf seconds\n", elapsed_time);



    // // Print some results
    // printf("Mesh values at some points:\n");
    // for (int i = 0; i < 5; i++) {
    //     printf("meshValue[%d] = %lf\n", i, meshValue[i]);
    // }

    // Free memory
    free(meshValue);
    free(points);

    return 0;
}
