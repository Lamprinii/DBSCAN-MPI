#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

struct Point {
    double x, y, z;
    double distance;
};

int comparePoints(const void *a, const void *b) {
    return ((struct Point*)a)->distance - ((struct Point*)b)->distance;
}

int main(int argc, char *argv[]) {
	int		epsilon, minpts, num_points, error, i, rand_seed;
	double		low, high, x, y, z;
	struct timeval	rand_init_timestamp;
	char		filename[256], output[256];
	FILE		*file;

	double tStart, tStop; //measure time


    	

    	
  	printf("Enter value for epsilon: ");
	scanf("%d", &epsilon);

	if (epsilon < 1) {
		printf("Invalid value for \"epsilon\".\n");
		exit(1);
	}

	printf("Enter value for minimum number of points that consist a cluster: ");
	scanf("%d", &minpts);

	if (minpts < 1) {
		printf("Invalid value for \"minpts\".\n");
		exit(1);
	}

	printf("Enter number of points to be clustered: ");
	scanf("%d", &num_points);

	if (num_points < 1) {
		printf("Invalid value for \"num_points\".\n");
		exit(1);
	}

	printf("Enter lowest value for point coordinates : ");
	scanf("%lf", &low);

	printf("Enter highest value for point coordinates: ");
	scanf("%lf", &high);

	printf("Enter output file name: ");
	scanf("%s", filename);

	struct Point *points = malloc(num_points * sizeof(struct Point));

	if (points == NULL) {
        	printf("Memory allocation error.\n");
       	 	exit(1);
    	}


	file = fopen(filename, "w");
	if (file == NULL) {
		printf("Could not open file \"%s\".\n", filename);
		exit(1);
	}

	error = gettimeofday(&rand_init_timestamp, NULL);
        if (error != 0) {
                printf("gettimeofday() returned an error.\n");
                exit(1);
        }

        rand_seed = rand_init_timestamp.tv_sec * 1000000 + rand_init_timestamp.tv_usec;

	srand48(rand_seed);

	sprintf(output, "%d %d %d\n", epsilon, minpts, num_points);
	fwrite(output, sizeof(char), strlen(output), file);


    for (i = 0; i < num_points - 1; i++) {
        x = low + (high - low) * drand48();
        y = low + (high - low) * drand48();
        z = low + (high - low) * drand48();
        //sprintf(output, "%.2f %.2f %.2f\n", x, y, z);
        //fwrite(output, sizeof(char), strlen(output), file);

        points[i].x = x;
        points[i].y = y;
        points[i].z = z;
        points[i].distance = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

    }

    x = low + (high - low) * drand48();
    y = low + (high - low) * drand48();
    z = low + (high - low) * drand48();
    //sprintf(output, "%.2f %.2f %.2f", x, y, z);
    //fwrite(output, sizeof(char), strlen(output), file);

    points[num_points - 1].x = x;
    points[num_points - 1].y = y;
    points[num_points - 1].z = z;
    points[num_points - 1].distance = sqrt(x*x + y*y + z*z);


    

    // Sort the points based on distance
    qsort(points, num_points, sizeof(struct Point), comparePoints);

    tStop = clock();
    double t = (tStop - tStart) / CLOCKS_PER_SEC;
    
    // Now 'points' array is sorted based on distance, you can use it as needed

	// Write the sorted points to the file
	for (i = 0; i < num_points; i++) {
	    sprintf(output, "%.2f %.2f %.2f\n", points[i].x, points[i].y, points[i].z);
	    fwrite(output, sizeof(char), strlen(output), file);
	}

	printf("Time taken: %f\n", t);

	fclose(file);



	    free(points);

	    return 0;
	}
