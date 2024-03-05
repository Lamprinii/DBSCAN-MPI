#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <string.h>

/******************************************************************************/

#define UNCLASSIFIED -1
#define NOISE -2

#define CORE_POINT 1
#define NOT_CORE_POINT 0

#define SUCCESS 0
#define FAILURE -3

/******************************************************************************/

typedef struct point_s point_t;
struct point_s {
    double x, y, z;
    int cluster_id;
};

typedef struct node_s node_t;
struct node_s {
    unsigned int index;
    node_t *next;
};

typedef struct epsilon_neighbours_s epsilon_neighbours_t;
struct epsilon_neighbours_s {
    unsigned int num_members;
    node_t *head, *tail;
};

typedef struct {
    unsigned int index;
    int cluster_id;
} IndexClusterPair;


/******************************************************************************/

node_t *create_node(unsigned int index);
int append_at_end(unsigned int index, epsilon_neighbours_t *en);
epsilon_neighbours_t *get_epsilon_neighbours(unsigned int index, point_t *points,unsigned int num_points, double epsilon, double (*dist)(point_t *a, point_t *b),unsigned int start_index,unsigned int end_index);
void print_epsilon_neighbours(point_t *points, epsilon_neighbours_t *en);
void destroy_epsilon_neighbours(epsilon_neighbours_t *en);
void dbscan(point_t *points,unsigned int num_points, double epsilon, unsigned int minpts, double (*dist)(point_t *a, point_t *b), unsigned int start_index,unsigned int end_index,int start,int end, int my_rank,
            int p,int *local_cluster_id,unsigned int overlap);
int expand(unsigned int index, unsigned int cluster_id, point_t *points, unsigned int num_points,
           double epsilon, unsigned int minpts, double (*dist)(point_t *a, point_t *b), unsigned int start_index,unsigned int end_index,int start,int end);
int spread(unsigned int index, epsilon_neighbours_t *seeds, unsigned int cluster_id, point_t *points, unsigned int num_points,
           double epsilon, unsigned int minpts, double (*dist)(point_t *a, point_t *b), unsigned int start_index,unsigned int end_index,int start,int end);
double euclidean_dist(point_t *a, point_t *b);
double adjacent_intensity_dist(point_t *a, point_t *b);
unsigned int parse_input(FILE *file, point_t **points, double *epsilon, unsigned int *minpts,int my_rank);
void print_points(point_t *points, unsigned int num_points);
void change_cluster_id(point_t *points, unsigned int num_points, int my_rank, int start, int end, int *all_cluster_ids);
void change_overlap(point_t *points, unsigned int num_points, int my_rank, int p,int start, int end, int start_index, int end_index, unsigned int overlap);


/******************************************************************************/

node_t *create_node(unsigned int index)
{
    node_t *n = (node_t *)calloc(1, sizeof(node_t));
    if (n == NULL)
        perror("Failed to allocate node.");
    else {
        n->index = index;
        n->next = NULL;
    }
    return n;
}

/******************************************************************************/

int append_at_end(unsigned int index, epsilon_neighbours_t *en)
{
    node_t *n = create_node(index);
    if (n == NULL) {
        free(en);
        return FAILURE;
    }
    if (en->head == NULL) {
        en->head = n;
        en->tail = n;
    } else {
        en->tail->next = n;
        en->tail = n;
    }
    ++(en->num_members);
    return SUCCESS;
}


/******************************************************************************/

epsilon_neighbours_t *get_epsilon_neighbours(unsigned int index, point_t *points,unsigned int num_points, double epsilon, double (*dist)(point_t *a, point_t *b),unsigned int start_index,unsigned int end_index)
{
    epsilon_neighbours_t *en = (epsilon_neighbours_t *)calloc(1, sizeof(epsilon_neighbours_t));
    if (en == NULL) {
        perror("Failed to allocate epsilon neighbours.");
        return en;
    }

    
    for (int i = start_index; i < end_index; ++i) {
        if (i == index)
            continue;
        if (dist(&points[index], &points[i]) > epsilon)
            continue;
        else {
            if (append_at_end(i, en) == FAILURE) {
                destroy_epsilon_neighbours(en);
                en = NULL;
                break;
            }
        }
    }
    return en;
}

/******************************************************************************/


void print_epsilon_neighbours(point_t *points, epsilon_neighbours_t *en)
{
    if (en) {
        node_t *h = en->head;
        while (h) {
            printf("(%lfm, %lf, %lf)\n", points[h->index].x, points[h->index].y, points[h->index].z);
            h = h->next;
        }
    }
}

/******************************************************************************/

void destroy_epsilon_neighbours(epsilon_neighbours_t *en)
{
    if (en) {
        node_t *t, *h = en->head;
        while (h) {
            t = h->next;
            free(h);
            h = t;
        }
        free(en);
    }
}

/******************************************************************************/

// void dbscan(point_t *points,unsigned int num_points, double epsilon, unsigned int minpts, double (*dist)(point_t *a, point_t *b), unsigned int start_index,unsigned int end_index,int start,int end,
//                     int my_rank,int p,int *local_cluster_id)
// {
//     unsigned int i,cluster_id = 0;

//     for (i = start_index; i <= end_index; ++i) {
//         if (points[i].cluster_id == UNCLASSIFIED) {
//             if (expand(i, cluster_id, points, num_points, epsilon, minpts, dist,start_index,end_index,start,end) == CORE_POINT){
//                 ++cluster_id;
//             }
//         }    
//     }

//     *local_cluster_id = cluster_id;
// }

void dbscan(point_t *points,unsigned int num_points, double epsilon, unsigned int minpts, double (*dist)(point_t *a, point_t *b), unsigned int start_index,unsigned int end_index,int start,int end,
                    int my_rank,int p,int *local_cluster_id, unsigned int overlap)
{
    unsigned int i,cluster_id = 0;

    for (i = start_index; i <= end_index; ++i) {
        if (points[i].cluster_id == UNCLASSIFIED) {
            if (expand(i, cluster_id, points, num_points, epsilon, minpts, dist,start_index,end_index,start,end) == CORE_POINT){
                ++cluster_id;
            }
        }    
    }

    *local_cluster_id = cluster_id;


}


/******************************************************************************/

// int expand(unsigned int index, unsigned int cluster_id, point_t *points, unsigned int num_points,
//            double epsilon, unsigned int minpts, double (*dist)(point_t *a, point_t *b), unsigned int start_index,unsigned int end_index,int start,int end)
// {
//     int return_value = NOT_CORE_POINT;
//     epsilon_neighbours_t *seeds = get_epsilon_neighbours(index, points, num_points, epsilon, dist,start_index,end_index);
//     if (seeds == NULL)
//         return FAILURE;
//     if(start <= index && index <= end){
//         if (seeds->num_members < minpts)
//             points[index].cluster_id = NOISE;
//         else {
//             points[index].cluster_id = cluster_id;
//             node_t *h = seeds->head;
//             while (h) {
//                 points[h->index].cluster_id = cluster_id;
//                 h = h->next;
//             }

//             h = seeds->head;
//             while (h) {
//                 spread(h->index, seeds, cluster_id, points, num_points, epsilon, minpts, dist,start_index,end_index,start,end);
//                 h = h->next;
//             }

//             return_value = CORE_POINT;
//         }
//     }
//     destroy_epsilon_neighbours(seeds);
//     return return_value;
// }



int expand(unsigned int index, unsigned int cluster_id, point_t *points, unsigned int num_points,
           double epsilon, unsigned int minpts, double (*dist)(point_t *a, point_t *b), unsigned int start_index,unsigned int end_index,int start,int end)
{
    int return_value = NOT_CORE_POINT;
    epsilon_neighbours_t *seeds = get_epsilon_neighbours(index, points, num_points, epsilon, dist,start_index,end_index);
    if (seeds == NULL)
        return FAILURE;
    
    if (seeds->num_members < minpts)
        points[index].cluster_id = NOISE;
    else {
        points[index].cluster_id = cluster_id;
        node_t *h = seeds->head;
        while (h) {
            points[h->index].cluster_id = cluster_id;
            h = h->next;
        }

        h = seeds->head;
        while (h) {
               spread(h->index, seeds, cluster_id, points, num_points, epsilon, minpts, dist,start_index,end_index,start,end);
               h = h->next;
        }
        if(start <= index && index <= end){
            return_value = CORE_POINT;
        }  
    }
    
    
    destroy_epsilon_neighbours(seeds);
    return return_value;
}



/******************************************************************************/

// int spread(unsigned int index, epsilon_neighbours_t *seeds, unsigned int cluster_id, point_t *points, unsigned int num_points,
//            double epsilon, unsigned int minpts, double (*dist)(point_t *a, point_t *b), unsigned int start_index,unsigned int end_index,int start,int end)
// {
//     epsilon_neighbours_t *spread = get_epsilon_neighbours(index, points, num_points, epsilon, dist,start_index,end_index);
//     if (spread == NULL)
//         return FAILURE;
//     if (start <= index && index <= end) {
//         if (spread->num_members >= minpts) {
//             node_t *n = spread->head;
//             point_t *d;
//             while (n) {
//                 d = &points[n->index];
//                 if (d->cluster_id == NOISE ||
//                     d->cluster_id == UNCLASSIFIED) {
//                     if (d->cluster_id == UNCLASSIFIED) {
//                         if (append_at_end(n->index, seeds) == FAILURE) {
//                             destroy_epsilon_neighbours(spread);
//                             return FAILURE;
//                         }
//                     }
//                     d->cluster_id = cluster_id;
//                 }
//                 n = n->next;
//             }
//         }
//     }

//     destroy_epsilon_neighbours(spread);
//     return SUCCESS;
// }

int spread(unsigned int index, epsilon_neighbours_t *seeds, unsigned int cluster_id, point_t *points, unsigned int num_points,
           double epsilon, unsigned int minpts, double (*dist)(point_t *a, point_t *b), unsigned int start_index,unsigned int end_index,int start,int end)
{
    epsilon_neighbours_t *spread = get_epsilon_neighbours(index, points, num_points, epsilon, dist,start_index,end_index);
    if (spread == NULL)
        return FAILURE;
    if (spread->num_members >= minpts) {
        node_t *n = spread->head;
        point_t *d;
        while (n) {
            d = &points[n->index];
            if (d->cluster_id == NOISE ||
                d->cluster_id == UNCLASSIFIED) {
                if (d->cluster_id == UNCLASSIFIED) {
                    if (append_at_end(n->index, seeds) == FAILURE) {
                        destroy_epsilon_neighbours(spread);
                        return FAILURE;
                    }
                }
                d->cluster_id = cluster_id;
            }
            n = n->next;
        }
    }

    destroy_epsilon_neighbours(spread);
    return SUCCESS;
}
/******************************************************************************/


double euclidean_dist(point_t *a, point_t *b)
{
    return sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2) + pow(a->z - b->z, 2));
}

/******************************************************************************/

unsigned int parse_input(FILE *file, point_t **points, double *epsilon, unsigned int *minpts,int my_rank)
{
    unsigned int num_points, i = 0;

    fscanf(file, "%lf %u %u\n", epsilon, minpts, &num_points);
    point_t *p = (point_t *)calloc(num_points, sizeof(point_t));
    if (p == NULL) {
        perror("Failed to allocate points array");
        return 0;
    }
    while (i < num_points) {
        fscanf(file, "%lf %lf %lf\n", &(p[i].x), &(p[i].y), &(p[i].z));
        p[i].cluster_id = UNCLASSIFIED;
        ++i;
    }
    *points = p;
    return num_points;
}

/******************************************************************************/

void print_points(point_t *points, unsigned int num_points)
{
    unsigned int i = 0;
    printf("Number of points: %u\n"
        " x     y     z     cluster_id\n"
        "-----------------------------\n"
        , num_points);
    while (i < num_points) {
          printf("%5.2lf %5.2lf %5.2lf: %d\n", points[i].x, points[i].y, points[i].z, points[i].cluster_id);
          ++i;
    }
}

/******************************************************************************/

void change_cluster_id(point_t *points, unsigned int num_points, int my_rank, int start, int end, int *all_cluster_ids)
{
    unsigned int i = 0;

    for(i=start ; i <= end ; i++) {
          if(points[i].cluster_id >= 0){
            if(my_rank != 0){
                points[i].cluster_id = points[i].cluster_id + all_cluster_ids[my_rank-1] - my_rank;;
            }
        }
    }
}

/******************************************************************************/

void change_overlap(point_t *points, unsigned int num_points, int my_rank, int p,int start, int end, int start_index, int end_index, unsigned int overlap)
{
    IndexClusterPair* pairs = malloc( (2*overlap) * sizeof(IndexClusterPair));
    int k=0;

    if(my_rank != 0){
        for(int j = start_index; j < start; ++j){
            if(points[j].cluster_id != -2){
                pairs[k].index = j;
                pairs[k].cluster_id = points[j].cluster_id;
                ++k;
            }
        }
    }
    if(my_rank != p-1){
        for(int j = end+1; j <= end_index; ++j){
            if(points[j].cluster_id != -2){
                pairs[k].index = j;
                pairs[k].cluster_id = points[j].cluster_id;
                ++k;
            }
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    //---------------Create an MPI_Datatype for the pairs array--------------------//

    MPI_Datatype MPI_Pairs_type;
    int nblocks = 2;
    int blocklen[2] = {1, 1};
    MPI_Datatype oldtypes[2] = {MPI_UNSIGNED,MPI_INT};
    MPI_Aint displ[2];

    // Get offsets of the struct members
    MPI_Get_address(&(((IndexClusterPair*)0)->index), &displ[0]);
    MPI_Get_address(&(((IndexClusterPair*)0)->cluster_id), &displ[1]);

    displ[1] -= displ[0];
    displ[0]  = 0;


    MPI_Type_create_struct(nblocks, blocklen, displ, oldtypes, &MPI_Pairs_type);
    MPI_Type_commit(&MPI_Pairs_type);


    //------------------------------------------------------------------------------//
    //-----------Gather the pairs array back to processor 0------------------------//

    IndexClusterPair *all_pairs;
    int *recvcounts;
    int *displs;


    recvcounts = (int *)malloc(p * sizeof(int));
        displs = (int *)malloc(p * sizeof(int));

    MPI_Gather(&k, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calculate recvcounts and displs for MPI_Gatherv 
    MPI_Bcast(recvcounts,p,MPI_INT,0,MPI_COMM_WORLD);

    int total_pairs = 0;
    for (int i = 0; i < p; ++i) {
        displs[i] = total_pairs;
        total_pairs += recvcounts[i];
    }
    

    //Allocate memory for final_points array
    all_pairs = malloc(total_pairs * sizeof(IndexClusterPair));


    MPI_Barrier(MPI_COMM_WORLD);

    //Gather each processor's points array in the range that it was processing.
    MPI_Gatherv(pairs, recvcounts[my_rank], MPI_Pairs_type, all_pairs, recvcounts, displs, MPI_Pairs_type, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);


    //Broadcast the all_pairs array to every process
    MPI_Bcast(all_pairs,total_pairs,MPI_Pairs_type,0,MPI_COMM_WORLD);

    //Make the changes in the recv

    for(int i=0; i<total_pairs;i++){
        int index = all_pairs[i].index;
        if(points[index].cluster_id == -2){
            points[index].cluster_id = all_pairs[i].cluster_id;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

/******************************************************************************/






int main(int argc, char *argv[]) {
    double tStart, tStop;
    point_t *points;
    point_t *points_copied; //the copied points array that gets broadcasted to the other processes
    double epsilon;
    unsigned int minpts;
    unsigned int num_points;
    
    unsigned int index_start, index_end; //start and end index with overlaping regions
    int start,end; //start and end index without overlaping regions
    int my_rank, p;
   
    //p = number of processors

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if(my_rank == 0 ){
        num_points = parse_input(stdin, &points, &epsilon, &minpts,my_rank);
    }
    

    MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&minpts, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_points, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&p, 1, MPI_INT, 0, MPI_COMM_WORLD);


    MPI_Barrier(MPI_COMM_WORLD);


    //-----------------------------------------------------------------------//


    
    MPI_Barrier(MPI_COMM_WORLD);

    

    unsigned int num = num_points / p;
    unsigned int overlap = num_points * 0.10; //10% of the points will be in the overlapping region


    //We need at least the same number of points as the processes.
    if (num_points < p) {
        if (my_rank == 0) {
            printf("You have to use at least %d total points\n", p);
        }  
        free(points);
        MPI_Finalize();
        exit(0);
    }


    start = (my_rank* num);
    end = start + num - 1;
    

    if((start-overlap) <= 0 || my_rank == 0){
        index_start = start;
    }
    else{
        index_start = start-overlap;
    }

    if((end+overlap) >= num_points || my_rank == p-1){
        end=num_points;
        index_end = end;
    }
    else{
        index_end = end + overlap;
    }


    //------------------------------------------------------------------------------//
    //---------------Create an MPI_Datatype for the points array--------------------//

    MPI_Datatype MPI_Points_type;
    int nblocks = 2;
    int blocklen[2] = {3, 1};
    MPI_Datatype oldtypes[2] = {MPI_DOUBLE,MPI_INT};
    MPI_Aint displ[2];
    point_t pt;
    MPI_Get_address(&(pt.x), &displ[0]);
    MPI_Get_address(&(pt.cluster_id), &displ[1]);

    displ[1] -= displ[0];
    displ[0]  = 0;



    MPI_Type_create_struct (nblocks, blocklen, displ, oldtypes, &MPI_Points_type);
    MPI_Type_commit(&MPI_Points_type);

    //------------------------------------------------------------------------------//
    //-----------Send the array of points from process 0 to all processes-----------//

    MPI_Barrier(MPI_COMM_WORLD);

    //Allocate memory for the points array for all processes except of process 0.
    if(my_rank != 0){
        points = (point_t *)malloc(num_points * sizeof(point_t));
    }
    
    MPI_Bcast(points, num_points, MPI_Points_type, 0, MPI_COMM_WORLD);

    //------------------------------------------------------------------------------//
    //---------------dbscan algorithm to every processor----------------------------//

    // Allocate memory for local_cluster_id
    int local_cluster_id;

    if (num_points > p) {
        if(my_rank == 0 ){
            tStart = clock();
        }

        dbscan(points, num_points, epsilon, minpts, euclidean_dist,index_start,index_end,start,end,my_rank,p,&local_cluster_id,overlap);
    }

    //Barrier to wait for all of the processes to finish their dbscan algorithms
    MPI_Barrier(MPI_COMM_WORLD);

    //Allocate memory for the array of the sums of cluster_id
    int *all_cluster_ids = (int *)malloc(p * sizeof(int));
    

    //Gather the local_cluster_ids of each processes to the all_cluster_ids array on process 0 (root)
    MPI_Gather(&local_cluster_id, 1, MPI_INT, all_cluster_ids, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Broadcast the cumulative sums back to all processes
    MPI_Bcast(all_cluster_ids, p, MPI_INT, 0, MPI_COMM_WORLD);



    //------------------------------------------------------------------------------//
    //---------------------------Rewrite the cluster_ids----------------------------//

    int increment = 0;
    for(int i=0;i<p ; i++){
        all_cluster_ids[i] = all_cluster_ids[i] + increment;
        increment = all_cluster_ids[i];
    }

    change_cluster_id(points,num_points,my_rank,start,end,all_cluster_ids);

    MPI_Barrier(MPI_COMM_WORLD);

    //------------------------------------------------------------------------------//
    //----------------Change border points in overlap region-----------------------//


    change_overlap(points, num_points, my_rank, p,start, end, index_start, index_end, overlap);

    MPI_Barrier(MPI_COMM_WORLD);

    //------------------------------------------------------------------------------//
    //-----------Gather the points array back to processor 0------------------------//


    point_t *all_points;
    int *recvcounts;
    int *displs;
    point_t *final_points;


    recvcounts = (int *)malloc(p * sizeof(int));
        displs = (int *)malloc(p * sizeof(int));

    // Calculate recvcounts and displs for MPI_Gatherv 
    //Last procesor has num_points - start points always


    for(int i=0;i<p ; i++){
        recvcounts[i] = (i == p - 1 ) ? (num_points - i*num) : num ;
    }

    
    for (int i = 0; i < p; i++) {
        displs[i] =  i* num ;
    }

    //Allocate memory for final_points array
    final_points = (point_t *)malloc(num_points * sizeof(point_t));


    MPI_Barrier(MPI_COMM_WORLD);

    //Gather each processor's points array in the range that it was processing.
    MPI_Gatherv(&points[displs[my_rank]], recvcounts[my_rank], MPI_Points_type, final_points, recvcounts, displs, MPI_Points_type, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);




    //------------------------------------------------------------------------------//
    //---------------------Print results--------------------------------------------//

    if (my_rank == 0) {

        //TODO BELOW
        //Gather the points array and update it on the golbal_custer_ids array
        //MPI_Gatherv(points,num,)


        tStop = clock();
        double t = (tStop-tStart)/CLOCKS_PER_SEC;
        printf("t = %f\n", t);
        printf("Epsilon: %lf\n", epsilon);
        printf("Minimum points: %u\n", minpts);
        print_points(final_points, num_points);
    }

    
    // Free memory
    if (my_rank == 0) {
        free(final_points);
        free(recvcounts);
        free(displs);
    }

    free(all_cluster_ids);
    free(points);
    free(points_copied);
    MPI_Finalize();

    return 0;
}

/******************************************************************************/

