#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include <float.h>


#define NUMLETTERS 26

//Linked List to store each row of data
struct record {
    float* arr;
    struct record* next;
};
typedef struct record record;

//converts row into float array 
float* convert(char* c, int col)
{  
	float *a=(float*)malloc(col*sizeof(float));
    char * token = strtok(c, ",");
    int i = 0;
    while( token != NULL ) {
        if(i != 0 && i != 1)
            a[i-2] = atof(token);
        i++;
        token = strtok(NULL, ",");
    }
    return a;
}

//Find Maximum in each column
float* findAnswer(float* data, int row, int col){
    float* res = (float*)malloc(sizeof(float) * col);
    for(int i = 0; i < col; i++){
        res[i] = FLT_MAX;
    }
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            if(res[j] > data[i*col+j])
                res[j] = data[i*col + j];
        }
    }
    return res;
}
 

int main(int argc, char *argv[] ) {

    MPI_Init(&argc, &argv);
    int my_rank, comm_size;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int rows = 0, col = 0;
    
    float *arr;
    if(my_rank == 0){
        FILE *my_file;
        char* str = (char*)malloc(sizeof(char) * 1000);
        my_file = fopen(argv[1], "r");
        int got = fscanf(my_file, "%s", str);

        //Calculating Number of Columns
        for(int i = 0; i < strlen(str); i++) {  
            if(str[i] == ','){
                col++;
            }     
        }

        col -= 1;
        printf("\n");
        free(str);
        char* s = (char*)malloc(sizeof(char) * 1000);

        //Reading one Row
        got = fscanf(my_file, "%s", s);
        record* head = NULL;
        while(got == 1){
            rows++;
            //Converting row to float array
            float *a=convert(s, col);
            //Saving float array to linkedlist
            if(!head){
                head = (record*)malloc(sizeof(record));
                head->arr = a;
                head->next = NULL;
            }
            else{
                record* temp = (record*)malloc(sizeof(record));
                temp->arr = a;
                temp->next = head;
                head = temp;
            }
            //Reading New Row
            got = fscanf(my_file, "%s", s);
        }
        fclose(my_file);
        
        printf("\nColumns = %d, Rows = %d\n", col, rows);

        arr = (float *)malloc(rows * col * sizeof(float));
        
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < col; j++){
                arr[i*col+j] = head->arr[j];
            }
            record* temp = head;
            head = head->next;
            free(temp);
        }
        
        printf("\nLinked List to Array Converted\n");
    }


    //MPI WORK Started
    clock_t t;
    t = clock();

    //BCAST Rows and Columns to other process from Rank-0
    MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&col, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //Calculating Number of Rows for each Rank
    int load_per_node = rows/comm_size;
    int load_on_first = rows - load_per_node * (comm_size-1);


    if(my_rank == 0){
        printf("Load on 0th rank = %d\nLoad on other = %d\ncomm_size = %d\n", load_on_first, load_per_node, comm_size);
    }

    //Stores for each process
    float *data;

    //Sending part of data to other processess
    if(my_rank == 0){
        MPI_Request request;
        for(int i = 1; i < comm_size; i++){
            MPI_Isend(&arr[load_per_node*(i-1)*col], load_per_node*col, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &request);
        }
        if(comm_size != 1)
            MPI_Wait(&request, MPI_STATUS_IGNORE);
        data  = (float *)malloc(load_on_first * col * sizeof(float));
        for(int i = load_per_node*(comm_size-1)*col; i < rows*col; i++){
                data[i-load_per_node*(comm_size-1)*col] = arr[i];
        }
    }
    else{
        //Data Receiving
        int received;
        data  = (float *)malloc(load_per_node * col * sizeof(float));
        MPI_Recv(data, load_per_node*col, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    
        }

    
    //Calculating Min value in each row
    float* result;
    if(my_rank == 0){
        result = findAnswer(data, load_on_first, col);
    }
    else{
        result = findAnswer(data, load_per_node, col);
    }
    
  

    //For Storing Minimum Value of each columns
    float* finalAnswer;
    if(my_rank == 0){
        finalAnswer = (float *)malloc(col * sizeof(float));
    } 
    
    //Reducing minimum value of each column from each processs
    MPI_Reduce(result, finalAnswer, col, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);


    //Minimum of Minimum in each column
    float allMin = FLT_MAX; 
    if(my_rank == 0){
        for(int j = 0; j < col; j++){
            if(allMin > finalAnswer[j])
                allMin = finalAnswer[j];
        }
    }

    //Calculating Time till here
    t = clock() - t;

    float time_taken = ((float)t)/CLOCKS_PER_SEC;
    float max_time;

    //Recucing Max time taken 
    MPI_Reduce(&time_taken, &max_time, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);


    if(my_rank == 0){
        //Writing in File
        FILE *filePointer ;
        filePointer = fopen("output.txt", "w") ;
        for(int i = 0; i < col; i++){
            if(i == col-1)
                fprintf(filePointer,"%0.2f", finalAnswer[i]);
            else
                fprintf(filePointer,"%0.2f,", finalAnswer[i]);
        }
        fprintf(filePointer, "\n%0.2f", allMin);
        fprintf(filePointer, "\n%f", max_time);
        
        fclose(filePointer) ;

        //Writing for plotting purpose
        FILE *FP ;
        FP = fopen("Temp_output.txt", "a") ;
        fprintf(FP, "%f\n",max_time);
        fclose(FP) ;

    }

    printf("\nRank %d Completed\n", my_rank);

    MPI_Finalize();
 
    return EXIT_SUCCESS;
} 


//Reference
//Geeksforgeeks.com
//stackoverflow.com
