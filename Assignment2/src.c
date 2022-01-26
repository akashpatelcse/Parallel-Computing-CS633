#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <time.h>

//Structure for every Rank
typedef struct proData{
    int dataSize;
    MPI_Comm Rack_comm;
    MPI_Comm Node_comm;
    MPI_Comm Node_leaders_comm;
    MPI_Comm Rack_leaders_comm;
    int Rack_colour;
    int Node_colour;
    int Node_Leader_colour;  
    int Rack_Leader_colour;     
    int world_size, rank_in_world;
    int rack_size, rank_in_rack;
    int node_size, rank_in_node;
    int real_Node_Leader;
    int real_Rack_Leader;
    int ProcessorNum;
    int root;      //Rank in World
} processData;


//Rack Details
int rack1[] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 31};                     //15
int rack2[] = {13, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32};         //16
int rack3[] = {33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46};                     //13
int rack4[] = {45, 47, 48, 49, 50, 51, 52, 53, 54, 56, 58, 59, 60, 61};                 //14
int rack5[] = {62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78};     //17
int rack6[] = {79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92};                 //14

//Checks if the node is in Rack or not
int inArray(int node, int rack[], int Rack_Size){
    for(int i=0; i<Rack_Size; i++)
        if(rack[i] == node)
            return 1;
    return 0;
}

//Gives Color to each node According to Rack in which it belong to.
int giveColor(int NodeNumber){
    int color = 0;
    if(inArray(NodeNumber, rack1, 15))
        return 1;
    if(inArray(NodeNumber, rack2, 16))
        return 2;
    if(inArray(NodeNumber, rack3, 13))
        return 3;
    if(inArray(NodeNumber, rack4, 14))
        return 4;
    if(inArray(NodeNumber, rack5, 17))
        return 5;
    if(inArray(NodeNumber, rack6, 14))
        return 6;
    return color;
}

//Converts System name to integer ex: CSEWS32 -> 32
int char2int (char *arr, size_t n)
{    
    int ProcessorNum = 0;
    for(int i = 0; i < n; i++){
        int num = arr[i];
        if(num <= 57 && num >= 48){
            if(ProcessorNum != 0)
                ProcessorNum *= 10;
            ProcessorNum += num-48;
        }
    }
    return ProcessorNum;
}


//Optimal BCAST
void Opti_Bcast(processData *pro){
  //Buffer to be send
    double buffer[pro->dataSize];
    for(int i = 0; i < pro->dataSize; i++){
      buffer[i] = pro->rank_in_world;
    }
    //Send data from root to Rank-0 Process.
    if(pro->root == pro->rank_in_world){
        MPI_Send(buffer,pro->dataSize,MPI_DOUBLE,0, 0, MPI_COMM_WORLD);
    }
    //Rank-0 Process will receive the data
    if(0 == pro->rank_in_world){
        MPI_Recv(buffer,pro->dataSize,MPI_DOUBLE, pro->root , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    //Now Rank-0 Process will send data to all Rack leaders
    if(pro->real_Rack_Leader == 1)
        MPI_Bcast(&buffer, pro->dataSize, MPI_DOUBLE, 0, pro->Rack_leaders_comm);
    //Rack Leaders sends data to Node Leaders
    if(pro->real_Node_Leader == 1)
        MPI_Bcast(&buffer, pro->dataSize, MPI_DOUBLE, 0, pro->Node_leaders_comm);
    //Node Leaders will send data to processess inside it.
    MPI_Bcast(&buffer, pro->dataSize, MPI_DOUBLE, 0, pro->Node_comm);
}


//Optimal Reduce
void Opti_Reduce(processData *pro){
    //We will use buffer and result in round robin nature to save the memory.
    double buffer[pro->dataSize];
    double result[pro->dataSize];
    //Generating Random Numbers to Buffer
    for(int i = 0; i < pro->dataSize; i++){
      buffer[i] = rand();
    }
    //Every process is sending data to its respective Node Leader
    MPI_Reduce(buffer, result, pro->dataSize, MPI_DOUBLE, MPI_SUM, 0, pro->Node_comm);
    //Node Leader will send data to Rack Leader.
    if(pro->real_Node_Leader)
        MPI_Reduce(result, buffer, pro->dataSize, MPI_DOUBLE, MPI_SUM, 0, pro->Node_leaders_comm);
    //Rack Leader Sends DAta to Root Process.
    if(pro->real_Rack_Leader)
        MPI_Reduce(buffer, result, pro->dataSize, MPI_DOUBLE, MPI_SUM, 0, pro->Rack_leaders_comm);

}


//Optimal Gather
void Opti_Gather(processData *pro, int ppn){
    //Every Process Data
    double my_value[pro->dataSize];
    //Stores Node Level Data
    double stores_node_level[pro->dataSize * pro->node_size];

    for(int i = 0; i < pro->dataSize; i++){
      my_value[i] = rand(); 
    }

    //Startng Node Level Collection
    if(pro->rank_in_node == 0)
    {
        MPI_Gather(&my_value, pro->dataSize, MPI_DOUBLE, stores_node_level, pro->dataSize, MPI_DOUBLE, 0, pro->Node_comm);
    }
    else
    {
        MPI_Gather(&my_value, pro->dataSize, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, pro->Node_comm);
    }

    

    if(pro->real_Node_Leader){
       //Starting Rack Level Collection
        double stores_Rack_level[pro->dataSize * pro->rack_size];
        if(pro->rank_in_rack == 0){
            MPI_Gather(&stores_node_level, pro->dataSize * pro->node_size, MPI_DOUBLE, stores_Rack_level, pro->dataSize * pro->node_size, MPI_DOUBLE, 0, pro->Node_leaders_comm);
        }
        else
        {
            MPI_Gather(&stores_node_level,  pro->dataSize * pro->node_size, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, pro->Node_leaders_comm);
        }

        if(pro->real_Rack_Leader){
            //Stating World Level Collection
            int total_racks;
            MPI_Comm_size(pro->Rack_leaders_comm, &total_racks);
            int allRacksData[pro->dataSize * pro->world_size * 8];
            if(pro->rank_in_world == 0){
                int t = 0;
                int willRecv;
                for(int i = 1; i < total_racks - 1; i++){
                  //First the rank0 will ask to each Rack leader that how much he is going to send then he accepts the data.
                    MPI_Recv(&willRecv, 1, MPI_INT, i, 0, pro->Rack_leaders_comm, MPI_STATUS_IGNORE);
                    MPI_Recv(&allRacksData[t], willRecv, MPI_DOUBLE, i, 45, pro->Rack_leaders_comm, MPI_STATUS_IGNORE);
                    t += willRecv;
                }
                //Sending the complete data to root process from Rank0 process
                MPI_Send(allRacksData, t, MPI_DOUBLE, pro->root, 0, MPI_COMM_WORLD);
                
            }
            else
            {
                int willsend = pro->dataSize * pro->rack_size;
                //Rack Leaders first informs the Rank0 that how much they are going to send, after informing they send the actual data.
                MPI_Send(&willsend, 1, MPI_INT, 0, 0, pro->Rack_leaders_comm);
                MPI_Send(stores_Rack_level, willsend, MPI_DOUBLE, 0, 45, pro->Rack_leaders_comm);
            } 
        }
    }
    //Root process Receiving the complete data from Rank 0 Process.
    if(pro->rank_in_world == pro->root){
        int allRacksData[pro->dataSize * pro->world_size * 8];
        MPI_Recv(allRacksData,pro->dataSize * pro->world_size * 8 , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
       
}  




void Opti_AllToAllv(processData *pro){   
   
    //Code Contains Bug
    // Define the buffer containing the values to send
    /* x
    double buffer_send[dataSize];
    int buffer_send_length = dataSize;
    for(int i = 0; i < dataSize; i++){
        buffer_send[i] = my_rank*10 + i; 
    }
 
    // Define my counts for sending (how many integers do I send to each process?)
    int send_matrix[pro->node_size][pro->node_size];
    if(my_rank == 0){
        for(int i=0; i<pro->node_size ; i++){
            int temp = 0;
            for(int j = 0; j < pro->node_size; j++){
                send_matrix[i][j] = rand()%(dataSize-temp+1);
                temp += send_matrix[i][j];
            }
        }
    }

    MPI_Bcast(&send_matrix, size * size, MPI_INT, 0, MPI_COMM_WORLD);
    int disp_send[size];

    int temp = 0;
    for(int i=0; i<size ; i++){
        if(send_matrix[my_rank][i] == 0)
            disp_send[i] = -1;
        else{
            
            disp_send[i] = temp;
            temp += send_matrix[my_rank][i];
        } 
    }

    int counts_send[size], counts_recv[size];
    temp = 0;
    for(int j = 0; j < size; j++){
        counts_send[j] = send_matrix[my_rank][j];
        counts_recv[j] = send_matrix[j][my_rank];
    }

    int disp_recv[size];

    temp = 0;
    for(int i=0; i<size ; i++){
        if(send_matrix[i][my_rank] == 0)
            disp_recv[i] = -1;
        else{
            disp_recv[i] = temp;
            temp += send_matrix[i][my_rank];
        } 
    }
    int total_data_to_recev = 0;
    for(int j = 0; j < size; j++){
        total_data_to_recev += send_matrix[my_rank][j];
    }
    double count_send[size];
    for(int j = 0; j < size; j++){
        count_send[j] = send_matrix[my_rank][j];
    }

    double recv_buffer[total_data_to_recev];
    
    MPI_Alltoallv(buffer_send, counts_send, disp_send, MPI_DOUBLE, recv_buffer, counts_recv, disp_recv, MPI_DOUBLE, MPI_COMM_WORLD);
    
    printf("Done: %d\n", my_rank);
    */

}


//In this Function first we create each Differnt Communicators and then we call each opertion according to need
void Optimized_Version(int dataSize, int Operation){
    MPI_Comm Rack_comm;                 
    MPI_Comm Node_comm;
    MPI_Comm Node_leaders_comm;
    MPI_Comm Rack_leaders_comm;

    int root;
    int root_node = 0;
    int root_rack = 0;

    int Rack_colour = 0;
    int Node_colour = 0;

    int Node_Leader_colour = 0;     //1 for processess Whose Rank is 0 in Node_comm
    int Rack_Leader_colour = 0;     //1 for processess Whose Rank is 0 in Rack_comm

    int world_size, rank_in_world;
    int rack_size, rank_in_rack;
    int node_size, rank_in_node;

    int real_Node_Leader = 0;
    int real_Rack_Leader = 0;

    char* RackName;
    int ProcessorNum;
    int len;                    //length of processor name

    //Getting Node name eg. CSEWS35
    char hostname[MPI_MAX_PROCESSOR_NAME];

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_in_world);
    if(rank_in_world == 0){
        root = rand()%world_size;
    }
    MPI_Bcast(&root, 1, MPI_INT, 0, MPI_COMM_WORLD);            //Broadcasting Root rank to all other rank
    
    MPI_Get_processor_name(hostname, &len);
    
    //Get Processor number
    ProcessorNum = char2int(hostname, len);
 
    // Determine the Rack_colour.
    Rack_colour = giveColor(ProcessorNum);
    int key;
    key = rank_in_world;
    if(Rack_colour == 1){
        RackName = "Rack1";
    }
    else if(Rack_colour == 2){
        RackName = "Rack2";
    }
    else if(Rack_colour == 3){
        RackName = "Rack3";
    }
    else if(Rack_colour == 4){
        RackName = "Rack4";
    }
    else if(Rack_colour == 5){
        RackName = "Rack5";
    }
    else if(Rack_colour == 6){
        RackName = "Rack6";
    }
    else{
        RackName = "Rack0";
    }
    
    // Split de global communicator make new Communicator of every rack
    MPI_Comm_split(MPI_COMM_WORLD, Rack_colour, key, &Rack_comm);
    
    MPI_Comm_rank(Rack_comm, &rank_in_rack);
    
    MPI_Comm_size(Rack_comm, &rack_size);

    if(rank_in_world == root){
        root_node = ProcessorNum;
        root_rack = Rack_colour;
    }
    //Broadcasting Root rack and root node to all other rank
    MPI_Bcast(&root_rack, 1, MPI_INT, root, MPI_COMM_WORLD);            
    MPI_Bcast(&root_node, 1, MPI_INT, root, MPI_COMM_WORLD);

//------------------------------------------------------------------------------------------------------------------------------------------
                //Starting creating Node Communicators for each node

    Node_colour = ProcessorNum;
    MPI_Comm_split(MPI_COMM_WORLD, ProcessorNum, key, &Node_comm);
    MPI_Comm_size(Node_comm, &node_size);
    MPI_Comm_rank(Node_comm, &rank_in_node);

//--------------------------------------------------------------------------------------------------------------------------------------------
            //Connecting Different Node Leaders and Rack Leaders
    if(rank_in_node == 0){                  //Will create circle having process which have 0 rank in node and are in same rack
        Node_Leader_colour = Rack_colour;
        real_Node_Leader = 1;
    }   
    if(rank_in_rack == 0){               //Will create circle having process which have 0 rank in node and are in same rack
        Rack_Leader_colour = 1;
        real_Rack_Leader = 1;
    }

    //Creating Node Leaders and Rack Leaders
    MPI_Comm_split(MPI_COMM_WORLD, Node_Leader_colour, key, &Node_leaders_comm);
    MPI_Comm_split(MPI_COMM_WORLD, Rack_Leader_colour, key, &Rack_leaders_comm);

    //Creating object to send data to Function easily.
    processData pro;
    pro.Rack_comm = Rack_comm;
    pro.Node_comm = Node_comm;
    pro.Node_leaders_comm = Node_leaders_comm;
    pro.Rack_leaders_comm = Rack_leaders_comm;
    pro.Rack_colour = Rack_colour;
    pro.Node_colour = Node_colour;
    pro.Node_Leader_colour = Node_Leader_colour;  
    pro.Rack_Leader_colour = Rack_Leader_colour;     
    pro.world_size = world_size;
    pro.rank_in_world = rank_in_world;;
    pro.rack_size = rack_size;
    pro.rank_in_rack = rank_in_rack;
    pro.node_size = node_size;
    pro.rank_in_node = rank_in_node;
    pro.real_Node_Leader = real_Node_Leader;
    pro.real_Rack_Leader = real_Rack_Leader;
    pro.ProcessorNum = ProcessorNum;
    pro.dataSize = dataSize;
    pro.root = root;

    //Finding p and ppn from World_comm size
    int p, ppn;
    if(world_size == 4){
      p = 4; ppn = 1;
    }
    else if(world_size == 16){
      p = 16; ppn = 1;
    }
    else if(world_size == 32){
      p = 4; ppn = 8;
    }
    else{
      p = 16; ppn = 8;
    }
    clock_t t;
    t = clock();

    FILE *fr;
    if(Operation == 1){                 //Bcasting
        for(int i = 0; i < 5; i++)
            Opti_Bcast(&pro);
        MPI_Barrier(MPI_COMM_WORLD);
        t = clock() - t;
        fr = fopen("Data_Bcast.txt", "a");
    }
    else if(Operation == 2){            //Reducing
        for(int i = 0; i < 5; i++)
            Opti_Reduce(&pro);
        MPI_Barrier(MPI_COMM_WORLD);
        t = clock() - t;
        fr = fopen("Data_Reduce.txt", "a");
    }
    else if(Operation == 3){            //Gathering
        for(int i = 0; i < 1; i++)
            Opti_Gather(&pro, ppn);
        MPI_Barrier(MPI_COMM_WORLD);
        t = clock() - t;
        fr = fopen("Data_Gather.txt", "a");
    }
    else if(Operation == 4){            //Alltoallv
        for(int i = 0; i < 1; i++)
            Opti_AllToAllv(&pro);
        MPI_Barrier(MPI_COMM_WORLD);
        t = clock() - t;
        fr = fopen("Data_Alltoallv.txt", "a");
    }

    //Writing in File
    if(rank_in_world == 0){
        int mode = 1;
        double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
        fprintf(fr, "%d, %d, %d, Optimized, %lf \n",dataSize, p, ppn, time_taken/5);
        printf("\t\tOptimal Version took %f seconds to execute %d \n", time_taken, dataSize);
    }
    fclose(fr);
}


//Normal Bcast
void Bcast(int dataSize){
    int my_rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int broadcast_root;
    if(my_rank == 0)
        broadcast_root =  rand()%size;
    MPI_Bcast(&broadcast_root, 1, MPI_INT,0, MPI_COMM_WORLD);

    double buffer[dataSize];
    for(int i = 0; i < dataSize; i++){
      buffer[i] = rand();
    }
    MPI_Bcast(&buffer, dataSize, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
}


//Normal Reduce
void Reduce(int dataSize){
    int size = 0, my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double buffer[dataSize], result[dataSize];
    int root_rank;  
    if(my_rank == 0)
        root_rank =  rand()%size;
    MPI_Bcast(&root_rank, 1, MPI_INT,0, MPI_COMM_WORLD);
      
    for(int i = 0; i < dataSize; i++){
      buffer[i] = i;
    }
    // Each MPI process sends its rank to reduction, root MPI process collects the result
    MPI_Reduce(buffer, result, dataSize, MPI_DOUBLE, MPI_SUM, root_rank, MPI_COMM_WORLD);    
}


//Normal Gather
void Gather(int dataSize){
    int size, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int root_rank;  
    if(my_rank == 0)
        root_rank =  rand()%size;
    MPI_Bcast(&root_rank, 1, MPI_INT,0, MPI_COMM_WORLD);
 
    double my_value[dataSize];
    for(int i = 0; i < dataSize; i++){
      my_value[i] = my_rank; 
    }

    if(my_rank == root_rank)
    {
        double buffer[dataSize * size];
        MPI_Gather(&my_value, dataSize, MPI_DOUBLE, buffer, dataSize, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Gather(&my_value, dataSize, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
    }
}


//Normal alltoallv
void AllToAllv(int dataSize){
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // Get my rank
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // Define the buffer containing the values to send
    double buffer_send[dataSize];
    int buffer_send_length = dataSize;
    for(int i = 0; i < dataSize; i++){
        buffer_send[i] = my_rank*10 + i; 
    }
 
    // Define my counts for sending (how many integers do I send to each process?)
    int send_matrix[size][size];
    if(my_rank == 0){
        for(int i=0; i<size ; i++){
            int temp = 0;
            for(int j = 0; j < size; j++){
                send_matrix[i][j] = rand()%(dataSize-temp+1);
                temp += send_matrix[i][j];
            }
        }
    }
    //Bcasting Sending and Receiving data of each process to each process 
    MPI_Bcast(&send_matrix, size * size, MPI_INT, 0, MPI_COMM_WORLD);
    int disp_send[size];

    int temp = 0;
    for(int i=0; i<size ; i++){
        if(send_matrix[my_rank][i] == 0)
            disp_send[i] = -1;
        else{
            
            disp_send[i] = temp;
            temp += send_matrix[my_rank][i];
        } 
    }

    int counts_send[size], counts_recv[size];
    temp = 0;
    for(int j = 0; j < size; j++){
        counts_send[j] = send_matrix[my_rank][j];
        counts_recv[j] = send_matrix[j][my_rank];
    }
    int disp_recv[size];
    temp = 0;
    for(int i=0; i<size ; i++){
        if(send_matrix[i][my_rank] == 0)
            disp_recv[i] = -1;
        else{
            disp_recv[i] = temp;
            temp += send_matrix[i][my_rank];
        } 
    }
    int total_data_to_recev = 0;
    for(int j = 0; j < size; j++){
        total_data_to_recev += send_matrix[my_rank][j];
    }
    double count_send[size];
    for(int j = 0; j < size; j++){
        count_send[j] = send_matrix[my_rank][j];
    }

    double recv_buffer[total_data_to_recev];
    
    MPI_Alltoallv(buffer_send, counts_send, disp_send, MPI_DOUBLE, recv_buffer, counts_recv, disp_recv, MPI_DOUBLE, MPI_COMM_WORLD);
 }



int main(int argc, char* argv[])
  {
    srand(time(0));
    int dataSize = atoi(argv[1]);
    int Operation = atoi(argv[2]);

    MPI_Init(&argc, &argv);
    
    int r, siz, p, ppn;
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    MPI_Comm_size(MPI_COMM_WORLD, &siz);

    if(siz == 4){
      p = 4; ppn = 1;
    }
    else if(siz == 16){
      p = 16; ppn = 1;
    }
    else if(siz == 32){
      p = 4; ppn = 8;
    }
    else{
      p = 16; ppn = 8;
    }

    clock_t t;
    t = clock();

    FILE *fr;
    if(Operation == 1){
      for(int i = 0; i < 5; i++)
        Bcast(dataSize);
      MPI_Barrier(MPI_COMM_WORLD);
      t = clock() - t;
      fr = fopen("Data_Bcast.txt", "a");
    }

    if(Operation == 2){
      for(int i = 0; i < 5; i++)
        Reduce(dataSize);
      MPI_Barrier(MPI_COMM_WORLD);
      t = clock() - t;
      fr = fopen("Data_Reduce.txt", "a");
    }
    
    if(Operation == 3){
      for(int i = 0; i < 5; i++)
        Gather(dataSize);
      MPI_Barrier(MPI_COMM_WORLD);
      t = clock() - t;
      fr = fopen("Data_Gather.txt", "a");
    }

    if(Operation == 4){
      for(int i = 0; i < 5; i++)
        AllToAllv(dataSize);
      MPI_Barrier(MPI_COMM_WORLD);
      t = clock() - t;
      fr = fopen("Data_Alltoallv.txt", "a");
    }

    if(r == 0){
        int mode = 0;
        double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
        fprintf(fr, "%d, %d, %d, Standard, %lf \n",dataSize, p, ppn, time_taken/5);
        printf("\t\tNormal Version took %f seconds to execute %d \n", time_taken, dataSize);
    }
    fclose(fr);
    
    printf("\nStarting Optimized with :  %d, %d, %d \n",dataSize, p, ppn);
    Optimized_Version(dataSize, Operation);
    
    
    MPI_Finalize();

    return EXIT_SUCCESS;    
}

