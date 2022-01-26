#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>

typedef struct data{
    int myrank;
    int size;
    int row;
    int col;
    double* arr;
    double* newarr;
    double* left;
    double* right;
    double* top;
    double* bottom;
    int sq;
}processData;

int findSqrt(int x) 
{  
    if (x == 0 || x == 1) 
    return x; 

    int i = 1, r = 1; 
    while (r <= x) 
    { 
      i++; 
      r = i * i; 
    } 
    return i - 1; 
} 


void PACK_SendRightReceiveLeft(processData *pro){
    MPI_Status stat;
    if(pro->myrank%pro->sq == 0){  //Left side processess, will only send
        double buffer[pro->col];
        int position = 0;
        for(int i = 0; i < pro->col ; i++){
           MPI_Pack( &pro->arr[i * pro->row + pro->col - 1], 1, MPI_DOUBLE, buffer, sizeof(double) * pro->row, &position, MPI_COMM_WORLD);
        }
        MPI_Send(buffer, position, MPI_PACKED, pro->myrank + 1, pro->myrank, MPI_COMM_WORLD );
    }

    else if((pro->myrank+1)%pro->sq == 0){  //Right
        int position = 0;
        double buff[pro->col];
        MPI_Recv(buff, pro->row * sizeof(double), MPI_PACKED, pro->myrank - 1, pro->myrank - 1, MPI_COMM_WORLD, &stat);
        for(int i = 0; i < pro->col; i++){
            MPI_Unpack(buff, pro->row, &position, &pro->left[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
        }
        
    }

    else{      //Middle
       double buffer[pro->col], buff[pro->col];;
        int position = 0;
        for(int i = 0; i < pro->col ; i++){
           MPI_Pack( &pro->arr[i * pro->row + pro->col - 1], 1, MPI_DOUBLE, buffer, sizeof(double) * pro->row, &position, MPI_COMM_WORLD);
        }
        MPI_Sendrecv(buffer, position, MPI_PACKED, pro->myrank + 1, pro->myrank, 
                buff, pro->row * sizeof(double), MPI_PACKED, pro->myrank - 1, pro->myrank - 1, MPI_COMM_WORLD, &stat);
       

        position = 0;
        
        for(int i = 0; i < pro->col; i++){
            MPI_Unpack(buff, pro->row, &position, &pro->left[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
        }
        
        
    }
}


void PACK_SendLeftReceiveRight(processData *pro){
    MPI_Status stat;
    if(pro->myrank%pro->sq == 0){  //Left
        int position = 0;
        double buff[pro->row];
        MPI_Recv(buff, pro->row * sizeof(double), MPI_PACKED, pro->myrank + 1, pro->myrank + 1, MPI_COMM_WORLD, &stat);
        for(int i = 0; i < pro->col; i++){
            MPI_Unpack(buff, pro->row, &position, &pro->right[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
        }
        
    }

    else if((pro->myrank+1)%pro->sq == 0){  //Right
        double buffer[pro->row];
        int position = 0;
        for(int i = 0; i < pro->row ; i++){
           MPI_Pack( &pro->arr[i*pro->row + 0], 1, MPI_DOUBLE, buffer, sizeof(double) * pro->row, &position, MPI_COMM_WORLD);
        }
        MPI_Send(buffer, position, MPI_PACKED, pro->myrank - 1, pro->myrank, MPI_COMM_WORLD );
    }

    else{      //Middle

        double buffer[pro->row], buff[pro->row];;
        int position = 0;
        for(int i = 0; i < pro->row ; i++){
           MPI_Pack( &pro->arr[i*pro->row + 0], 1, MPI_DOUBLE, buffer, sizeof(double) * pro->row, &position, MPI_COMM_WORLD);
        }
        MPI_Sendrecv(buffer, position, MPI_PACKED, pro->myrank - 1, pro->myrank,
            buff, pro->row * sizeof(double), MPI_PACKED, pro->myrank + 1, pro->myrank + 1, MPI_COMM_WORLD, &stat);

        position = 0;
        for(int i = 0; i < pro->col; i++){
            MPI_Unpack(buff, pro->row, &position, &pro->right[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
        }
        

    }
}



void PACK_SendTopReceiveBottom(processData *pro){                   //Send Top Receive Bottom
    MPI_Status stat; 
    if(pro->myrank < pro->sq){     //Top Process
        int position = 0;
        double buff[pro->col];
        MPI_Recv(buff, pro->row * sizeof(double), MPI_PACKED, pro->myrank + pro->sq , pro->myrank + pro->sq , MPI_COMM_WORLD, &stat);
        for(int i = 0; i < pro->col; i++){
            MPI_Unpack(buff, pro->row, &position, &pro->bottom[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
        }
    }

    else if(pro->myrank >= pro->sq*(pro->sq-1)){ //Bottom Process
        double buffer[pro->col];
        int position = 0;
        for(int i = 0; i < pro->col ; i++){
           MPI_Pack( &pro->arr[0 + i], 1, MPI_DOUBLE, buffer, sizeof(double) * pro->col, &position, MPI_COMM_WORLD);
        }
        MPI_Send(buffer, position, MPI_PACKED, pro->myrank - pro->sq, pro->myrank, MPI_COMM_WORLD );
    }

    else{   //Middle Process
        double buffer[pro->col], buff[pro->col];
        int position = 0;
        for(int i = 0; i < pro->col ; i++){
           MPI_Pack( &pro->arr[0+i], 1, MPI_DOUBLE, buffer, sizeof(double) * pro->col, &position, MPI_COMM_WORLD);
        }

        MPI_Sendrecv(buffer, position, MPI_PACKED, pro->myrank - pro->sq, pro->myrank, 
                    buff, pro->row * sizeof(double), MPI_PACKED, pro->myrank + pro->sq , pro->myrank + pro->sq , MPI_COMM_WORLD, &stat);

        position = 0;
        for(int i = 0; i < pro->col; i++){
            MPI_Unpack(buff, pro->row, &position, &pro->bottom[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
        }
        
    }
}


void PACK_SendBottomReceiveTop(processData *pro){              //Send Bottom Receive Top
    MPI_Status stat;
    if(pro->myrank < pro->sq){             //Top Process
        double buffer[pro->col];
        int position = 0;
        for(int i = 0; i < pro->col ; i++){
           MPI_Pack( &pro->arr[pro->row- 1 * pro->row + i], 1, MPI_DOUBLE, buffer, sizeof(double) * pro->col, &position, MPI_COMM_WORLD);
        }
        MPI_Send(buffer, position, MPI_PACKED, pro->myrank + pro->sq, pro->myrank, MPI_COMM_WORLD );
    }

    else if(pro->myrank >= pro->sq*(pro->sq-1)){ //Bottom Process
        int position = 0;
        double buff[pro->col];
        MPI_Recv(buff, pro->row * sizeof(double), MPI_PACKED, pro->myrank - pro->sq , pro->myrank - pro->sq , MPI_COMM_WORLD, &stat);
        for(int i = 0; i < pro->col; i++){
            MPI_Unpack(buff, pro->row, &position, &pro->top[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
        }
    }

    else{   //Middle Process
        double buffer[pro->col], buff[pro->col];
        int position = 0;
        for(int i = 0; i < pro->col ; i++){
           MPI_Pack( &pro->arr[pro->row- 1 * pro->row + i], 1, MPI_DOUBLE, buffer, sizeof(double) * pro->col, &position, MPI_COMM_WORLD);
        }
        MPI_Sendrecv(buffer, position, MPI_PACKED, pro->myrank + pro->sq, pro->myrank,
                buff, pro->row * sizeof(double), MPI_PACKED, pro->myrank - pro->sq , pro->myrank - pro->sq , MPI_COMM_WORLD, &stat);

        position = 0;
        for(int i = 0; i < pro->col; i++){
            MPI_Unpack(buff, pro->row, &position, &pro->top[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
        }
        

    }
}


int TagGenerator(int r, int c){         //Tag Generator Because there will be multiple transaction between same processess in Normal Mode
	int res;                            //It Stores Row and Column value in one integer
	res = r;
	res <<= 8;
	res += c;
	return res;
} 


void NORMAL_SendRightReceiveLeft(processData *pro){
    MPI_Status stat;    
    if(pro->myrank%pro->sq == 0){  //Left Side Processess
        for(int i = 0; i < pro->row; i++){
            double send = pro->arr[i* pro->row + pro->col-1];
            int tag_send = TagGenerator(i,0);
            MPI_Send(&send, 1, MPI_DOUBLE, pro->myrank + 1, tag_send, MPI_COMM_WORLD);           
        }           
    }
    else if((pro->myrank+1)%pro->sq == 0){  //Right Side Processess
        for(int i = 0; i < pro->row ; i++){
            int tag_recv = TagGenerator(i, 0);
            MPI_Recv(&pro->left[i], 1, MPI_DOUBLE, pro->myrank - 1, tag_recv, MPI_COMM_WORLD, &stat);            
        }
    }
    else{      //Middle
        pro->left = (double *)malloc(pro->col * sizeof(double));
        for(int i = 0; i < pro->row ; i++){
            double send = pro->arr[i* pro->row + 0];
            int tag_send = TagGenerator(i, 0);            //send to
            int tag_recv = TagGenerator(i, 0);
            MPI_Sendrecv(&send, 1, MPI_DOUBLE, pro->myrank + 1, tag_send, 
                &pro->left[i], 1, MPI_DOUBLE, pro->myrank - 1 , tag_recv, MPI_COMM_WORLD, &stat);
        }
    }
}

void NORMAL_SendLeftReceiveRight(processData *pro){
    MPI_Status stat;
    if(pro->myrank%pro->sq == 0){  //Left side Processess
        for(int i = 0; i < pro->row ; i++){
            int tag_recv = TagGenerator(i, pro->col-1);
            MPI_Recv(&pro->right[i], 1, MPI_DOUBLE, pro->myrank + 1, tag_recv, MPI_COMM_WORLD, &stat);
        }    
    }
    else if((pro->myrank+1)%pro->sq == 0){  //Right side Processess
        for(int i = 0; i < pro->row; i++){
            double send = pro->arr[i* pro->row + 0];
            int tag_send = TagGenerator(i,pro->col-1);
            MPI_Send(&send, 1, MPI_DOUBLE, pro->myrank - 1, tag_send, MPI_COMM_WORLD);           
        }
    }
    else{      //Middle
        pro->right = (double *)malloc(pro->col * sizeof(double));
        for(int i = 0; i < pro->row ; i++){
            double send = pro->arr[i* pro->row + 0];
            int tag_send = TagGenerator(i, pro->col-1);           
            int tag_recv = TagGenerator(i, pro->col-1);
            MPI_Sendrecv(&send, 1, MPI_DOUBLE, pro->myrank - 1, tag_send, 
                &pro->right[i], 1, MPI_DOUBLE, pro->myrank + 1 , tag_recv, MPI_COMM_WORLD, &stat);
            
        }
    }
}



void NORMAL_SendTopReceiveBottom(processData *pro){                   //Send Top Receive Bottom
    MPI_Status stat[pro->size]; 
    
    if(pro->myrank < pro->sq){     //Top Process
        
        for(int i = 0; i < pro->col ; i++){
            int tag_recv = TagGenerator(pro->row-1, i);
            MPI_Recv(&pro->bottom[i], 1, MPI_DOUBLE, pro->myrank + pro->sq, tag_recv, MPI_COMM_WORLD, &stat[pro->myrank]);
        }
    }
    else if(pro->myrank >= pro->sq*(pro->sq-1)){ //Bottom Process
        for(int i = 0; i < pro->col; i++){
            double send = pro->arr[0* pro->row + i];
            int tag_send = TagGenerator(pro->row-1,i);
            MPI_Send(&send, 1, MPI_DOUBLE, pro->myrank - pro->sq, tag_send, MPI_COMM_WORLD);
            
        }
    }
    else{   //Middle Process
        
        for(int i = 0; i < pro->col ; i++){
            double send = pro->arr[0* pro->row + i];
            int tag_send = TagGenerator(pro->row-1,i);              //send to
            int tag_recv = TagGenerator(pro->row-1, i);
            MPI_Sendrecv(&send, 1, MPI_DOUBLE, pro->myrank - pro->sq, tag_send, 
                &pro->bottom[i], 1, MPI_DOUBLE, pro->myrank + pro->sq , tag_recv, MPI_COMM_WORLD, &stat[pro->myrank]);
        }
    }
}


void NORMAL_SendBottomReceiveTop(processData *pro){              //Send Bottom Receive Top
    MPI_Status stat[pro->size];
  
    if(pro->myrank < pro->sq){                            //Top Process
        for(int i = 0; i < pro->col; i++){
            double send = pro->arr[pro->row-1* pro->row + i];
            int tag_send = TagGenerator(0,i);
            MPI_Send(&send, 1, MPI_DOUBLE, pro->myrank + pro->sq, tag_send, MPI_COMM_WORLD);
        }
    }
    else if(pro->myrank >= pro->sq*(pro->sq-1)){ //Bottom Process
        
        for(int i = 0; i < pro->col ; i++){
            int tag_recv = TagGenerator(0, i);
            MPI_Recv(&pro->top[i], 1, MPI_DOUBLE, pro->myrank - pro->sq, tag_recv, MPI_COMM_WORLD, &stat[pro->myrank]);
        }
    }
    else{   //Middle Process
        
        for(int i = 0; i < pro->col ; i++){
            double send = pro->arr[pro->row-1* pro->row + i];
            int tag_send = TagGenerator(0,i);              
            int tag_recv = TagGenerator(0, i);
            MPI_Sendrecv(&send, 1, MPI_DOUBLE, pro->myrank + pro->sq, tag_send, 
                &pro->top[i], 1, MPI_DOUBLE, pro->myrank - pro->sq , tag_recv, MPI_COMM_WORLD, &stat[pro->myrank]);
        }
    }
}


void VECTOR_SendRightReceiveLeft(processData *pro){
    MPI_Status stat;
    MPI_Datatype sendvac, recvvac;
    MPI_Type_vector(pro->row, 1, pro->row-1, MPI_DOUBLE, &sendvac);
    MPI_Type_vector(pro->row, 1, 1, MPI_DOUBLE, &recvvac);
    MPI_Type_commit(&sendvac);
    MPI_Type_commit(&recvvac);

    if(pro->myrank%pro->sq == 0){  //Left Side Processess ,so only send
        int tag = pro->myrank;
        MPI_Send(&pro->arr[pro->row-1], 1, sendvac ,pro->myrank + 1, tag , MPI_COMM_WORLD);                   
    }
    else if((pro->myrank+1)%pro->sq == 0){  //Right Side Processess
        int tag = pro->myrank-1;
        MPI_Recv(&pro->left[0], 1, recvvac, pro->myrank - 1, tag, MPI_COMM_WORLD, &stat);
    }
    else{      //Middle
        int tag = pro->myrank;
        int tag1 = pro->myrank-1;
        MPI_Sendrecv(&pro->arr[pro->row-1], 1, sendvac, pro->myrank + 1, tag ,&pro->left[0], 1, recvvac, pro->myrank - 1, tag1, MPI_COMM_WORLD, &stat);
    }
}

void VECTOR_SendLeftReceiveRight(processData *pro){
    MPI_Status stat;
    MPI_Datatype sendvac, recvvac;
    MPI_Type_vector(pro->row, 1, pro->row-1, MPI_DOUBLE, &sendvac);
    MPI_Type_vector(pro->row, 1, 1, MPI_DOUBLE, &recvvac);
    MPI_Type_commit(&sendvac);
    MPI_Type_commit(&recvvac);
    if(pro->myrank%pro->sq == 0){  //Left
        int tag = pro->myrank+1;
        MPI_Recv(&pro->right[0], 1, recvvac, pro->myrank + 1, tag, MPI_COMM_WORLD, &stat);
    }
    else if((pro->myrank+1)%pro->sq == 0){  //Right
        int tag = pro->myrank;
        MPI_Send(&pro->arr[0], 1, sendvac ,pro->myrank - 1, tag , MPI_COMM_WORLD);   
    }
    else{      //Middle
        int tag = pro->myrank;
        int tag2 = pro->myrank+1;
        MPI_Sendrecv(&pro->arr[0], 1, sendvac ,pro->myrank - 1, tag , &pro->right[0], 1, recvvac, pro->myrank + 1, tag2, MPI_COMM_WORLD, &stat);
    }
}

void VECTOR_SendTopReceiveBottom(processData *pro){                   //Send Top Receive Bottom
    MPI_Status stat;
    MPI_Datatype sendvac, recvvac;
    MPI_Type_vector(pro->row, 1, 1, MPI_DOUBLE, &sendvac);
    MPI_Type_vector(pro->row, 1, 1, MPI_DOUBLE, &recvvac);
    MPI_Type_commit(&sendvac);
    MPI_Type_commit(&recvvac);
    if(pro->myrank < pro->sq){     //Top Process
        int tag = pro->myrank + pro->sq;
        MPI_Recv(&pro->bottom[0], 1, recvvac, pro->myrank + pro->sq , tag, MPI_COMM_WORLD, &stat);          
    }
    else if(pro->myrank >= pro->sq*(pro->sq-1)){ //Bottom Process
        int tag = pro->myrank;
        MPI_Send(&pro->arr[0], 1, sendvac ,pro->myrank - pro->sq, tag , MPI_COMM_WORLD);   
    }
    else{   //Middle Process
        int tag = pro->myrank;
        int tag2 = pro->myrank + pro->sq;
        MPI_Sendrecv(&pro->arr[0], 1, sendvac ,pro->myrank - pro->sq, tag , &pro->bottom[0], 1, recvvac, pro->myrank + pro->sq , tag2, MPI_COMM_WORLD, &stat);
    }
}


void VECTOR_SendBottomReceiveTop(processData *pro){              //Send Bottom Receive Top
    MPI_Status stat;
    MPI_Datatype sendvac, recvvac;
    MPI_Type_vector(pro->row, 1, 1, MPI_DOUBLE, &sendvac);
    MPI_Type_vector(pro->row, 1, 1, MPI_DOUBLE, &recvvac);
    MPI_Type_commit(&sendvac);
    MPI_Type_commit(&recvvac);
    if(pro->myrank < pro->sq){                            //Top Process
        int tag = pro->myrank;
        MPI_Send(&pro->arr[pro->row * (pro->row - 1)], 1, sendvac ,pro->myrank + pro->sq, tag , MPI_COMM_WORLD);   
    }
    else if(pro->myrank >= pro->sq*(pro->sq-1)){ //Bottom Process
        int tag = pro->myrank - pro->sq;
        MPI_Recv(&pro->top[0], 1, recvvac, pro->myrank - pro->sq , tag, MPI_COMM_WORLD, &stat);
    }
    else{   //Middle Process
        int tag = pro->myrank;
        int tag2 = pro->myrank - pro->sq;
        MPI_Sendrecv(&pro->arr[pro->row * (pro->row - 1)], 1, sendvac ,pro->myrank + pro->sq, tag , &pro->top[0], 1, recvvac, pro->myrank - pro->sq , tag2, MPI_COMM_WORLD, &stat);
    }
}

//Corner Case of Corner Processess
void TopLeftCorner(processData *pro, int i, int j){
    pro->newarr[0] = (pro->arr[0*pro->row + 1]+pro->arr[1 * pro->row + 0])/2;
}

void TopRightCorner(processData *pro, int i, int j){
    pro->newarr[0 + pro->col - 1] = (pro->arr[pro->col-2]+pro->arr[1 * pro->row + pro->col - 1])/2;
}

void BottomLeftCorner(processData *pro, int i, int j){
    pro->newarr[(pro->row - 1) * pro->row +0] = (pro->arr[(pro->row -1) * pro->row + 1] + pro->arr[(pro->row-2) * pro->row + 1])/2;
}

void BottomRightCorner(processData *pro, int i, int j){
    pro->newarr[(pro->row-1) * pro->row + (pro->col-1)] = (pro->arr[(pro->row-2) * pro->row + (pro->col-1)] + pro->arr[(pro->row-1) * pro->row + (pro->col-2)])/2;
}



//Elements in boundary But are adjacent to different processess
void TopInnerRight(processData *pro, int i, int j){
    pro->newarr[(i) * pro->row + (j)] = (pro->arr[(i+1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j-1)] + pro->right[i])/3;
}

void BottomInnerRight(processData *pro, int i, int j){
    pro->newarr[(i) * pro->row + (j)] = (pro->arr[(i-1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j-1)] + pro->right[i])/3;
}

void TopInnerLeft(processData *pro, int i, int j){
    pro->newarr[(i) * pro->row + (j)] = (pro->arr[(i+1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j+1)] + pro->left[i])/3;
}

void BottomInnerLeft(processData *pro, int i, int j){
    pro->newarr[(i) * pro->row + (j)] = (pro->arr[(i-1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j+1)] + pro->left[i])/3;
}

void TopLeft(processData *pro, int i, int j){
    pro->newarr[(i) * pro->row + (j)] = (pro->arr[(i+1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j+1)] + pro->top[j])/3;
}

void BottomLeft(processData *pro, int i, int j){
    pro->newarr[(i) * pro->row + (j)] = (pro->arr[(i-1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j+1)] + pro->bottom[j])/3;
}

void TopRight(processData *pro, int i, int j){
    pro->newarr[(i) * pro->row + (j)] = (pro->arr[(i+1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j-1)] + pro->top[j])/3;
}

void BottomRight(processData *pro, int i, int j){
    pro->newarr[(i) * pro->row + (j)] = (pro->arr[(i-1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j-1)] + pro->bottom[j])/3;
}



//Boundary Middle Elements of processess
void LeftBoundary(processData *pro, int i, int j){
    pro->newarr[(i) * pro->row + (j)] = (pro->arr[(i-1) * pro->row + (j)] + pro->arr[(i+1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j+1)] )/3; 
}

void RightBoundary(processData *pro, int i, int j){
    pro->newarr[(i) * pro->row + (j)] = (pro->arr[(i-1) * pro->row + (j)] + pro->arr[(i+1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j-1)] )/3;
}

void TopBoundary(processData *pro, int i, int j){
    pro->newarr[(i) * pro->row + (j)] = (pro->arr[(i) * pro->row + (j)] + pro->arr[(i) * pro->row + (j)] + pro->arr[(i) * pro->row + (j)] )/3;
}

void BottomBoundary(processData *pro, int i, int j){
    pro->newarr[(i) * pro->row + (j)] = (pro->arr[(i) * pro->row + (j)] + pro->arr[(i) * pro->row + (j)] + pro->arr[(i) * pro->row + (j)] )/3;
}


//Inner Corners
void InnerBottomLeft(processData *pro, int i, int j){
    pro->newarr[(i)*pro->row+(j)]=(pro->arr[(i-1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j+1)] + pro->left[i] + pro->bottom[j])/4;
}

void InnerBottomRight(processData *pro, int i, int j){
    pro->newarr[(i)*pro->row+(j)]=(pro->arr[(i-1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j-1)] + pro->right[i] + pro->bottom[j])/4;
}

void InnerTopLeft(processData *pro, int i, int j){
    pro->newarr[(i)*pro->row+(j)]=(pro->arr[(i+1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j+1)] + pro->left[i] + pro->top[j])/4;
}

void InnerTopRight(processData *pro, int i, int j){
    pro->newarr[i * pro->row + j]=(pro->arr[(i+1) * pro->row + (j)] + pro->arr[i* pro->row + (j-1)] + pro->right[i] + pro->top[j])/4;
}

//Inner Boundaries
void InnerTop(processData *pro, int i, int j){
    pro->newarr[(i)*pro->row+(j)]=(pro->arr[(i+1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j-1)] + pro->arr[(i) * pro->row + (j+1)] + pro->top[j])/4;
}

void InnerBottom(processData *pro, int i, int j){
    pro->newarr[(i)*pro->row+(j)]=(pro->arr[(i-1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j+1)] + pro->arr[(i) * pro->row + (j-1)] + pro->bottom[j])/4;
}

void InnerLeft(processData *pro, int i, int j){
    pro->newarr[(i)*pro->row+(j)]=(pro->arr[(i+1) * pro->row + (j)] + pro->arr[(i-1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j+1)] + pro->left[i])/4;
}

void InnerRight(processData *pro, int i, int j){
    pro->newarr[(i)*pro->row+(j)]=(pro->arr[(i+1) * pro->row + (j)] + pro->arr[(i-1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j-1)] + pro->right[i])/4;
}

//INNER
void MiddleMiddle(processData *pro, int i, int j){
    pro->newarr[(i)*pro->row+(j)]=(pro->arr[(i+1) * pro->row + (j)] + pro->arr[(i-1) * pro->row + (j)] + pro->arr[(i) * pro->row + (j-1)] + pro->arr[(i) * pro->row + (j+1)])/4;
}


void LSP(processData *pro){    
    if(pro->myrank == 0){ //TopLeftCorner Process
     TopLeftCorner(pro, 0, 0);
     for(int i = 1; i < pro->col - 1; i++){
       TopBoundary(pro, 0, i);
     }
     TopInnerRight(pro, 0, pro->col-1);  
    }
    else{
      TopLeft(pro, 0, 0);
      for(int i = 1; i < pro->col - 1; i++){
        InnerTop(pro, 0, i);
      }
      InnerTopRight(pro, 0, pro->col - 1);
    }
    for(int i = 1; i < pro->row - 1 ; i++){
        for(int j = 1; j < pro->col ; j++){
            if(j == pro->col - 1) 
                InnerRight(pro, i, j);
            else
                MiddleMiddle(pro, i , j);
        }
    }
    for(int i = 1; i < pro->row - 1; i++){
      LeftBoundary(pro, i, 0);
    }  
    if(pro->myrank == pro->sq*(pro->sq-1)){ //BottomLeftProcess
      BottomLeftCorner(pro, pro->row-1 , 0);
      for(int i = 1; i < pro->col - 1; i++){
        BottomBoundary(pro, pro->row - 1, i);
      }
      BottomInnerRight(pro, pro->row-1, pro->col-1);
    }
    else{
      BottomLeft(pro, pro->row-1, 0);
      for(int i = 1; i < pro->col - 1; i++){
        InnerBottom(pro, pro->row - 1, i);
      }
      InnerBottomRight(pro, pro->row-1, pro->col-1);
    }
}


void RSP(processData *pro){
   if(pro->myrank == pro->sq - 1){   //TOP RIGHT CORNER PROCESS
    TopRightCorner(pro, 0 , pro->col - 1);
    for(int i = 1; i < pro->col - 1; i++){
       TopBoundary(pro, 0, i);
    }
    TopInnerLeft(pro, 0,0);
  }
  else{
    TopRight(pro, 0 , pro->col - 1);
    for(int i = 1; i < pro->col - 1; i++){
       InnerTop(pro, 0, i);
    }
    InnerTopLeft(pro, 0,0);
  }
  for(int i = 1; i < pro->row - 1 ; i++){
      for(int j = 0; j < pro->col - 1 ; j++){
          if(j == 0)
              InnerLeft(pro, i, j);
          else
              MiddleMiddle(pro, i , j);
      }
  }
  for(int i = 1; i < pro->row - 1; i++){
    RightBoundary(pro, i, pro->col-1);
  } 
  if(pro->myrank == pro->size - 1){  //Bottom Right Corner Preocess
    BottomRightCorner(pro, pro->row-1 , pro->col-1);
    for(int i = 1; i < pro->col - 1; i++){
      BottomBoundary(pro, pro->row - 1, i);
    }
    BottomInnerLeft(pro, pro->row-1, 0);
  }
  else{
    BottomRight(pro, pro->row-1, pro->col-1);
    for(int i = 1; i < pro->col - 1; i++){
      InnerBottom(pro, pro->row - 1, i);
    }
    InnerBottomLeft(pro, pro->row-1, pro->col-1);
  }
}

void TM(processData *pro){
    for(int i = 1; i < pro->row - 1 ; i++){
        for(int j = 1; j < pro->col - 1 ; j++){
                MiddleMiddle(pro, i, j);
        }
    }
    TopInnerLeft(pro, 0,0);
    TopInnerRight(pro, 0, pro->col-1);
    for(int i = 1; i < pro->col-1; i++){
      TopBoundary(pro, 0, i);
      InnerBottom(pro, pro->row-1, i); 
    }
    for(int i = 1; i < pro->row -1; i++){
      InnerLeft(pro, i, 0);
      InnerRight(pro, i, pro->col-1);
    }
    InnerBottomLeft(pro, pro->row - 1, 0);
    InnerBottomRight(pro, pro->row, pro->col);
}

void LM(processData *pro){
    for(int i = 1; i < pro->row - 1 ; i++){
        for(int j = 1; j < pro->col - 1 ; j++){
                MiddleMiddle(pro, i, j);
        }
    }
    InnerTopLeft(pro, 0,0);
    InnerTopRight(pro, 0, pro->col-1);
    for(int i = 1; i < pro->col-1; i++){
      InnerTop(pro, 0, i);
      BottomBoundary(pro, pro->row-1, i); 
    }
    for(int i = 1; i < pro->row -1; i++){
      InnerLeft(pro, i, 0);
      InnerRight(pro, i, pro->col-1);
    }
    BottomInnerLeft(pro, pro->row - 1, 0);
    BottomInnerRight(pro, pro->row, pro->col);
}

void Middle(processData *pro){
    for(int i = 1; i < pro->row - 1 ; i++){
        for(int j = 1; j < pro->col - 1 ; j++){
                MiddleMiddle(pro, i, j);
        }
    }
    for(int i = 1; i < pro->col-1; i++){
      InnerTop(pro, 0, i);
      InnerBottom(pro, pro->row-1, i); 
    }
    for(int i = 1; i < pro->row -1; i++){
      InnerLeft(pro, i, 0);
      InnerRight(pro, i, pro->col-1);
    }
    InnerTopLeft(pro, 0,0);
    InnerTopRight(pro, 0, pro->col - 1);
    InnerBottomLeft(pro, pro->row - 1, 0);
    InnerBottomRight(pro, pro->row - 1, pro->col - 1);
}

void calculateProcess(processData *pro){
    if((pro->myrank) % pro->sq == 0){     //Left side Processess
    LSP(pro);
  }
  else if((pro->myrank + 1) % pro->sq == 0){  //Right side Processess
    RSP(pro);
  }
  else if(pro->myrank < pro->sq){     //Top Middle
    TM(pro);
  }
  else if(pro->myrank > pro->sq*( pro->sq - 1 )){ //Lower Middle
    LM(pro);
  }
  else{
      Middle(pro);     
  }
pro->arr = pro->newarr;
pro->newarr = NULL;
}


int main( int argc, char *argv[])
{
    char *temp = argv[1], *temp2 = argv[2];
    int N = atoi(temp);
    int choice = atoi(temp2);           //{1 = Normal, 2 = Vector, 3 = Pack Transfer}
    int myrank, size, position=0, count;
    int row = N, col = N;             //row and column in one process
    
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank) ;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    processData pro;
    pro.myrank = myrank; pro.size = size; pro.row = row; pro.col = col; pro.sq = findSqrt(size);
    pro.arr = (double *)malloc(row* col * sizeof(double *)); 
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col ; j++)
            pro.arr[i * pro.row + j] = (double)rand();
    }

    pro.top = (double *)malloc(pro.col * sizeof(double));
    pro.bottom = (double *)malloc(pro.col * sizeof(double));
    pro.right = (double *)malloc(pro.col * sizeof(double));
    pro.left = (double *)malloc(pro.col * sizeof(double));
     
    double starttime, endtime;
    starttime = MPI_Wtime();
    if(choice == 1){                        //Normal Transfer
        for(int i = 0; i < 50; i++){
            NORMAL_SendTopReceiveBottom(&pro);
            NORMAL_SendBottomReceiveTop(&pro);
            NORMAL_SendLeftReceiveRight(&pro);
            NORMAL_SendRightReceiveLeft(&pro);
            pro.newarr = (double *)malloc(row* col * sizeof(double *));
            calculateProcess(&pro);
        }
    }
    else if(choice == 2){                  //Vector Mode
        for(int i = 0; i < 50; i++){
            VECTOR_SendTopReceiveBottom(&pro);
            VECTOR_SendBottomReceiveTop(&pro);
            VECTOR_SendLeftReceiveRight(&pro);
            VECTOR_SendRightReceiveLeft(&pro);
            pro.newarr = (double *)malloc(row* col * sizeof(double *));
            calculateProcess(&pro);
        }
    }
    else if(choice == 3){                 //Pack Mode
        for(int i = 0; i < 50; i++){
            PACK_SendTopReceiveBottom(&pro);
            PACK_SendBottomReceiveTop(&pro);
            PACK_SendLeftReceiveRight(&pro);
            PACK_SendRightReceiveLeft(&pro);
            pro.newarr = (double *)malloc(row* col * sizeof(double *));
            calculateProcess(&pro);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    endtime   = MPI_Wtime();
    MPI_Finalize();
    
    if(myrank == 0)
        printf("P = %d, N = %d , Type = %d, Time = %lf \n",size, N, choice, endtime-starttime);
    if(myrank == 0){
        FILE *fr;
        char s[40];
        int n;
        if(choice == 1){
            fr = fopen("Data_NormalVersion.txt", "a");
        }
        else if(choice == 3){
            fr = fopen("Data_PackVersion.txt", "a");
        }
        else{
            fr = fopen("Data_VectorVersion.txt", "a");
        }
                
        fprintf(fr, "%d, %d, %lf \n",size,N, endtime-starttime);
        fclose(fr);
    }
    return 0;
}
