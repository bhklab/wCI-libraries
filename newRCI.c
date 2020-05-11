/*
 compute wci in O(nlogn) using merge sort idea
 in this code we use merge sort 3 times
 by Farnoosh Khodakarami
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>

// Macros here are careful to not cause double evaluation (so that i++ only increments once).
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })


#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })


struct Res{
    double rCI;
    double valid_pairs;
};

typedef struct Res Res;


struct Items{
    double value; // value
    long int link_info; // link to the info array (this can be a polong inter as well)
    long int move; // the number of steps which is required to move.
    long int move_forward;
    int state; // the state of the item becuase we have main item && axilary item so if its 1 it is main item otherwise it is axilary item
    long int move_distance;
    long int type; //
};

typedef struct Items Items;

struct Items_info {
    double value; // value of the info which is main value
    long int link_main; // it can be a polong inter (link to the main value)
    long int link_axilary_left; // it can be a polong inter (link to the axilary value main-distance)
    long int link_axilary_main; // it can be a polong inter (link to the main value)
    long int link_axilary_right; // it can be a polong inter (link to the axilary value main+distance)
    long int link_axilary_lr; // it can be a polong inter (link to the axilary value main -distance, place in the right of the main )
    long int link_axilary_rl; // it can be a polong inter (link to the axilary value main+ distance, place in the left of the main) 
};

typedef struct Items_info Items_info;


void set_value_item(Items *x, double value_,long int link_info_,long int move_, int state_, long int type_ ){
        (*x).value = value_;
        (*x).link_info = link_info_; //it can be a polong inter
        (*x).move = move_;
        (*x).state = state_;
        (*x).type = type_;
        (*x).move_forward = 0;
        (*x).move_distance = 0;

    }

void set_value_info(Items_info *x, double value_,long int link_main_,long int link_axilary_left_,long int link_axilary_right_){
        (*x).value = value_;
        (*x).link_main = link_main_;
        (*x).link_axilary_left = link_axilary_left_;
        (*x).link_axilary_right = link_axilary_right_;
    }

 void set_axilary_info(Items_info *x,long int link_axilary_lr_,long int link_axilary_rl_){
        (*x).link_axilary_lr = link_axilary_lr_;
        (*x).link_axilary_rl = link_axilary_rl_;
    }



void merge(Items *arr, long int l, long int m, long int r)
{
    long int i, j, k;
    long int n1 = m - l + 1;
    long int n2 =  r - m;
    
    /* create temp arrays */
    //Items L[n1], R[n2];
    // TODO: Write allocation of memory here. 
    Items *L = malloc(n1 * sizeof(Items));
    Items *R = malloc(n2 * sizeof(Items));
    
    /* Copy data to temp arrays L[] && R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1+ j];
    
    /* Merge the temp arrays back long into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2)
    {   if(L[i].value == R[j].value){
        if(L[i].type <= R[j].type)
        {
            arr[k] = L[i];
            i++;
            
        }else{
            arr[k] = R[j];
            j++;
        }
    }
    else if (L[i].value < R[j].value)
    {
        arr[k] = L[i];
        i++;
    }else //if(L[i].value > R[j].value)
    {
        arr[k] = R[j];
        j++;
    }
        k++;
    }
    
    /* Copy the remaining elements of L[], if there
     are any */
    while (i < n1)
    {
        arr[k] = L[i];
        i++;
        k++;
    }
    
    /* Copy the remaining elements of R[], if there
     are any */
    while (j < n2)
    {
        arr[k] = R[j];
        j++;
        k++;
    }
    // TODO: replace with free statements

    free(R);
    free(L);
}

/* l is for left index and r is right index of the
 sub-array of arr to be sorted */
void mergeSort(Items *arr, long int l, long int r)
{
    if (l < r)
    {
        // Same as (l+r)/2, but avoids overflow for
        // large l && h
        long int m = l+(r-l)/2;
        
        // Sort first && second halves
        mergeSort(arr, l, m);
        mergeSort(arr, m+1, r);
        merge(arr, l, m, r);
    }
}

void Iterative_mergesort(Items *arr,long int n)
{
    long int size;
    long int l = 0;
    for (size=1; size<=n-1; size= 2*size)
    {
        for(l=0;l<n-1;l += 2*size)
        {
            long int m=min(l+size-1,n-1);
            long int r = min(l+2*size-1,n-1);
            merge(arr,l,m,r);
        }
        
    }
}


void merge_CI(Items *arr, long int l, long int m, long int r)
{
    long int i, j, k;
    long int n1 = m - l + 1;
    long int n2 =  r - m;
    
    /* create temp arrays */
    //Items L[n1], R[n2];

    // TODO: write allocation of memory here
    // printf("Got to merge_CI function\n");

    // printf("%d\n", n2);

    Items *L = malloc(n1 * sizeof(Items));
    Items *R = malloc(n2 * sizeof(Items));
    long int *steps = malloc((n1+1) * sizeof(long int));


    steps[n1] = 0;
    /* Copy data to temp arrays L[] && R[] */
    for (i = n1-1; i >= 0; i--){
        L[i] = arr[l + i];
        // printf("%d\n", L[i].state);
        if(L[i].state)
            steps[i] = steps[i+1] +1;
        else
            steps[i] = steps[i+1];
        
    }
    
    //long int move_forward = 0;
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1+ j];
    
    /* Merge the temp arrays back long into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2)
    {   if(L[i].value == R[j].value){
        if(L[i].type<= R[j].type )
        {
            arr[k] = L[i];
            i++;
        }
        else{
            arr[k] = R[j];
            arr[k].move = arr[k].move + steps[i];
            j++;
        }
    }else if (L[i].value < R[j].value)
    {
        arr[k] = L[i];
        i++;
    }else
    {
        arr[k] = R[j];
        arr[k].move = arr[k].move + steps[i];
        j++;
    }
        k++;
    }
    
    /* Copy the remaining elements of L[], if there
     are any */
    while (i < n1)
    {
        arr[k] = L[i];
        i++;
        k++;
    }
    
    /* Copy the remaining elements of R[], if there
     are any */
    while (j < n2)
    {
        arr[k] = R[j];
        j++;
        k++;
    }

    // TODO: replace with free statements
    free(steps);
    free(R);
    free(L);
}




/* l is for left index && r is right index of the
 sub-array of arr to be sorted */
void mergeSort_CI(Items *arr, long int l, long int r)
{
    if (l < r)
    {
        // Same as (l+r)/2, but avoids overflow for
        // large l && h
        long int m = l+(r-l)/2;
        
        // Sort first && second halves
        mergeSort_CI(arr, l, m);
        mergeSort_CI(arr, m+1, r);
        
        merge_CI(arr, l, m, r);
    }
}

void Iterative_mergesort_CI(Items *arr,long int n)
{
    long int size;
    long int l = 0;
    for (size=1; size<=n-1; size= 2*size)
    {
        for(l=0;l<n-1;l += 2*size)
        {
            long int m=min(l+size-1,n-1);
            long int r = min(l+2*size-1,n-1);
            merge_CI(arr,l,m,r);
        }
        
    }
}

void print_item(Items x){
  printf("#### Item Start: ####\n");
  printf("value: %e\n", x.value);
  printf("link_info: %d\n", x.link_info);
  printf("move: %d\n", x.move);
  printf("move_forward: %d\n", x.move_forward);
  printf("state: %d\n", x.state);
  printf("move_distance: %d\n", x.move_distance);
  printf("type: %d\n", x.type);
  printf("#### Item End ####\n\n");

}


Res rci(double *arr_X,double *arr_Y,double dist_x,double dist_y,long int n){
    
    Items *x = malloc(3*n * sizeof(Items));
    Items *y = malloc(6*n * sizeof(Items));
    Items *z = malloc(6*n * sizeof(Items));

    Items_info *x_info = malloc(n * sizeof(Items_info));
    Items_info *y_info = malloc(n * sizeof(Items_info));
    Items_info *z_info = malloc(n * sizeof(Items_info));

    // printf("Allocated arrays\n");

    for(long int i=0;i<n;i++){
        set_value_item(&x[3*i], arr_X[i]-dist_x,i,0,0,-1);
        set_value_item(&x[3*i+1], arr_X[i],i,0,1,0);
        set_value_item(&x[3*i+2], arr_X[i]+dist_x,i,0,0,+1);
        set_value_info(&x_info[i], arr_X[i],3*i+1,3*i,3*i+2);
        
    }

    // printf("Address of first element input data: %p\n", &arr_X[0]);
    // printf("Address of first Items data: %p\n", &(x[0].value));

    // printf("Initialized x arrays\n");

    /* sort the first array*/
    //mergeSort(x,0,3*n-1);
    Iterative_mergesort(x, 3*n);


    // printf("mergesort done\n");

    /* put the information in the info*/
    for(long int i=0;i<3*n;i++){
      // printf("\n%d\n", i);
      // printf("%d", x[i].type);
      // print_item(x[i]);
        switch (x[i].type) {
            case 0:
                x_info[x[i].link_info].link_main = i;
                break;
            case -1:
                x_info[x[i].link_info].link_axilary_left = i;
                break;
            case 1:
                x_info[x[i].link_info].link_axilary_right = i;
                break;
        }
        
    }
    // printf("xinfo initialized\n");

    /* based on the first array put the second array long into the write place */
    for(long int i=0;i<n;i++){
        set_value_item(&y[2*x_info[i].link_axilary_left], arr_Y[i]+dist_y,i,0,0,2);
        set_value_item(&y[2*x_info[i].link_axilary_left+1], arr_Y[i]-dist_y,i,0,0,-1);
        set_value_item(&y[2*x_info[i].link_main], arr_Y[i],i,0,0,0);
        set_value_item(&y[2*x_info[i].link_main+1], arr_Y[i],i,0,1,0);
        set_value_item(&y[2*x_info[i].link_axilary_right+1], arr_Y[i]-dist_y,i,0,0,-2);
        set_value_item(&y[2*x_info[i].link_axilary_right], arr_Y[i]+dist_y,i,0,0,1);
    }

    // printf("y initialized\n");

    // for(int i = 0; i < 6*n; i ++){
    //     print_item(y[i]);
    // }


    /* sort the second array && find the steps which each item needs to move */
    //mergeSort_CI(y, 0,6*n-1);
    Iterative_mergesort_CI(y, 6*n);
    // printf("mergesort CI y\n");


    for(long int i=0;i<6*n;i++){
        if(y[i].state){
            y_info[y[i].link_info].link_main = i;
        }else{
            switch ( y[i].type ) {
                case 2 :
                    y_info[y[i].link_info].link_axilary_lr =i;
                    
                    break;
                case -2 :
                    y_info[y[i].link_info].link_axilary_rl =i;
                    
                    break;
                case -1 :
                    y_info[y[i].link_info].link_axilary_left=i;
                    break;
                case 1 :
                    
                    y_info[y[i].link_info].link_axilary_right=i;
                    break;
                    
            }
        }
    }

    /*display(y, 6*n-1);*/
    // find the number of all steps which items moved
    long int steps = 0;
    for(long int i = 0;i<6*n;i++)
        if(!y[i].state && y[i].type == 2){
            steps = y[i].move + steps;
        }
    
    // to consider the items which come && placed inside  the long interval we consider z
    for(long int i=0;i<n;i++){
        set_value_item(&z[2*x_info[i].link_axilary_left], -1 * arr_Y[i]+dist_y,i,0,0,2);
        set_value_item(&z[2*x_info[i].link_axilary_left+1], -1* arr_Y[i]-dist_y,i,0,0,-1);
        set_value_item(&z[2*x_info[i].link_main], -1*arr_Y[i],i,0,0,0);
        set_value_item(&z[2*x_info[i].link_main+1], -1*arr_Y[i],i,0,1,0);
        set_value_item(&z[2*x_info[i].link_axilary_right+1], -1*arr_Y[i]-dist_y,i,0,0,-2);
        set_value_item(&z[2*x_info[i].link_axilary_right], -1*arr_Y[i]+dist_y,i,0,0,1);
        
    }

    //mergeSort_CI(z, 0, 6*n-1);
    Iterative_mergesort_CI(z,6*n);
    
    for(long int i=0;i<6*n;i++){
        if(z[i].state){
            z_info[z[i].link_info].link_main = i;
        }else{
            switch ( z[i].type ) {
                case 2 :
                    z_info[z[i].link_info].link_axilary_lr =i;
                    
                    break;
                case -2 :
                    z_info[z[i].link_info].link_axilary_rl =i;
                    
                    break;
                case -1 :
                    z_info[z[i].link_info].link_axilary_left=i;
                    break;
                case 1 :
                    
                    z_info[z[i].link_info].link_axilary_right=i;
                    break;
                    
            }
        }
    }
    
    long int *y_location = malloc((6*n)*sizeof(long int));

    if(y[0].state)
        y_location[0] = 1;
    else
        y_location[0] = 0;
    for(long int i=1;i<6*n;i++){
        if(!y[i].state)
            y_location[i] =y_location[i-1] ;
        else
            y_location[i]=y_location[i-1]+1;
    }
    
    long int *x_location = malloc((3*n) * sizeof(long int));

    if(!x[0].state)
        x_location[0] =0;
    else
        x_location[0]=1;
    for(long int i=1;i<3*n;i++){
        if(!x[i].state)
            x_location[i] =x_location[i-1] ;
        else
            x_location[i]=x_location[i-1]+1;
    }
    
    long int move_mines;
    long int move_plus;
    long int invalid_pairs = 0;
    long int incommon;
    // to consider the valid pairs we need to remove the pairs which are in long intervals && add those who removed twice
    for(long int i=0;i<n;i++){
        incommon = 0;
        if(y[y_info[i].link_axilary_left].move > 0){
            move_mines = y[y_info[i].link_axilary_rl].move - y[y_info[i].link_axilary_left].move;
            move_plus = y[y_info[i].link_axilary_right].move - y[y_info[i].link_axilary_lr].move;
            incommon = move_mines - move_plus;
            
        }else if(z[z_info[i].link_axilary_lr].move > 0 ){
            move_plus = ( z[z_info[i].link_axilary_lr].move - z[z_info[i].link_axilary_right].move);
            move_mines = ( z[z_info[i].link_axilary_left].move - z[z_info[i].link_axilary_rl].move);
            incommon = move_plus - move_mines;
        }else if(y[y_info[i].link_axilary_left].move == 0  &&  z[z_info[i].link_axilary_left].move ==0){
            move_mines = y[y_info[i].link_axilary_rl].move;
            move_plus = y[y_info[i].link_axilary_right].move;
            incommon = (move_mines - move_plus);
        }
        invalid_pairs = ( (y_location[y_info[i].link_axilary_right] - y_location[y_info[i].link_axilary_left] -1) + (x_location[x_info[i].link_axilary_right] - x_location[x_info[i].link_axilary_left] -1 ) - incommon) + invalid_pairs+1;
    }

    double valid_pairs = ((double)n*(n-1))/2;
    valid_pairs = ((long)valid_pairs - ((double)invalid_pairs)/2);
    double rci = (valid_pairs-steps) / valid_pairs;


    // TODO: rewrite as memfree
    // TODO: DONE?

    free(x_location);
    free(y_location);

    free(x_info);
    free(y_info);
    free(z_info);


    free(x);
    free(y);
    free(z);


    Res res;
    res.rCI = rci;
    res.valid_pairs = valid_pairs;

    return res;

}


SEXP rCI(SEXP pin_x,
                  SEXP pin_y,
                  SEXP pdelta_x,
                  SEXP pdelta_y,
                  SEXP pn,
                  SEXP pdiscard_x_ties, 
                  SEXP pdiscard_y_ties){

  long N = *INTEGER(pn);
  // printf("N: %d\\n", N);
  // TODO: Problem here?
  // int discard_x_ties = *INTEGER(pdiscard_x_ties);
  // int discard_y_ties = *INTEGER(pdiscard_y_ties);

  double delta_x = *REAL(pdelta_x);
  double delta_y = *REAL(pdelta_y);

  double *in_x = REAL(PROTECT(allocVector(REALSXP, N)));
  double *in_y = REAL(PROTECT(allocVector(REALSXP, N)));

  SEXP pout_rci = PROTECT(allocVector(REALSXP,1));
  SEXP pout_valid = PROTECT(allocVector(REALSXP,1));

  double *out_rci = REAL(pout_rci);
  double *out_valid = REAL(pout_valid);

  memcpy(in_x, REAL(pin_x), N * sizeof(*in_x));
  memcpy(in_y, REAL(pin_y), N * sizeof(*in_y));



  SEXP out = PROTECT(allocVector(VECSXP, 2));

  // printf("First x element is: %f\n", in_x[0]);
  Res res = rci(in_x, in_y, delta_x, delta_y, N);
  out_rci[0] = res.rCI;
  out_valid[0] = res.valid_pairs;

  // TODO: these should be REALSXP pointers
  SET_VECTOR_ELT(out, 0, pout_rci);
  SET_VECTOR_ELT(out, 1, pout_valid);

  UNPROTECT(5);

  return out;

}

static const R_CallMethodDef callMethods[]  = {
  {"rCI", (DL_FUNC) &rCI, 7},
  {NULL, NULL, 0}
};


