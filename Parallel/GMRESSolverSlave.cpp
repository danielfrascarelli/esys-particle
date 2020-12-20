/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2014 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0          //
//                                                         //
/////////////////////////////////////////////////////////////



#include "GMRESSolverSlave.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

using namespace std;



GMRESSolverSlave::GMRESSolverSlave(
  MPI_Comm* comm,
  TML_Comm* tml_comm,
  vector<double> B0,
  vector<int> row,
  vector<int> col,
  vector<double> v,
  int number,
  int nz
)
{
  m_worker_comm=comm;
  m_tml_global_comm=tml_comm;
  MPI_Comm_rank(*m_worker_comm, &id);
  MPI_Comm_size(*m_worker_comm, &nproc);

  Number=number;
  nrowsI=Number/nproc;
  B.resize(Number);
  for(int n=0;n<Number;n++){
    B[n]=B0[n];
  }
  nnz=nz;
  irow.resize(nnz); icol.resize(nnz);var.resize(nnz);
  for(int n=0;n<nnz;n++){
    irow[n]=row[n];
    icol[n]=col[n];
    var[n]=v[n];
  }
}

GMRESSolverSlave::~GMRESSolverSlave()
{
  vector<int> zero;
  irow.swap(zero);icol.swap(zero);ia.swap(zero);ja.swap(zero);
  vector<double> zero_db;
  var.swap(zero_db);A.swap(zero_db);B.swap(zero_db);b.swap(zero_db);
}



void GMRESSolverSlave::LocalMatrixSolving()
{
  /* local vars */
  int precond_flag = 0;//0 for normal GMRES and 1 for preconditioned GMRES algorithm
  double gmres_res = 1.e-5; //gmres residual tolerance
  int i = 0, j = 0, k, ss,  k_end = 200, mm = 10, k_sum=0;//k_end=number of GMRES iterations; m=number of restarts;

  LocalMatrixSetting();

  //localize
  int **receive = NULL;
  int *cntowner= NULL;
  int **to_be_sent= NULL;
  int *cntsent = NULL;
  localize(&receive, &cntowner, &to_be_sent, &cntsent);

  int n_total = max_int_array(ja, my_nnz)+1;
  double norm_r0 = 0.;

  double *u = (double *)calloc( nrowsI , sizeof(double) );
  double **v = (double **)malloc( 1 * sizeof(double *) );
  double **h = (double **)malloc( 1 * sizeof(double *) );
  double *x0 = (double *)calloc( n_total , sizeof(double) );
  double *r0 = (double *)calloc( nrowsI , sizeof(double) );
  double *g = (double *)calloc( 1 , sizeof(double) );
  double *c = NULL, *s = NULL, *resi = NULL;
  double delta = 0., gamma = 0.;

  vector<double> x_local;
  x_local.resize(nrowsI);

  //pick the diagonal entries for preconditioning (if it is preconditioned)
  double *diag = (double *) calloc( nrowsI , sizeof(double) );
  double *yk = (double *) calloc( n_total , sizeof(double) );
  if( precond_flag ) //if it is preconditioned
  {
    for (i = 0 ; i < nrowsI; i++ )
      for (j = ia[i] ; j < ia[i+1]; j++)
        if(ja[j] == i )
          diag[i] = A[j];
  }

  //init x0
  for(i = 0; i < nrowsI; i++)
    x0[i] = 1.0;

  //restarting vars
  int rst = 0,  completed = 0, total_itr = 1;

  //restarting loop
  for(rst = 0; rst < mm; rst++){
    gather_sync(receive, cntowner, to_be_sent, cntsent, x0);
    //matrix vector product
    sparse_matmult(x0, r0);
    //updating the final r0
    for(i = 0; i < nrowsI; i++)
      r0[i] = b[i] - r0[i];

    // normalization ...
    dotproduct(r0, r0, &norm_r0); //dot product
    norm_r0 = sqrt(norm_r0); //root
    for(i = 0; i < nrowsI; i++)
      r0[i] /= norm_r0; //normalizing

    // initial allocation of v
    h[0] = (double *)calloc(1 , sizeof(double) );
    v[0] = (double *)calloc(n_total , sizeof(double) );
    for(i = 0; i < nrowsI; i++)
      v[0][i] = r0[i]; //initializing ...

    g[0] = norm_r0; //initial g

    //gmres loop
    for(k = 0; k < k_end; k++)
    {
      //reallocating vars to have enough space for the next items
      g = (double *) realloc (g, (k+2) * sizeof(double) );
      g[k+1] = 0.;
      c = (double *) realloc (c, (k+1) * sizeof(double) );
      s = (double *) realloc (s, (k+1) * sizeof(double) );
      resi = (double *) realloc (resi, (k+1) * sizeof(double) );
      v = (double **)realloc (v, (k+2) * sizeof(double *) );
      v[k+1] = (double *)calloc(n_total , sizeof(double) );
      h = (double **)realloc (h , (k+2) * sizeof(double *) );
      h[k+1] = (double *)calloc(1 , sizeof(double) );

      for (j = 0; j <= (k+1); j++)
        h[j] = (double *)realloc( h[j], (k+1) *sizeof(double) );

      if( precond_flag )
      {
        for (i = 0 ; i < nrowsI; i++ )
          if( diag[i] )
            yk[i] = v[k][i] / diag[i];
          else
            yk[i] = v[k][i];
        //gather/sync
        gather_sync(receive, cntowner, to_be_sent, cntsent, yk);
        //matrix vector product
        sparse_matmult(yk, u);
      }
      else
      {
        //gather/sync
        gather_sync(receive, cntowner, to_be_sent, cntsent, v[k]);
        //matrix vector product
        sparse_matmult(v[k], u);
      }

      for (j = 0 ; j <= k; j++)
      {
        dotproduct(v[j], u, &h[j][k]); //dot product
        for( ss = 0; ss < nrowsI; ss++)
          u[ss] -= h[j][k] * v[j][ss];
      }
      dotproduct(u, u, &h[k+1][k]); //dot product
      h[k+1][k] = sqrt(h[k+1][k]); //norm
      //updating v[k+1]
      for(ss = 0; ss < nrowsI; ss++)
        v[k+1][ss] = u[ss] / h[k+1][k];

      for (j = 0 ; j < k; j++)
      {
        delta = h[j][k];
        h[j][k] = c[j] * delta + s[j] * h[j+1][k];
        h[j+1][k] = -s[j] * delta + c[j] * h[j+1][k];
      }
      gamma = sqrt(h[k][k] * h[k][k] + h[k+1][k]*h[k+1][k]);
      c[k] = h[k][k] / gamma;
      s[k] = h[k+1][k] / gamma;
      h[k][k] = gamma;
      h[k+1][k] = 0.;
      delta = g[k];
      g[k] = c[k] * delta + s[k] * g[k+1];
      g[k+1] = -s[k] * delta + c[k] * g[k+1];
      resi[k] = fabs(g[k+1]);

      if (resi[k] <= gmres_res)
      {
        completed = 1; //set the completed flag
        break; //terminate gmres loop ASAP!
      }

      if (rst==mm-1 && k==k_end-1 && resi[k] > gmres_res && id==0){
        cerr << "Fails to converge, please try a larger number of GMRES iterations or restarts \n";
        exit(0);
      };

    } //end of the main gmres loop

    k--; //reducing k to emulate the inside the loop effect
    k_sum+=k;

    //compute alpha
    double *alpha = (double *)calloc(k+1 , sizeof(double) );
    //solve backward
    for (j = k ; j >= 0; j--)
    {
      alpha[j] = g[j]/h[j][j];
      for (ss = (j+1) ; ss <= k; ss++)
        alpha[j] -= (h[j][ss]/h[j][j] * alpha[ss]);
    }
    //compute zk
    double *zk = (double *)calloc(nrowsI , sizeof(double) );
    for (j = 0 ; j <= k; j++)
      for (ss = 0 ; ss < nrowsI; ss++)
        zk[ss] += alpha[j] * v[j][ss];

    if(precond_flag )
      for (i = 0 ; i < nrowsI; i++ )
        if (diag[i] )
          zk[i] /= diag[i];

    //compute solution
    for (ss = 0 ; ss < nrowsI; ss++){
      x_local[ss] = x0[ss] + zk[ss];
    }

    //preparing for restart operation ...
    for(i = 0; i < nrowsI; i++)
	  x0[i] = x_local[i]; //put last x_local into x0

    //clean ups
    free(alpha);
    free(zk);
    for (j = 0 ; j <= (k+1); j++)
    {
      free(v[j]);
      free(h[j]);
    }

    if(completed)
      break; //terminate restart loop ASAP!

  } //end of restarting loop
  //Sending local results for gather by master
  m_tml_global_comm->send_gather(x_local,0);

  free(u);
  free(x0);
  free(r0);
  free(g);
  free(c);
  free(diag);
  free(yk);
  free(s);
  free(resi);
  free(v);
  free(h);
  free(cntowner);
  free(cntsent);
  for (j = 0 ; j < nproc; j++)
  {
      free(receive[j]);
      free(to_be_sent[j]);
  }
  free(receive);
  free(to_be_sent);
}



void GMRESSolverSlave::LocalMatrixSetting()
{
  //Setting local CRS matrix
  int i, i_start, i_end, DI;

  DI = Number/nproc;
  int *col_at_row = (int *)calloc( DI , sizeof(int) );
  i_start = id * DI;
  i_end = (id+1) * DI;

  int irow2,icol2;
  for(i = 0; i < nnz; i++) {
    irow2=irow[i]-1; //converting to zero-based
    if ( (i_start <= irow2) && (irow2 < i_end) ) {
      col_at_row[(irow2-i_start)]++;
    };
  }

  //find my total number of non-zero items (my_nnz) and also update ia[]
  ia.resize(DI+1);
  my_nnz=0;
  for(i= 0; i < DI; i++){
    ia[i] = my_nnz;
    my_nnz += col_at_row[i];
  }
  ia[i] = my_nnz;

  //allocate A one big chunck of doubles for local CRS array for storing data
  A.resize(my_nnz);
  //allocate (*ja)
  ja.resize(my_nnz);

  //fill data matrix A[], ia[] and ja[]
  for(i= 0; i < DI; i++)
    col_at_row[i] = 0; //reset col_at_row[]

  for( i = 0; i < nnz; i++){
    irow2=irow[i]-1;
    icol2=icol[i]-1;//converting to zero-based
    if ((i_start <= irow2) && (irow2 < i_end) ){
      irow2 -= i_start;
      A[ia[irow2] + col_at_row[irow2]] = var[i];
      ja[ia[irow2] + col_at_row[irow2]] = icol2;
      col_at_row[irow2]++;
    }
  }

  //setting local rhs
  b.resize(nrowsI);
  for(i = 0; i < Number; i++) {
    if ( (i_start <= i) && (i < i_end) )
      b[(i-i_start)] = B[i]; //put it in local [b]
  }

  //clean-up
  free(col_at_row);
}


//finds the maximum of the given integer array with size "n"
int GMRESSolverSlave::max_int_array(vector<int> input, int n)
{
  //locals
  int i = 0;
  int max = input[0];
  for( i = 0; i < n ; i++)
    if( input[i] > max)
      max = input[i];
  //returning the maximum value
  return max;
}


void GMRESSolverSlave::localize(int ***receive, int **cntowner, int ***to_be_sent, int **cntsent)
{
  // locals
  int i, j, indx, owner;
  int phn = nrowsI;
  int *tag = (int *)calloc(Number , sizeof(int) );
  (*receive) = (int **)malloc(nproc * sizeof(int *));
  (*to_be_sent) = (int **)malloc(nproc * sizeof(int *));
  int **expect = (int **)malloc(nproc * sizeof(int *));
  (*cntowner) = (int *)calloc(nproc , sizeof(int));
  (*cntsent) = (int *)calloc(nproc , sizeof(int));

  int *tmp_receive = NULL, *tmp_expect = NULL;
  int msg_tag1 = 40, msg_tag2 = 60;
  MPI_Status status1[nproc];
  MPI_Request request1[nproc];
  MPI_Status status2[nproc];
  MPI_Request request2[nproc];

  //initializing to NULL (safe)
  for ( i = 0; i < nproc; i++){
    (*receive)[i] = NULL;
    expect[i] = NULL;
    (*to_be_sent)[i] = NULL;
  }

  // main loop for hacking ja[]
  for ( i = 0; i < nrowsI; i++) { // i-local always starts from zero
    for ( indx = ia[i]; indx < ia[i+1]; indx++) //ia[] always starts from zero for this process
    {
      j=ja[indx]; //find the global column (zero-based)
      owner = (int)(j/nrowsI); //and the owner rank
      if (owner == id)
        ja[indx] -= (id*nrowsI); // NOTE: my_rank*nlocal=vect_start_indx for this process
      else
      {
        if( !tag[j] ) //is it already localized? if no do it!
        {
	  //updating receive list
  	  tmp_receive = (int *) realloc ((*receive)[owner],  ((*cntowner)[owner]+1)*sizeof(int));
   	  if(tmp_receive == NULL) {
	    cout<<"can't realloc receive list on processor "<<id<<endl;
	    fflush(stdout);
	    exit(0);
	  } else
	    (*receive)[owner] = tmp_receive;

	  (*receive)[owner][(*cntowner)[owner]] = phn;
	  //updating expect list
	  tmp_expect = (int *) realloc (expect[owner],  ((*cntowner)[owner]+1) * sizeof(int));
	  if(tmp_expect == NULL) {
            cout<<"can't realloc expect list on processor "<<id<<endl;
	    fflush(stdout);
	    exit(0);
	  } else
	    expect[owner] = tmp_expect;

	  expect[owner][(*cntowner)[owner]] = j;
	  // updating countors
	  tag[j] = phn;
	  ((*cntowner)[owner])++;
	  phn++;
	}
        ja[indx] = tag[j];
      }
    }
  }

  // sending expect list to other processes
  int rq_inc1 = 0, rq_inc2 = 0; //increaments for requests
  for( owner = 0; owner < nproc; owner++) //loop over process which we want to send expect list
    if ( owner == id) //dont send to myself
      continue;
    else
    {
      //sending the numbers of expected values to that process
      MPI_Isend(&(*cntowner)[owner], 1, MPI_INT, owner, msg_tag1, *m_worker_comm, request1+rq_inc1++);
      //sending expected values to that process
      MPI_Isend(expect[owner], (*cntowner)[owner], MPI_INT, owner, msg_tag2, *m_worker_comm, request2+ rq_inc2++);
    }

  // receiving to_be_sent list from other processes
  int st_inc1 = 0, st_inc2 = 0; //increaments for status
  for( owner = 0; owner < nproc; owner++) //loop over process which we want to receive from them.
    if ( owner == id)
      continue;
    else
    {
      //receiving the numbers of expected values for that process
      MPI_Recv(&(*cntsent)[owner], 1, MPI_INT, owner, msg_tag1, *m_worker_comm, status1+st_inc1++ );
      //once get the size of data, reallocating enough space
      tmp_expect = (int *)realloc ((*to_be_sent)[owner], (*cntsent)[owner]*sizeof(int) );
      if(tmp_expect == NULL){
        printf("\ncan't realloc to_be_sent list!\n");
        exit(0);
      }
      else {
        (*to_be_sent)[owner] = tmp_expect;
      }
      // reciving data and putting them in place
      MPI_Recv((*to_be_sent)[owner],(*cntsent)[owner], MPI_INT, owner, msg_tag2, *m_worker_comm, status2+st_inc2++ );
      for ( i = 0; i < (*cntsent)[owner]; i++) //localizing data
        (*to_be_sent)[owner][i] -= (id*nrowsI);
    }

  // wait until all send and receives ar complete
  MPI_Waitall(nproc-1, request1, status1);
  MPI_Waitall(nproc-1, request2, status2);

  // clean-up
  for( i = 0; i < nproc; i++)
    if( i != id)
      free(expect[i]);
  free(expect);
  free(tag);
}



// for gathering data in one location in the local process and synchronization among processes
void GMRESSolverSlave::gather_sync(int **receive, int *cntowner, int **to_be_sent, int *cntsent, double *qk)
{
  int i = 0;
  int owner = 0;
  int msg_tag = 40;
  MPI_Status status[nproc];
  MPI_Request request[nproc];
  double **buff_snd = (double **)malloc(nproc * sizeof(double *));
  double **buff_rcv = (double **)malloc(nproc * sizeof(double *));

  // sending doubles to other processes
  int rq_inc = 0; //increament for requests
  for(owner = 0; owner < nproc; owner++) //loop over process which we want to send
    if (owner == id) //dont send to myself
      continue;
    else
    {
      // first pack doubles into the send buffer
      buff_snd[owner] = (double *)malloc(cntsent[owner] * sizeof(double));
      for( i = 0; i < cntsent[owner]; i++)
	    buff_snd[owner][i] = qk[to_be_sent[owner][i]];

      //sending expected values to that process
      MPI_Isend(buff_snd[owner], cntsent[owner], MPI_DOUBLE, owner, msg_tag, *m_worker_comm, request+rq_inc++);
    }

  // receiving doubles from other processes and putting them in the right places
  int st_inc = 0; //increament for status
  for( owner = 0; owner < nproc; owner++) //loop over process which we want to receive from them.
    if ( owner == id)
      continue;
    else
    {
      // first allocate receive buffer
      buff_rcv[owner] = (double *)malloc(cntowner[owner] * sizeof(double));
      // reciving data and putting them in receive buffer
      MPI_Recv(buff_rcv[owner], cntowner[owner], MPI_DOUBLE, owner, msg_tag, *m_worker_comm, status+st_inc++ );
      for ( i = 0; i < cntowner[owner]; i++) //putting in right place indicated by ja[]
        qk[receive[owner][i]] = buff_rcv[owner][i];
    }

  // wait until all send and receives ar complete
  MPI_Waitall(nproc-1, request, status);

  //clean - ups
  for( i = 0; i < nproc ; i++)
    if( i != id) {
      free(buff_snd[i]);
      free(buff_rcv[i]);
    }
    free(buff_rcv);
    free(buff_snd);
}


void GMRESSolverSlave::sparse_matmult(double *q, double *u)
{
  //locals
  int r, j;
  // reset u
  for( r = 0 ; r < nrowsI; r++)
    u[r] = 0.;

  for( r = 0 ; r < nrowsI; r++)
    for ( j = ia[r]; j < ia[r+1]; j++) {
     u[r] += (A[j] * q[ja[j]]);
    }
}


void GMRESSolverSlave::dotproduct(double *a, double *b, double *c)
{
  //locals
  int i = 0;
  double local_dot = 0.;
  //perform local dot product
  for( i = 0; i < nrowsI; i++)
    local_dot += (a[i] * b[i]);

  //reduceall to single values
  MPI_Allreduce(&local_dot,c,1,MPI_DOUBLE,MPI_SUM,*m_worker_comm);
}



