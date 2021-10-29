//#define PY_SSIZE_T_CLEAN
//#include <Python.h>

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <time.h> 

#include <pthread.h>

//////////////////
#define VREDUCE1
//////////////////
#define VREDUCE2


//#define COMBIDX


// Extract hom cycles
#define HOM_CYCLES

// Store V adaptively to extract homology basis
#define ADAPTIVE_V_STORAGE

// minimize lengths of birth cycles

//#define MINIMIZE_BIRTH_CYCLES

#define STORE_LENGTHS_CYCLES

//#define ADD_0PERS_CYCLES


//#define MINIMIZE_HOM_CYCLES

/////////////////////////////////////////////////
// min of max dist. reduction
/////////////////////////////////////////////////

#define DISTMAT_MINMAX

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

//#define POINTCLOUD_MINMAX

#define RECORD_V_USAGE

/////////////////
// SAVING MACROS
/////////////////

//#define SAVEV

//#define SAVEPD

//#define PRINT

/////////////////
// DEBUG MACROS
/////////////////
//#define DEBUGCOMBIDX

//#define VDEBUG

//#define DEBUGPIVOTS

//#define COH1DEBUG


//#define COMB_IDX(a, b)((a) > (b) ? 
//                            self->g_edges_comb_idx[(EDGE_ID)((a*(a-1))/2 + b)] : \
//                                        ( (a) < (b) ? self->g_edges_comb_idx[(EDGE_ID)((b*(b-1))/2 + a)] \
//                                              : self->g_n_valid_edges))
//




#define COMB_IDX0(a,b) ( a > b ? (EDGE_ID)((a*(a-1))/2 + b) : (EDGE_ID)((b*(b-1))/2 + a) )

#define COMB_IDX(a,b) ( a == b ? self->g_n_valid_edges : self->g_edges_comb_idx[COMB_IDX0(a, b)])


typedef unsigned long long int BIGINT;

typedef unsigned int EDGE_ID;
typedef unsigned int VERT_ID;
typedef double  PAR;


typedef struct{
    VERT_ID verts[2];
    PAR length;
}F1;

typedef struct{

  VERT_ID neighbor;
  EDGE_ID order;

}Neighbors;

typedef struct{
    
    EDGE_ID key1;
    EDGE_ID key2;

}simplex;

typedef struct{

    EDGE_ID col_idx;
    EDGE_ID o_ab;
    
}H0_pivots;


typedef struct{

  EDGE_ID key2;
  EDGE_ID col_idx;
  EDGE_ID bndry;

}H1_cohom_pivots;

typedef struct{

  EDGE_ID key2;
  EDGE_ID col_idx;
  simplex bndry;

}H2_cohom_pivots;

typedef struct{
    
  int a_ptr;
  int b_ptr;
  EDGE_ID o_ab;
  simplex low;

}coboundary_H1;

typedef struct{

  // The simplex
  simplex triangle;
  // Note: triangle.key1 is o_ab
  // Note: triangle.key2 is c

  VERT_ID a_ptr;
  VERT_ID b_ptr;
  VERT_ID c_ptr;


  // The low of the simplex
  simplex low;
  //key1 is 0: ab, 1: ad, 2: bd, 3: cd
  int vertex;

  //vertex is -1 should mean empty. But have not been consistent.
  //low.key1 = n_valid_edges also means empty.

}coboundary_H2;

typedef struct{
      
      int original;

      EDGE_ID len;
      EDGE_ID max_len;
      
      int flag_first;
      int flag_reduce;
      int flag_red_w_complex;
      int flag_append_to_complex;
      int flag_empty;

      EDGE_ID pivot;
      simplex triangle;

      EDGE_ID R_col_idx;
      EDGE_ID reduce_with_len;
      EDGE_ID* trivial_boundary;

      
}boundary_H1_ws;


typedef struct{
      
      int original;

      EDGE_ID len;
      EDGE_ID max_len;

      int flag_first;
      int flag_reduce;
      int flag_red_w_complex;
      int flag_append_to_complex;
      int flag_empty;

      simplex pivot;
      simplex tetrahedron;

      EDGE_ID R_col_idx;
      EDGE_ID reduce_with_len;
      simplex* trivial_boundary;
      
}boundary_H2_ws;


typedef struct{
      
      int original;

      EDGE_ID cob;

      EDGE_ID len;
      EDGE_ID max_len;
      
      int flag_red_w_complex;
      int flag_append_to_complex;
      int flag_non_empty;

      EDGE_ID pivot;

      
}boundary_H0_ws;

typedef struct{
    
    EDGE_ID max_len;
    EDGE_ID last;
    EDGE_ID* o_ab;
}edges_list;

typedef struct{
    
    EDGE_ID k2;
    EDGE_ID o_ab;
    EDGE_ID a_ptr;
    EDGE_ID b_ptr;
    int flag_next;
}implicit_keys2;


typedef struct{
    
    EDGE_ID k1;

    int flag_empty;

    implicit_keys2* keys2;
    EDGE_ID max_len;
    EDGE_ID last;

}implicit_keys1;


typedef struct{
    
    EDGE_ID k1_ptr;
    EDGE_ID k2_ptr;

    EDGE_ID edge;

    implicit_keys1* keys1;
    edges_list v_edges;

    EDGE_ID max_len;
    EDGE_ID last;

    simplex pivot;


    int flag_first;
    int flag_red_w_complex;
    int flag_red_w_trivial;

    EDGE_ID reduce_w_bndry;
    EDGE_ID V_col_idx;

    int flag_append_to_complex;
    int flag_non_empty;


}coboundary_H1_ws;



typedef struct{
    
    EDGE_ID max_len;
    EDGE_ID last;
    simplex* o_abc;

}triangles_list;

typedef struct{
    
    EDGE_ID k2;
    simplex o_abc;

    VERT_ID a_ptr;
    VERT_ID b_ptr;
    VERT_ID c_ptr;
    int vertex;

    int flag_next;

}coH2_implicit_keys2;


typedef struct{
    
    EDGE_ID k1;

    int flag_empty;

    coH2_implicit_keys2* keys2;
    EDGE_ID max_len;
    EDGE_ID last;

}coH2_implicit_keys1;


typedef struct{
    
    EDGE_ID k1_ptr;
    EDGE_ID k2_ptr;

    simplex triangle;

    coH2_implicit_keys1* keys1;
    triangles_list v_triangles;

    EDGE_ID max_len;
    EDGE_ID last;

    simplex pivot;

    int flag_first;
    int flag_red_w_complex;
    int flag_red_w_trivial;
    simplex reduce_w_bndry;
    EDGE_ID V_col_idx;

    int flag_append_to_complex;
    int flag_non_empty;


}coboundary_H2_ws;


typedef struct{
    
    EDGE_ID len;
    EDGE_ID max_len;
    EDGE_ID* VV;
}hom1_birth;

typedef struct{
    
    EDGE_ID* RR;
    int original;
    EDGE_ID len;
    EDGE_ID max_len;
    
}R_struct;


typedef struct{
    
    EDGE_ID birth_edge;
    EDGE_ID death_triangle_key1;
    EDGE_ID R_col_idx;
    
}homH1_pers;


typedef struct{
    
    EDGE_ID len;
    EDGE_ID max_len;
    simplex* VV;
}hom2_birth;

typedef struct{
    
    EDGE_ID* RR;
    int original;
    EDGE_ID len;
    EDGE_ID max_len;
    
}R_struct_H2;


typedef struct{
    
    simplex birth_simplex;
    EDGE_ID death_edge;
    EDGE_ID R_col_idx;
    
}homH2_pers;

typedef struct{
    
    // o_ab
    EDGE_ID key1;
    // o_cd
    EDGE_ID key2;

    EDGE_ID o_ac;
    EDGE_ID o_ad;
    EDGE_ID o_bd;
    EDGE_ID o_bc;

}H2_preprocess;

typedef struct{

  EDGE_ID key2;
  EDGE_ID col_idx;
  simplex tetrahedron;

}H2_pivots;

typedef struct{
    
  EDGE_ID coface;

  EDGE_ID V_usage;

#ifdef RECORD_V_USAGE
  EDGE_ID V_depth;
#endif

  int V_stored;

  EDGE_ID V_len;
  EDGE_ID* VV;

}V_H0;

typedef struct{

  simplex coface;

  EDGE_ID V_usage;

#ifdef RECORD_V_USAGE
  EDGE_ID V_depth;
#endif

  int V_stored;

  EDGE_ID V_len;
  simplex* VV;

}V_H1;


typedef struct{
    
    EDGE_ID cycid;
    EDGE_ID Lidx;

    EDGE_ID V_len;
    EDGE_ID* VV;

    //EDGE_ID ops_len;
    //EDGE_ID* ops;

}min_update_V;


typedef struct{
    
    EDGE_ID cycid;
    EDGE_ID Lidx;

    EDGE_ID V_len;
    simplex* VV;

    //EDGE_ID ops_len;
    //EDGE_ID* ops;

}min_update_V_H2;


//#ifdef MINIMIZE_BIRTH_CYCLES

typedef struct{

    EDGE_ID* boundary;
    EDGE_ID len;

    EDGE_ID redw;
    EDGE_ID diff;

    //EDGE_ID* ops;
    //EDGE_ID ops_len;

    PAR perspair[2];
    EDGE_ID Lidx;

    PAR updated_birth;

    //EDGE_ID* in_cycles;
    //EDGE_ID in_cycles_len;
    //EDGE_ID in_cycles_max_len;

}cyc_info;


typedef struct{
    
    EDGE_ID cyc;
    int flag;

}update_in_cyc;

typedef struct{
    EDGE_ID cj;
    EDGE_ID diff;
}cyc_in_cyc;



typedef struct{

    simplex* boundary;
    EDGE_ID len;

    EDGE_ID redw;
    EDGE_ID diff;

    PAR perspair[2];
    EDGE_ID Lidx;

    PAR updated_birth;

}cyc_info_H2;


//#endif




int compare_neighbors_vertex(Neighbors s1, Neighbors s2){
    
      if (s1.neighbor < s2.neighbor) return -1;
      else if (s1.neighbor > s2.neighbor) return 1;
      else return 0;

 
}

int compare_neighbors_order(Neighbors s1, Neighbors s2){
    
      if (s1.order < s2.order) return -1;
      else if (s1.order > s2.order) return 1;
      else return 0;

 
}

// This is for tim sort
int compare_simplex(simplex s1, simplex s2){
    
    if (s1.key1 > s2.key1) return 1;
    if (s1.key1 < s2.key1) return -1;
    else{

        if (s1.key2 > s2.key2) return 1;
        else if (s1.key2 < s2.key2) return -1;
        else return 0;

    }

      
}



// This is for tim sort
int compare_cob_H1(coboundary_H1 s1, coboundary_H1 s2){
    
    if (s1.o_ab > s2.o_ab) return 1;
    else if (s1.o_ab < s2.o_ab) return -1;
    else return 0;

      
}

int compare_coboundary_H2(coboundary_H2 s1, coboundary_H2 s2){
    

      if (s1.triangle.key1 < s2.triangle.key1) return -1;
      else if (s1.triangle.key1 > s2.triangle.key1) return 1;
      else{
            if (s1.triangle.key2 < s2.triangle.key2) return -1;
            else if (s1.triangle.key2 > s2.triangle.key2) return 1;
            else return 0;
      }


 
}

int compare_EDGE_ID ( EDGE_ID x, EDGE_ID y){
    
    if (x < y) return -1;
    else if (x > y) return 1;
    else return 0;
}

#define SORT_NAME sorter
//#define SORT_TYPE int64_t
#define SORT_TYPE Neighbors
#define SORT_CMP(x, y) compare_neighbors_vertex((x), (y))
/* You can redefine the comparison operator.
   The default is
#define SORT_CMP(x, y)  ((x) < (y) ? -1 : ((x) == (y) ? 0 : 1))
   but the one below is often faster for integer types.
*/
//#define SORT_CMP(x, y) ((x->key1) < (y->key1) ? -1 : ((x->key1) == (y->key1) ? ((x->key2 < y->key2 ? -1 : ((x->key2) == (y->key2) ? 0: 1))) : 1))

//#define SORT_CMP(x, y,) compare_partial_cob_dec((x), (y))
#include "sort.h"


#undef SORT_NAME
#undef SORT_TYPE
#undef SORT_CMP

#define SORT_NAME sorter2
//#define SORT_TYPE int64_t
#define SORT_TYPE Neighbors
#define SORT_CMP(x, y) compare_neighbors_order((x), (y))
/* You can redefine the comparison operator.
   The default is
#define SORT_CMP(x, y)  ((x) < (y) ? -1 : ((x) == (y) ? 0 : 1))
   but the one below is often faster for integer types.
*/
//#define SORT_CMP(x, y) ((x->key1) < (y->key1) ? -1 : ((x->key1) == (y->key1) ? ((x->key2 < y->key2 ? -1 : ((x->key2) == (y->key2) ? 0: 1))) : 1))

//#define SORT_CMP(x, y,) compare_partial_cob_dec((x), (y))
#include "sort2.h"


#undef SORT_NAME
#undef SORT_TYPE
#undef SORT_CMP

#define SORT_NAME sorter3
//#define SORT_TYPE int64_t
#define SORT_TYPE EDGE_ID
#define SORT_CMP(x, y) compare_EDGE_ID(x, y)
/* You can redefine the comparison operator.
   The default is
#define SORT_CMP(x, y)  ((x) < (y) ? -1 : ((x) == (y) ? 0 : 1))
   but the one below is often faster for integer types.
*/
//#define SORT_CMP(x, y) ((x->key1) < (y->key1) ? -1 : ((x->key1) == (y->key1) ? ((x->key2 < y->key2 ? -1 : ((x->key2) == (y->key2) ? 0: 1))) : 1))

//#define SORT_CMP(x, y,) compare_partial_cob_dec((x), (y))
#include "sort3.h"


#undef SORT_NAME
#undef SORT_TYPE
#undef SORT_CMP

#define SORT_NAME sorter4
//#define SORT_TYPE int64_t
#define SORT_TYPE simplex
#define SORT_CMP(x, y) compare_simplex((x), (y))
//#define SORT_CMP(x, y) compare_EDGE_ID(x, y)
/* You can redefine the comparison operator.
   The default is
#define SORT_CMP(x, y)  ((x) < (y) ? -1 : ((x) == (y) ? 0 : 1))
   but the one below is often faster for integer types.
*/
//#define SORT_CMP(x, y) ((x->key1) < (y->key1) ? -1 : ((x->key1) == (y->key1) ? ((x->key2 < y->key2 ? -1 : ((x->key2) == (y->key2) ? 0: 1))) : 1))

//#define SORT_CMP(x, y,) compare_partial_cob_dec((x), (y))
#include "sort4.h"

#undef SORT_NAME
#undef SORT_TYPE
#undef SORT_CMP

#define SORT_NAME sorter5
//#define SORT_TYPE int64_t
#define SORT_TYPE coboundary_H1
#define SORT_CMP(x, y) compare_cob_H1((x), (y))
//#define SORT_CMP(x, y) compare_EDGE_ID(x, y)
/* You can redefine the comparison operator.
   The default is
#define SORT_CMP(x, y)  ((x) < (y) ? -1 : ((x) == (y) ? 0 : 1))
   but the one below is often faster for integer types.
*/
//#define SORT_CMP(x, y) ((x->key1) < (y->key1) ? -1 : ((x->key1) == (y->key1) ? ((x->key2 < y->key2 ? -1 : ((x->key2) == (y->key2) ? 0: 1))) : 1))

//#define SORT_CMP(x, y,) compare_partial_cob_dec((x), (y))
#include "sort5.h"


#undef SORT_NAME
#undef SORT_TYPE
#undef SORT_CMP

#define SORT_NAME sorter6
//#define SORT_TYPE int64_t
#define SORT_TYPE coboundary_H2
#define SORT_CMP(x, y) compare_coboundary_H2((x), (y))
//#define SORT_CMP(x, y) compare_EDGE_ID(x, y)
/* You can redefine the comparison operator.
   The default is
#define SORT_CMP(x, y)  ((x) < (y) ? -1 : ((x) == (y) ? 0 : 1))
   but the one below is often faster for integer types.
*/
//#define SORT_CMP(x, y) ((x->key1) < (y->key1) ? -1 : ((x->key1) == (y->key1) ? ((x->key2 < y->key2 ? -1 : ((x->key2) == (y->key2) ? 0: 1))) : 1))

//#define SORT_CMP(x, y,) compare_partial_cob_dec((x), (y))
#include "sort6.h"




#undef SORT_NAME
#undef SORT_TYPE
#undef SORT_CMP

#define SORT_NAME sorter7
//#define SORT_TYPE int64_t
#define SORT_TYPE implicit_keys2

int compare_implicit_keys2(implicit_keys2 s1, implicit_keys2 s2){
    
      if (s1.k2 < s2.k2) return -1;
      else if (s1.k2 > s2.k2) return 1;
      else{

            if (s1.o_ab < s2.o_ab) return -1;
            else if (s1.o_ab > s2.o_ab) return 1;
            else return 0;

      }

 
}

#define SORT_CMP(x, y) compare_implicit_keys2((x), (y))
//#define SORT_CMP(x, y) compare_EDGE_ID(x, y)
/* You can redefine the comparison operator.
   The default is
#define SORT_CMP(x, y)  ((x) < (y) ? -1 : ((x) == (y) ? 0 : 1))
   but the one below is often faster for integer types.
*/
//#define SORT_CMP(x, y) ((x->key1) < (y->key1) ? -1 : ((x->key1) == (y->key1) ? ((x->key2 < y->key2 ? -1 : ((x->key2) == (y->key2) ? 0: 1))) : 1))

//#define SORT_CMP(x, y,) compare_partial_cob_dec((x), (y))
#include "sort7.h"




#undef SORT_NAME
#undef SORT_TYPE
#undef SORT_CMP

#define SORT_NAME sorter8
//#define SORT_TYPE int64_t
#define SORT_TYPE EDGE_ID

int compare_edges(EDGE_ID s1, EDGE_ID s2){
    
      if (s1 < s2) return -1;
      else if (s1 > s2) return 1;
      else return 0;

 
}

#define SORT_CMP(x, y) compare_edges((x), (y))
//#define SORT_CMP(x, y) compare_EDGE_ID(x, y)
/* You can redefine the comparison operator.
   The default is
#define SORT_CMP(x, y)  ((x) < (y) ? -1 : ((x) == (y) ? 0 : 1))
   but the one below is often faster for integer types.
*/
//#define SORT_CMP(x, y) ((x->key1) < (y->key1) ? -1 : ((x->key1) == (y->key1) ? ((x->key2 < y->key2 ? -1 : ((x->key2) == (y->key2) ? 0: 1))) : 1))

//#define SORT_CMP(x, y,) compare_partial_cob_dec((x), (y))
#include "sort8.h"


#undef SORT_NAME
#undef SORT_TYPE
#undef SORT_CMP

#define SORT_NAME sorter9
//#define SORT_TYPE int64_t
#define SORT_TYPE coH2_implicit_keys2

int coH2_compare_implicit_keys2(coH2_implicit_keys2 s1, coH2_implicit_keys2 s2){
    
      if (s1.k2 < s2.k2) return -1;
      else if (s1.k2 > s2.k2) return 1;
      else{

            if (s1.o_abc.key1 < s2.o_abc.key1) return -1;
            else if (s1.o_abc.key1 > s2.o_abc.key1) return 1;
            else{

                if (s1.o_abc.key2 < s2.o_abc.key2) return -1;
                else if (s1.o_abc.key2 > s2.o_abc.key2) return 1;
                else return 0;
                  
            }

      }

 
}

#define SORT_CMP(x, y) coH2_compare_implicit_keys2((x), (y))
//#define SORT_CMP(x, y) compare_EDGE_ID(x, y)
/* You can redefine the comparison operator.
   The default is
#define SORT_CMP(x, y)  ((x) < (y) ? -1 : ((x) == (y) ? 0 : 1))
   but the one below is often faster for integer types.
*/
//#define SORT_CMP(x, y) ((x->key1) < (y->key1) ? -1 : ((x->key1) == (y->key1) ? ((x->key2 < y->key2 ? -1 : ((x->key2) == (y->key2) ? 0: 1))) : 1))

//#define SORT_CMP(x, y,) compare_partial_cob_dec((x), (y))
#include "sort9.h"




#undef SORT_NAME
#undef SORT_TYPE
#undef SORT_CMP

#define SORT_NAME sorter10
//#define SORT_TYPE int64_t
#define SORT_TYPE H2_preprocess

int compare_tetra_key2(H2_preprocess s1, H2_preprocess s2){
    
     if (s1.key2 < s2.key2) return -1;
     else if (s1.key2 > s2.key2) return 1;
     else return 0;
                  

 
}

#define SORT_CMP(x, y) compare_tetra_key2((x), (y))
//#define SORT_CMP(x, y) compare_EDGE_ID(x, y)
/* You can redefine the comparison operator.
   The default is
#define SORT_CMP(x, y)  ((x) < (y) ? -1 : ((x) == (y) ? 0 : 1))
   but the one below is often faster for integer types.
*/
//#define SORT_CMP(x, y) ((x->key1) < (y->key1) ? -1 : ((x->key1) == (y->key1) ? ((x->key2 < y->key2 ? -1 : ((x->key2) == (y->key2) ? 0: 1))) : 1))

//#define SORT_CMP(x, y,) compare_partial_cob_dec((x), (y))
#include "sort10.h"






int compare_implicit(implicit_keys2 s1, coboundary_H1 phi){
    
      if (s1.k2 < phi.low.key2) return 0;
      else if (s1.k2 > phi.low.key2) return 1;
      else{

            if (s1.o_ab < phi.o_ab) return 0;
            else return 1;
          

      }

 
}


int coH2_compare_implicit(coH2_implicit_keys2 s1, coboundary_H2 phi){
    
      if (s1.k2 < phi.low.key2) return 0;
      else if (s1.k2 > phi.low.key2) return 1;
      else{

            if (s1.o_abc.key1 < phi.triangle.key1) return 0;
            else if (s1.o_abc.key1 > phi.triangle.key1) return 1;
            else{

                if (s1.o_abc.key2 < phi.triangle.key2) return 0;
                else return 1;

            }

      }

 
}


typedef struct{

    char* filename;


#ifdef SAVEPD
    char* g_H0_pers_file;
    char* g_H1_pers_file;
    char* g_H2_pers_file;
#endif

#ifdef SAVEV
    char* g_coH1_V_file;
    char* g_coH2_V_file;
#endif


    int g_suppress_output;

//#ifdef MINIMIZE_BIRTH_CYCLES
    char* g_minimal_V_H0_file;
    //char* g_minimal_V_H0_in_cycles_file;

//#ifdef STORE_LENGTHS_CYCLES
    char* g_V_H0_birthcyc_lens_file;
    char* g_minimal_V_H0_birthcyc_lens_file;

    char* g_birth_subset_points_file_H0;

    char* g_V_H1_birthcyc_lens_file;
    char* g_minimal_V_H1_birthcyc_lens_file;

//#endif

    char* g_minimal_V_H1_file;
//#endif

#ifdef MINIMIZE_HOM_CYCLES
    char* g_minimal_V_hom_H1_file;
    char* g_minimal_V_hom_H2_file;

#ifdef STORE_LENGTHS_CYCLES
    char* g_minimal_V_H1_homcyc_lens;
    char* g_minimal_V_H2_homcyc_lens;
#endif

#endif

#ifdef RECORD_V_USAGE
    char* g_V_H0_usage_file;
    char* g_V_H1_usage_file;
#endif
    
    //char* g_file_prefix;
    //const char* g_H1_boundaries;
    //const char* g_H1_indices;
    //const char* g_H2_boundaries;


    char* g_source;
    char* g_target;
    

    int g_cpu_count;

    int g_dim_lim;

    int g_compute_cycles;

    int g_reduce_cyc_lengths;

    int g_filetype;


    PAR g_thresh;

    VERT_ID g_n_vert;

    EDGE_ID g_n_valid_edges;

    BIGINT g_n_all_simp;

    EDGE_ID* g_edges_list;

    PAR* g_edge_parameter;

#ifdef COMBIDX
    // Combinatorial idx
    EDGE_ID g_n_edges;
    EDGE_ID* g_edges_comb_idx;

#endif


    // Neighbor data structures
    Neighbors** g_Neighbors;
    Neighbors** g_Neighbors_e;
    VERT_ID* g_Neigh_len;
    EDGE_ID g_max_neighbors;



    // WORKSPACE PARAMETERS
    int g_workspace_size;
    int g_ws_pre_alloc;
    int g_ws_counter;


    ////////////////////////////////////
    // Parallel job allocation
    ////////////////////////////////////
    int* g_jobs;

    int g_sleeping_threads;
    int g_processed_threads;

    int g_thread_id;

    int g_delete_threads;

    pthread_t *g_threads;

    pthread_mutex_t g_thread_lock;
    pthread_cond_t g_start_boss;
    pthread_cond_t g_start_workers;


    ////////////////////////////////////
    // H0 Structures
    ////////////////////////////////////

    // Pivots for H0
    // i is pivot of A[i]
    EDGE_ID* g_pivots_H0;

    // STORE R for H0
    EDGE_ID* g_R_sparse_H0;
    EDGE_ID g_R_sparse_ptr_H0;
    EDGE_ID g_R_sparse_max_H0;

    // Mapping of R columns to sparse R linear
    EDGE_ID* g_R_col_indices_H0;
    EDGE_ID g_R_col_indices_max_H0;
    EDGE_ID g_R_col_indices_ptr_H0;


    // Store pivot for H0
    EDGE_ID* g_edges_with_pivots_H0;

    // H0 WORKSPACE STRUCTURES
    EDGE_ID** g_R_ws_H0; 

    boundary_H0_ws* g_R_ws_H0_info;

    ////////////////////////////////////
    // cohomology H1 structures
    ////////////////////////////////////
   
    EDGE_ID g_this_edge;

    coboundary_H1* g_coH1_all_lows;
    
    // V SPARSE
    
    EDGE_ID* g_V_sparse_H1;
    EDGE_ID g_V_sparse_max;
    EDGE_ID g_V_sparse_ptr;
    EDGE_ID g_V_sparse_beg_ptr;
    EDGE_ID g_V_sparse_end_ptr;

    EDGE_ID* g_V_col_indices;
    EDGE_ID g_V_col_indices_max;
    EDGE_ID g_V_col_indices_ptr;

    // V workspace
    coboundary_H1_ws* g_V_ws_H1;


    // PIVOTS OF H1 COHOMOLOGY
    H1_cohom_pivots** g_H1_cohom_pivots;
    EDGE_ID* g_H1_cohom_pivots_len;
    EDGE_ID* g_H1_cohom_pivots_max_len;

    // Pers pairs
    EDGE_ID g_H1_pers_pairs_max_len;
    EDGE_ID g_H1_pers_pairs_len;
    PAR* g_H1_pers_pairs;
 
    
    ////////////////////////////////////
    // cohomology H2 structures
    ////////////////////////////////////
    
    simplex* g_V_sparse_H2;

    // WORKSPACE STRUCTURES
    int g_cohom_ws_size;
    coboundary_H2_ws* g_V_ws_H2; 

    
    // NEW PIVOTS OF H2 COHOMOLOGY
    H2_cohom_pivots** g_H2_cohom_pivots;
    EDGE_ID* g_H2_cohom_pivots_len;
    EDGE_ID* g_H2_cohom_pivots_max_len;


    // Pers pairs
    EDGE_ID g_H2_pers_pairs_max_len;
    EDGE_ID g_H2_pers_pairs_len;
    PAR* g_H2_pers_pairs;

    ////////////////////////////////////
    // Timers
    ////////////////////////////////////
    double g_timer_H2_low;
    double g_timer_H2_next;
    double g_timer_H2_greater;

    struct timespec g_start_wall_clock;
    struct timespec g_finish_wall_clock;

    double g_timer_process_input;
    double g_timer_neigh;
    double g_timer_H0;
    double g_timer_coH1;
    double g_timer_coH2;


    double g_timer_computeH1;
    double g_timer_computeH2;

    double g_timer_H1cycles;
    double g_timer_H2cycles;

    double g_timer_minimize_H1cycles;
    double g_timer_minimize_H2cycles;

    double g_timer_minimize_H1_homcycles;
    double g_timer_minimize_H2_homcycles;

    double g_timer_coH2_serial;
    double g_timer_coH2_parallel;

    BIGINT g_n_H1_birth_cycles;
    BIGINT g_n_H2_birth_cycles;

    BIGINT g_n_H0_stored_V;
    BIGINT g_n_H1_stored_V;

    // Temporary
    int g_p_flag;
    EDGE_ID g_counter;

    ////////////////////////////////////
    // homology H1 structures
    ////////////////////////////////////

    char* g_homH1_cycles_file;

    EDGE_ID g_R_max_len_H1;
    EDGE_ID g_R_len_H1;
    EDGE_ID* g_R_H1;


    EDGE_ID g_R_col_idx_max_len_H1;
    EDGE_ID* g_R_col_idx_H1;

    EDGE_ID g_R_col_idx_H1_ptr;


    EDGE_ID* g_pivots_H1;
    //

    EDGE_ID** g_workspace_H1;
    boundary_H1_ws* g_workspace_H1_info;

    ////////////////////////////////////
    // For H1 birth cycles
    ////////////////////////////////////
#ifdef HOM_CYCLES
    V_H0* g_H0_pivot_of;
#endif


    R_struct g_temp_R_birth_cycles;
    hom1_birth g_temp_V_primary;
    EDGE_ID g_homH1_pers_len;
    EDGE_ID g_homH1_pers_max_len;
    homH1_pers* g_homH1_pers;

    EDGE_ID* g_H1_undead;
    EDGE_ID g_H1_undead_ptr;
    EDGE_ID g_H1_undead_max;

    EDGE_ID g_depth;



    ////////////////////////////////////
    // homology H2 structures
    ////////////////////////////////////

    char* g_homH2_cycles_file;

    EDGE_ID g_R_max_len_H2;
    EDGE_ID g_R_len_H2;
    simplex* g_R_H2;

    EDGE_ID g_R_col_idx_max_len_H2;
    EDGE_ID* g_R_col_idx_H2;

    EDGE_ID g_R_col_idx_H2_ptr;

    H2_pivots** g_H2_pivots;
    EDGE_ID* g_H2_pivots_len;
    EDGE_ID* g_H2_pivots_max_len;


    simplex** g_workspace_H2;
    boundary_H2_ws* g_workspace_H2_info;




    ////////////////////////////////////
    // For H2 birth cycles
    ////////////////////////////////////
    
#ifdef HOM_CYCLES
    V_H1* g_H1_pivot_of;
#endif

    R_struct_H2 g_temp_R_H2_birth_cycles;
    hom2_birth g_temp_V_H2_primary;

    // Pers pairs info
    EDGE_ID g_homH2_pers_len;
    EDGE_ID g_homH2_pers_max_len;
    homH2_pers* g_homH2_pers;

    simplex* g_H2_undead;
    EDGE_ID g_H2_undead_ptr;
    EDGE_ID g_H2_undead_max;

    int g_extract_cycles;

#ifdef ADAPTIVE_V_STORAGE
    EDGE_ID g_cycle_usage_thresh;
    EDGE_ID g_cycle_depth_thresh;


    EDGE_ID g_store_V_for_len;
    EDGE_ID g_store_V_for_max_len;
    EDGE_ID* g_store_V_for;

    simplex* g_store_V_voids_for;

#endif

//#ifdef MINIMIZE_BIRTH_CYCLES

    EDGE_ID g_global_minimizer;

    EDGE_ID g_all_V_stored_num;
    EDGE_ID g_all_V_stored_max_num;

    cyc_info* g_all_V_H0_stored;


    //EDGE_ID** g_all_V_H0_stored;

    cyc_info_H2* g_all_V_H1_stored;

    EDGE_ID** g_edges_in_cycles;
    EDGE_ID* g_edges_in_cycles_len;

//#endif

#ifdef MINIMIZE_HOM_CYCLES

    cyc_info* g_all_V_hom_H1_stored;

    EDGE_ID g_all_V_hom_stored_num;
    
    // LEGACY
    EDGE_ID g_all_V_hom_stored_max_num;
    EDGE_ID* g_all_V_hom_stored_len;
    EDGE_ID** g_all_V_hom_H1_stored;
    simplex** g_all_V_hom_H2_stored;

#endif



    int g_new_debug;
    int g_new_debug2;
    EDGE_ID g_debug_edge;
    simplex g_debug_triangle;


    // Cycle minimization birth threshold
    //PAR g_cycle_min_birth_thresh;
    
} filtration;



int simplex1_check(VERT_ID, VERT_ID, PAR, PAR);

int simplex2_check(VERT_ID, VERT_ID, VERT_ID);

int simplex3_check(VERT_ID, VERT_ID, VERT_ID, VERT_ID);


// MERGE SORT ALGORITHM
void mergeSort(PAR* , EDGE_ID* , EDGE_ID , EDGE_ID ) ;
void merge(PAR* , EDGE_ID* , EDGE_ID , EDGE_ID , EDGE_ID ) ;

// MERGE SORT ALGORITHM 2
void mergeSort_V_H0(EDGE_ID* , EDGE_ID** , EDGE_ID*, EDGE_ID*, EDGE_ID , EDGE_ID ) ;
void merge_V_H0(EDGE_ID* , EDGE_ID** , EDGE_ID*, EDGE_ID*, EDGE_ID , EDGE_ID , EDGE_ID ) ;

// MERGE SORT ALGORITHM 3
void mergeSort_V_H1(EDGE_ID* , simplex** , EDGE_ID , EDGE_ID ) ;
void merge_V_H1(EDGE_ID* , simplex** , EDGE_ID , EDGE_ID , EDGE_ID ) ;
//
// MERGE SORT ALGORITHM 4
void mergeSort_update_V(min_update_V* , EDGE_ID , EDGE_ID ) ;
void merge_update_V(min_update_V* , EDGE_ID , EDGE_ID , EDGE_ID ) ;

// MERGE SORT ALGORITHM 5
void mergeSort_update_V_byLidx(min_update_V* , EDGE_ID , EDGE_ID ) ;
void merge_update_V_byLidx(min_update_V* , EDGE_ID , EDGE_ID , EDGE_ID ) ;

// MERGE SORT ALGORITHM 6: Sort Lcycid, Llen, Lupdated by Llen
void mergeSort_Llen(EDGE_ID*, EDGE_ID*,  EDGE_ID*, EDGE_ID, EDGE_ID ) ;
void merge_Llen(EDGE_ID*, EDGE_ID*, EDGE_ID*, EDGE_ID, EDGE_ID, EDGE_ID ) ;


// MERGE SORT ALGORITHM 7: Sort 
void mergeSort_temp_par(PAR*, EDGE_ID*, EDGE_ID, EDGE_ID ) ;
void merge_temp_par(PAR*, EDGE_ID*, EDGE_ID, EDGE_ID, EDGE_ID ) ;

// MERGE SORT ALGORITHM 8: Sort 
void mergeSort_incycleslen(EDGE_ID*, cyc_info*, EDGE_ID, EDGE_ID ) ;
void merge_incycleslen(EDGE_ID*, cyc_info*, EDGE_ID, EDGE_ID, EDGE_ID ) ;

// MERGE SORT ALGORITHM 9: Sort 
void mergeSort_edges_in_cycles(EDGE_ID*, cyc_info*, EDGE_ID, EDGE_ID ) ;
void merge_edges_in_cycles(EDGE_ID*, cyc_info*, EDGE_ID, EDGE_ID, EDGE_ID ) ;

// MERGE SORT ALGORITHM 10: Sort 
void mergeSort_edges_in_cycles_bycycid(EDGE_ID*, EDGE_ID, EDGE_ID ) ;
void merge_edges_in_cycles_bycycid(EDGE_ID*, EDGE_ID, EDGE_ID, EDGE_ID ) ;



//////////////////////////////////////////////////////////////
//       NEIGHBOR CREATION AND SEARCH ALGORITHMS
//////////////////////////////////////////////////////////////

void update_neighbors_new(filtration* , VERT_ID , VERT_ID , EDGE_ID);

VERT_ID search_Neighbors(filtration* , VERT_ID , VERT_ID , VERT_ID , VERT_ID);
VERT_ID search_Neighbors_e(filtration* , VERT_ID , EDGE_ID , VERT_ID , VERT_ID, EDGE_ID);

VERT_ID bin_search_min_geq_Ne(Neighbors* , VERT_ID, VERT_ID, VERT_ID, EDGE_ID);
VERT_ID bin_search_min_geq_N(Neighbors* , VERT_ID, VERT_ID, VERT_ID, EDGE_ID);

EDGE_ID bin_search_cycle_ops(EDGE_ID*, EDGE_ID, EDGE_ID, EDGE_ID, EDGE_ID);

EDGE_ID bin_search_cyc_in_cyc(cyc_in_cyc* , EDGE_ID , EDGE_ID , EDGE_ID , EDGE_ID );

//////////////////////////////////////////////////////////////



// H0 HOMOLOGY FUNCTIONS

// Parallel homology reduction H0
//main reduction
void reduce_ws_H0(filtration* );
//reduction with complex
void* reduce_with_complex_H0(void* );
//reduction with self
void reduce_with_self_H0(filtration* );
//Update R
void update_R_H0(filtration* , int );


void allocate_jobs(filtration*, int);



BIGINT compute_num_simplices(filtration* );

// H1 cohomology functions


void update_V_coH1 (filtration*, int);

void find_H1_cohom_next (filtration* , coboundary_H1* );
void find_H1_cohom_low(filtration* , coboundary_H1* );
void find_H1_cohom_greater(filtration* , coboundary_H1* , simplex* );




void insert_in_implicit_v(filtration* , int, coboundary_H1*, int);
void print_v_implicit(filtration*, int );

void reduce_ws_coH1(filtration* );
void reduce_with_self_coH1(filtration* );
void* reduce_with_complex_coH1(void* );
void reduce_hash_table_coH1(filtration*, int );

void coH2_insert_in_implicit_v(filtration*, int , coboundary_H2* , int );
void reduce_hash_table_coH2(filtration*, int );
void coH2_print_v_implicit(filtration*, int );


void find_H2_cohom_next (filtration* , coboundary_H2* );
void find_H2_cohom_low(filtration* , coboundary_H2* );
void find_H2_cohom_greater(filtration* , coboundary_H2* , simplex* );

int H2_case1 (filtration*, coboundary_H2*);
void H2_case2 (filtration*, coboundary_H2*);


EDGE_ID search_H1_cohom_pivots(H1_cohom_pivots* , EDGE_ID , EDGE_ID , EDGE_ID , EDGE_ID ); 
EDGE_ID search_H2_cohom_pivots(H2_cohom_pivots* , EDGE_ID , EDGE_ID , EDGE_ID , EDGE_ID ); 


void H2_reduce (filtration*, coboundary_H2*, EDGE_ID, int);
void update_V_coH2(filtration* , int );
void add_coH2_pivot(filtration*, simplex, simplex, EDGE_ID);

void reduce_ws_coH2(filtration* );
void* reduce_with_complex_coH2(void* );
void reduce_with_self_coH2(filtration* );



// H1 HOMOLOGY FUNCTIONS
//main reduction
void reduce_ws_H1(filtration* );
//reduction with complex
void* reduce_with_complex_H1(void* );
//reduction with self
void reduce_with_self_H1(filtration* );
//Update R
void update_R_H1(filtration* , int );


EDGE_ID bin_search_max_less_V(EDGE_ID* , EDGE_ID , EDGE_ID , EDGE_ID , EDGE_ID );

EDGE_ID bin_search_min_greater_updated_V_byLidx(EDGE_ID* , EDGE_ID , EDGE_ID , EDGE_ID , EDGE_ID );

EDGE_ID find_first_diff_H0(EDGE_ID* \
                          , EDGE_ID \
                          , EDGE_ID** );

// H1 cycles
void compute_H1_homology_cycles(filtration* );
void get_birth_cycle(filtration*, EDGE_ID);
void find_V_recursively_edges(filtration*, EDGE_ID, EDGE_ID);

void shuffle_cyc(cyc_info*, EDGE_ID);

//#ifdef MINIMIZE_BIRTH_CYCLES
void minimize_birth_cycles_H0(filtration*, EDGE_ID**, EDGE_ID*, EDGE_ID, char*);

void minimize_birth_cycles_H0_v2(filtration*, EDGE_ID**, EDGE_ID*, EDGE_ID, char*);

void minimize_birth_cycles_H1(filtration*);

void minimize_birth_cycles_H1_v2(filtration* \
                              , cyc_info_H2* \
                              , EDGE_ID \
                              , char* \
                              , char* \
                              , char* \
                              );



void minimize_birth_cycles_H0_v3(filtration* \
                              , cyc_info* \
                              , EDGE_ID \
                              , char* \
                              , char* \
                              , char* \
                              , char* \
                              );

void minimize_birth_cycles_H0_v4(filtration* \
                              , cyc_info* \
                              , EDGE_ID \
                              , char* \
                              , char* \
                              );

void minimal_CASE1(EDGE_ID , cyc_info* , EDGE_ID* , EDGE_ID*, EDGE_ID );

void minimal_CASE2(filtration*, EDGE_ID , cyc_info* , EDGE_ID* , EDGE_ID* \
                  , EDGE_ID* , EDGE_ID );

void update_diff(filtration* , EDGE_ID , EDGE_ID* , int \
                , EDGE_ID* , cyc_info* , EDGE_ID);

void find_first_diff(filtration* , EDGE_ID , EDGE_ID*\
                , EDGE_ID* , cyc_info* , EDGE_ID);

//#endif

void store_V_H0(filtration* );
void reduce_temp_V_H0(filtration* );

// H2 HOMOLOGY FUNCTIONS
//main reduction
void reduce_ws_H2(filtration* );
//reduction with complex
void* reduce_with_complex_H2(void* );
//reduction with self
void reduce_with_self_H2(filtration* );
//Update R
void update_R_H2(filtration* , int );
//add_pivot 
void add_H2_pivot (filtration* , simplex , simplex , EDGE_ID );
//search pivot
EDGE_ID search_H2_pivots(H2_pivots* , EDGE_ID , EDGE_ID , EDGE_ID , EDGE_ID );

// H2 cycles
void compute_boundary_triangle(filtration* , simplex , EDGE_ID* );
void compute_boundary_tetra(filtration* , simplex , simplex* );
void compute_H2_homology_cycles(filtration* );
void get_birth_void(filtration*, simplex);
void find_V_recursively_triangles(filtration*, simplex, EDGE_ID);


void store_V_H1(filtration* );
void reduce_temp_V_H1(filtration* );

// DEALLOCATE
void deallocator(filtration*);



int main(int argc, char* argv[]){
//static PyObject *compute_PH(PyObject *self2, PyObject *args){

     struct timespec start_wall_clock, finish_wall_clock;
     clock_gettime(CLOCK_MONOTONIC, &start_wall_clock);

    // Filetype = 0 : Distance matrix
    // Filetype = 1 : Locations
    // Filetype = 2 : Edge list with edge length

     //printf("%ld", (long)getpid());

     filtration* self;
     self = (filtration*)malloc(sizeof(filtration));

     self->g_new_debug = 0;

     //////////////////////////////////////////////////////
     // Set testing timers and test counters
     //////////////////////////////////////////////////////
     self->g_counter = 0;

     self->g_timer_H2_low = 0;
     self->g_timer_H2_next = 0;
     self->g_timer_H2_greater = 0;

     self->g_timer_coH2_serial = 0;
     self->g_timer_coH2_parallel = 0;
     //////////////////////////////////////////////////////
     //////////////////////////////////////////////////////

     //if (!PyArg_ParseTuple(args, "sdiisiiii"\
     //                               , &(self->g_source), &(self->g_thresh)\
     //                               , &(self->g_filetype), &(self->g_cpu_count)\
     //                               , &(self->g_target), &(self->g_dim_lim)\
     //                               , &(self->g_compute_cycles), &(self->g_reduce_cyc_lengths)\
     //                               , &(self->g_suppress_output)\
     //                               )){
     //    printf("\nERROR in parse args");
     //    return NULL;
     //                                //, &(self->g_cycle_min_birth_thresh)\

     //}

     //int file_len = strlen(self->g_target) + 100;
     
     int file_len = strlen(argv[5]) + 100;


     self->g_thresh = atof(argv[2]);

     self->g_filetype = atoi(argv[3]);

     self->g_cpu_count = atoi(argv[4]);

     self->g_dim_lim = atoi(argv[6]);
      
     self->g_compute_cycles = atoi(argv[7]);

     self->g_reduce_cyc_lengths = atoi(argv[8]);

     self->g_suppress_output = atoi(argv[9]);



     char* duplicate = (char*)malloc(file_len*sizeof(char));

     //strcpy(duplicate, self->g_target);

     strcpy(duplicate, argv[1]);
     self->filename = strdup(duplicate);
     
     strcpy(duplicate, argv[5]);
     self->g_target = strdup(duplicate);


     strcpy(duplicate, self->g_target);
     strcat(duplicate, "homH1_cycles.txt");
     self->g_homH1_cycles_file = strdup(duplicate);


     strcpy(duplicate, self->g_target);
     strcat(duplicate, "homH2_cycles.txt");
     self->g_homH2_cycles_file = strdup(duplicate);
     

//#ifdef MINIMIZE_BIRTH_CYCLES
     strcpy(duplicate, self->g_target);
     strcat(duplicate, "minimal_V_birth_H1.txt");
     self->g_minimal_V_H0_file = strdup(duplicate);

     strcpy(duplicate, self->g_target);
     strcat(duplicate, "birth_subsets_H1.txt");
     self->g_birth_subset_points_file_H0 = strdup(duplicate);

     //strcpy(duplicate, argv[5]);
     //strcat(duplicate, "minimal_V_birth_H1_in_cycles.txt");
     //self->g_minimal_V_H0_in_cycles_file = strdup(duplicate);


     strcpy(duplicate, self->g_target);
     strcat(duplicate, "minimal_V_birth_H2.txt");
     self->g_minimal_V_H1_file = strdup(duplicate);

//#ifdef STORE_LENGTHS_CYCLES


     strcpy(duplicate, self->g_target);
     strcat(duplicate, "V_birth_len_H1.txt");
     self->g_V_H0_birthcyc_lens_file = strdup(duplicate);

     strcpy(duplicate, self->g_target);
     strcat(duplicate, "minimal_V_birth_len_H1.txt");
     self->g_minimal_V_H0_birthcyc_lens_file = strdup(duplicate);


     strcpy(duplicate, self->g_target);
     strcat(duplicate, "V_birth_len_H2.txt");
     self->g_V_H1_birthcyc_lens_file = strdup(duplicate);

     strcpy(duplicate, self->g_target);
     strcat(duplicate, "minimal_V_birth_len_H2.txt");
     self->g_minimal_V_H1_birthcyc_lens_file = strdup(duplicate);

  
//#endif




//#endif









#ifdef MINIMIZE_HOM_CYCLES
     if (self->g_reduce_cyc_lengths){

          strcpy(duplicate, self->g_target);
          strcat(duplicate, "minimal_V_hom_H1.txt");
          self->g_minimal_V_hom_H1_file = strdup(duplicate);

          strcpy(duplicate, self->g_target);
          strcat(duplicate, "minimal_V_hom_H2.txt");
          self->g_minimal_V_hom_H2_file = strdup(duplicate);

     }



#endif


#ifdef RECORD_V_USAGE
     strcpy(duplicate, self->g_target);
     strcat(duplicate, "V_H0_usage.txt");
     self->g_V_H0_usage_file = strdup(duplicate);

     strcpy(duplicate, self->g_target);
     strcat(duplicate, "V_H1_usage.txt");
     self->g_V_H1_usage_file = strdup(duplicate);
#endif


#if defined(SAVEPD) || defined(SAVEV)

     strcpy(duplicate, self->g_target);
     strcat(duplicate, "H0_pers_data.txt");

     self->g_H0_pers_file = strdup(duplicate);


     if (self->g_dim_lim > 0){

          strcpy(duplicate, self->g_target);
          strcat(duplicate, "H1_pers_data.txt");

          self->g_H1_pers_file = strdup(duplicate);

#ifdef SAVEV
          strcpy(duplicate, self->g_target);
          strcat(duplicate, "coH1_V_data.txt");
          self->g_coH1_V_file = strdup(duplicate);
#endif

          if (self->g_dim_lim > 1){

               strcpy(duplicate, self->g_target);
               strcat(duplicate, "H2_pers_data.txt");
               self->g_H2_pers_file = strdup(duplicate);

#ifdef SAVEV
               strcpy(duplicate, self->g_target);
               strcat(duplicate, "coH2_V_data.txt");
               self->g_coH2_V_file = strdup(duplicate);
#endif

          }

     }

#endif

     free(duplicate);


     //self->g_extract_cycles = atoi(argv[6]);

     //self->g_cycle_birth_limit = atof(argv[7]);

#ifdef ADAPTIVE_V_STORAGE
     self->g_cycle_usage_thresh = 2;
     self->g_cycle_depth_thresh = 1;
#endif
   
     
     omp_set_num_threads(self->g_cpu_count);


     FILE *fp = fopen(self->filename, "r");  

     if (fp == NULL){
          perror("Unable to open file!");
          exit(1);
     }

     char* line = NULL;
     size_t len = 0;
     char* dist;
     PAR dist_d;
     char* end;

     getline(&line, &len, fp);
     dist = strtok(line, " ,");

     fclose(fp);

     fp = fopen(self->filename, "r");  

     int prealloc = 100000;

     VERT_ID row = 0;
     VERT_ID col = 0;

     self->g_edges_list = (EDGE_ID*)malloc(2*prealloc*sizeof(EDGE_ID));
     self->g_edge_parameter = (PAR*)malloc(prealloc*sizeof(PAR));

     self->g_n_valid_edges = 0;

     if (self->g_filetype == 0){

#ifdef DISTMAT_MINMAX
        PAR dist_min_max, dist_max;
        dist_min_max = INFINITY; 
        // this is a distance matrix
        while(getline(&line, &len, fp) != -1) {

            //col = 0;
            
            dist = strtok(line, " ,");
            dist_max = 0;
            while(dist != NULL){
              dist_d = strtod(dist, &end);
              dist = strtok(NULL, ",");

              if (dist_d > dist_max){
                  dist_max = dist_d;
              }

              //col += 1;
            }


            if (dist_max < dist_min_max){
                dist_min_max = dist_max;
            }
            //row += 1;
        }
        self->g_thresh = dist_min_max;
	
        rewind(fp);

#endif

        EDGE_ID edge_list_ptr = 0;

        // this is a distance matrix
        while(getline(&line, &len, fp) != -1) {

            col = 0;
            
            dist = strtok(line, " ,");
            while(dist != NULL){
              dist_d = strtod(dist, &end);
              //if (dist_d != 0) dist_d = 1/dist_d;
              dist = strtok(NULL, ",");
              if (col > row){
                   
                   //if (simplex1_check(row, col, dist_d, self->g_thresh)){
                   if (dist_d < self->g_thresh){

                         //self->g_edges_list[self->g_n_valid_edges] = (EDGE_ID*)malloc(2*sizeof(EDGE_ID));

                         // Note that g_edges_list is sorted, row < col
                         //self->g_edges_list[self->g_n_valid_edges][0] = row;
                         //self->g_edges_list[self->g_n_valid_edges][1] = col;

                         self->g_edges_list[edge_list_ptr++] = row;
                         self->g_edges_list[edge_list_ptr++] = col;


		                     // parameter
                         self->g_edge_parameter[self->g_n_valid_edges] = dist_d;

                         self->g_n_valid_edges += 1;
                         if (self->g_n_valid_edges == prealloc){
                               
                               prealloc += 100000;
                               self->g_edges_list = (EDGE_ID*)realloc(self->g_edges_list, 2*prealloc*sizeof(EDGE_ID));
                               self->g_edge_parameter = (PAR*)realloc(self->g_edge_parameter, prealloc*sizeof(PAR));

                         }
                   }

              }
              col += 1;
            }
            row += 1;
        }

	      self->g_n_vert = row;

     }
     else if (self->g_filetype == 1){

       self->g_thresh = self->g_thresh * self->g_thresh;

       //Locations information
       
       if (!self->g_suppress_output){
          printf("extracting edges");
       }
       int dim_space = 0;

       while(getline(&line, &len, fp) != -1) {
          dist = strtok(line, " ,");
          while(dist != NULL){
            dist_d = strtod(dist, &end);
            dist = strtok(NULL, ",");
            dim_space++;
          }
          break;
       }

       rewind(fp);

       PAR** locations;
       locations = (PAR**)malloc(sizeof(PAR*));


        while(getline(&line, &len, fp) != -1) {

            col = 0;
            if (!self->g_suppress_output){
                printf("\rrow %d", row);
            }
            
            locations = (PAR**)realloc(locations, (row+1)*sizeof(PAR*));
            locations[row] = (PAR*)malloc(dim_space*sizeof(PAR));

            dist = strtok(line, " ,");
            while(dist != NULL){

              dist_d = strtod(dist, &end);
              //if (dist_d != 0) dist_d = 1/dist_d;
              dist = strtok(NULL, ",");

              locations[row][col++] = dist_d; 

            }

            row++;

        }

        PAR diff;

#ifdef POINTCLOUD_MINMAX

        if (!self->g_suppress_output){
            printf("\nthresh is %lf", self->g_thresh);
        }
        PAR dist_min_max, dist_max;
        dist_min_max = INFINITY; 
        for (int i = 0; i < row; i++){

              //if (i%1000 == 0)
              if (!self->g_suppress_output){
                  printf("\n%d", i);
              }

              dist_max = 0;
              for (int j = 0; j < row; j++){
                  for (int k = 0; k < dim_space; k++){
                      diff = locations[i][k] - locations[j][k];
                      dist_d += diff*diff;                        
                  }
                  if (dist_d > dist_max){
                      dist_max = dist_d;
                  }
              }
              if (dist_max < dist_min_max){
                  dist_min_max = dist_max;
              }

        }

        if (self->g_thresh > dist_min_max){
            self->g_thresh = dist_min_max;
        }

        if (!self->g_suppress_output){
            printf("\nupdated thresh is %lf", self->g_thresh);
        }

#endif

        EDGE_ID edge_list_ptr = 0;

        if (!self->g_suppress_output){
        printf("\n");
        }
        for (int i = 0; i < row-1; i++){

              if (!self->g_suppress_output){
                  printf("Done %f percent edges %d\r", (float)i/(float)(row-1), self->g_n_valid_edges);
              }

              for (int j = i+1; j < row; j++){
                    
                  dist_d = 0;
                  for (int k = 0; k < dim_space; k++){
                          
                        diff = locations[i][k] - locations[j][k];
                        dist_d += diff*diff;                        

                  }

                  //dist_d = sqrt(dist_d);

                  //if (simplex1_check(i, j, dist_d, self->g_thresh)){
                  if (dist_d < self->g_thresh){

                        //self->g_edges_list[self->g_n_valid_edges] = (EDGE_ID*)malloc(2*sizeof(EDGE_ID));

                        // Note that g_edges_list is sorted, row < col
                        //self->g_edges_list[self->g_n_valid_edges][0] = i;
                        //self->g_edges_list[self->g_n_valid_edges][1] = j;

                        self->g_edges_list[edge_list_ptr++] = i;
                        self->g_edges_list[edge_list_ptr++] = j;


		                    // parameter
                        self->g_edge_parameter[self->g_n_valid_edges] = dist_d;

                        self->g_n_valid_edges += 1;
                        if (self->g_n_valid_edges == prealloc){
                              
                              prealloc += 100000;
                              self->g_edges_list = (EDGE_ID*)realloc(self->g_edges_list, 2*prealloc*sizeof(EDGE_ID));
                              self->g_edge_parameter = (PAR*)realloc(self->g_edge_parameter, prealloc*sizeof(PAR));

                        }
                  }
                    
              }
              
        }


        for (int i = 0; i < row; i++) free(locations[i]);
        free(locations);

	      self->g_n_vert = row;

        if (!self->g_suppress_output){
            printf("\nExtracted edges\n");
        }

     }
     else if (self->g_filetype == 2){
          
        //List of edges with lengths
        //Format is v1, v2, length
        if (!self->g_suppress_output){
            printf("extracting edges");
        }

        int n_edges = 0;

        row = 0;

        //while(getline(&line, &len, fp) != -1)
        //      n_edges++;

        //printf("\nnumber of edges %d", n_edges);
        //rewind(fp);

        int vv1, vv2;
	      int max_v = 0;

        EDGE_ID edge_list_ptr = 0;

        while(getline(&line, &len, fp) != -1) {

            col = 0;
            
            //self->g_edges_list[self->g_n_valid_edges] = (EDGE_ID*)malloc(2*sizeof(EDGE_ID));

            dist = strtok(line, " ,");
            while(dist != NULL){

                if (col == 0){

                    vv1 = atoi(dist);
		    if (vv1 > max_v)
			max_v = vv1;

                }
                else if (col == 1){

                    vv2 = atoi(dist);
		    if (vv2 > max_v)
			max_v = vv2;
                    
                }
                else if (col == 2){


                    PAR edge_length = strtod(dist, &end);
                    //if (simplex1_check(vv1, vv2, edge_length, self->g_thresh)){
                    if (edge_length < self->g_thresh){
                        self->g_edge_parameter[self->g_n_valid_edges] = edge_length;

                        if (vv1 < vv2){

                            //self->g_edges_list[self->g_n_valid_edges][0] = vv1;
                            //self->g_edges_list[self->g_n_valid_edges][1] = vv2;
                            self->g_edges_list[edge_list_ptr++] = vv1;
                            self->g_edges_list[edge_list_ptr++] = vv2;

                        }
                        else {

                            //self->g_edges_list[self->g_n_valid_edges][0] = vv2;
                            //self->g_edges_list[self->g_n_valid_edges][1] = vv1;
                            self->g_edges_list[edge_list_ptr++] = vv2;
                            self->g_edges_list[edge_list_ptr++] = vv1;

                        }

            	          self->g_n_valid_edges++;

            	          if (self->g_n_valid_edges == prealloc){
            	                
            	                prealloc += 100000;
            	                self->g_edges_list = (EDGE_ID*)realloc(self->g_edges_list, 2*prealloc*sizeof(EDGE_ID));
            	                self->g_edge_parameter = (PAR*)realloc(self->g_edge_parameter, prealloc*sizeof(PAR));

            	          }

                    }

                }

                dist = strtok(NULL, ",");
                col++;

            }

            row++;

        }

     	self->g_n_vert = max_v+1;
        if (!self->g_suppress_output){
            printf("\nExtracted edges\n");
        }

     }



     self->g_edges_list = (EDGE_ID*)realloc(self->g_edges_list, 2*self->g_n_valid_edges*sizeof(EDGE_ID));

     self->g_edge_parameter = (PAR*)realloc(self->g_edge_parameter, self->g_n_valid_edges*sizeof(PAR));
 
     fclose(fp);
     free(line);    

     if (!self->g_suppress_output){
        printf("\nNumber of vertices %d", self->g_n_vert);
     }

     mergeSort(self->g_edge_parameter, self->g_edges_list, 0, self->g_n_valid_edges-1);

     if (!self->g_suppress_output){
        printf("\nSorted %d edges\n", self->g_n_valid_edges);
     }
     //exit(0);


     clock_gettime(CLOCK_MONOTONIC, &finish_wall_clock);
     self->g_timer_process_input = (finish_wall_clock.tv_sec - start_wall_clock.tv_sec);
     self->g_timer_process_input += (finish_wall_clock.tv_nsec - start_wall_clock.tv_nsec) / 1000000000.0;


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
//                           STEP 1
//                  Generate Neighbor matrices
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

     //printf("\nPress key to start...");
     //getchar();

     clock_gettime(CLOCK_MONOTONIC, &start_wall_clock);
  
     // Initiate the Neighbor data structures
     self->g_Neighbors_e = (Neighbors**)malloc(self->g_n_vert*sizeof(Neighbors*));
     self->g_Neighbors = (Neighbors**)malloc(self->g_n_vert*sizeof(Neighbors*));
     self->g_Neigh_len = (VERT_ID*)calloc(self->g_n_vert, sizeof(VERT_ID));

     self->g_pivots_H0 = (EDGE_ID*)calloc(self->g_n_vert, sizeof(EDGE_ID));


     EDGE_ID* n_neigh = (EDGE_ID*)calloc(self->g_n_vert, sizeof(EDGE_ID));

     VERT_ID vv;

     for (EDGE_ID i = 0; i < self->g_n_valid_edges; i++){

          n_neigh[self->g_edges_list[2*i]]++;
          n_neigh[self->g_edges_list[(2*i)+1]]++;
          
          
     }

     for (VERT_ID i = 0; i < self->g_n_vert; i++){
          
          self->g_Neighbors_e[i] = (Neighbors*)malloc(n_neigh[i]*sizeof(Neighbors));
          self->g_Neighbors[i] = (Neighbors*)malloc(n_neigh[i]*sizeof(Neighbors));

     }

     free(n_neigh);
     

 
     if (!self->g_suppress_output){
        printf("\nCreating neighbors...");
     }
    
     //double time_create_neigh = omp_get_wtime();
     self->g_max_neighbors = 0;

     for (EDGE_ID i = 0; i < self->g_n_valid_edges; i++){
        
          VERT_ID v1 = self->g_edges_list[2*i];
          VERT_ID v2 = self->g_edges_list[(2*i)+1];

          len = self->g_Neigh_len[v1];

          self->g_Neighbors[v1][len].order = i;
          self->g_Neighbors[v1][len].neighbor = v2;

          self->g_Neighbors_e[v1][len].order = i;
          self->g_Neighbors_e[v1][len].neighbor = v2;

          self->g_Neigh_len[v1]++;

          len = self->g_Neigh_len[v2];

          self->g_Neighbors[v2][len].order = i;
          self->g_Neighbors[v2][len].neighbor = v1;

          self->g_Neighbors_e[v2][len].order = i;
          self->g_Neighbors_e[v2][len].neighbor = v1;

          self->g_Neigh_len[v2]++;


     }

     //printf("Time taken %f", omp_get_wtime() - time_create_neigh);

     if (!self->g_suppress_output){
        printf("\nSorting neighbors...");
     }
     //double time_sort_neigh = omp_get_wtime();

     #pragma omp parallel for schedule(static) shared(self)
     for (EDGE_ID i = 0; i < self->g_n_vert; i++){

          //self->g_Neighbors[i] = (Neighbors*)realloc(self->g_Neighbors[i],\
          //                                      self->g_Neigh_len[i]*sizeof(Neighbors));

          //self->g_Neighbors_e[i] = (Neighbors*)realloc(self->g_Neighbors_e[i],\
          //                                      self->g_Neigh_len[i]*sizeof(Neighbors));

          if (self->g_Neigh_len[i] > 1){

              sorter_tim_sort(self->g_Neighbors[i], self->g_Neigh_len[i]);
              sorter2_tim_sort(self->g_Neighbors_e[i], self->g_Neigh_len[i]);

          }
          
     }


#ifdef COMBIDX
     self->g_n_edges = (EDGE_ID)((self->g_n_vert) * (self->g_n_vert-1))/2;

     self->g_edges_comb_idx = (EDGE_ID*)malloc(self->g_n_edges*sizeof(EDGE_ID));

     for (EDGE_ID mm = 0; mm < self->g_n_edges; mm++){

          self->g_edges_comb_idx[mm] = self->g_n_valid_edges;
          
     }

     for (EDGE_ID mm = 0; mm < self->g_n_valid_edges; mm++){

          
          EDGE_ID idx = COMB_IDX0(self->g_edges_list[2*mm], self->g_edges_list[2*mm+1]);

          self->g_edges_comb_idx[idx] = mm;
#ifdef DEBUGCOMBIDX

          VERT_ID idx2 = search_Neighbors(self\
                                          , self->g_edges_list[2*mm]\
                                          , self->g_edges_list[2*mm+1]\
                                          , 0\
                                          , self->g_Neigh_len[self->g_edges_list[2*mm]] - 1);
          if (idx2 == self->g_n_vert){

              if (idx != self->g_n_valid_edges){
                  printf("\nERRRRROR 0");
                  getchar();
              }

          }
          else{
              if (self->g_Neighbors[self->g_edges_list[2*mm]][idx2].order != self->g_edges_comb_idx[idx]){
                  printf("\nERRRRROR 1");
                  getchar();

              }
          }

#endif

     }
#endif

     clock_gettime(CLOCK_MONOTONIC, &finish_wall_clock);
     self->g_timer_neigh = (finish_wall_clock.tv_sec - start_wall_clock.tv_sec);
     self->g_timer_neigh += (finish_wall_clock.tv_nsec - start_wall_clock.tv_nsec) / 1000000000.0;

     //printf("Time taken %f", omp_get_wtime() - time_sort_neigh);

     //compute_num_simplices(self);
     //exit(1);

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
//            STEP H0.1: Reduce the edges using column method
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

     if (!self->g_suppress_output){
     printf("\n\n---------------");
     printf("\nComputing H0...");
     printf("\n---------------\n");
     }

     clock_gettime(CLOCK_MONOTONIC, &start_wall_clock);


     // R Sparse
     self->g_R_sparse_max_H0 = 1000;
     self->g_R_sparse_H0 = (EDGE_ID*)malloc(self->g_R_sparse_max_H0*sizeof(EDGE_ID));
     self->g_R_sparse_ptr_H0 = 0;


     // R sparse col mapping
     self->g_R_col_indices_max_H0 = 100;
     self->g_R_col_indices_H0 = (EDGE_ID*)malloc(self->g_R_col_indices_max_H0*sizeof(EDGE_ID));
     self->g_R_col_indices_ptr_H0 = 1;

     
     // Note which edges have pivots in H0
     self->g_edges_with_pivots_H0 = \
                            (EDGE_ID*)calloc(self->g_n_valid_edges, sizeof(EDGE_ID));


//#ifdef HOM_CYCLES
     // For birth cycles
     if (self->g_compute_cycles){
          self->g_H0_pivot_of = (V_H0*)malloc(self->g_n_vert*sizeof(V_H0));
     }
//#endif


     /////////////
     // WORKSPACE
     /////////////
     self->g_ws_pre_alloc = 100;
     self->g_workspace_size = 1000;

     // H0 workspace structures
     self->g_R_ws_H0 = \
                     (EDGE_ID**)malloc(self->g_workspace_size*sizeof(EDGE_ID*));

     // H0 workspace info
     self->g_R_ws_H0_info = (boundary_H0_ws*)malloc(self->g_workspace_size*sizeof(boundary_H0_ws));


     // Initialize ws counter
     self->g_ws_counter = 0;

     for (int ws_counter = 0; ws_counter < self->g_workspace_size; ws_counter++){

         self->g_R_ws_H0_info[ws_counter].max_len = self->g_ws_pre_alloc;

         self->g_R_ws_H0[ws_counter] = (EDGE_ID*)malloc(2*self->g_R_ws_H0_info[ws_counter].max_len*sizeof(EDGE_ID));
         
     }


     ////////////////////////////////
     // Allocate jobs for parallel H0
     ////////////////////////////////
     
     self->g_jobs = (int*)malloc((self->g_cpu_count + 1)*sizeof(int));

     allocate_jobs(self, self->g_workspace_size);

     int rtn;

     self->g_threads = (pthread_t *)malloc(self->g_cpu_count*sizeof(pthread_t));

     if ((rtn = pthread_mutex_init(&(self->g_thread_lock), NULL)) !=0)
        fprintf(stderr, "pthread_mutex_init %s", strerror(rtn)), exit(-1);

     if ((rtn = pthread_cond_init(&(self->g_start_boss), NULL)) !=0)
        fprintf(stderr, "pthread_cond_init %s", strerror(rtn)), exit(-1);

     if ((rtn = pthread_cond_init(&(self->g_start_workers), NULL)) !=0)
        fprintf(stderr, "pthread_cond_init %s", strerror(rtn)), exit(-1);


     // Initialize thread creation
     self->g_thread_id = 0;
     self->g_sleeping_threads = 0;
     self->g_delete_threads = 0;

     for (int i = 0; i < self->g_cpu_count; i++){

        if ((rtn = pthread_create( \
                                &(self->g_threads[i]) \
                                , NULL \
                                , reduce_with_complex_H0 \
                                , (void*)self)!= 0))
          fprintf(stderr, "pthread_create %d", rtn), exit(-1);
      
     }

     // Wait for threads to be initialized
     pthread_mutex_lock(&(self->g_thread_lock));

     while(self->g_sleeping_threads != self->g_cpu_count){
        
          pthread_cond_wait(&(self->g_start_boss) \
                          , &(self->g_thread_lock));

     }

     ////////////////////////////////

     ////////////////////////////////
     // Main H0 Homology loop
     ////////////////////////////////
     
     for (EDGE_ID i = 0; i < self->g_n_valid_edges; i++){

               //printf("Percentage %f\r", (float)i/(float)self->g_n_valid_edges);
               //

               //if (i%10000 == 0){
               //    printf("\rProcessing edge %d", i);
               //}

               ////////////////////
               // Append to workspace_H0
               ////////////////////
                //self->g_ws_simplices_H0[self->g_ws_counter] = i;
                
                boundary_H0_ws* this_ws = self->g_R_ws_H0_info + self->g_ws_counter;

                // coboundary
                this_ws->cob = i;

                // Initially, the original is at 0
                this_ws->original = 0;

                // Length
                this_ws->len = 2;

                // Non empty
                this_ws->flag_non_empty = 1;
                
                // Recall: edge_list has v_max at 1 and v_min at 0
                self->g_R_ws_H0[self->g_ws_counter][0] = self->g_edges_list[2*i];
                self->g_R_ws_H0[self->g_ws_counter][1] = self->g_edges_list[2*i+1];
                
                // Pivot
                this_ws->pivot = self->g_edges_list[2*i+1];
                     
                self->g_ws_counter += 1;

                if (self->g_ws_counter == self->g_workspace_size){

                     reduce_ws_H0(self);
                     
                
                }



     }


     // Reduction of final batch
     while (self->g_ws_counter){

          // Allocate the last batch of size g_ws_counter
          allocate_jobs(self, self->g_ws_counter);

          reduce_ws_H0(self);

     }


     self->g_R_sparse_H0 = (EDGE_ID*)realloc( \
                                    self->g_R_sparse_H0\
                                  , (self->g_R_sparse_ptr_H0+1)*sizeof(EDGE_ID));


     self->g_R_col_indices_H0 = (EDGE_ID*)realloc( \
                                   self->g_R_col_indices_H0 \
                                  , (self->g_R_col_indices_ptr_H0+1)*sizeof(EDGE_ID));
                         
     
     /////////////////////////
     // Cancel the threads
     /////////////////////////

     self->g_delete_threads = 1;

     pthread_cond_broadcast(&(self->g_start_workers));

     pthread_mutex_unlock(&(self->g_thread_lock));

     for (int i = 0; i < self->g_cpu_count; i++){

        pthread_join(self->g_threads[i], NULL);
      
     }

     free(self->g_threads);
     free(self->g_jobs);

     //////////////////////////////////////////////////
     // Clear H0 parallel workspace
     //////////////////////////////////////////////////

     for (int ws_counter = 0; ws_counter < self->g_workspace_size; ws_counter++){
          
          free(self->g_R_ws_H0[ws_counter]);

     }

     free(self->g_R_ws_H0);
     free(self->g_R_ws_H0_info);
    
     /////////////////////////
     // Write H0 deaths to file
     /////////////////////////

     //// BINARY FILE
     //FILE* fp2 = fopen("H0_pers_pairs.bin", "wb");
     //fwrite(self->g_H0_pers_pairs, sizeof(PAR),self->g_H0_pers_pairs_len, fp2);
     //fclose(fp2);

#ifdef SAVEPD
     // TEXT FILE
     FILE* fp2 = fopen(self->g_H0_pers_file, "w");

     if (self->g_filetype == 1){

            for (EDGE_ID it = 0; it < self->g_n_valid_edges; it++){
                  
                  if (self->g_edges_with_pivots_H0[it]){
                     fprintf(fp2, "%.12lf,", sqrt(self->g_edge_parameter[it]));
                  }
                  
            }

     }
     else{

              for (EDGE_ID it = 0; it < self->g_n_valid_edges; it++){
                    
                    if (self->g_edges_with_pivots_H0[it]){
                       fprintf(fp2, "%.12lf,", self->g_edge_parameter[it]);
                    }
                    
              }
     }

     fclose(fp2);

#endif

#ifdef PRINT
 
     if (!self->g_suppress_output){
     printf("\nPers pairs in dim 0");
     }
     if (self->g_filetype == 1){

            for (EDGE_ID it = 0; it < self->g_n_valid_edges; it++){
                  
                  if (self->g_edges_with_pivots_H0[it]){
                     printf("\n%.12lf,", sqrt(self->g_edge_parameter[it]));
                  }
                  
            }

     }
     else{

              for (EDGE_ID it = 0; it < self->g_n_valid_edges; it++){
                    
                    if (self->g_edges_with_pivots_H0[it]){
                       printf("\n%.12lf,", self->g_edge_parameter[it]);
                    }
                    
              }
     }
      
#endif


     clock_gettime(CLOCK_MONOTONIC, &finish_wall_clock);
     self->g_timer_H0 = (finish_wall_clock.tv_sec - start_wall_clock.tv_sec);
     self->g_timer_H0 += (finish_wall_clock.tv_nsec - start_wall_clock.tv_nsec) / 1000000000.0;
     


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
//            STEP coH1.1: Find cohomology now for the edges
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

     if (!self->g_suppress_output){
     printf("\n\n-----------------");
     printf("\nComputing coH1...");
     printf("\n-----------------\n");
     }
     //double time_compute_coH1 = omp_get_wtime();
     clock_gettime(CLOCK_MONOTONIC, &start_wall_clock);

     // V sparse 
     EDGE_ID pre_alloc = 1000;

     self->g_V_sparse_max = pre_alloc;
     self->g_V_sparse_H1 = (EDGE_ID*)malloc(self->g_V_sparse_max*sizeof(EDGE_ID));
     self->g_V_sparse_ptr = 1;
     self->g_V_sparse_beg_ptr = 1;
     self->g_V_sparse_end_ptr = 1;

     self->g_V_col_indices_max = pre_alloc;
     self->g_V_col_indices = (EDGE_ID*)malloc(self->g_V_col_indices_max*sizeof(EDGE_ID));
     self->g_V_col_indices_ptr = 1;


     ////////////////////////////////////////////////////
     // INITIALIZE WORKSPACE
     ////////////////////////////////////////////////////
     self->g_cohom_ws_size = 100;
     self->g_V_ws_H1 = (coboundary_H1_ws*)malloc(self->g_cohom_ws_size*sizeof(coboundary_H1_ws));
     for (EDGE_ID mm = 0; mm < self->g_cohom_ws_size; mm++){
          
          self->g_V_ws_H1[mm].max_len = 10;
          self->g_V_ws_H1[mm].last = 0;
          self->g_V_ws_H1[mm].keys1 = (implicit_keys1*)malloc(self->g_V_ws_H1[mm].max_len*sizeof(implicit_keys1));

          for (EDGE_ID nn = 0; nn < self->g_V_ws_H1[mm].max_len; nn++){
                
                self->g_V_ws_H1[mm].keys1[nn].max_len = 10;
                self->g_V_ws_H1[mm].keys1[nn].last = 0;
                self->g_V_ws_H1[mm].keys1[nn].flag_empty = 1;

                self->g_V_ws_H1[mm].keys1[nn].keys2 =\
                                        (implicit_keys2*)malloc(self->g_V_ws_H1[mm].keys1[nn].max_len*sizeof(implicit_keys2));

          }

          self->g_V_ws_H1[mm].v_edges.max_len = 100;
          self->g_V_ws_H1[mm].v_edges.last = 0;
          self->g_V_ws_H1[mm].v_edges.o_ab = (EDGE_ID*)malloc(self->g_V_ws_H1[mm].v_edges.max_len*sizeof(EDGE_ID));
          

     }

     ////////////////////////////////////////////////////


     // H1 pivots 

     self->g_H1_cohom_pivots = (H1_cohom_pivots**)malloc(self->g_n_valid_edges*sizeof(H1_cohom_pivots*));

     self->g_H1_cohom_pivots_len = (EDGE_ID*)calloc(self->g_n_valid_edges, sizeof(EDGE_ID));
     self->g_H1_cohom_pivots_max_len = (EDGE_ID*)calloc(self->g_n_valid_edges, sizeof(EDGE_ID));


     // H1 Pers pairs
     self->g_H1_pers_pairs_max_len = 1000;
     self->g_H1_pers_pairs_len = 0;
     self->g_H1_pers_pairs = (PAR*)malloc(self->g_H1_pers_pairs_max_len*sizeof(PAR));

//#ifdef HOM_CYCLES
     if (self->g_compute_cycles){
          self->g_H1_undead_ptr = 0;
          self->g_H1_undead_max = 10;
          self->g_H1_undead = (EDGE_ID*)malloc(self->g_H1_undead_max*sizeof(EDGE_ID));
     }
//#endif


     int new_debug = 0;

     self->g_coH1_all_lows = (coboundary_H1*)malloc(self->g_n_valid_edges*sizeof(coboundary_H1));


     ////////////////////////////////////////////////////////////////
     // Allocate jobs/threads for parallel coH1
     ////////////////////////////////////////////////////////////////
     
     self->g_jobs = (int*)malloc((self->g_cpu_count + 1)*sizeof(int));

     allocate_jobs(self, self->g_cohom_ws_size);

     self->g_threads = (pthread_t *)malloc(self->g_cpu_count*sizeof(pthread_t));

     if ((rtn = pthread_mutex_init(&(self->g_thread_lock), NULL)) !=0)
        fprintf(stderr, "pthread_mutex_init %s", strerror(rtn)), exit(-1);

     if ((rtn = pthread_cond_init(&(self->g_start_boss), NULL)) !=0)
        fprintf(stderr, "pthread_cond_init %s", strerror(rtn)), exit(-1);

     if ((rtn = pthread_cond_init(&(self->g_start_workers), NULL)) !=0)
        fprintf(stderr, "pthread_cond_init %s", strerror(rtn)), exit(-1);


     // Initialize thread creation
     self->g_thread_id = 0;
     self->g_sleeping_threads = 0;
     self->g_delete_threads = 0;

     for (int i = 0; i < self->g_cpu_count; i++){

        if ((rtn = pthread_create( \
                                &(self->g_threads[i]) \
                                , NULL \
                                , reduce_with_complex_coH1 \
                                , (void*)self)!= 0))
          fprintf(stderr, "pthread_create %d", rtn), exit(-1);
      
     }

     // Wait for threads to be initialized
     pthread_mutex_lock(&(self->g_thread_lock));

     while(self->g_sleeping_threads != self->g_cpu_count){
        
          pthread_cond_wait(&(self->g_start_boss) \
                          , &(self->g_thread_lock));

     }

    ////////////////////////////////


     #pragma omp parallel for schedule(static) shared(self)
     for (EDGE_ID mm = 0; mm < self->g_n_valid_edges; mm++) {

          self->g_coH1_all_lows[mm].o_ab = mm; 
          find_H1_cohom_low(self, &(self->g_coH1_all_lows[mm]));
          // Need to find a_ptr and b_ptr if first low.key1 > e
          if (self->g_coH1_all_lows[mm].low.key1 > self->g_coH1_all_lows[mm].o_ab){

              VERT_ID a = self->g_edges_list[2*self->g_coH1_all_lows[mm].o_ab];
              VERT_ID b = self->g_edges_list[2*self->g_coH1_all_lows[mm].o_ab+1];

              self->g_coH1_all_lows[mm].a_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[a], 0, self->g_Neigh_len[a]-1\
                                                      , self->g_coH1_all_lows[mm].low.key1, self->g_Neigh_len[a]);
              
              self->g_coH1_all_lows[mm].b_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[b], 0, self->g_Neigh_len[b]-1\
                                                      , self->g_coH1_all_lows[mm].low.key1, self->g_Neigh_len[b]);
          }

     }



     self->g_this_edge = self->g_n_valid_edges;

     ///////////////////////////////////////////////////
     // MAIN coH1 loop 
     ///////////////////////////////////////////////////

     //getchar();

    
     
     self->g_new_debug = 0;
     self->g_new_debug2 = 0;

     self->g_ws_counter = 0;

     //self->g_debug_edge = self->g_n_valid_edges;
     self->g_debug_edge = 307605;

     while(self->g_this_edge){
          
          self->g_this_edge--;

          //if (self->g_this_edge%10000 == 0){
          //  printf("\nProcessing edge %d", self->g_this_edge);
          //}

          ///////////////////////////////////////////////////
          // CLEARING ALGORITHM
          // Does this edge have a pivot?
          ///////////////////////////////////////////////////
          

          if (self->g_edges_with_pivots_H0[self->g_this_edge]){
            //This edge has a pivot in H0. So, skip it. Continue;
            //skip++;
            
#ifdef COH1DEBUG
              if (self->g_this_edge == self->g_debug_edge ){
                printf("\nskipping because cleared. so, cannot have anything relevant");
              }
#endif
            continue;
          }
          ///////////////////////////////////////////////////
          ///////////////////////////////////////////////////


          if (self->g_coH1_all_lows[self->g_this_edge].low.key1 == self->g_n_valid_edges){
            // This edge has no coboundary
            //if (self->g_new_debug){
            //  printf("\nno cob, skipping");
            //}

              // Add this as undead in H1
            
              if (self->g_H1_pers_pairs_len+2 == self->g_H1_pers_pairs_max_len){
                    self->g_H1_pers_pairs_max_len += 1000;
                    self->g_H1_pers_pairs = (PAR*)realloc(self->g_H1_pers_pairs\
                                            , self->g_H1_pers_pairs_max_len*sizeof(PAR));
              
              }
              self->g_H1_pers_pairs[self->g_H1_pers_pairs_len++] = \
                                                   self->g_edge_parameter[self->g_this_edge];

              self->g_H1_pers_pairs[self->g_H1_pers_pairs_len++] = -1;


//#ifdef HOM_CYCLES
              if (self->g_compute_cycles){
                self->g_H1_undead[self->g_H1_undead_ptr++] = self->g_this_edge;
                if (self->g_H1_undead_ptr == self->g_H1_undead_max){
                    self->g_H1_undead_max += 100;
                    self->g_H1_undead = (EDGE_ID*)realloc(self->g_H1_undead\
                                                          , self->g_H1_undead_max*sizeof(EDGE_ID));
                }
              }
//#endif



#ifdef COH1DEBUG
              if (self->g_this_edge == self->g_debug_edge ){
                printf("\nskipping because has no cob");
              }
#endif
            continue;
          }


          // This is a trivial pair
          if (self->g_coH1_all_lows[self->g_this_edge].low.key1 == self->g_this_edge){

#ifdef COH1DEBUG
              if (self->g_this_edge == self->g_debug_edge ){
                printf("\nis a trivial pers pair. dont have to add pivot.");
              }
#endif
            // I DO NOT KNOW WHY I HAVE THIS HERE
            //self->g_edges_with_pivots_H0[self->g_this_edge]  = 10;

            continue;

          }

#ifdef COH1DEBUG
          if (self->g_this_edge == self->g_debug_edge ){
              printf("\nhave to start reduction");
              self->g_new_debug = 1;
              self->g_new_debug2 = 1;
          }
#endif
          coboundary_H1_ws* this_ws = self->g_V_ws_H1 + self->g_ws_counter;


          this_ws->edge = self->g_this_edge;

          this_ws->pivot = self->g_coH1_all_lows[this_ws->edge].low;

          this_ws->flag_first = 1;
          this_ws->flag_red_w_complex = 0;
          this_ws->flag_red_w_trivial = 0;
          this_ws->flag_append_to_complex = 0;

          this_ws->flag_non_empty = 1;

          
          // FIRST ENTRY IN hash-table

          this_ws->k1_ptr = 0;
          this_ws->k2_ptr = 0;
          this_ws->last = 1;

          this_ws->keys1[0].last = 1;
          this_ws->keys1[0].flag_empty = 0;
          this_ws->keys1[0].k1 = this_ws->pivot.key1;

          this_ws->keys1[0].keys2[0].k2 = this_ws->pivot.key2;
          this_ws->keys1[0].keys2[0].o_ab = self->g_coH1_all_lows[this_ws->edge].o_ab;
          this_ws->keys1[0].keys2[0].a_ptr = self->g_coH1_all_lows[this_ws->edge].a_ptr;
          this_ws->keys1[0].keys2[0].b_ptr = self->g_coH1_all_lows[this_ws->edge].b_ptr;
          this_ws->keys1[0].keys2[0].flag_next = 1;

          this_ws->v_edges.last = 0;

          self->g_ws_counter++;

          if (self->g_ws_counter == self->g_cohom_ws_size){
                reduce_ws_coH1(self);
          }


    }

    while(self->g_ws_counter){

        allocate_jobs(self, self->g_ws_counter);
        reduce_ws_coH1(self);
        
    }


     /////////////////////////
     // Cancel the threads used in getting next during reduction
     /////////////////////////

     self->g_delete_threads = 1;

     pthread_cond_broadcast(&(self->g_start_workers));

     pthread_mutex_unlock(&(self->g_thread_lock));

     for (int i = 0; i < self->g_cpu_count; i++){

        pthread_join(self->g_threads[i], NULL);
      
     }

     free(self->g_threads);
     free(self->g_jobs);

     if (!self->g_suppress_output){
     printf("\nsparse V coH1 length %d", self->g_V_sparse_ptr);
     }


    self->g_H1_pers_pairs = (PAR*)realloc(self->g_H1_pers_pairs, self->g_H1_pers_pairs_len*sizeof(PAR));

    //// BINARY FILE
    //FILE* fp2 = fopen("H1_pers_pairs.bin", "wb");
    //fwrite(self->g_H1_pers_pairs, sizeof(PAR),self->g_H1_pers_pairs_len, fp2);
    //fclose(fp2);

#ifdef SAVEPD
    // TEXT FILE
    fp2 = fopen(self->g_H1_pers_file, "w");

    PAR ddeath;

    if (self->g_filetype == 1){

        for (EDGE_ID it = 0; it < self->g_H1_pers_pairs_len; it+=2){

              if (self->g_H1_pers_pairs[it+1] == -1){
                      ddeath = -1;
              }
              else{
                      ddeath = sqrt(self->g_H1_pers_pairs[it+1]);
              }
              
              fprintf(fp2, "%0.12lf, %0.12lf\n", sqrt(self->g_H1_pers_pairs[it]), ddeath);
              
        }

    }
    else{

        for (EDGE_ID it = 0; it < self->g_H1_pers_pairs_len; it+=2){
              
              fprintf(fp2, "%0.12lf, %0.12lf\n", self->g_H1_pers_pairs[it], self->g_H1_pers_pairs[it+1]);
              
        }

    }

    fclose(fp2);
#endif


#ifdef PRINT
    // TEXT FILE

     if (!self->g_suppress_output){
     printf("\nPers pairs in dim 1");
     }
    if (self->g_filetype == 1){

        for (EDGE_ID it = 0; it < self->g_H1_pers_pairs_len; it+=2){
              
              printf("\n%0.12lf, %0.12lf", sqrt(self->g_H1_pers_pairs[it]), sqrt(self->g_H1_pers_pairs[it+1]));
              
        }

    }
    else{

        for (EDGE_ID it = 0; it < self->g_H1_pers_pairs_len; it+=2){
              
              printf("\n%0.12lf, %0.12lf", self->g_H1_pers_pairs[it], self->g_H1_pers_pairs[it+1]);
              
        }

    }

#endif



#ifdef SAVEV
    // TEXT FILE
    fp2 = fopen(self->g_coH1_V_file, "w");

        for (EDGE_ID it = 1; it < self->g_V_sparse_max; it++){
              
              fprintf(fp2, "%d\n", self->g_V_sparse_H1[it]);

    }

    fclose(fp2);
#endif


    // FREE coH1 Workspace
    for (int bb = 0; bb < self->g_cohom_ws_size; bb++){
         
          for (int mm = 0; mm < self->g_V_ws_H1[bb].max_len; mm++){
                
                free(self->g_V_ws_H1[bb].keys1[mm].keys2);
                
          }

          free(self->g_V_ws_H1[bb].keys1);
          free(self->g_V_ws_H1[bb].v_edges.o_ab);

    }

    free(self->g_V_ws_H1);

    // FREE V_sparse
    free(self->g_V_sparse_H1);



    clock_gettime(CLOCK_MONOTONIC, &finish_wall_clock);
    self->g_timer_coH1 = (finish_wall_clock.tv_sec - start_wall_clock.tv_sec);
    self->g_timer_coH1 += (finish_wall_clock.tv_nsec - start_wall_clock.tv_nsec) / 1000000000.0;
    //printf("Time: %lf\n", elapsed_wall_clock);

    /////////////////////////////////////////
    // Computing H1 homology and birth cycles
    /////////////////////////////////////////
    clock_gettime(CLOCK_MONOTONIC, &start_wall_clock);

    self->g_timer_computeH1 = 0;
    self->g_timer_H1cycles = 0;
    self->g_timer_minimize_H1cycles = 0;

//#ifdef HOM_CYCLES
    if (self->g_compute_cycles) compute_H1_homology_cycles(self);
//#endif


    if (self->g_dim_lim == 1){

      if (!self->g_suppress_output){

          printf("\nTime to process input : %lf" , self->g_timer_process_input);
          printf("\nTime to create neigh: %lf" , self->g_timer_neigh);
          printf("\nTime to compute H0: %lf"   , self->g_timer_H0);
          printf("\nTime to compute coH1: %lf" , self->g_timer_coH1);
//#ifdef HOM_CYCLES
          if (self->g_compute_cycles){

              printf("\nTime to compute H1: %lf" , self->g_timer_computeH1);
              printf("\nTime to compute %llu H1 cycles: %lf" , self->g_n_H1_birth_cycles, self->g_timer_H1cycles);
              printf("\nStored V_H0 %llu" , self->g_n_H0_stored_V);

          }
//#endif

//#ifdef MINIMIZE_BIRTH_CYCLES
          if (self->g_reduce_cyc_lengths) printf("\nTime to minimize H1 birth cycles: %lf" , self->g_timer_minimize_H1cycles);
//#endif

//#ifdef MINIMIZE_HOM_CYCLES
//          printf("\nTime to minimize H1 hom cycles: %lf" , self->g_timer_minimize_H1_homcycles);
//#endif

          printf("\nTotal time taken: %lf", \
                                          self->g_timer_process_input\
                                        + self->g_timer_neigh\
                                        + self->g_timer_H0\
                                        + self->g_timer_coH1\
                                        + self->g_timer_computeH1\
                                        + self->g_timer_H1cycles\
                                        + self->g_timer_minimize_H1cycles\
                                        + self->g_timer_minimize_H1_homcycles\
                                        );
      }

      deallocator(self);

      if (!self->g_suppress_output){
          printf("\nQuitting after coH1");
      }

      //Py_RETURN_NONE;

      
      exit(0);
    }

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
//            STEP coH2.1: Find cohomology now for the triangles
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


     clock_gettime(CLOCK_MONOTONIC, &start_wall_clock);

     if (!self->g_suppress_output){
     printf("\n\n--------------");
     printf("\nComputing coH2");
     printf("\n--------------");
     }


    // sparse V coH2
    self->g_V_sparse_max = 100000;
    self->g_V_sparse_H2 = (simplex*)malloc(self->g_V_sparse_max*sizeof(simplex));
    self->g_V_sparse_ptr = 1;
    self->g_V_sparse_beg_ptr = 1;
    self->g_V_sparse_end_ptr = 1;

    self->g_V_col_indices_ptr = 1;



     ////////////////////////////////////////////////////
     // INITIALIZE WORKSPACE
     ////////////////////////////////////////////////////
     self->g_cohom_ws_size = 100;
     self->g_V_ws_H2 = (coboundary_H2_ws*)malloc(self->g_cohom_ws_size*sizeof(coboundary_H2_ws));
     for (EDGE_ID mm = 0; mm < self->g_cohom_ws_size; mm++){
          
          self->g_V_ws_H2[mm].max_len = 10;
          self->g_V_ws_H2[mm].last = 0;
          self->g_V_ws_H2[mm].keys1 = (coH2_implicit_keys1*)malloc(self->g_V_ws_H2[mm].max_len*sizeof(coH2_implicit_keys1));

          for (EDGE_ID nn = 0; nn < self->g_V_ws_H2[mm].max_len; nn++){
                
                self->g_V_ws_H2[mm].keys1[nn].max_len = 10;
                self->g_V_ws_H2[mm].keys1[nn].last = 0;
                self->g_V_ws_H2[mm].keys1[nn].flag_empty = 1;

                self->g_V_ws_H2[mm].keys1[nn].keys2 =\
                                        (coH2_implicit_keys2*)malloc(self->g_V_ws_H2[mm].keys1[nn].max_len*sizeof(coH2_implicit_keys2));

          }

          self->g_V_ws_H2[mm].v_triangles.max_len = 100;
          self->g_V_ws_H2[mm].v_triangles.last = 0;
          self->g_V_ws_H2[mm].v_triangles.o_abc = (simplex*)malloc(self->g_V_ws_H2[mm].v_triangles.max_len*sizeof(simplex));
          

     }
     ////////////////////////////////////////////////////



    // PIVOTS
    self->g_H2_cohom_pivots = (H2_cohom_pivots**)malloc(self->g_n_valid_edges*sizeof(H2_cohom_pivots*));

    self->g_H2_cohom_pivots_len = (EDGE_ID*)calloc(self->g_n_valid_edges, sizeof(EDGE_ID));
    self->g_H2_cohom_pivots_max_len = (EDGE_ID*)calloc(self->g_n_valid_edges, sizeof(EDGE_ID));


     // H1 Pers pairs
     self->g_H2_pers_pairs_max_len = 1000;
     self->g_H2_pers_pairs_len = 0;
     self->g_H2_pers_pairs = (PAR*)malloc(self->g_H2_pers_pairs_max_len*sizeof(PAR));



//#ifdef HOM_CYCLES
     if (self->g_compute_cycles){
          self->g_H2_undead_ptr = 0;
          self->g_H2_undead_max = 10;
          self->g_H2_undead = (simplex*)malloc(self->g_H2_undead_max*sizeof(simplex));
     }
//#endif


     ////////////////////////////////////////////////////////////////
     // Allocate jobs/threads for parallel coH2
     ////////////////////////////////////////////////////////////////
     
     self->g_jobs = (int*)malloc((self->g_cpu_count + 1)*sizeof(int));

     allocate_jobs(self, self->g_cohom_ws_size);

     self->g_threads = (pthread_t *)malloc(self->g_cpu_count*sizeof(pthread_t));

     if ((rtn = pthread_mutex_init(&(self->g_thread_lock), NULL)) !=0)
        fprintf(stderr, "pthread_mutex_init %s", strerror(rtn)), exit(-1);

     if ((rtn = pthread_cond_init(&(self->g_start_boss), NULL)) !=0)
        fprintf(stderr, "pthread_cond_init %s", strerror(rtn)), exit(-1);

     if ((rtn = pthread_cond_init(&(self->g_start_workers), NULL)) !=0)
        fprintf(stderr, "pthread_cond_init %s", strerror(rtn)), exit(-1);


     // Initialize thread creation
     self->g_thread_id = 0;
     self->g_sleeping_threads = 0;
     self->g_delete_threads = 0;

     for (int i = 0; i < self->g_cpu_count; i++){

        if ((rtn = pthread_create( \
                                &(self->g_threads[i]) \
                                , NULL \
                                , reduce_with_complex_coH2 \
                                , (void*)self)!= 0))
          fprintf(stderr, "pthread_create %d", rtn), exit(-1);
      
     }

     // Wait for threads to be initialized
     pthread_mutex_lock(&(self->g_thread_lock));

     while(self->g_sleeping_threads != self->g_cpu_count){
        
          pthread_cond_wait(&(self->g_start_boss) \
                          , &(self->g_thread_lock));

     }

    ////////////////////////////////



    // BUFFER
    EDGE_ID buffer_len = 1000000;
    int buffer_ptr = 0;
    coboundary_H2* coH2_lows_buffer = (coboundary_H2*)malloc(buffer_len*sizeof(coboundary_H2));
    


    ///////////////////////////////////////////////////
    // MAIN coH2 loop 
    ///////////////////////////////////////////////////
    //
    
    //EDGE_ID i = self->g_n_valid_edges;

    coboundary_H2* temp_temp_triangles = (coboundary_H2*)malloc(self->g_n_vert*sizeof(coboundary_H2));
    VERT_ID temp_temp_len = 0;

    coboundary_H2 temp_triangle;


    self->g_debug_triangle.key1 = self->g_n_valid_edges ;
    self->g_debug_triangle.key2 = self->g_n_valid_edges ;

    //self->g_debug_triangle.key1 = 46494 ;
    //self->g_debug_triangle.key2 = 269 ;

    VERT_ID a, b, c, a_ptr, b_ptr;

    EDGE_ID ac, bc, has_pivot;

    EDGE_ID ab = self->g_n_valid_edges;

    while(ab){

          ab--;

          if (!self->g_suppress_output){
            if (self->g_n_valid_edges > 1000){
                if (ab % (self->g_n_valid_edges/100) == 0){
                  printf("\rProcessing coH2 for edge %d out of %d", ab, self->g_n_valid_edges);
                }
            }
            else{
                printf("\rProcessing coH2 for edge %d out of %d", ab, self->g_n_valid_edges);
            }
          }

          a = self->g_edges_list[2*ab];
          b = self->g_edges_list[2*ab+1];

          // Find the faces which are created when this edge is formed
          // That means, the o_max will be ab

          //ab = i;

          a_ptr = 0;
          b_ptr = 0;


          while ((a_ptr < self->g_Neigh_len[a])\
                && (b_ptr < self->g_Neigh_len[b])){


                if (self->g_Neighbors[a][a_ptr].neighbor < self->g_Neighbors[b][b_ptr].neighbor)
                {
                      a_ptr++;
                }

                else if (self->g_Neighbors[a][a_ptr].neighbor > self->g_Neighbors[b][b_ptr].neighbor)
                {
                      b_ptr++;
                }

                else{

                      c = self->g_Neighbors[a][a_ptr].neighbor;

                      //if (!simplex2_check(a, b, c)) continue;

                      ac = self->g_Neighbors[a][a_ptr++].order;
                      bc = self->g_Neighbors[b][b_ptr++].order;

                      if ((ac > ab) \
                          || (bc > ab))
                        continue;

                      temp_triangle.triangle.key1 = ab;
                      temp_triangle.triangle.key2 = (EDGE_ID)c;

                      ///////////////////////////////////////////////////
                      // CLEARING ALGORITHM
                      // Does this triangle have a pivot in coH1?
                      ///////////////////////////////////////////////////
                      // Check whether the triangle is pivot of a trivial pair in coH1
                      if ((self->g_coH1_all_lows[temp_triangle.triangle.key1].low.key1 == temp_triangle.triangle.key1)\
                          &&(self->g_coH1_all_lows[temp_triangle.triangle.key1].low.key2 == temp_triangle.triangle.key2)){
                            //printf("\nSkipping");
                            continue;
                      }

                      // Check whether the triangle is a pivot in coH1
                      if (self->g_H1_cohom_pivots_len[temp_triangle.triangle.key1]){

                          has_pivot = search_H1_cohom_pivots(self->g_H1_cohom_pivots[temp_triangle.triangle.key1]\
                                        , 0 \
                                        , self->g_H1_cohom_pivots_len[temp_triangle.triangle.key1] - 1\
                                        , temp_triangle.triangle.key2 \
                                        , self->g_n_valid_edges);
                          
                          if (has_pivot != self->g_n_valid_edges){
                            //This triangle has a pivot in H1. So, skip it. Continue;
                                //printf("\nSkipping");
                                //getchar();
                                continue;
                          }

                      }
                      // END OF CLEARING ALGORITHM
                      ///////////////////////////////////////////////////

                      temp_temp_triangles[temp_temp_len++] = temp_triangle;

                }

          }
          

          while (temp_temp_len > 0){
              
                temp_temp_len--;
                coH2_lows_buffer[buffer_ptr++].triangle = temp_temp_triangles[temp_temp_len].triangle;

                if (buffer_ptr == buffer_len){

                      #pragma omp parallel for schedule(static) \
                                               shared(self, coH2_lows_buffer)
                      for (EDGE_ID mm = 0; mm < buffer_len; mm++) {
                           find_H2_cohom_low(self, &(coH2_lows_buffer[mm]));
                      }


                      EDGE_ID mm = 0;
                      while (mm < buffer_len){

                           if (coH2_lows_buffer[mm].vertex == -1){
                                //If it has empty cob, then it is undead cycle
                                      
                                //printf("\n Adding undead for H2");

                                if (self->g_H2_pers_pairs_len+2 == self->g_H2_pers_pairs_max_len){
                                      self->g_H2_pers_pairs_max_len += 1000;
                                      self->g_H2_pers_pairs = (PAR*)realloc(self->g_H2_pers_pairs\
                                                                    , self->g_H2_pers_pairs_max_len*sizeof(PAR));
                                
                                }
                                self->g_H2_pers_pairs[self->g_H2_pers_pairs_len++] =\
                                                                      self->g_edge_parameter[coH2_lows_buffer[mm].triangle.key1];
                                self->g_H2_pers_pairs[self->g_H2_pers_pairs_len++] = -1;

//#ifdef HOM_CYCLES
                                if (self->g_compute_cycles){
                                      self->g_H2_undead[self->g_H2_undead_ptr++] = coH2_lows_buffer[mm].triangle;
                                      if (self->g_H2_undead_ptr == self->g_H2_undead_max){
                                          self->g_H2_undead_max += 100;
                                          self->g_H2_undead = (simplex*)realloc(self->g_H2_undead\
                                                                                , self->g_H2_undead_max*sizeof(simplex));
                                      }
                                }
//#endif


                                mm++;
                                continue;
                           }

                           // Is this is a trivial pair?
                           if ((coH2_lows_buffer[mm].low.key1 == coH2_lows_buffer[mm].triangle.key1)\
                               &&(self->g_edges_list[2*coH2_lows_buffer[mm].low.key2+1]== coH2_lows_buffer[mm].triangle.key2)){

                                mm++;
                                continue;
                           }


                            coboundary_H2_ws* this_ws = self->g_V_ws_H2 + self->g_ws_counter;

                            this_ws->triangle = coH2_lows_buffer[mm].triangle;

                            this_ws->pivot = coH2_lows_buffer[mm].low;

                            this_ws->flag_first = 1;
                            this_ws->flag_red_w_complex = 0;
                            this_ws->flag_red_w_trivial = 0;
                            this_ws->flag_append_to_complex = 0;

                            this_ws->flag_non_empty = 1;

                            this_ws->k1_ptr = 0;
                            this_ws->k2_ptr = 0;
                            this_ws->last = 1;

                            this_ws->keys1[0].last = 1;
                            this_ws->keys1[0].flag_empty = 0;
                            this_ws->keys1[0].k1 = this_ws->pivot.key1;

                            this_ws->keys1[0].keys2[0].k2 = this_ws->pivot.key2;
                            this_ws->keys1[0].keys2[0].o_abc  = coH2_lows_buffer[mm].triangle;
                            this_ws->keys1[0].keys2[0].a_ptr  = coH2_lows_buffer[mm].a_ptr;
                            this_ws->keys1[0].keys2[0].b_ptr  = coH2_lows_buffer[mm].b_ptr;
                            this_ws->keys1[0].keys2[0].c_ptr  = coH2_lows_buffer[mm].c_ptr;
                            this_ws->keys1[0].keys2[0].vertex = coH2_lows_buffer[mm].vertex;
                            this_ws->keys1[0].keys2[0].flag_next = 1;

                            this_ws->v_triangles.last = 0;

                            self->g_ws_counter++;

                            if (self->g_ws_counter == self->g_cohom_ws_size){
                                  reduce_ws_coH2(self);
                            }

                            mm++;
                            
                      }
                      
                      
                      buffer_ptr = 0;

                }



          }

    }




    #pragma omp parallel for schedule(static) \
                             shared(self, coH2_lows_buffer)
    for (EDGE_ID mm = 0; mm < buffer_ptr; mm++) {

          find_H2_cohom_low(self, &(coH2_lows_buffer[mm]));

    }


    //for (EDGE_ID mm = 0; mm < buffer_ptr; mm++) {
    EDGE_ID mm = 0;
    while (mm < buffer_ptr){

          if (coH2_lows_buffer[mm].vertex == -1){
               mm++;
               continue;
          }
          
          // Is this is a trivial pair?
          if ((coH2_lows_buffer[mm].low.key1 == coH2_lows_buffer[mm].triangle.key1)\
              &&(self->g_edges_list[2*coH2_lows_buffer[mm].low.key2+1]== coH2_lows_buffer[mm].triangle.key2)){
               
               mm++;
               continue;
          }
          
          
          coboundary_H2_ws* this_ws = self->g_V_ws_H2 + self->g_ws_counter;
          
          this_ws->triangle = coH2_lows_buffer[mm].triangle;
          
          this_ws->pivot = coH2_lows_buffer[mm].low;
          
          this_ws->flag_first = 1;
          this_ws->flag_red_w_complex = 0;
          this_ws->flag_red_w_trivial = 0;
          this_ws->flag_append_to_complex = 0;
          
          this_ws->flag_non_empty = 1;
          
          this_ws->k1_ptr = 0;
          this_ws->k2_ptr = 0;
          this_ws->last = 1;
          
          this_ws->keys1[0].last = 1;
          this_ws->keys1[0].flag_empty = 0;
          this_ws->keys1[0].k1 = this_ws->pivot.key1;
          
          this_ws->keys1[0].keys2[0].k2 = this_ws->pivot.key2;
          this_ws->keys1[0].keys2[0].o_abc  = coH2_lows_buffer[mm].triangle;
          this_ws->keys1[0].keys2[0].a_ptr  = coH2_lows_buffer[mm].a_ptr;
          this_ws->keys1[0].keys2[0].b_ptr  = coH2_lows_buffer[mm].b_ptr;
          this_ws->keys1[0].keys2[0].c_ptr  = coH2_lows_buffer[mm].c_ptr;
          this_ws->keys1[0].keys2[0].vertex = coH2_lows_buffer[mm].vertex;
          this_ws->keys1[0].keys2[0].flag_next = 1;
          
          this_ws->v_triangles.last = 0;
          
          self->g_ws_counter++;
          
          if (self->g_ws_counter == self->g_cohom_ws_size){
                reduce_ws_coH2(self);
          }
          
          mm++;
          
    }


    while(self->g_ws_counter){

        allocate_jobs(self, self->g_ws_counter);
        reduce_ws_coH2(self);
        
    }




    //////////////////////////////////////////////////
    // Cancel the threads used in getting next during reduction
    //////////////////////////////////////////////////

    self->g_delete_threads = 1;

    pthread_cond_broadcast(&(self->g_start_workers));

    pthread_mutex_unlock(&(self->g_thread_lock));

    for (int i = 0; i < self->g_cpu_count; i++){

       pthread_join(self->g_threads[i], NULL);
     
    }

    free(self->g_threads);
    free(self->g_jobs);
    //printf("\nTime taken to compute coH1 %f", omp_get_wtime() - time_compute_coH1);

     if (!self->g_suppress_output){
     printf("\nsparse V coH2 length %d", self->g_V_sparse_ptr);
     }

    //////////////////////////////////////////////////
    // FREE coH2 Workspace
    //////////////////////////////////////////////////
    for (int bb = 0; bb < self->g_cohom_ws_size; bb++){
         
          for (int mm = 0; mm < self->g_V_ws_H2[bb].max_len; mm++){
                
                free(self->g_V_ws_H2[bb].keys1[mm].keys2);
                
          }

          free(self->g_V_ws_H2[bb].keys1);
          free(self->g_V_ws_H2[bb].v_triangles.o_abc);

    }

    free(self->g_V_ws_H2);

    // FREE temp_temp_triangles
    free(temp_temp_triangles);


    self->g_H2_pers_pairs = (PAR*)realloc(self->g_H2_pers_pairs, self->g_H2_pers_pairs_len*sizeof(PAR));

    // BINARY FILE
    //fp2 = fopen("H2_pers_pairs.bin", "wb");
    //fwrite(self->g_H2_pers_pairs, sizeof(PAR),self->g_H2_pers_pairs_len, fp2);
    //fclose(fp2);


#ifdef SAVEPD

    // TEXT FILE
    fp2 = fopen(self->g_H2_pers_file, "w");

    if (self->g_filetype == 1){

        for (EDGE_ID it = 0; it < self->g_H2_pers_pairs_len; it+=2){
              
              if (self->g_H2_pers_pairs[it+1] == -1){
                      ddeath = -1;
              }
              else{
                      ddeath = sqrt(self->g_H2_pers_pairs[it+1]);
              }

              fprintf(fp2, "%0.12lf, %0.12lf\n", sqrt(self->g_H2_pers_pairs[it]), ddeath);
              
        }

    }
    else{

        for (EDGE_ID it = 0; it < self->g_H2_pers_pairs_len; it+=2){
              
              fprintf(fp2, "%0.12lf, %0.12lf\n", self->g_H2_pers_pairs[it], self->g_H2_pers_pairs[it+1]);
              
        }

    }

    fclose(fp2);

#endif


#ifdef PRINT

    // TEXT FILE

     if (!self->g_suppress_output){
      printf("\nPers pairs in dim 2");
     }

    if (self->g_filetype == 1){

        for (EDGE_ID it = 0; it < self->g_H2_pers_pairs_len; it+=2){
              
              printf("\n%0.12lf, %0.12lf", sqrt(self->g_H2_pers_pairs[it]), sqrt(self->g_H2_pers_pairs[it+1]));
              
        }

    }
    else{

        for (EDGE_ID it = 0; it < self->g_H2_pers_pairs_len; it+=2){
              
              printf("\n%0.12lf, %0.12lf", self->g_H2_pers_pairs[it], self->g_H2_pers_pairs[it+1]);
              
        }

    }


#endif


#ifdef SAVEV
    // TEXT FILE
    fp2 = fopen(self->g_coH2_V_file, "w");

        for (EDGE_ID it = 1; it < self->g_V_sparse_max; it++){
              
              fprintf(fp2, "%d\n", self->g_V_sparse_H2[it]);

    }

    fclose(fp2);
#endif
    
    // FREE V_sparse
    free(self->g_V_sparse_H2);

    clock_gettime(CLOCK_MONOTONIC, &finish_wall_clock);
    self->g_timer_coH2 = (finish_wall_clock.tv_sec - start_wall_clock.tv_sec);
    self->g_timer_coH2 += (finish_wall_clock.tv_nsec - start_wall_clock.tv_nsec) / 1000000000.0;




    /////////////////////////////////////////
    // Computing H2 homology and birth cycles
    /////////////////////////////////////////

    self->g_timer_computeH2 = 0;
    self->g_timer_H2cycles = 0;
    self->g_timer_minimize_H2cycles = 0;

//#ifdef HOM_CYCLES
    if (self->g_compute_cycles) compute_H2_homology_cycles(self);
//#endif


     if (!self->g_suppress_output){
    printf("\nTime to process input : %lf" , self->g_timer_process_input);
    printf("\nTime to create neigh: %lf" , self->g_timer_neigh);
    printf("\nTime to compute H0: %lf"   , self->g_timer_H0);
    printf("\nTime to compute coH1: %lf" , self->g_timer_coH1);
    printf("\nTime to compute coH2: %lf" , self->g_timer_coH2);
     }
//#ifdef HOM_CYCLES
    if (self->g_compute_cycles){

     if (!self->g_suppress_output){
        printf("\nTime to compute H1: %lf" , self->g_timer_computeH1);
        printf("\nTime to compute %llu H1cycles: %lf" , self->g_n_H1_birth_cycles, self->g_timer_H1cycles);
        printf("\nStored V_H0 %llu" , self->g_n_H0_stored_V);
    
    //#ifdef MINIMIZE_BIRTH_CYCLES
        if (self->g_reduce_cyc_lengths){
            printf("\nTime to minimize H1 birth cycles: %lf" , self->g_timer_minimize_H1cycles);
        }
    //#endif
    
    //#ifdef MINIMIZE_HOM_CYCLES
    //      printf("\nTime to minimize H1 hom cycles: %lf" , self->g_timer_minimize_H1_homcycles);
    //#endif
    
        printf("\nTime to compute H2: %lf" , self->g_timer_computeH2);
        printf("\nTime to compute %llu H2cycles: %lf" , self->g_n_H2_birth_cycles, self->g_timer_H2cycles);
        printf("\nStored V_H1 %llu" , self->g_n_H1_stored_V);
    //#ifdef MINIMIZE_BIRTH_CYCLES
        if (self->g_reduce_cyc_lengths){
            printf("\nTime to minimize H2 cycles: %lf" , self->g_timer_minimize_H2cycles);
        }
     }
    //#endif
    //#ifdef MINIMIZE_HOM_CYCLES
        //printf("\nTime to minimize H2 hom cycles: %lf" , self->g_timer_minimize_H2_homcycles);
    //#endif
    }
//#endif
    //printf("Time to compute coH2 serial: %lf\n" , self->g_timer_coH2_serial);
    //printf("Time to compute coH2 parallel: %lf\n" , self->g_timer_coH2_parallel);
    if (!self->g_suppress_output){
    printf("\nTotal time taken: %lf",\
                                    self->g_timer_process_input\
                                    + self->g_timer_neigh\
                                    + self->g_timer_H0\
                                    + self->g_timer_coH1\
                                    + self->g_timer_coH2\
                                    + self->g_timer_computeH1\
                                    + self->g_timer_computeH2\
                                    + self->g_timer_H1cycles\
                                    + self->g_timer_H2cycles\
                                    + self->g_timer_minimize_H1cycles\
                                    + self->g_timer_minimize_H2cycles\
                                    + self->g_timer_minimize_H1_homcycles\
                                    + self->g_timer_minimize_H2_homcycles\
                                    );
     }
    //printf("Time in H2_low: %lf\n", self->g_timer_H2_low);
    //printf("Time in H2_greater: %lf\n", self->g_timer_H2_greater);
    //printf("Time in H2_next: %lf\n", self->g_timer_H2_next);



    deallocator(self);

     if (!self->g_suppress_output){
    printf("\nQuitting after coH2");
     }

    //Py_RETURN_NONE;
    


}

VERT_ID search_Neighbors(filtration* self, VERT_ID v1, VERT_ID v2, VERT_ID l, VERT_ID r){

    if (r >= l) { 
        VERT_ID mid = l + (r - l) / 2; 
  
        // If the element is present at the middle 
        // itself 
        if (self->g_Neighbors[v1][mid].neighbor == v2) 
            return mid; 
  
        // If element is smaller than mid, then 
        // it can only be present in left subarray 
        if (self->g_Neighbors[v1][mid].neighbor > v2){ 
            if (!mid) return self->g_n_vert;
            return search_Neighbors(self, v1, v2, l, mid - 1); 
        }
  
        // Else the element can only be present 
        // in right subarray 
        return search_Neighbors(self, v1, v2, mid + 1, r); 
    } 
  
    // We reach here when element is not 
    // present in array 
    return self->g_n_vert; 

}


VERT_ID search_Neighbors_e(filtration* self, VERT_ID v1, EDGE_ID order, VERT_ID l, VERT_ID r, EDGE_ID len){
  
    int mid = l + (r - l) / 2; 

    if (self->g_Neighbors_e[v1][mid].order < order){
        
        if (mid < len-1)
          if (self->g_Neighbors_e[v1][mid+1].order > order)
            return mid+1;

        return search_Neighbors_e(self, v1, order, mid+1, r, len); 
          
    }
    else if (self->g_Neighbors_e[v1][mid].order > order){

        if (!mid) return 0;
        return search_Neighbors_e(self, v1, order, l, mid-1, len); 

    }
    else{
        
        return mid+1;

    }

}

//////////////////////////////////////////////////////////
// MERGING ALGORITHMS
//////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////
// END OF MERGING ALGORITHMS
//////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////
// CUSTOM BOUNDARIES
//////////////////////////////////////////////////////////
int simplex1_check(VERT_ID v1, VERT_ID v2, PAR dist, PAR thresh){

  if (dist == -1){
    return 0;
  }
  //if (dist == 0){
  //  return 0;
  //}
  
  if (dist > thresh){
    return 0;
  }

  return 1;

}

int simplex2_check(VERT_ID v1, VERT_ID v2, VERT_ID v3){

  return 1;

}

int simplex3_check(VERT_ID v1, VERT_ID v2, VERT_ID v3, VERT_ID v4){

  return 1;

}
//////////////////////////////////////////////////////////



int compare_simplices(simplex* s1, simplex* s2){
    
  // returns 1 if s1 > s2
  // returns 0 if s1 < s2
  // returns -1 if s1 = s2
      
    if (s1->key1 > s2->key1) return 1;
    else if (s1->key1 < s2->key1) return 0;
    else{
      if (s1->key2 > s2->key2) return 1;
      else if (s1->key2 < s2->key2) return 0;
      else return -1;
    }
      
}

int compare_simplices_keys(EDGE_ID key11, EDGE_ID key12 \
                         , EDGE_ID key21, EDGE_ID key22 \
                          ){
    
  // returns 1 if s1 > s2
  // returns 0 if s1 < s2
  // returns -1 if s1 = s2
      
    if (key11 > key21) return 1;
    else if (key11 < key21) return 0;
    else{
      if (key12 > key22) return 1;
      else if (key12 < key22) return 0;
      else return -1;
    }
      
}



void find_H1_cohom_low(filtration* self, coboundary_H1* V_info){


    VERT_ID a = self->g_edges_list[2*V_info->o_ab];
    VERT_ID b = self->g_edges_list[2*V_info->o_ab+1];

    V_info->a_ptr = 0;
    V_info->b_ptr = 0;

    EDGE_ID o_min = self->g_n_valid_edges;
    EDGE_ID v_min = 0;

    EDGE_ID o_s, v_s;


    while(1){

      if ((V_info->a_ptr < self->g_Neigh_len[a]) \
          && (V_info->b_ptr < self->g_Neigh_len[b])){


            if (self->g_Neighbors[a][V_info->a_ptr].neighbor < self->g_Neighbors[b][V_info->b_ptr].neighbor){

                V_info->a_ptr++;

            }
            else if (self->g_Neighbors[a][V_info->a_ptr].neighbor > self->g_Neighbors[b][V_info->b_ptr].neighbor){

                V_info->b_ptr++;

            }
            else{

               EDGE_ID ac = self->g_Neighbors[a][V_info->a_ptr].order;

               EDGE_ID bc = self->g_Neighbors[b][V_info->b_ptr].order;

               EDGE_ID c = self->g_Neighbors[a][V_info->a_ptr].neighbor;

               //printf("\n a, b, c: %d, %d, %d", a, b, c);
               //printf("\n ab, ac, bc: %d, %d, %d", V_info->o_ab, ac, bc);

               o_s = ac;
               v_s = b;

               if (bc > ac){
                  o_s = bc;
                  v_s = a;
               }

               if (o_s < V_info->o_ab){
                    
                    V_info->low.key1 = V_info->o_ab;
                    V_info->low.key2 = c;
                    return;
               }

               // Reached here means o_s > o_ab

               //printf("\n o_s, v_s, : %d, %d", o_s, v_s);
               if (o_s < o_min){
                 o_min = o_s;
                 v_min = v_s;
               }
               else if ((o_s == o_min) && (v_s < v_min)){
                 v_min = v_s;
               }
               //printf("\n o_min, v_min, : %d, %d", o_min, v_min);

               V_info->a_ptr++;
               V_info->b_ptr++;

            }

      }
      else{

          //printf("\nReturning lowest of > a(%d)b(%d) (%d, %d)", o_s, v_s);
          V_info->low.key1 = o_min;
          V_info->low.key2 = v_min;
          return;
        
      }

    }

}








// A recursive binary search function. It returns 
// location of x in given array arr[l..r] is present, 
// otherwise -1 
EDGE_ID search_H1_cohom_pivots(H1_cohom_pivots* arr, EDGE_ID l, EDGE_ID r, EDGE_ID key2, EDGE_ID max) 
{ 
    if (r >= l) { 
        EDGE_ID mid = l + (r - l) / 2; 

        if (arr[mid].key2 == key2) 
            return mid; 
  
        // If element is smaller than mid, then 
        // it can only be present in left subarray 
        if (arr[mid].key2 > key2) 
        {

          /// PRECAUTIONARY: CAN REMOVE LATER
            if (!mid){
                return max; 
                printf("\nMID 0 WILL GIVE ERROR FOR UNSIGNED NEXT");
                getchar();
            }
          ///////////////////
            return search_H1_cohom_pivots(arr, l, mid - 1, key2, max); 
        }
  
        // Else the element can only be present 
        // in right subarray 
        return search_H1_cohom_pivots(arr, mid + 1, r, key2, max); 
    } 
  
    // We reach here when element is not 
    // present in array 
    //printf("\nNOT FOUND");
    return max; 
} 


//// A recursive binary search function. It returns 
//// location of x in given array arr[l..r] is present, 
//// otherwise -1 
//EDGE_ID bin_search_min_update_V(min_update_V* arr, EDGE_ID l, EDGE_ID r, EDGE_ID mm, EDGE_ID max) 
//{ 
//    if (r >= l) { 
//        EDGE_ID mid = l + (r - l) / 2; 
//
//        if (arr[mid].mm == mm) 
//            return mid; 
//  
//        // If element is smaller than mid, then 
//        // it can only be present in left subarray 
//        if (arr[mid].mm > mm) 
//        {
//
//          /// PRECAUTIONARY: CAN REMOVE LATER
//            if (!mid){
//                return max; 
//                printf("\nMID 0 WILL GIVE ERROR FOR UNSIGNED NEXT");
//                getchar();
//            }
//          ///////////////////
//            return bin_search_min_update_V(arr, l, mid - 1, mm, max); 
//        }
//  
//        // Else the element can only be present 
//        // in right subarray 
//        return bin_search_min_update_V(arr, mid + 1, r, mm, max); 
//    } 
//  
//    // We reach here when element is not 
//    // present in array 
//    //printf("\nNOT FOUND");
//    return max; 
//} 


// A recursive binary search function. It returns 
// location of x in given array arr[l..r] is present, 
// otherwise -1 
EDGE_ID bin_search_cycle_ops(EDGE_ID* arr, EDGE_ID l, EDGE_ID r, EDGE_ID mm, EDGE_ID max) 
{ 
    if (r >= l) { 
        EDGE_ID mid = l + (r - l) / 2; 

        if (arr[mid] == mm) 
            return mid; 
  
        // If element is smaller than mid, then 
        // it can only be present in left subarray 
        if (arr[mid] > mm) 
        {

          /// PRECAUTIONARY: CAN REMOVE LATER
            if (!mid){
                return max; 
                printf("\nMID 0 WILL GIVE ERROR FOR UNSIGNED NEXT");
                getchar();
            }
          ///////////////////
            return bin_search_cycle_ops(arr, l, mid - 1, mm, max); 
        }
  
        // Else the element can only be present 
        // in right subarray 
        return bin_search_cycle_ops(arr, mid + 1, r, mm, max); 
    } 
  
    // We reach here when element is not 
    // present in array 
    //printf("\nNOT FOUND");
    return max; 
} 


// A recursive binary search function. It returns 
// location of x in given array arr[l..r] is present, 
// otherwise -1 
EDGE_ID bin_search_cyc_in_cyc(cyc_in_cyc* arr, EDGE_ID l, EDGE_ID r, EDGE_ID mm, EDGE_ID max) 
{ 
    if (r >= l) { 
        EDGE_ID mid = l + (r - l) / 2; 

        if (arr[mid].cj == mm) 
            return mid; 
  
        // If element is smaller than mid, then 
        // it can only be present in left subarray 
        if (arr[mid].cj > mm) 
        {

          /// PRECAUTIONARY: CAN REMOVE LATER
            if (!mid){
                return max; 
                printf("\nMID 0 WILL GIVE ERROR FOR UNSIGNED NEXT");
                getchar();
            }
          ///////////////////
            return bin_search_cyc_in_cyc(arr, l, mid - 1, mm, max); 
        }
  
        // Else the element can only be present 
        // in right subarray 
        return bin_search_cyc_in_cyc(arr, mid + 1, r, mm, max); 
    } 
  
    // We reach here when element is not 
    // present in array 
    //printf("\nNOT FOUND");
    return max; 
} 



// returns greater than or equal to
void find_H1_cohom_greater(filtration* self, coboundary_H1* V_info, simplex* pivot){


    //EDGE_ID o_min = self->g_n_valid_edges;
    //EDGE_ID v_min = 0;


    if (pivot->key1 < V_info->o_ab){
        
          // Find first low of o_ab
          find_H1_cohom_low(self, V_info);

          // If it has a low
          if (V_info->low.key1 < self->g_n_valid_edges){

                // Need to find a_ptr and b_ptr if first low.key1 > e
                if (V_info->low.key1 > V_info->o_ab){

                    VERT_ID a = self->g_edges_list[2*V_info->o_ab];
                    VERT_ID b = self->g_edges_list[2*V_info->o_ab+1];

                    V_info->a_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[a], 0, self->g_Neigh_len[a]-1\
                                                         , V_info->low.key1, self->g_Neigh_len[a]);

                    V_info->b_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[b], 0, self->g_Neigh_len[b]-1\
                                                            , V_info->low.key1, self->g_Neigh_len[b]);

                }

                
          }

          return;

          
    }
    else if (pivot->key1 == V_info->o_ab){

        VERT_ID a = self->g_edges_list[2*V_info->o_ab];
        VERT_ID b = self->g_edges_list[2*V_info->o_ab+1];

        V_info->a_ptr = 0;
        V_info->b_ptr = 0;

        EDGE_ID o_min = self->g_n_valid_edges;
        EDGE_ID v_min;

        EDGE_ID o_s;
        EDGE_ID v_s;


        V_info->a_ptr = bin_search_min_geq_N(self->g_Neighbors[a], 0, self->g_Neigh_len[a]-1\
                                                , pivot->key2, self->g_Neigh_len[a]);
        
        V_info->b_ptr = bin_search_min_geq_N(self->g_Neighbors[b], 0, self->g_Neigh_len[b]-1\
                                                             , pivot->key2, self->g_Neigh_len[b]);

        if ((self->g_Neighbors[a][V_info->a_ptr].neighbor == pivot->key2) \
            && (self->g_Neighbors[b][V_info->b_ptr].neighbor == pivot->key2)){

            V_info->low.key1 = V_info->o_ab;
            V_info->low.key2 = pivot->key2;
            return;
            
        }

        
        while(1) {
              
            if ((V_info->a_ptr < self->g_Neigh_len[a]) \
                && (V_info->b_ptr < self->g_Neigh_len[b])){

                if (self->g_Neighbors[a][V_info->a_ptr].neighbor < self->g_Neighbors[b][V_info->b_ptr].neighbor) {
                      
                      V_info->a_ptr++;

                }
                else if (self->g_Neighbors[a][V_info->a_ptr].neighbor > self->g_Neighbors[b][V_info->b_ptr].neighbor) {
                      
                      V_info->b_ptr++;

                }
                else {

                      EDGE_ID o_ac = self->g_Neighbors[a][V_info->a_ptr].order;
                      EDGE_ID o_bc = self->g_Neighbors[b][V_info->b_ptr].order;
                      EDGE_ID c = self->g_Neighbors[b][V_info->b_ptr].neighbor;

                      o_s = o_ac;
                      v_s = b;
                      if (o_bc > o_s){
                        o_s = o_bc;
                        v_s = a;
                      }

                      if (o_s < V_info->o_ab){

                          if ((c > pivot->key2) || (c == pivot->key2)) {
                                V_info->low.key1 = V_info->o_ab;
                                V_info->low.key2 = c;
                                return;
                          }
                          
                      }
                      else{

                          if (o_s < o_min){
                            o_min = o_s;
                            v_min = v_s;
                          }
                          else if (o_s == o_min){
                            if (v_s < v_min){
                              v_min = v_s;
                            }
                          }

                      }

                      V_info->a_ptr++;
                      V_info->b_ptr++;
                }

            }
            else{
                
                  if (o_min != self->g_n_valid_edges){

                      V_info->a_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[a], 0, self->g_Neigh_len[a]-1\
                                                            , o_min, self->g_Neigh_len[a]);
                      V_info->b_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[b], 0, self->g_Neigh_len[b]-1\
                                                            , o_min, self->g_Neigh_len[b]);

                  }

                  V_info->low.key1 = o_min;
                  V_info->low.key2 = v_min;
                  return;
            }

        }
        

    }
    else {

        VERT_ID a = self->g_edges_list[2*V_info->o_ab];
        VERT_ID b = self->g_edges_list[2*V_info->o_ab+1];

        V_info->a_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[a], 0, self->g_Neigh_len[a]-1\
                                              , pivot->key1, self->g_Neigh_len[a]);
        V_info->b_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[b], 0, self->g_Neigh_len[b]-1\
                                              , pivot->key1, self->g_Neigh_len[b]);

        //printf("\na_len b_len %d, %d", self->g_Neigh_len[a], self->g_Neigh_len[b]);

        while (1) {
              
            //printf("\nptrs are %d, %d", V_info->a_ptr, V_info->b_ptr);
            //getchar();
            if ((V_info->a_ptr < self->g_Neigh_len[a]) \
                && (V_info->b_ptr < self->g_Neigh_len[b])){

                  if (self->g_Neighbors_e[a][V_info->a_ptr].order < self->g_Neighbors_e[b][V_info->b_ptr].order){
                        
                        EDGE_ID o_ac = self->g_Neighbors_e[a][V_info->a_ptr].order;
                        VERT_ID c = self->g_Neighbors_e[a][V_info->a_ptr].neighbor;
                        VERT_ID idx = search_Neighbors(self, b, c, 0, self->g_Neigh_len[b]-1);
                        if (idx < self->g_n_vert) {
                              
                            EDGE_ID o_bc = self->g_Neighbors[b][idx].order;
                            if (o_bc < o_ac){

                                  if ((o_ac > pivot->key1)\
                                    ||((o_ac == pivot->key1) && (b > pivot->key2))\
                                    ||((o_ac == pivot->key1) && (b == pivot->key2))){

                                      V_info->low.key1 = o_ac;
                                      V_info->low.key2 = b;
                                      return;

                                  }

                            }

                        }
                        //else{
                        V_info->a_ptr++;
                        //}

                  }
                  else {
                        
                        EDGE_ID o_bc = self->g_Neighbors_e[b][V_info->b_ptr].order;
                        VERT_ID c = self->g_Neighbors_e[b][V_info->b_ptr].neighbor;

#ifdef COMBIDX
                        EDGE_ID o_ac = COMB_IDX(a, c);
                        if (o_ac != self->g_n_valid_edges){
#else

                        VERT_ID idx = search_Neighbors(self, a, c, 0, self->g_Neigh_len[a]-1);
                        if (idx < self->g_n_vert) {
                            EDGE_ID o_ac = self->g_Neighbors[a][idx].order;
#endif

                            if (o_ac < o_bc){
                                  if ((o_bc > pivot->key1)\
                                    ||((o_bc == pivot->key1) && (a > pivot->key2))\
                                    ||((o_bc == pivot->key1) && (a == pivot->key2))){

                                      V_info->low.key1 = o_bc;
                                      V_info->low.key2 = a;
                                      return;
                                  }

                            }

                        }



                        //else{
                        V_info->b_ptr++;
                        //}

                  }

            }
            else if (V_info->a_ptr < self->g_Neigh_len[a]){
                  
                //Here b_ptr has reached end. So, o_bc should be less than o_ac

                EDGE_ID o_ac = self->g_Neighbors_e[a][V_info->a_ptr].order;
                VERT_ID c = self->g_Neighbors_e[a][V_info->a_ptr].neighbor;

#ifdef COMBIDX
                EDGE_ID o_bc = COMB_IDX(b, c);
                if (o_bc != self->g_n_valid_edges){

#else
                VERT_ID idx = search_Neighbors(self, b, c, 0, self->g_Neigh_len[b]-1);
                if (idx < self->g_n_vert) {
#endif

                    // check if errors 
                    //if (o_bc > o_ac){
                    //    printf("\nError check obc_oac");
                    //    getchar();
                    //}

                    if ((o_ac > pivot->key1)\
                      ||((o_ac == pivot->key1) && (b > pivot->key2))\
                      ||((o_ac == pivot->key1) && (b == pivot->key2))){
                        V_info->low.key1 = o_ac;
                        V_info->low.key2 = b;
                        return;
                    }

                }

                //else{
                V_info->a_ptr++;
                //}
                
                  
            }
            else if (V_info->b_ptr < self->g_Neigh_len[b]){
                  
                //Here b_ptr has reached end. So, o_ac should be less than o_bc

                EDGE_ID o_bc = self->g_Neighbors_e[b][V_info->b_ptr].order;
                VERT_ID c = self->g_Neighbors_e[b][V_info->b_ptr].neighbor;

#ifdef COMBIDX
                EDGE_ID o_ac = COMB_IDX(a, c);
                if (o_ac != self->g_n_valid_edges){

#else
                VERT_ID idx = search_Neighbors(self, a, c, 0, self->g_Neigh_len[a]-1);
                if (idx < self->g_n_vert) {
#endif
                    // check if errors 
                    //if (o_ac > o_bc){
                    //    printf("\nError check obc_oac");
                    //    getchar();
                    //}

                    if ((o_bc > pivot->key1)\
                      ||((o_bc == pivot->key1) && (a > pivot->key2))\
                      ||((o_bc == pivot->key1) && (a == pivot->key2))){
                        V_info->low.key1 = o_bc;
                        V_info->low.key2 = a;
                        return;
                    }

                }

                //else{
                V_info->b_ptr++;
                //}
                
                  
            }
            else{
              break;
            }

              
        }

        V_info->low.key1 = self->g_n_valid_edges;
        return;

    }

}

void find_H1_cohom_next (filtration* self, coboundary_H1* V_info){


      VERT_ID a = self->g_edges_list[2*V_info->o_ab];
      VERT_ID b = self->g_edges_list[2*V_info->o_ab+1];


      //printf("\nLOW of %d is (%d, %d)", V_info->o_ab, V_info->low.key1, V_info->low.key2);
      //

      //printf("\nCurrent c, ac, c, bc %d, %d, %d, %d", self->g_Neighbors[a][V_info->a_ptr].neighbor\
      //                                , self->g_Neighbors[a][V_info->a_ptr].order\
      //                                , self->g_Neighbors[b][V_info->b_ptr].neighbor\
      //                                , self->g_Neighbors[b][V_info->b_ptr].order\
      //                                );

      if (V_info->o_ab == V_info->low.key1){

          V_info->a_ptr++;
          V_info->b_ptr++;
          //printf("\naptr, max, bptr, max %d, %d, %d, %d", V_info->a_ptr\
          //                                              , self->g_Neigh_len[a]\
          //                                              , V_info->b_ptr\
          //                                              , self->g_Neigh_len[b]);
          //  
          //printf("\nstarting loop1");
          //EDGE_ID o_min = self->g_n_valid_edges;
          //EDGE_ID v_min;

          while (1){
                
                //printf("\naptr, max, bptr, max %d, %d, %d, %d", V_info->a_ptr\
                //                                        , self->g_Neigh_len[a]\
                //                                        , V_info->b_ptr\
                //                                        , self->g_Neigh_len[b]);
                
                if ((V_info->a_ptr < self->g_Neigh_len[a]) \
                    && (V_info->b_ptr < self->g_Neigh_len[b])){

                    if (self->g_Neighbors[a][V_info->a_ptr].neighbor < self->g_Neighbors[b][V_info->b_ptr].neighbor){
                          
                          V_info->a_ptr++;
                    }
                    else if (self->g_Neighbors[a][V_info->a_ptr].neighbor > self->g_Neighbors[b][V_info->b_ptr].neighbor){
                          
                          V_info->b_ptr++;
                    }
                    else{


                          EDGE_ID o_ac = self->g_Neighbors[a][V_info->a_ptr].order;
                          EDGE_ID o_bc = self->g_Neighbors[b][V_info->b_ptr].order;
                          EDGE_ID c = self->g_Neighbors[b][V_info->b_ptr].neighbor;
                          //printf("\nINSIDE FOUND NEXT COMMON %d", c);
                          EDGE_ID o_s = o_ac;
                          EDGE_ID v_s = b;
                          if (o_bc > o_ac) {
                              o_s = o_bc;
                              v_s = a;
                          }
                          
                          if (o_s < V_info->low.key1) {
                              
                                V_info->low.key2 = c;
                                //printf("\nReturn 1 ");
                                return;

                          }

                          //else if (o_s < o_min){
                          //      
                          //      o_min = o_s;
                          //      v_min = v_s;
                          //      
                          //}
                          //else if ((o_s == o_min) && (v_s < v_min)){

                          //      v_min = v_s;

                          //}

                          V_info->a_ptr++;
                          V_info->b_ptr++;


                    }
                    

                }
                else {
                    break;
                }

                
          }
          
          //printf("\nended loop1");

          //printf("\nsearching ptrs");

          V_info->a_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[a], 0, self->g_Neigh_len[a] - 1\
                                                , V_info->o_ab, self->g_Neigh_len[a]);
          V_info->b_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[b], 0, self->g_Neigh_len[b] - 1\
                                                , V_info->o_ab, self->g_Neigh_len[b]);


          //if (o_min < self->g_n_valid_edges){
          //      
          //      V_info->a_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[a], 0, self->g_Neigh_len[a] - 1\
          //                                            , o_min, self->g_Neigh_len[a]);
          //      V_info->b_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[b], 0, self->g_Neigh_len[b] - 1\
          //                                            , o_min, self->g_Neigh_len[b]);


          //}

          //V_info->low.key1 = o_min;
          //V_info->low.key2 = v_min;

          //printf("\nReturn 3 ");
          
          //return;
            
      }

      // Here 
      //V_info->low.key1 > V_info->o_ab


      //printf("\norders are ac: %d, bc: %d", self->g_Neighbors_e[a][V_info->a_ptr].order\
      //                          , self->g_Neighbors_e[b][V_info->b_ptr].order);

     if ((V_info->a_ptr < self->g_Neigh_len[a]) \
         && (V_info->b_ptr < self->g_Neigh_len[b])){

            if ((self->g_Neighbors_e[a][V_info->a_ptr].order < self->g_Neighbors_e[b][V_info->b_ptr].order)){

                 V_info->a_ptr++;

            }
            else{
                 V_info->b_ptr++;
            }

     }
     else{

            if (V_info->a_ptr < self->g_Neigh_len[a]){
                        V_info->a_ptr++;
            }
            else{
                        V_info->b_ptr++;
            }

     }


      //printf("\nstarting loop2");
      while (1){
            
            //printf("\naptr, max, bptr, max %d, %d, %d, %d", V_info->a_ptr\
            //                                        , self->g_Neigh_len[a]\
            //                                        , V_info->b_ptr\
            //                                        , self->g_Neigh_len[b]);
            if ((V_info->a_ptr < self->g_Neigh_len[a]) \
                && (V_info->b_ptr < self->g_Neigh_len[b])){
                
                 if (self->g_Neighbors_e[a][V_info->a_ptr].order < self->g_Neighbors_e[b][V_info->b_ptr].order) {

                      EDGE_ID o_ac = self->g_Neighbors_e[a][V_info->a_ptr].order;
                      EDGE_ID c = self->g_Neighbors_e[a][V_info->a_ptr].neighbor;

#ifdef COMBIDX
                      EDGE_ID o_bc = COMB_IDX(b, c);
                      if (o_bc != self->g_n_valid_edges){

#else
                      VERT_ID idx = search_Neighbors(self, b, c, 0, self->g_Neigh_len[b]-1);
                      if (idx < self->g_n_vert) {
                          EDGE_ID o_bc = self->g_Neighbors[b][idx].order;
#endif


                          if (o_bc < o_ac){
                                V_info->low.key1 = o_ac;
                                V_info->low.key2 = b;
                                //printf("\nReturn 4 ");
                                return;

                          }

                      }

                      V_info->a_ptr++;
                 }
                 else{
                      EDGE_ID o_bc = self->g_Neighbors_e[b][V_info->b_ptr].order;
                      EDGE_ID c = self->g_Neighbors_e[b][V_info->b_ptr].neighbor;
#ifdef COMBIDX
                      EDGE_ID o_ac = COMB_IDX(a, c);
                      if (o_ac != self->g_n_valid_edges){
#else
                      VERT_ID idx = search_Neighbors(self, a, c, 0, self->g_Neigh_len[a]-1);
                      if (idx < self->g_n_vert) {
                          EDGE_ID o_ac = self->g_Neighbors[a][idx].order;
#endif


                          if (o_ac < o_bc){
                                V_info->low.key1 = o_bc;
                                V_info->low.key2 = a;
                                //printf("\nReturn 5 ");
                                return;

                          }

                      }

                      V_info->b_ptr++;

                 }


            }
            else if (V_info->a_ptr < self->g_Neigh_len[a]){

                      EDGE_ID o_ac = self->g_Neighbors_e[a][V_info->a_ptr].order;
                      EDGE_ID c = self->g_Neighbors_e[a][V_info->a_ptr].neighbor;

#ifdef COMBIDX
                      EDGE_ID o_bc = COMB_IDX(b, c);
                      if (o_bc != self->g_n_valid_edges){
#else
                      VERT_ID idx = search_Neighbors(self, b, c, 0, self->g_Neigh_len[b]-1);
                      if (idx < self->g_n_vert) {
#endif
                            
                          //EDGE_ID o_bc = self->g_Neighbors[b][idx];
                          // SHOULD HAVE THE NEED TO CHECK
                          //if (o_bc < o_ac){
                          V_info->low.key1 = o_ac;
                          V_info->low.key2 = b;
                          //printf("\nReturn 5 ");
                          return;

                          //}

                      }

                      V_info->a_ptr++;
                  
            }
            else if (V_info->b_ptr < self->g_Neigh_len[b]){

                      EDGE_ID o_bc = self->g_Neighbors_e[b][V_info->b_ptr].order;
                      EDGE_ID c = self->g_Neighbors_e[b][V_info->b_ptr].neighbor;

#ifdef COMBIDX
                      EDGE_ID o_ac = COMB_IDX(a, c);
                      if (o_ac != self->g_n_valid_edges){
#else
                      VERT_ID idx = search_Neighbors(self, a, c, 0, self->g_Neigh_len[a]-1);
                      if (idx < self->g_n_vert) {
#endif
                            
                          //EDGE_ID o_bc = self->g_Neighbors[b][idx];
                          // SHOULD HAVE THE NEED TO CHECK
                          //if (o_bc < o_ac){
                          V_info->low.key1 = o_bc;
                          V_info->low.key2 = a;
                          //printf("\nReturn 6 ");
                          return;

                          //}

                      }

                      V_info->b_ptr++;
                  
            }
            else{
              break;
            }

      }
      //printf("\nended loop2");

      V_info->low.key1 = self->g_n_valid_edges;
      //printf("\nReturn 7 ");
      return;

}



EDGE_ID bin_search_min_geq_Ne(Neighbors* arr, VERT_ID l, VERT_ID r, VERT_ID x, EDGE_ID MAX){

    if (arr[r].order < x){
      return MAX;
    }

    if (arr[l].order > x){
      return l;
    }

    VERT_ID mid = l + (r-l)/2;

    if (arr[mid].order < x){
        
        l = mid + 1;
        if ((arr[l].order > x) || (arr[l].order == x)){
          return l;
        }
        bin_search_min_geq_Ne(arr, l, r, x, MAX);

    }
    else{
        r = mid;
        if (arr[r].order == x) return r;
        bin_search_min_geq_Ne(arr, l , r, x, MAX);

    }
     
}

EDGE_ID bin_search_min_geq_N(Neighbors* arr, VERT_ID l, VERT_ID r, VERT_ID x, EDGE_ID MAX){

    if (arr[r].neighbor < x){
      return MAX;
    }

    if (arr[l].neighbor > x){
      return l;
    }

    VERT_ID mid = l + (r-l)/2;

    if (arr[mid].neighbor < x){
        
        l = mid + 1;
        if ((arr[l].neighbor > x) || (arr[l].neighbor == x)){
          return l;
        }
        bin_search_min_geq_N(arr, l, r, x, MAX);

    }
    else{
        r = mid;
        if (arr[r].neighbor == x) return r;
        bin_search_min_geq_N(arr, l , r, x, MAX);

    }
     
}



void find_H2_cohom_low (filtration* self, coboundary_H2* V_info){


      V_info->c_ptr = 0;

      //int flag = H2_case1 (self, V_info);
      if (H2_case1(self, V_info)){
        return;
      }

      
      VERT_ID a = self->g_edges_list[2*V_info->triangle.key1];
      VERT_ID b = self->g_edges_list[2*V_info->triangle.key1+1];


      V_info->a_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[a], 0, self->g_Neigh_len[a]-1\
                                              , V_info->triangle.key1, self->g_Neigh_len[a]);
      
      V_info->b_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[b], 0, self->g_Neigh_len[b]-1\
                                              , V_info->triangle.key1, self->g_Neigh_len[b]);


      V_info->a_ptr++;
      V_info->b_ptr++;

      H2_case2(self, V_info);


}



void find_H2_cohom_next (filtration* self, coboundary_H2* V_info){
    
      //clock_gettime(CLOCK_MONOTONIC, &(self->g_start_wall_clock));

      //printf("\nfinding H2 next");
      
      EDGE_ID o_ab = V_info->triangle.key1;
      //VERT_ID c = V_info->triangle.key2;
      
      VERT_ID a = self->g_edges_list[2*o_ab];
      VERT_ID b = self->g_edges_list[2*o_ab+1];

      //EDGE_ID o_ad, o_bd;
      //VERT_ID idxa, idxb, d;

      int flag = 0;

      //if (V_info->low.key1 == o_ab){ 
      if (V_info->vertex == 0){ 

            V_info->c_ptr++;

            if (H2_case1(self, V_info)){
              return;
            }

            // Here means that we did not return and 
            // not a_ptr and b_ptr are at o_ab
            // So, both need to be incremented

            V_info->a_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[a], 0, self->g_Neigh_len[a]-1\
                                                 , V_info->triangle.key1, self->g_Neigh_len[a]);

            V_info->b_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[b], 0, self->g_Neigh_len[b]-1\
                                                    , V_info->triangle.key1, self->g_Neigh_len[b]);
            V_info->a_ptr++;
            V_info->b_ptr++;
            flag = 1;

      }

      if (!flag){

          if (V_info->vertex == 1){
                V_info->a_ptr++;
          }
          else if (V_info->vertex == 2){
                V_info->b_ptr++;
          }
          else if (V_info->vertex == 3){
                V_info->c_ptr++;
          }
            
            
      }


      H2_case2(self, V_info);


}


void find_H2_cohom_greater (filtration* self, coboundary_H2* V_info, simplex* pivot){

    //if (self->g_p_flag){
    //    printf("\nfinding H2 greater than (%d, %d) for (%d, %d)", pivot->key1, pivot->key2
    //                                                            , V_info->triangle.key1\
    //                                                            , V_info->triangle.key2);
    //    getchar();
    //}

    //EDGE_ID o_ab;
    VERT_ID c, a, b;


    if (pivot->key1 < V_info->triangle.key1){
        
          // Find first low of o_ab
          find_H2_cohom_low(self, V_info);
          return;

    }
    else if (pivot->key1 == V_info->triangle.key1){

          //o_ab = V_info->triangle.key1;
          
          VERT_ID c = V_info->triangle.key2;

          //a = self->g_edges_list[o_ab][0];
          //b = self->g_edges_list[o_ab][1];

          V_info->c_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[c], 0, self->g_Neigh_len[c]-1\
                                                  , pivot->key2, self->g_Neigh_len[c]);

          if (self->g_Neighbors_e[c][V_info->c_ptr].order == pivot->key2){
              
                V_info->low = *pivot;
                V_info->vertex= 0;
                return;

          }

          if (H2_case1(self, V_info)){
                return;
          }
          

          // Here means that we did not return and 
          // not a_ptr and b_ptr are at o_ab
          // So, both need to be incremented

          a = self->g_edges_list[2*V_info->triangle.key1];
          b = self->g_edges_list[2*V_info->triangle.key1+1];

          V_info->a_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[a], 0, self->g_Neigh_len[a]-1\
                                               , V_info->triangle.key1, self->g_Neigh_len[a]);

          V_info->b_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[b], 0, self->g_Neigh_len[b]-1\
                                                  , V_info->triangle.key1, self->g_Neigh_len[b]);

          V_info->a_ptr++;
          V_info->b_ptr++;
                

    }
    else{


          c = V_info->triangle.key2;

          a = self->g_edges_list[2*V_info->triangle.key1];
          b = self->g_edges_list[2*V_info->triangle.key1+1];

          V_info->a_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[a], 0, self->g_Neigh_len[a]-1\
                                               , pivot->key1, self->g_Neigh_len[a]);

          V_info->b_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[b], 0, self->g_Neigh_len[b]-1\
                                                  , pivot->key1, self->g_Neigh_len[b]);

          V_info->c_ptr = bin_search_min_geq_Ne(self->g_Neighbors_e[c], 0, self->g_Neigh_len[c]-1\
                                                  , pivot->key1, self->g_Neigh_len[c]);


    }

    while (1){

          H2_case2(self, V_info);

          if (((V_info->low.key1 == pivot->key1) && (V_info->low.key2 > pivot->key2))\
              ||(V_info->low.key1 > pivot->key1) || (V_info->low.key1 == self->g_n_valid_edges)\
              ||((V_info->low.key1 == pivot->key1) && (V_info->low.key2 == pivot->key2)) ){
              break;
          }

          if (V_info->vertex == 1){
                
                V_info->a_ptr++;
          }
          else if (V_info->vertex == 2){
                
                V_info->b_ptr++;
          }
          else if (V_info->vertex == 3){
                
                V_info->c_ptr++;
          }

    }


}



// A recursive binary search function. It returns 
// location of x in given array arr[l..r] is present, 
// otherwise -1 
EDGE_ID search_H2_cohom_pivots(H2_cohom_pivots* arr, EDGE_ID l, EDGE_ID r, EDGE_ID key2, EDGE_ID max) 
{ 
    if (r >= l) { 
        EDGE_ID mid = l + (r - l) / 2; 

        if (arr[mid].key2 == key2) 
            return mid; 
  
        // If element is smaller than mid, then 
        // it can only be present in left subarray 
        if (arr[mid].key2 > key2) 
        {

          /// PRECAUTIONARY: CAN REMOVE LATER
            if (!mid){
                return max; 
                printf("\nMID 0 WILL GIVE ERROR FOR UNSIGNED NEXT");
                getchar();
            }
          ///////////////////
            return search_H2_cohom_pivots(arr, l, mid - 1, key2, max); 
        }
  
        // Else the element can only be present 
        // in right subarray 
        return search_H2_cohom_pivots(arr, mid + 1, r, key2, max); 
    } 
  
    // We reach here when element is not 
    // present in array 
    //printf("\nNOT FOUND");
    return max; 
} 


// Reduces with complex in parallel

void* reduce_with_complex_H0(void* arg){
      
      filtration* self = arg;

      pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, 0);

      pthread_mutex_lock(&(self->g_thread_lock));

      int tid = ++self->g_thread_id;



      for (;;){

          self->g_sleeping_threads++;
          
          if (self->g_sleeping_threads == self->g_cpu_count)
              pthread_cond_signal(&(self->g_start_boss));

          pthread_cond_wait(&(self->g_start_workers), &(self->g_thread_lock));

          if (self->g_delete_threads){
            //printf("\nexiting from thread %d", tid);
            pthread_mutex_unlock(&(self->g_thread_lock));
            pthread_exit(NULL);
          }

          self->g_sleeping_threads--;

          pthread_mutex_unlock(&(self->g_thread_lock));

          for (int ws_counter = self->g_jobs[tid - 1]; ws_counter < self->g_jobs[tid]; ws_counter++){


              boundary_H0_ws* this_ws = self->g_R_ws_H0_info + ws_counter;

              EDGE_ID* orig = self->g_R_ws_H0[ws_counter] + this_ws->original*this_ws->max_len;

              this_ws->flag_red_w_complex = 0;
              this_ws->flag_append_to_complex = 1;

              EDGE_ID idx = self->g_pivots_H0[this_ws->pivot];

              while(idx){

                    //reduced_col = self->g_pivots[self->g_dim_now][idx].red_col;
                    //reduced_col = self->g_pivots_H0[idx];


                    EDGE_ID red_start_idx = self->g_R_col_indices_H0[idx];

                    EDGE_ID red_finish_idx = self->g_R_col_indices_H0[idx+1];

                    EDGE_ID red_len = red_finish_idx - red_start_idx;

                    if ((this_ws->len + red_len) > this_ws->max_len){
                          
                        if (this_ws->original){
                              
                                for (EDGE_ID it=0; it < this_ws->len; it++){

                                      self->g_R_ws_H0[ws_counter][it] = \
                                                        self->g_R_ws_H0[ws_counter][it + this_ws->max_len];

                                }

                                this_ws->original = 0;

                        }

                        this_ws->max_len = this_ws->len + red_len + 100;

                        pthread_mutex_lock(&(self->g_thread_lock));

                        self->g_R_ws_H0[ws_counter] = (EDGE_ID*)realloc(self->g_R_ws_H0[ws_counter]\
                                                                          , 2*this_ws->max_len*sizeof(EDGE_ID));

                        pthread_mutex_unlock(&(self->g_thread_lock));

                        orig = self->g_R_ws_H0[ws_counter];


                          
                    }


                    EDGE_ID* scratch = self->g_R_ws_H0[ws_counter] + (1-this_ws->original)*this_ws->max_len;



                    EDGE_ID orig_ptr = 0;
                    EDGE_ID red_ptr = red_start_idx;
                    EDGE_ID scratch_ptr = 0;


                    while ((orig_ptr < this_ws->len) && (red_ptr < red_finish_idx)){

                        if (orig[orig_ptr] < self->g_R_sparse_H0[red_ptr]){

                              scratch[scratch_ptr++] = orig[orig_ptr++];

                        }
                        else if (orig[orig_ptr] > self->g_R_sparse_H0[red_ptr]){

                              scratch[scratch_ptr++] = self->g_R_sparse_H0[red_ptr++];

                        }
                        else{
                              orig_ptr++;
                              red_ptr++;
                        }

                    }

                    while (orig_ptr < this_ws->len){

                        scratch[scratch_ptr++] = orig[orig_ptr++];

                    }

                    while (red_ptr < red_finish_idx){
                          
                        scratch[scratch_ptr++] = self->g_R_sparse_H0[red_ptr++];
                          
                    }


                    this_ws->len = scratch_ptr;


                    if (!this_ws->len){

                        //idx = self->g_n_reduced_simplex[self->g_dim_now];
                        //idx = -1;
                        break;

                    }
                    else{
                    

                        this_ws->original = 1 - this_ws->original;

                        orig = self->g_R_ws_H0[ws_counter] + this_ws->original*this_ws->max_len;

                        this_ws->pivot = orig[this_ws->len-1];
                        
                        idx = self->g_pivots_H0[this_ws->pivot];

                    }
                    

              }


          }

          pthread_mutex_lock(&(self->g_thread_lock));

          self->g_processed_threads++;

            
      }


}


void allocate_jobs(filtration* self, int ws_size){
      
     int x = (int)(ws_size/self->g_cpu_count);
     int y = (ws_size % self->g_cpu_count);

     self->g_jobs[0] = 0;

     for (int i = 1; i < self->g_cpu_count+1; i++){
          
          if (i < y + 1)
            self->g_jobs[i] = self->g_jobs[i-1] + x + 1;
          else
            self->g_jobs[i] = self->g_jobs[i-1] + x;

     }
      
}

void reduce_ws_H0(filtration* self){


      //if (self->g_n_reduced_simplex_H0 > 0){

        
            self->g_processed_threads = 0;


            pthread_cond_broadcast(&(self->g_start_workers));


            while (self->g_processed_threads != self->g_cpu_count){
                  
                  pthread_cond_wait(&(self->g_start_boss) \
                                  ,&(self->g_thread_lock));
            }



      //}


      reduce_with_self_H0( \
                            self \
                            );

      int count_valid = 0;

      for (int ws_counter=0; ws_counter < self->g_ws_counter; ws_counter++){

            if (!self->g_R_ws_H0_info[ws_counter].len){continue;}

            if (self->g_R_ws_H0_info[ws_counter].flag_append_to_complex){

                  update_R_H0(self \
                            , ws_counter
                            );
                  continue;
                    
            }


            // Swap R
            EDGE_ID* temp = self->g_R_ws_H0[count_valid];
            self->g_R_ws_H0[count_valid] = self->g_R_ws_H0[ws_counter];
            self->g_R_ws_H0[ws_counter] = temp;

            // Swap R info
            boundary_H0_ws temp2 = self->g_R_ws_H0_info[count_valid];
            self->g_R_ws_H0_info[count_valid] = self->g_R_ws_H0_info[ws_counter];
            self->g_R_ws_H0_info[ws_counter] = temp2;


            // At this point, this has to be a non-zero column
            self->g_R_ws_H0_info[count_valid].flag_non_empty = 1;

            count_valid += 1;

      }

      self->g_ws_counter = count_valid;

      //if (dim)
      //  self->g_H0_MAX = self->g_n_reduced_simplex[dim];

}


void reduce_with_self_H0( \
                      filtration* self \
                      ){


    int m;
    EDGE_ID orig_ptr, scratch_ptr, m_ptr, idx;
    EDGE_ID *orig, *scratch, *original_m;
    

    for (int ws_counter=0; ws_counter < self->g_ws_counter; ws_counter++){

        boundary_H0_ws* this_ws = self->g_R_ws_H0_info + ws_counter;

        // If the simplex has already been reduced to 0
      // then continue
        if (!this_ws->len){ 
          this_ws->flag_append_to_complex = 0;
          continue;

        }



        m = 0;
        while (m < ws_counter){

            boundary_H0_ws* m_ws = self->g_R_ws_H0_info + m;

            if (!m_ws->len){
                m++;
                continue;
            }

            if (m_ws->pivot > this_ws->pivot){
                  
                  if (m_ws->flag_red_w_complex){
                        
                        this_ws->flag_append_to_complex = 0;
                        break;
                  }
                  m++;
                  continue;
                  
            }


            if (m_ws->pivot < this_ws->pivot){
                    m++;
                    continue;
            }

            if (!m_ws->flag_append_to_complex){
                    m++;
                    continue;
            }

            orig = self->g_R_ws_H0[ws_counter] + this_ws->original*this_ws->max_len;

            if ((this_ws->len + m_ws->len) > this_ws->max_len){

                if (this_ws->original){
                      
                    for (EDGE_ID it = 0; it < this_ws->len; it++){
                          
                          orig[it] = orig[it + this_ws->max_len];
                    }
                       
                    this_ws->original = 0;

                }

                this_ws->max_len = this_ws->len + m_ws->len + 100;

                self->g_R_ws_H0[ws_counter] = (EDGE_ID*)realloc(self->g_R_ws_H0[ws_counter]\
                                                        , 2*this_ws->max_len*sizeof(EDGE_ID));

                orig = self->g_R_ws_H0[ws_counter];

                
            }

            scratch = self->g_R_ws_H0[ws_counter] + (1-this_ws->original)*this_ws->max_len;

            original_m = self->g_R_ws_H0[m] + m_ws->original*m_ws->max_len;

            
            // Store the result in scratch

            orig_ptr = 0;
            scratch_ptr = 0;
            m_ptr = 0;

            while ((orig_ptr < this_ws->len) && (m_ptr < m_ws->len)){
                  
                  if (orig[orig_ptr] < original_m[m_ptr]){
                        
                        scratch[scratch_ptr++]  = orig[orig_ptr++];
                  }
                  else if (orig[orig_ptr] > original_m[m_ptr]){
                        
                        scratch[scratch_ptr++]  = original_m[m_ptr++];
                  }
                  else{
                        orig_ptr++;
                        m_ptr++;
                  }


            }

            while (orig_ptr < this_ws->len){
                
                scratch[scratch_ptr++]  = orig[orig_ptr++];

            }

            while (m_ptr < m_ws->len){
                
                scratch[scratch_ptr++]  = original_m[m_ptr++];

            }


            this_ws->len = scratch_ptr;

            if (!scratch_ptr){
                  
                  this_ws->flag_append_to_complex = 0;
                  break;
                  
            }


            this_ws->pivot = scratch[scratch_ptr - 1];

            this_ws->original = 1 - this_ws->original;


            //if (self->g_n_reduced_simplex_H0){

                idx = self->g_pivots_H0[this_ws->pivot];
                // If the pivot is in red complex, then this has to be reduced w/ complex
                //if (idx != self->g_n_reduced_simplex[self->g_dim_now]){
                if (idx){
                      
                      this_ws->flag_red_w_complex = 1;
                      this_ws->flag_append_to_complex = 0;
                      break;

                }

            //}

            m = 0;

        }//End of m loop

    }

}//End of red_ws_w_self_single


void update_R_H0(filtration* self, int ws_counter){


      boundary_H0_ws* this_ws = self->g_R_ws_H0_info + ws_counter;

      EDGE_ID* orig = self->g_R_ws_H0[ws_counter] + this_ws->original*this_ws->max_len;


      // Check space for R Sparse
      if ((this_ws->len + self->g_R_sparse_ptr_H0) > self->g_R_sparse_max_H0 ){

          self->g_R_sparse_max_H0 = this_ws->len + self->g_R_sparse_ptr_H0 + 1000;

          self->g_R_sparse_H0 = (EDGE_ID*)realloc(self->g_R_sparse_H0\
                                                , self->g_R_sparse_max_H0*sizeof(EDGE_ID));
            
      }

      // Check space for R col indices
      if ((self->g_R_col_indices_ptr_H0 + 3) > self->g_R_col_indices_max_H0){
            
            self->g_R_col_indices_max_H0 += 100;
            self->g_R_col_indices_H0 = (EDGE_ID*)realloc(self->g_R_col_indices_H0\
                                                       , self->g_R_col_indices_max_H0*sizeof(EDGE_ID));
      }



      self->g_pivots_H0[this_ws->pivot] = self->g_R_col_indices_ptr_H0;

      self->g_R_col_indices_H0[self->g_R_col_indices_ptr_H0++] = self->g_R_sparse_ptr_H0;


      for (EDGE_ID j=0; j < this_ws->len; j++){

          self->g_R_sparse_H0[self->g_R_sparse_ptr_H0++] = orig[j];

      }

      self->g_R_col_indices_H0[self->g_R_col_indices_ptr_H0++] = self->g_R_sparse_ptr_H0;


      // Update edges with pivots for H0 to be used in clearing algo
      self->g_edges_with_pivots_H0[this_ws->cob] = 1;
      

//#ifdef HOM_CYCLES      
      // Update vertex in pivot to edge mapping
      if (self->g_compute_cycles){
          self->g_H0_pivot_of[this_ws->pivot].coface = this_ws->cob;
      }
//#endif


}



//////////////////////////////////////////////////////////
// MERGE SORT ALGORITHMS
//////////////////////////////////////////////////////////

// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge(PAR* arr, EDGE_ID* aux, EDGE_ID l, EDGE_ID m, EDGE_ID r) 
{ 
    EDGE_ID i, j, k; 
    EDGE_ID n1 = m - l + 1; 
    EDGE_ID n2 =  r - m; 
    //printf("\nn1, n2: %u, %u", n1, n2);
  
    /* create temp arrays */
    PAR *L, *R;
    L = (PAR*)malloc(n1*sizeof(PAR));
    R = (PAR*)malloc(n2*sizeof(PAR));

    /* create temp arrays */
    EDGE_ID* L_aux;
    EDGE_ID* R_aux;
    L_aux = (EDGE_ID*)malloc(2*n1*sizeof(EDGE_ID));
    R_aux = (EDGE_ID*)malloc(2*n2*sizeof(EDGE_ID));

    
    //int L[n1], R[n2]; 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){ 
        L[i] = arr[l + i]; 
        L_aux[2*i] = aux[2*(l + i)]; 
        L_aux[2*i+1] = aux[2*(l + i) + 1]; 
    }

    for (j = 0; j < n2; j++) {
        R[j] = arr[m + 1+ j]; 
        R_aux[2*j] = aux[2*(m + 1+ j)]; 
        R_aux[2*j+1] = aux[2*(m + 1+ j)+1]; 
    }
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 


    while (i < n1 && j < n2) 
    { 

          if (L[i] <= R[j]) 
          { 
              arr[k] = L[i]; 
	            aux[2*k] = L_aux[2*i];
	            aux[2*k+1] = L_aux[2*i+1];
              i++; 
          } 
          else
          { 
              arr[k] = R[j]; 
              aux[2*k] = R_aux[2*j]; 
              aux[2*k+1] = R_aux[2*j+1]; 
              j++; 
          } 

          k++;

    }

  
    /* Copy the remaining elements of L[], if there 
       are any */
    while (i < n1) 
    { 
        arr[k] = L[i]; 
        aux[2*k] = L_aux[2*i]; 
        aux[2*k+1] = L_aux[2*i+1]; 
        i++; 
        k++; 
    } 
  
    /* Copy the remaining elements of R[], if there 
       are any */
    while (j < n2) 
    { 
        arr[k] = R[j]; 
        aux[2*k] = R_aux[2*j]; 
        aux[2*k+1] = R_aux[2*j+1]; 
        j++; 
        k++; 
    } 

    free(L);
    free(R);
    free(L_aux);
    free(R_aux);
} 
  
/* l is for left index and r is right index of the 
   sub-array of arr to be sorted */
void mergeSort(PAR* arr, EDGE_ID* aux, EDGE_ID l, EDGE_ID r) 
{ 
    if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        EDGE_ID m = l+(r-l)/2; 
  
        // Sort first and second halves 
        mergeSort(arr, aux, l, m); 
        mergeSort(arr, aux, m+1, r); 
  
        merge(arr, aux, l, m, r); 
    } 
} 


//////////////////////////////////////////////////////////
// MERGE SORT ALGORITHMS FOR V_H0
//////////////////////////////////////////////////////////

// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge_V_H0(EDGE_ID* arr, EDGE_ID** aux, EDGE_ID* aux2, EDGE_ID* aux3, EDGE_ID l, EDGE_ID m, EDGE_ID r) 
{ 
    EDGE_ID i, j, k; 
    EDGE_ID n1 = m - l + 1; 
    EDGE_ID n2 =  r - m; 
    //printf("\nn1, n2: %u, %u", n1, n2);
  
    /* create temp arrays */
    EDGE_ID *L, *R;
    L = (EDGE_ID*)malloc(n1*sizeof(EDGE_ID));
    R = (EDGE_ID*)malloc(n2*sizeof(EDGE_ID));

    /* create temp arrays */
    EDGE_ID** L_aux;
    EDGE_ID** R_aux;
    L_aux = (EDGE_ID**)malloc(n1*sizeof(EDGE_ID*));
    R_aux = (EDGE_ID**)malloc(n2*sizeof(EDGE_ID*));

    /* create temp arrays */
    EDGE_ID* L_aux2;
    EDGE_ID* R_aux2;
    L_aux2 = (EDGE_ID*)malloc(n1*sizeof(EDGE_ID));
    R_aux2 = (EDGE_ID*)malloc(n2*sizeof(EDGE_ID));

    /* create temp arrays */
    EDGE_ID* L_aux3;
    EDGE_ID* R_aux3;
    L_aux3 = (EDGE_ID*)malloc(n1*sizeof(EDGE_ID));
    R_aux3 = (EDGE_ID*)malloc(n2*sizeof(EDGE_ID));

    //int L[n1], R[n2]; 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){ 
        L[i] = arr[l + i]; 
        L_aux[i] = aux[l + i]; 
        L_aux2[i] = aux2[l + i]; 
        L_aux3[i] = aux3[l + i]; 
    }

    for (j = 0; j < n2; j++) {
        R[j] = arr[m + 1+ j]; 
        R_aux[j] = aux[m + 1+ j]; 
        R_aux2[j] = aux2[m + 1+ j]; 
        R_aux3[j] = aux3[m + 1+ j]; 
    }
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 


    while (i < n1 && j < n2) 
    { 

          if (L[i] >= R[j]) 
          { 
              arr[k] = L[i]; 
	            aux[k] = L_aux[i];
	            aux2[k] = L_aux2[i];
	            aux3[k] = L_aux3[i];
              i++; 
          } 
          else
          { 
              arr[k] = R[j]; 
              aux[k] = R_aux[j]; 
              aux2[k] = R_aux2[j]; 
              aux3[k] = R_aux3[j]; 
              j++; 
          } 

          k++;

    }

  
    /* Copy the remaining elements of L[], if there 
       are any */
    while (i < n1) 
    { 
        arr[k] = L[i]; 
        aux[k] = L_aux[i]; 
        aux2[k] = L_aux2[i]; 
        aux3[k] = L_aux3[i]; 
        i++; 
        k++; 
    } 
  
    /* Copy the remaining elements of R[], if there 
       are any */
    while (j < n2) 
    { 
        arr[k] = R[j]; 
        aux[k] = R_aux[j]; 
        aux2[k] = R_aux2[j]; 
        aux3[k] = R_aux3[j]; 
        j++; 
        k++; 
    } 

    free(L);
    free(R);
    free(L_aux);
    free(R_aux);
    free(L_aux2);
    free(R_aux2);
    free(L_aux3);
    free(R_aux3);
} 
  
/* l is for left index and r is right index of the 
   sub-array of arr to be sorted */
void mergeSort_V_H0(EDGE_ID* arr, EDGE_ID** aux, EDGE_ID* aux2, EDGE_ID* aux3, EDGE_ID l, EDGE_ID r) 
{ 
    if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        EDGE_ID m = l+(r-l)/2; 
  
        // Sort first and second halves 
        mergeSort_V_H0(arr, aux, aux2, aux3, l, m); 
        mergeSort_V_H0(arr, aux, aux2, aux3, m+1, r); 
  
        merge_V_H0(arr, aux, aux2, aux3, l, m, r); 
    } 
} 


//////////////////////////////////////////////////////////
// MERGE SORT ALGORITHMS FOR V_H1
//////////////////////////////////////////////////////////

// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge_V_H1(EDGE_ID* arr, simplex** aux, EDGE_ID l, EDGE_ID m, EDGE_ID r) 
{ 
    EDGE_ID i, j, k; 
    EDGE_ID n1 = m - l + 1; 
    EDGE_ID n2 =  r - m; 
    //printf("\nn1, n2: %u, %u", n1, n2);
  
    /* create temp arrays */
    EDGE_ID *L, *R;
    L = (EDGE_ID*)malloc(n1*sizeof(EDGE_ID));
    R = (EDGE_ID*)malloc(n2*sizeof(EDGE_ID));

    /* create temp arrays */
    simplex** L_aux;
    simplex** R_aux;
    L_aux = (simplex**)malloc(n1*sizeof(simplex*));
    R_aux = (simplex**)malloc(n2*sizeof(simplex*));

    
    //int L[n1], R[n2]; 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){ 
        L[i] = arr[l + i]; 
        L_aux[i] = aux[l + i]; 
    }

    for (j = 0; j < n2; j++) {
        R[j] = arr[m + 1+ j]; 
        R_aux[j] = aux[m + 1+ j]; 
    }
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 


    while (i < n1 && j < n2) 
    { 

          if (L[i] >= R[j]) 
          { 
              arr[k] = L[i]; 
	            aux[k] = L_aux[i];
              i++; 
          } 
          else
          { 
              arr[k] = R[j]; 
              aux[k] = R_aux[j]; 
              j++; 
          } 

          k++;

    }

  
    /* Copy the remaining elements of L[], if there 
       are any */
    while (i < n1) 
    { 
        arr[k] = L[i]; 
        aux[k] = L_aux[i]; 
        i++; 
        k++; 
    } 
  
    /* Copy the remaining elements of R[], if there 
       are any */
    while (j < n2) 
    { 
        arr[k] = R[j]; 
        aux[k] = R_aux[j]; 
        j++; 
        k++; 
    } 

    free(L);
    free(R);
    free(L_aux);
    free(R_aux);
} 
  
/* l is for left index and r is right index of the 
   sub-array of arr to be sorted */
void mergeSort_V_H1(EDGE_ID* arr, simplex** aux, EDGE_ID l, EDGE_ID r) 
{ 
    if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        EDGE_ID m = l+(r-l)/2; 
  
        // Sort first and second halves 
        mergeSort_V_H1(arr, aux, l, m); 
        mergeSort_V_H1(arr, aux, m+1, r); 
  
        merge_V_H1(arr, aux, l, m, r); 
    } 

} 


//////////////////////////////////////////////////////////
// MERGE SORT ALGORITHMS FOR Llen
//////////////////////////////////////////////////////////

// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge_Llen(EDGE_ID* arr, EDGE_ID* aux, EDGE_ID* aux2, EDGE_ID l, EDGE_ID m, EDGE_ID r) 
{ 
    EDGE_ID i, j, k; 
    EDGE_ID n1 = m - l + 1; 
    EDGE_ID n2 =  r - m; 
    //printf("\nn1, n2: %u, %u", n1, n2);
  
    /* create temp arrays */
    EDGE_ID *L, *R;
    L = (EDGE_ID*)malloc(n1*sizeof(EDGE_ID));
    R = (EDGE_ID*)malloc(n2*sizeof(EDGE_ID));

    /* create temp arrays */
    EDGE_ID *L_aux, *R_aux;
    L_aux = (EDGE_ID*)malloc(n1*sizeof(EDGE_ID));
    R_aux = (EDGE_ID*)malloc(n2*sizeof(EDGE_ID));

    /* create temp arrays */
    EDGE_ID *L_aux2, *R_aux2;
    L_aux2 = (EDGE_ID*)malloc(n1*sizeof(EDGE_ID));
    R_aux2 = (EDGE_ID*)malloc(n2*sizeof(EDGE_ID));

    //int L[n1], R[n2]; 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){ 
        L[i] = arr[l + i]; 
        L_aux[i] = aux[l + i]; 
        L_aux2[i] = aux2[l + i]; 
    }

    for (j = 0; j < n2; j++) {
        R[j] = arr[m + 1+ j]; 
        R_aux[j] = aux[m + 1+ j]; 
        R_aux2[j] = aux2[m + 1+ j]; 
    }
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 


    while (i < n1 && j < n2) 
    { 

          if (L[i] >= R[j]) 
          { 
              arr[k] = L[i]; 
	            aux[k] = L_aux[i];
	            aux2[k] = L_aux2[i];
              i++; 
          } 
          else
          { 
              arr[k] = R[j]; 
              aux[k] = R_aux[j]; 
              aux2[k] = R_aux2[j]; 
              j++; 
          } 

          k++;

    }

  
    /* Copy the remaining elements of L[], if there 
       are any */
    while (i < n1) 
    { 
        arr[k] = L[i]; 
        aux[k] = L_aux[i]; 
        aux2[k] = L_aux2[i]; 
        i++; 
        k++; 
    } 
  
    /* Copy the remaining elements of R[], if there 
       are any */
    while (j < n2) 
    { 
        arr[k] = R[j]; 
        aux[k] = R_aux[j]; 
        aux2[k] = R_aux2[j]; 
        j++; 
        k++; 
    } 

    free(L);
    free(R);
    free(L_aux);
    free(R_aux);
    free(L_aux2);
    free(R_aux2);
} 
  
/* l is for left index and r is right index of the 
   sub-array of arr to be sorted */
void mergeSort_Llen(EDGE_ID* arr, EDGE_ID* aux, EDGE_ID* aux2, EDGE_ID l, EDGE_ID r) 
{ 
    if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        EDGE_ID m = l+(r-l)/2; 
  
        // Sort first and second halves 
        mergeSort_Llen(arr, aux, aux2, l, m); 
        mergeSort_Llen(arr, aux, aux2, m+1, r); 
  
        merge_Llen(arr, aux, aux2, l, m, r); 
    } 
} 


//////////////////////////////////////////////////////////
// MERGE SORT ALGORITHMS FOR temp_par
//////////////////////////////////////////////////////////

// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge_temp_par(PAR* arr, EDGE_ID* aux, EDGE_ID l, EDGE_ID m, EDGE_ID r) 
{ 
    EDGE_ID i, j, k; 
    EDGE_ID n1 = m - l + 1; 
    EDGE_ID n2 =  r - m; 
    //printf("\nn1, n2: %u, %u", n1, n2);
  
    /* create temp arrays */
    PAR *L, *R;
    L = (PAR*)malloc(n1*sizeof(PAR));
    R = (PAR*)malloc(n2*sizeof(PAR));

    /* create aux arrays */
    EDGE_ID *L_aux, *R_aux;
    L_aux = (EDGE_ID*)malloc(n1*sizeof(EDGE_ID));
    R_aux = (EDGE_ID*)malloc(n2*sizeof(EDGE_ID));
    //int L[n1], R[n2]; 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){ 
        L[i] = arr[l + i]; 
        L_aux[i] = aux[l + i]; 
    }

    for (j = 0; j < n2; j++) {
        R[j] = arr[m + 1+ j]; 
        R_aux[j] = aux[m + 1+ j]; 
    }
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 

    while (i < n1 && j < n2) 
    { 

          if (L[i] >= R[j]) 
          { 
              arr[k] = L[i]; 
	            aux[k] = L_aux[i];
              i++; 
          } 
          else
          { 
              arr[k] = R[j]; 
              aux[k] = R_aux[j]; 
              j++; 
          } 

          k++;

    }

  
    /* Copy the remaining elements of L[], if there 
       are any */
    while (i < n1) 
    { 
        arr[k] = L[i]; 
        aux[k] = L_aux[i]; 
        i++; 
        k++; 
    } 
  
    /* Copy the remaining elements of R[], if there 
       are any */
    while (j < n2) 
    { 
        arr[k] = R[j]; 
        aux[k] = R_aux[j]; 
        j++; 
        k++; 
    } 

    free(L);
    free(R);
    free(L_aux);
    free(R_aux);
} 
  
/* l is for left index and r is right index of the 
   sub-array of arr to be sorted */
void mergeSort_temp_par(PAR* arr, EDGE_ID* aux, EDGE_ID l, EDGE_ID r) 
{ 
    if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        EDGE_ID m = l+(r-l)/2; 
  
        // Sort first and second halves 
        mergeSort_temp_par(arr, aux, l, m); 
        mergeSort_temp_par(arr, aux, m+1, r); 
  
        merge_temp_par(arr, aux, l, m, r); 
    } 
} 


//////////////////////////////////////////////////////////
// MERGE SORT ALGORITHMS FOR in_cycles_len
//////////////////////////////////////////////////////////

// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge_incycleslen(EDGE_ID* arr, cyc_info* aux, EDGE_ID l, EDGE_ID m, EDGE_ID r) 
{ 
    EDGE_ID i, j, k; 
    EDGE_ID n1 = m - l + 1; 
    EDGE_ID n2 =  r - m; 
    //printf("\nn1, n2: %u, %u", n1, n2);
  
    /* create temp arrays */
    EDGE_ID *L, *R;
    L = (EDGE_ID*)malloc(n1*sizeof(EDGE_ID));
    R = (EDGE_ID*)malloc(n2*sizeof(EDGE_ID));

    //int L[n1], R[n2]; 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){ 
        L[i] = arr[l + i]; 
    }

    for (j = 0; j < n2; j++) {
        R[j] = arr[m + 1+ j]; 
    }
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 

    while (i < n1 && j < n2) 
    { 

          if (aux[L[i]].len <= aux[R[j]].len) 
          { 
              arr[k] = L[i]; 
              i++; 
          } 
          else
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

    free(L);
    free(R);
} 
  
/* l is for left index and r is right index of the 
   sub-array of arr to be sorted */
void mergeSort_incycleslen(EDGE_ID* arr, cyc_info* aux, EDGE_ID l, EDGE_ID r) 
{ 
    if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        EDGE_ID m = l+(r-l)/2; 
  
        // Sort first and second halves 
        mergeSort_incycleslen(arr, aux, l, m); 
        mergeSort_incycleslen(arr, aux, m+1, r); 
  
        merge_incycleslen(arr, aux, l, m, r); 
    } 
} 

//////////////////////////////////////////////////////////
// MERGE SORT ALGORITHMS FOR edges_in_cycles
//////////////////////////////////////////////////////////

// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge_edges_in_cycles(EDGE_ID* arr, cyc_info* aux, EDGE_ID l, EDGE_ID m, EDGE_ID r) 
{ 
    EDGE_ID i, j, k; 
    EDGE_ID n1 = m - l + 1; 
    EDGE_ID n2 =  r - m; 
    //printf("\nn1, n2: %u, %u", n1, n2);
  
    /* create temp arrays */
    EDGE_ID *L, *R;
    L = (EDGE_ID*)malloc(n1*sizeof(EDGE_ID));
    R = (EDGE_ID*)malloc(n2*sizeof(EDGE_ID));

    //int L[n1], R[n2]; 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){ 
        L[i] = arr[l + i]; 
    }

    for (j = 0; j < n2; j++) {
        R[j] = arr[m + 1+ j]; 
    }
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 

    while (i < n1 && j < n2) 
    { 

          if (aux[L[i]].len >= aux[R[j]].len) 
          { 
              arr[k] = L[i]; 
              i++; 
          } 
          else
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

    free(L);
    free(R);
} 
  
/* l is for left index and r is right index of the 
   sub-array of arr to be sorted */
void mergeSort_edges_in_cycles(EDGE_ID* arr, cyc_info* aux, EDGE_ID l, EDGE_ID r) 
{ 
    if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        EDGE_ID m = l+(r-l)/2; 
  
        // Sort first and second halves 
        mergeSort_edges_in_cycles(arr, aux, l, m); 
        mergeSort_edges_in_cycles(arr, aux, m+1, r); 
  
        merge_edges_in_cycles(arr, aux, l, m, r); 
    } 
} 






//////////////////////////////////////////////////////////
// MERGE SORT ALGORITHMS FOR edges_in_cycles by cycid
//////////////////////////////////////////////////////////


// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge_edges_in_cycles_bycycid(EDGE_ID* arr, EDGE_ID l, EDGE_ID m, EDGE_ID r) 
{ 
    EDGE_ID i, j, k; 
    EDGE_ID n1 = m - l + 1; 
    EDGE_ID n2 =  r - m; 
    //printf("\nn1, n2: %u, %u", n1, n2);
  
    /* create temp arrays */
    EDGE_ID *L, *R;
    L = (EDGE_ID*)malloc(n1*sizeof(EDGE_ID));
    R = (EDGE_ID*)malloc(n2*sizeof(EDGE_ID));

    //int L[n1], R[n2]; 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){ 
        L[i] = arr[l + i]; 
    }

    for (j = 0; j < n2; j++) {
        R[j] = arr[m + 1+ j]; 
    }
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 

    while (i < n1 && j < n2) 
    { 

          if (L[i] <= R[j]) 
          { 
              arr[k] = L[i]; 
              i++; 
          } 
          else
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

    free(L);
    free(R);
} 
  
/* l is for left index and r is right index of the 
   sub-array of arr to be sorted */
void mergeSort_edges_in_cycles_bycycid(EDGE_ID* arr, EDGE_ID l, EDGE_ID r) 
{ 
    if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        EDGE_ID m = l+(r-l)/2; 
  
        // Sort first and second halves 
        mergeSort_edges_in_cycles_bycycid(arr, l, m); 
        mergeSort_edges_in_cycles_bycycid(arr, m+1, r); 
  
        merge_edges_in_cycles_bycycid(arr, l, m, r); 
    } 
} 







#ifdef COMBIDX

      int H2_case1(filtration* self, coboundary_H2* V_info){
            
            //if (self->g_p_flag){
            //    printf("\nstarting H2 case 1");
            //    getchar();
            //}
      
            //EDGE_ID o_ab = V_info->triangle.key1;
            
            VERT_ID a = self->g_edges_list[2*V_info->triangle.key1];
            VERT_ID b = self->g_edges_list[2*V_info->triangle.key1+1];
            VERT_ID c = V_info->triangle.key2;
            
            VERT_ID idxa, idxb, idxc;
      
            while ((V_info->c_ptr < self->g_Neigh_len[c])\
                  && (self->g_Neighbors_e[c][V_info->c_ptr].order < V_info->triangle.key1)){
      
      
      
                  VERT_ID d = self->g_Neighbors_e[c][V_info->c_ptr].neighbor;
      
                  if ((d == a) || (d == b)){
                        V_info->c_ptr++;
                        continue;
                  }
      
      
                  if (COMB_IDX(a, d) > V_info->triangle.key1){
                        V_info->c_ptr++;
                        continue;
                  }
      
      
      
                  if (COMB_IDX(b, d) > V_info->triangle.key1){
                        V_info->c_ptr++;
                        continue;
                  }
      
      
                  V_info->low.key1 = V_info->triangle.key1;
                  V_info->low.key2 = self->g_Neighbors_e[c][V_info->c_ptr].order;
                  V_info->vertex = 0;
                  
                  return 1;
      
            }
      
            return 0;
      
            
            
      }
      
      void H2_case2 ( filtration* self, coboundary_H2* V_info){
            
            //if (self->g_p_flag){
            //    printf("\nstarting H2 case 2");
            //    getchar();
            //}
            VERT_ID idxa, idxb, idxc, idx;
            VERT_ID a, b, c;
            EDGE_ID o_ad, o_bd, o_cd;
      
            c = V_info->triangle.key2;
            
            a = self->g_edges_list[2*V_info->triangle.key1];
            b = self->g_edges_list[2*V_info->triangle.key1+1];
      
      
            
            while (1){
                  
                  
                EDGE_ID ep = self->g_n_valid_edges;
                VERT_ID d;
                int flag = -1;
      
                if (V_info->a_ptr < self->g_Neigh_len[a]){
                      
                      ep = self->g_Neighbors_e[a][V_info->a_ptr].order;
                      flag = 1;
      
                }
      
                if (V_info->b_ptr < self->g_Neigh_len[b]){
                      
                      if (self->g_Neighbors_e[b][V_info->b_ptr].order < ep){
      
                          ep = self->g_Neighbors_e[b][V_info->b_ptr].order;
                          flag = 2;
      
                      }
      
                }
      
                if (V_info->c_ptr < self->g_Neigh_len[c]){
                      
                      if (self->g_Neighbors_e[c][V_info->c_ptr].order < ep){
      
                          ep = self->g_Neighbors_e[c][V_info->c_ptr].order;
                          flag = 3;
      
                      }
      
                }
      
                if (flag == -1){
                      
                      V_info->low.key1 = ep;
                      V_info->vertex = -1;
                      return;
                }
                else if (flag == 1){
                      
                      d = self->g_Neighbors_e[a][V_info->a_ptr].neighbor;
      
                      if ((d == b) || (d == c)){
                          V_info->a_ptr++;
                          continue;
                      }
      
      
                      o_bd = COMB_IDX(b, d);
                      if (o_bd > ep){
                          V_info->a_ptr++;
                          continue;
                      }
      
                      o_cd = COMB_IDX(c, d);
                      if (o_cd > ep){
                          V_info->a_ptr++;
                          continue;
                      }
      
                      V_info->low.key1 = ep;
      
                      //o_bc = COMB_IDX(b, c);
                      V_info->low.key2 = COMB_IDX(b, c);
      
      
      
                      V_info->vertex = 1;
                      return;
      
      
      
      
                }
                else if (flag == 2){
                      
                      d = self->g_Neighbors_e[b][V_info->b_ptr].neighbor;
      
      
                      if ((d == a) || (d == c)){
                          V_info->b_ptr++;
                          continue;
                      }
      
      
                      o_ad = COMB_IDX(a, d);
                      if (o_ad > ep){
                          V_info->b_ptr++;
                          continue;
                      }
      
                      o_cd = COMB_IDX(c, d);
                      if (o_cd > ep){
                          V_info->b_ptr++;
                          continue;
                      }
      
      
                      V_info->low.key1 = ep;
      
                      //o_ac = COMB_IDX(a, c);
                      V_info->low.key2 = COMB_IDX(a, c);
      
                      V_info->vertex = 2;
                      return;
      
      
                }
                else if (flag == 3){
                      
                      d = self->g_Neighbors_e[c][V_info->c_ptr].neighbor;
      
      
      
                      if ((d == a) || (d == b)){
                          V_info->c_ptr++;
                          continue;
                      }
      
      
                      o_ad = COMB_IDX(a, d);
                      if (o_ad > ep){
                          V_info->c_ptr++;
                          continue;
                      }
      
                      o_bd = COMB_IDX(b, d);
                      if (o_bd > ep){
                          V_info->c_ptr++;
                          continue;
                      }
      
                      V_info->low.key1 = ep;
      
                      V_info->low.key2 = V_info->triangle.key1;
      
                      V_info->vertex = 3;
                      return;
      
                }
                  
            }
      
            //V_info->low.key1 = self->g_n_valid_edges;
            //V_info->vertex = -1;
            
            
      }







#else

      int H2_case1(filtration* self, coboundary_H2* V_info){
            
            //if (self->g_p_flag){
            //    printf("\nstarting H2 case 1");
            //    getchar();
            //}
      
            //EDGE_ID o_ab = V_info->triangle.key1;
            
            VERT_ID a = self->g_edges_list[2*V_info->triangle.key1];
            VERT_ID b = self->g_edges_list[2*V_info->triangle.key1+1];
            VERT_ID c = V_info->triangle.key2;
            
            VERT_ID idxa, idxb, idxc;
      
            while ((V_info->c_ptr < self->g_Neigh_len[c])\
                  && (self->g_Neighbors_e[c][V_info->c_ptr].order < V_info->triangle.key1)){
      
      
                  VERT_ID d = self->g_Neighbors_e[c][V_info->c_ptr].neighbor;
      
                  idxa = search_Neighbors(self, a, d, 0, self->g_Neigh_len[a] - 1);
      
                  if (idxa == self->g_n_vert){
                        V_info->c_ptr++;
                        continue;
                  }
      
                  if (self->g_Neighbors[a][idxa].order > V_info->triangle.key1){
                        V_info->c_ptr++;
                        continue;
                  }
      
                  idxb = search_Neighbors(self, b, d, 0, self->g_Neigh_len[b] - 1);
      
                  if (idxb == self->g_n_vert){
                        V_info->c_ptr++;
                        continue;
                  }
      
                  if (self->g_Neighbors[b][idxb].order > V_info->triangle.key1){
                        V_info->c_ptr++;
                        continue;
                  }
      
                  V_info->low.key1 = V_info->triangle.key1;
                  V_info->low.key2 = self->g_Neighbors_e[c][V_info->c_ptr].order;
                  V_info->vertex = 0;
                  return 1;
      
            }
      
            return 0;
            
            
            
      }
      
      void H2_case2 ( filtration* self, coboundary_H2* V_info){
            
            //if (self->g_p_flag){
            //    printf("\nstarting H2 case 2");
            //    getchar();
            //}
            VERT_ID idxa, idxb, idxc, idx;
            VERT_ID a, b, c;
            EDGE_ID o_ad, o_bd, o_cd;
      
            c = V_info->triangle.key2;
            
            a = self->g_edges_list[2*V_info->triangle.key1];
            b = self->g_edges_list[2*V_info->triangle.key1+1];
            
            while (1){
                  
                  
                EDGE_ID ep = self->g_n_valid_edges;
                VERT_ID d;
                int flag = -1;
      
                if (V_info->a_ptr < self->g_Neigh_len[a]){
                      
                      ep = self->g_Neighbors_e[a][V_info->a_ptr].order;
                      flag = 1;
      
                }
      
                if (V_info->b_ptr < self->g_Neigh_len[b]){
                      
                      if (self->g_Neighbors_e[b][V_info->b_ptr].order < ep){
      
                          ep = self->g_Neighbors_e[b][V_info->b_ptr].order;
                          flag = 2;
      
                      }
      
                }
      
                if (V_info->c_ptr < self->g_Neigh_len[c]){
                      
                      if (self->g_Neighbors_e[c][V_info->c_ptr].order < ep){
      
                          ep = self->g_Neighbors_e[c][V_info->c_ptr].order;
                          flag = 3;
      
                      }
      
                }
      
                if (flag == -1){
                      
                      V_info->low.key1 = ep;
                      V_info->vertex = -1;
                      return;
                }
                else if (flag == 1){
                      
                      d = self->g_Neighbors_e[a][V_info->a_ptr].neighbor;
                      idxb = search_Neighbors(self, b, d, 0, self->g_Neigh_len[b]-1);
      
                      if (idxb == self->g_n_vert){
                          V_info->a_ptr++;
                          continue;
                      }
      
                      o_bd = self->g_Neighbors[b][idxb].order;
      
                      if (o_bd > ep){
                          V_info->a_ptr++;
                          continue;
                      }
      
                      idxc = search_Neighbors(self, c, d, 0, self->g_Neigh_len[c]-1);
      
                      if (idxc == self->g_n_vert){
                          V_info->a_ptr++;
                          continue;
                      }
      
                      o_cd = self->g_Neighbors[c][idxc].order;
      
                      if (o_cd > ep){
                          V_info->a_ptr++;
                          continue;
                      }
      
                      V_info->low.key1 = ep;
      
                      idx = search_Neighbors(self, b, c, 0, self->g_Neigh_len[b]-1);
                      V_info->low.key2 = self->g_Neighbors[b][idx].order;
      
                      V_info->vertex = 1;
                      return;
      
      
      
      
                }
                else if (flag == 2){
                      
                      d = self->g_Neighbors_e[b][V_info->b_ptr].neighbor;
      
                      idxa = search_Neighbors(self, a, d, 0, self->g_Neigh_len[a]-1);
                      if (idxa == self->g_n_vert){
                          V_info->b_ptr++;
                          continue;
                      }
      
                      o_ad = self->g_Neighbors[a][idxa].order;
      
                      if (o_ad > ep){
                          V_info->b_ptr++;
                          continue;
                      }
      
                      idxc = search_Neighbors(self, c, d, 0, self->g_Neigh_len[c]-1);
      
                      if (idxc == self->g_n_vert){
                          V_info->b_ptr++;
                          continue;
                      }
      
                      o_cd = self->g_Neighbors[c][idxc].order;
      
                      if (o_cd > ep){
                          V_info->b_ptr++;
                          continue;
                      }
      
                      V_info->low.key1 = ep;
      
                      idx = search_Neighbors(self, a, c, 0, self->g_Neigh_len[a]-1);
                      V_info->low.key2 = self->g_Neighbors[a][idx].order;
      
                      V_info->vertex = 2;
                      return;
      
      
                }
                else if (flag == 3){
                      
                      d = self->g_Neighbors_e[c][V_info->c_ptr].neighbor;
      
                      idxb = search_Neighbors(self, b, d, 0, self->g_Neigh_len[b]-1);
      
                      if (idxb == self->g_n_vert){
                          V_info->c_ptr++;
                          continue;
                      }
      
                      o_bd = self->g_Neighbors[b][idxb].order;
      
                      if (o_bd > ep){
                          V_info->c_ptr++;
                          continue;
                      }
      
                      idxa = search_Neighbors(self, a, d, 0, self->g_Neigh_len[a]-1);
      
                      if (idxa == self->g_n_vert){
                          V_info->c_ptr++;
                          continue;
                      }
      
                      o_ad = self->g_Neighbors[a][idxa].order;
      
                      if (o_ad > ep){
                          V_info->c_ptr++;
                          continue;
                      }
      
                      V_info->low.key1 = ep;
      
                      //idx = search_Neighbors(self, a, c, 0, self->g_Neigh_len[a]-1);
                      V_info->low.key2 = V_info->triangle.key1;
      
                      V_info->vertex = 3;
                      return;
      
                }
                  
            }
      
            //V_info->low.key1 = self->g_n_valid_edges;
            //V_info->vertex = -1;
            
            
      }


#endif


void update_V_coH1(filtration* self, int ws_counter){

    EDGE_ID red_col = 0;

    coboundary_H1_ws* this_ws = self->g_V_ws_H1 + ws_counter;

    //if (self->g_new_debug){
    // 
    //    printf("\n ADDDDDING %d, %d, %d", this_ws->edge\
    //                     , this_ws->pivot.key1\
    //                     , this_ws->pivot.key2\
    //                     );
    //    getchar();

    //}


        //printf("\n%d, %d, %d", this_ws->edge\
        //                 , this_ws->pivot.key1\
        //                 , this_ws->pivot.key2\
        //                 );



    self->g_V_sparse_beg_ptr = self->g_V_sparse_ptr;

    if (this_ws->v_edges.last){

        if ((this_ws->v_edges.last + self->g_V_sparse_ptr) + 1 > self->g_V_sparse_max){

                    self->g_V_sparse_max = self->g_V_sparse_ptr + this_ws->v_edges.last + 10000;
                    self->g_V_sparse_H1 = (EDGE_ID*)realloc(self->g_V_sparse_H1\
                                                        , self->g_V_sparse_max*sizeof(EDGE_ID));

        }




        if (this_ws->v_edges.last > 1){

#ifdef VREDUCE1

            sorter8_tim_sort(this_ws->v_edges.o_ab, this_ws->v_edges.last);          
          
            int coeff = 1;

            for (EDGE_ID vv = 0; vv < this_ws->v_edges.last-1; vv++){
                  
                if (this_ws->v_edges.o_ab[vv] == this_ws->v_edges.o_ab[vv+1])
                {
                    coeff = 1 - coeff;
                }
                else{
                    if (coeff){
                        self->g_V_sparse_H1[self->g_V_sparse_ptr++] = this_ws->v_edges.o_ab[vv];
                    }

                    coeff = 1;
                }
                  
                  
            }

            if (coeff){
                 self->g_V_sparse_H1[self->g_V_sparse_ptr++] = this_ws->v_edges.o_ab[this_ws->v_edges.last-1];
            }

#else

            for (EDGE_ID vv = 0; vv < this_ws->v_edges.last; vv++){

                 self->g_V_sparse_H1[self->g_V_sparse_ptr++] = this_ws->v_edges.o_ab[vv];

            }

#endif


        }
        else if (this_ws->v_edges.last == 1){

            self->g_V_sparse_H1[self->g_V_sparse_ptr++] = this_ws->v_edges.o_ab[0];
            
        }


        //if (this_ws->edge == self->g_debug_edge){

        //    printf("\nAfter adding to V sparse ");
        //    for (EDGE_ID bb = self->g_V_sparse_beg_ptr; bb < self->g_V_sparse_ptr; bb++){

        //        printf("%d, ", self->g_V_sparse_H1[bb]);
        //    }
        //  getchar();
        //}
        // All have been added
        this_ws->v_edges.last = 0;


        if ((self->g_V_sparse_ptr - self->g_V_sparse_beg_ptr) > 0){

              red_col = self->g_V_col_indices_ptr;

              if (self->g_V_col_indices_ptr+1 == self->g_V_col_indices_max){
                    
                    self->g_V_col_indices_max += 1000;
                    self->g_V_col_indices = (EDGE_ID*)realloc(self->g_V_col_indices
                                                          , self->g_V_col_indices_max*sizeof(EDGE_ID));
              
              }


              self->g_V_col_indices[self->g_V_col_indices_ptr] = self->g_V_sparse_beg_ptr;
              self->g_V_col_indices[self->g_V_col_indices_ptr+1] = self->g_V_sparse_ptr;

              self->g_V_col_indices_ptr++;


        }

    }



#ifdef COH1DEBUG
        if (this_ws->edge == self->g_debug_edge){

            printf("\n%d, %d, %d, %d: "\
                                            , this_ws->pivot.key1\
                                            , this_ws->pivot.key2\
                                            , this_ws->edge\
                                            , self->g_V_sparse_ptr - self->g_V_sparse_beg_ptr\
                                            );
            for (EDGE_ID mm = self->g_V_sparse_beg_ptr; mm < self->g_V_sparse_ptr; mm++){
                printf("%d, ", self->g_V_sparse_H1[mm]);
            }
            getchar();
              
        }
#endif

#ifdef VDEBUG
        if (self->g_V_sparse_ptr - self->g_V_sparse_beg_ptr > 0){

            printf("\n%d, %d, %d, %d"\
                                            , this_ws->pivot.key1\
                                            , this_ws->pivot.key2\
                                            , this_ws->edge\
                                            , self->g_V_sparse_ptr - self->g_V_sparse_beg_ptr\
                                            );
            //getchar();

        }
#endif

    // ADDING THE LOW

    if (!self->g_H1_cohom_pivots_max_len[this_ws->pivot.key1]){
          
        self->g_H1_cohom_pivots_max_len[this_ws->pivot.key1] = 2;
        self->g_H1_cohom_pivots[this_ws->pivot.key1] = \
                               (H1_cohom_pivots*)malloc(self->g_H1_cohom_pivots_max_len[this_ws->pivot.key1]*sizeof(H1_cohom_pivots));

    }

    if (self->g_H1_cohom_pivots_len[this_ws->pivot.key1]\
                            == self->g_H1_cohom_pivots_max_len[this_ws->pivot.key1]){
          
          self->g_H1_cohom_pivots_max_len[this_ws->pivot.key1] += 5;
          self->g_H1_cohom_pivots[this_ws->pivot.key1] = (H1_cohom_pivots*)realloc( \
                          self->g_H1_cohom_pivots[this_ws->pivot.key1] \
                          , self->g_H1_cohom_pivots_max_len[this_ws->pivot.key1]*sizeof(H1_cohom_pivots));
          //self->g_cohom_ALL_pivots_len += 5;

    }
          


    EDGE_ID old_ptr = self->g_H1_cohom_pivots_len[this_ws->pivot.key1];
    EDGE_ID new_ptr = self->g_H1_cohom_pivots_len[this_ws->pivot.key1];

    while (old_ptr){
          
          old_ptr--;
          
          if (self->g_H1_cohom_pivots[this_ws->pivot.key1][old_ptr].key2 > this_ws->pivot.key2){

                self->g_H1_cohom_pivots[this_ws->pivot.key1][new_ptr--] =\
                                                       self->g_H1_cohom_pivots[this_ws->pivot.key1][old_ptr];
                continue;

          }
          break;

    }

    //printf("\nAdding pivot (%d, %d) for edge %d"\
    //                                          , this_ws->pivot.key1\
    //                                          , this_ws->pivot.key2\
    //                                          , this_ws->edge\
    //                                          );

    self->g_H1_cohom_pivots[this_ws->pivot.key1][new_ptr].key2 = this_ws->pivot.key2;
    self->g_H1_cohom_pivots[this_ws->pivot.key1][new_ptr].col_idx = red_col;
                            
    self->g_H1_cohom_pivots[this_ws->pivot.key1][new_ptr].bndry = this_ws->edge;

    self->g_H1_cohom_pivots_len[this_ws->pivot.key1]++;

    
    // PERS PAIRS
    // Add non-zero barcodes
        
    PAR birth = self->g_edge_parameter[this_ws->edge];
    PAR death = self->g_edge_parameter[this_ws->pivot.key1];
    if (birth != death){

           //printf("\nNon trivial pers pair (%f, %f)", birth, death);
           

#ifdef DEBUGPIVOTS
           printf("\nBirth, death (%lf, %lf)", birth, death);
           printf("\n%d at pair (%d, %d)", this_ws->edge\
                                                , this_ws->pivot.key1\
                                                , this_ws->pivot.key2);
           getchar();
          
#endif
           //if (birth > death){
           //  
           //}


          if (self->g_H1_pers_pairs_len+2 == self->g_H1_pers_pairs_max_len){
                self->g_H1_pers_pairs_max_len += 1000;
                self->g_H1_pers_pairs = (PAR*)realloc(self->g_H1_pers_pairs\
                                              , self->g_H1_pers_pairs_max_len*sizeof(PAR));
          
          }
          self->g_H1_pers_pairs[self->g_H1_pers_pairs_len++] = birth;
          self->g_H1_pers_pairs[self->g_H1_pers_pairs_len++] = death;
          
    }


}



void deallocator(filtration* self){

      struct timespec start_wall_clock;
      struct timespec finish_wall_clock;
      double timer;

      if (!self->g_suppress_output){
          clock_gettime(CLOCK_MONOTONIC, &start_wall_clock);
      }

      free(self->filename);

      //free(self->g_homH1_cycles_file);

      // Deallocate edges
      free(self->g_edge_parameter);
      free(self->g_edges_list);


      for (EDGE_ID mm = 0; mm < self->g_n_valid_edges; mm++){
          
            if (self->g_dim_lim > 0){

                if (self->g_H1_cohom_pivots_max_len[mm]){
                      
                      free(self->g_H1_cohom_pivots[mm]);

                }

                if (self->g_dim_lim > 1){

                      if (self->g_H2_cohom_pivots_max_len[mm]){
                            
                            free(self->g_H2_cohom_pivots[mm]);

                      }
                      
                }
                  
            }

      }



      // Deallocate Neighbors
      for (EDGE_ID mm = 0; mm < self->g_n_vert; mm++){
          
            if (self->g_Neigh_len[mm]){

                free(self->g_Neighbors[mm]);
                free(self->g_Neighbors_e[mm]);

            }

      }

      free(self->g_Neighbors);
      free(self->g_Neighbors_e);
      free(self->g_Neigh_len);



      // Deallocate R0
      free(self->g_pivots_H0);
      free(self->g_R_sparse_H0);
      free(self->g_R_col_indices_H0);
      free(self->g_edges_with_pivots_H0);

#ifdef SAVEPD
      free(self->g_H0_pers_file);
#endif

      if (self->g_dim_lim > 0){
          
            free(self->g_coH1_all_lows);
            free(self->g_H1_cohom_pivots);
            free(self->g_H1_cohom_pivots_len);
            free(self->g_H1_cohom_pivots_max_len);
#ifdef SAVEPD
            free(self->g_H1_pers_file);
#endif

#ifdef SAVEV
            free(self->g_coH1_V_file);
#endif
            free(self->g_H1_pers_pairs);

            free(self->g_V_col_indices);
            
            if (self->g_dim_lim > 1){
                  
                  free(self->g_H2_cohom_pivots);
                  free(self->g_H2_cohom_pivots_len);
                  free(self->g_H2_cohom_pivots_max_len);
#ifdef SAVEPD
                  free(self->g_H2_pers_file);
#endif

#ifdef SAVEV
                  free(self->g_coH2_V_file);
#endif
                  free(self->g_H2_pers_pairs);
                  
            }

            
      }

      if (self->g_compute_cycles){
            
            free(self->g_H1_undead);
            if (self->g_dim_lim > 1){
                free(self->g_H2_undead);
            }

      }



      free(self);
      


      if (!self->g_suppress_output){
        clock_gettime(CLOCK_MONOTONIC, &finish_wall_clock);
        timer = (finish_wall_clock.tv_sec - start_wall_clock.tv_sec);
        timer += (finish_wall_clock.tv_nsec - start_wall_clock.tv_nsec) / 1000000000.0;
        printf("\nTime taken to deallocate: %lf", timer);
      }

}





void insert_in_implicit_v(filtration* self, int ws_counter, coboundary_H1* phi, int flag_next){


      if (phi->low.key1 == self->g_n_valid_edges){
          return;
      }

      coboundary_H1_ws* this_ws = self->g_V_ws_H1 + ws_counter;
      

      if (phi->low.key1 == this_ws->keys1[this_ws->k1_ptr].k1){


            if (this_ws->keys1[this_ws->k1_ptr].last ==\
                                            this_ws->keys1[this_ws->k1_ptr].max_len){

                  self->g_V_ws_H1[ws_counter].keys1[self->g_V_ws_H1[ws_counter].k1_ptr].max_len += 10;
                  self->g_V_ws_H1[ws_counter].keys1[self->g_V_ws_H1[ws_counter].k1_ptr].keys2 = \
                                                        (implicit_keys2*)realloc\
                                                             (self->g_V_ws_H1[ws_counter].keys1[self->g_V_ws_H1[ws_counter].k1_ptr].keys2\
                                                            , self->g_V_ws_H1[ws_counter].keys1[self->g_V_ws_H1[ws_counter].k1_ptr].max_len\
                                                                  *sizeof(implicit_keys2));
                  
            }


            EDGE_ID mm = this_ws->keys1[this_ws->k1_ptr].last;

            int compare;


            while (1){


                //int compare = compare_implicit(this_ws->keys1[this_ws->k1_ptr].keys2[mm-1], *phi);


                if (this_ws->keys1[this_ws->k1_ptr].keys2[mm-1].k2 < phi->low.key2) compare = 0;
                else if (this_ws->keys1[this_ws->k1_ptr].keys2[mm-1].k2 > phi->low.key2) compare = 1;
                else{

                      if (this_ws->keys1[this_ws->k1_ptr].keys2[mm-1].o_ab < phi->o_ab) compare = 0;
                      else compare = 1;

                }

                  
                if (compare){

                      this_ws->keys1[this_ws->k1_ptr].keys2[mm] =\
                                                        this_ws->keys1[this_ws->k1_ptr].keys2[mm-1];
                }
                else{
                    
                      this_ws->keys1[this_ws->k1_ptr].keys2[mm].k2 = phi->low.key2; 
                  
                      this_ws->keys1[this_ws->k1_ptr].keys2[mm].k2 = phi->low.key2; 
                      this_ws->keys1[this_ws->k1_ptr].keys2[mm].o_ab = phi->o_ab; 
                      this_ws->keys1[this_ws->k1_ptr].keys2[mm].a_ptr = phi->a_ptr; 
                      this_ws->keys1[this_ws->k1_ptr].keys2[mm].b_ptr = phi->b_ptr; 
                      this_ws->keys1[this_ws->k1_ptr].keys2[mm].flag_next = flag_next; 

                      this_ws->keys1[this_ws->k1_ptr].last++;

                      return;

                }

                mm--;

                //// ERROR CHECKING, REMOVE LATER
                //if (!mm){
                //    printf("\nk2_ptr %d", v_implicit->k2_ptr);
                //    printf("\nADDING %d:(%d, %d) to ", phi->o_ab, phi->low.key1, phi->low.key2);
                //    print_v_implicit(self);
                //    exit(0);
                //    
                //}

                if (mm == this_ws->k2_ptr){

                      this_ws->keys1[this_ws->k1_ptr].keys2[mm].k2 = phi->low.key2; 
                      this_ws->keys1[this_ws->k1_ptr].keys2[mm].o_ab = phi->o_ab; 
                      this_ws->keys1[this_ws->k1_ptr].keys2[mm].a_ptr = phi->a_ptr; 
                      this_ws->keys1[this_ws->k1_ptr].keys2[mm].b_ptr = phi->b_ptr; 
                      this_ws->keys1[this_ws->k1_ptr].keys2[mm].flag_next = flag_next; 

                      this_ws->keys1[this_ws->k1_ptr].last++;


                      return;
                    

                }

            }
            
            
      }

      for (EDGE_ID mm = 0; mm < this_ws->last; mm++){

            
            if (this_ws->keys1[mm].k1 == phi->low.key1){
                
                  
                  //check_space_implicit_keys2(&(v_implicit->keys1[mm]));

                  if (this_ws->keys1[mm].last ==\
                                                  this_ws->keys1[mm].max_len){

                        self->g_V_ws_H1[ws_counter].keys1[mm].max_len += 10;
                        self->g_V_ws_H1[ws_counter].keys1[mm].keys2 = (implicit_keys2*)realloc\
                                                                 (self->g_V_ws_H1[ws_counter].keys1[mm].keys2\
                                                                , self->g_V_ws_H1[ws_counter].keys1[mm].max_len*sizeof(implicit_keys2));
                        
                  }

                  this_ws->keys1[mm].flag_empty = 0;

                  this_ws->keys1[mm].keys2[this_ws->keys1[mm].last].k2 = phi->low.key2;
                  this_ws->keys1[mm].keys2[this_ws->keys1[mm].last].o_ab = phi->o_ab;
                  this_ws->keys1[mm].keys2[this_ws->keys1[mm].last].a_ptr = phi->a_ptr;
                  this_ws->keys1[mm].keys2[this_ws->keys1[mm].last].b_ptr = phi->b_ptr;
                  this_ws->keys1[mm].keys2[this_ws->keys1[mm].last].flag_next = flag_next;

                  this_ws->keys1[mm].last++;

                  return;

                  
            }
            
            
      }

      //if (self->g_new_debug){
      //    printf("\nBefore inserting c4");
      //    print_v_implicit(self);
      //    getchar();

      //}
      //check_space_implicit_keys1(v_implicit);

      if (self->g_V_ws_H1[ws_counter].last == self->g_V_ws_H1[ws_counter].max_len){
            
            EDGE_ID mm = self->g_V_ws_H1[ws_counter].max_len;

            self->g_V_ws_H1[ws_counter].max_len += 10;
            self->g_V_ws_H1[ws_counter].keys1 = (implicit_keys1*)realloc(self->g_V_ws_H1[ws_counter].keys1\
                                                                  , self->g_V_ws_H1[ws_counter].max_len*sizeof(implicit_keys1));

            while (mm < self->g_V_ws_H1[ws_counter].max_len){
                  
                  this_ws->keys1[mm].flag_empty = 1;
                  this_ws->keys1[mm].max_len = 10;
                  this_ws->keys1[mm].last = 0;
                  self->g_V_ws_H1[ws_counter].keys1[mm].keys2 = (implicit_keys2*)malloc(10*sizeof(implicit_keys2));

                  mm++;
                
                  
            }


      }




      this_ws->keys1[this_ws->last].flag_empty = 0;
      this_ws->keys1[this_ws->last].k1 = phi->low.key1;
      this_ws->keys1[this_ws->last].keys2[0].k2 = phi->low.key2;
      this_ws->keys1[this_ws->last].keys2[0].o_ab = phi->o_ab;
      this_ws->keys1[this_ws->last].keys2[0].a_ptr = phi->a_ptr;
      this_ws->keys1[this_ws->last].keys2[0].b_ptr = phi->b_ptr;
      this_ws->keys1[this_ws->last].keys2[0].flag_next = flag_next;
      this_ws->keys1[this_ws->last].last = 1;

      this_ws->last++;


      return;



      
        
}



void print_v_implicit(filtration* self, int ws_counter){
      

    if (self->g_V_ws_H1[ws_counter].edge == self->g_debug_edge){


          EDGE_ID k1_ptr = 0;

          if (k1_ptr == self->g_V_ws_H1[ws_counter].last){
            printf("\nv implicit is empty");
            return;
          }

          //EDGE_ID k2_ptr = self->g_v_implicit.k2_ptr;
          EDGE_ID k2_ptr = 0;


          while (k1_ptr < self->g_V_ws_H1[ws_counter].last){

              //printf("\n%d, %d, last %d, flag_e %d ", k1_ptr\
              //                                      , self->g_v_implicit.keys1[k1_ptr].k1\
              //                                      , self->g_v_implicit.keys1[k1_ptr].last\
              //                                      , self->g_v_implicit.keys1[k1_ptr].flag_empty\
              //                                      );

              if (self->g_V_ws_H1[ws_counter].keys1[k1_ptr].flag_empty){
                printf("empty", k1_ptr, self->g_V_ws_H1[ws_counter].keys1[k1_ptr].k1);
                k1_ptr++;
                continue;
              }


              //if (k1_ptr == self->g_v_implicit.k1_ptr){
              //    printf("\nk1_ptr is %d, k2_ptr is %d", self->g_v_implicit.k1_ptr\
              //                                         , self->g_v_implicit.k2_ptr);
              //}
                
              //printf("\nentries in %d are %d", k1_ptr, self->g_v_implicit.keys1[k1_ptr].last);
              
              printf("\n");
              printf("idx %d, last %d:: ", k1_ptr, self->g_V_ws_H1[ws_counter].keys1[k1_ptr].last);

              while (k2_ptr < self->g_V_ws_H1[ws_counter].keys1[k1_ptr].last){
                  printf("%d:(%d, %d):%d,  ", self->g_V_ws_H1[ws_counter].keys1[k1_ptr].keys2[k2_ptr].o_ab\
                                         , self->g_V_ws_H1[ws_counter].keys1[k1_ptr].k1\
                                         , self->g_V_ws_H1[ws_counter].keys1[k1_ptr].keys2[k2_ptr].k2\
                                         , self->g_V_ws_H1[ws_counter].keys1[k1_ptr].keys2[k2_ptr].flag_next\
                                         );
                  k2_ptr++;
              }
              k2_ptr = 0;
              k1_ptr++;

                
          }

    }
      
}


void coH2_print_v_implicit(filtration* self, int ws_counter){
      


          printf("\nk1ptr is %d, k2ptr is %d", self->g_V_ws_H2[ws_counter].k1_ptr\
                                             , self->g_V_ws_H2[ws_counter].k2_ptr);

          EDGE_ID k1_ptr = 0;

          if (k1_ptr == self->g_V_ws_H2[ws_counter].last){
            printf("\nv implicit is empty");
            return;
          }

          //EDGE_ID k2_ptr = self->g_v_implicit.k2_ptr;
          EDGE_ID k2_ptr = 0;


          while (k1_ptr < self->g_V_ws_H2[ws_counter].last){

              //printf("\n%d, %d, last %d, flag_e %d ", k1_ptr\
              //                                      , self->g_v_implicit.keys1[k1_ptr].k1\
              //                                      , self->g_v_implicit.keys1[k1_ptr].last\
              //                                      , self->g_v_implicit.keys1[k1_ptr].flag_empty\
              //                                      );

              if (self->g_V_ws_H2[ws_counter].keys1[k1_ptr].flag_empty){
                printf("empty", k1_ptr, self->g_V_ws_H2[ws_counter].keys1[k1_ptr].k1);
                k1_ptr++;
                continue;
              }


              //if (k1_ptr == self->g_v_implicit.k1_ptr){
              //    printf("\nk1_ptr is %d, k2_ptr is %d", self->g_v_implicit.k1_ptr\
              //                                         , self->g_v_implicit.k2_ptr);
              //}
                
              //printf("\nentries in %d are %d", k1_ptr, self->g_v_implicit.keys1[k1_ptr].last);
              
              printf("\n");
              printf("idx %d, last %d:: ", k1_ptr, self->g_V_ws_H2[ws_counter].keys1[k1_ptr].last);

              while (k2_ptr < self->g_V_ws_H2[ws_counter].keys1[k1_ptr].last){
                  printf("(%d, %d):%d,  "\
                                         , self->g_V_ws_H2[ws_counter].keys1[k1_ptr].keys2[k2_ptr].o_abc.key1\
                                         , self->g_V_ws_H2[ws_counter].keys1[k1_ptr].keys2[k2_ptr].o_abc.key2\
                                         , self->g_V_ws_H2[ws_counter].keys1[k1_ptr].keys2[k2_ptr].flag_next\
                                         );
                  k2_ptr++;
              }
              k2_ptr = 0;
              k1_ptr++;

                                         //, self->g_V_ws_H2[ws_counter].keys1[k1_ptr].k1\
                                         //, self->g_V_ws_H2[ws_counter].keys1[k1_ptr].keys2[k2_ptr].k2\
                
          }

      
}


void* reduce_with_complex_coH1(void* arg){
        
      filtration* self = arg;


      pthread_mutex_lock(&(self->g_thread_lock));

      int tid = ++self->g_thread_id;


      pthread_mutex_unlock(&(self->g_thread_lock));

      for (;;){

          pthread_mutex_lock(&(self->g_thread_lock));

          self->g_sleeping_threads++;
          self->g_processed_threads++;

          if (self->g_sleeping_threads == self->g_cpu_count){

              pthread_cond_signal(&(self->g_start_boss));
          }


          pthread_cond_wait(&(self->g_start_workers), &(self->g_thread_lock));

          if (self->g_delete_threads){
            pthread_mutex_unlock(&(self->g_thread_lock));
            pthread_exit(NULL);
          }

          self->g_sleeping_threads--;

          pthread_mutex_unlock(&(self->g_thread_lock));


          for (int ws_counter = self->g_jobs[tid - 1]; ws_counter < self->g_jobs[tid]; ws_counter++){


              coboundary_H1_ws* this_ws = self->g_V_ws_H1 + ws_counter;


              if (!this_ws->flag_non_empty){
                  // We are sure that we will exit only if there is no reduction
                  // required with existing complex or with trivial pair
                  this_ws->flag_red_w_complex = 0;
                  this_ws->flag_red_w_trivial = 0;
                  this_ws->flag_append_to_complex = 0;
                  continue;
              }

              if (this_ws->flag_append_to_complex){
                  continue;
              }

              // If being processed for the first time...
              if (this_ws->flag_first){

                  this_ws->flag_first = 0;

                  if ((self->g_coH1_all_lows[this_ws->pivot.key1].low.key1 == this_ws->pivot.key1)\
                      && (self->g_coH1_all_lows[this_ws->pivot.key1].low.key2 == this_ws->pivot.key2)){


                        this_ws->reduce_w_bndry = this_ws->pivot.key1;
                        this_ws->V_col_idx = 0;

                        this_ws->flag_red_w_trivial = 1;


                  }
                  else{

                        // If this low is not a pivot
                        if (!self->g_H1_cohom_pivots_len[this_ws->pivot.key1]){
#ifdef COH1DEBUG
                              if (this_ws->edge == self->g_debug_edge ){
                                  printf("\n(%d, %d) pivot of %d is not a pivot 1"\
                                                          , this_ws->pivot.key1\
                                                          , this_ws->pivot.key2\
                                                          , this_ws->edge\
                                                          );
                              } 
#endif

                              this_ws->flag_append_to_complex = 1;
                              continue;

                        }
                        else{

                              EDGE_ID idx = search_H1_cohom_pivots(self->g_H1_cohom_pivots[this_ws->pivot.key1]\
                                                        , 0 \
                                                        , self->g_H1_cohom_pivots_len[this_ws->pivot.key1] - 1\
                                                        , this_ws->pivot.key2 \
                                                        , self->g_n_valid_edges);

                              // If this low is not a pivot
                              if (idx == self->g_n_valid_edges){

#ifdef COH1DEBUG
                                    if (this_ws->edge == self->g_debug_edge ){
                                        printf("\n(%d, %d) pivot of %d is not a pivot 2"\
                                                                , this_ws->pivot.key1\
                                                                , this_ws->pivot.key2\
                                                                , this_ws->edge\
                                                                );
                                    } 
#endif
                                    this_ws->flag_append_to_complex = 1;
                                    continue;


                              }
                              else{

                                  this_ws->flag_red_w_complex = 1;
                                  this_ws->reduce_w_bndry = self->g_H1_cohom_pivots[this_ws->pivot.key1][idx].bndry;
                                  this_ws->V_col_idx = self->g_H1_cohom_pivots[this_ws->pivot.key1][idx].col_idx;

                              }
                        
                        }

                  }

              }

              if ((!this_ws->flag_red_w_trivial) && (!this_ws->flag_red_w_complex)){
                  this_ws->flag_append_to_complex = 1;
                  continue;
              }


              // We know that parallel will end only when there are no more red. to be with trivial and complex
              this_ws->flag_red_w_complex = 0;
              this_ws->flag_red_w_trivial = 0;
              this_ws->flag_append_to_complex = 1;
                    

              while(1){

                    EDGE_ID check_len = this_ws->v_edges.last + 1;

                    if (this_ws->V_col_idx){

                          check_len += self->g_V_col_indices[this_ws->V_col_idx+1] -\
                                                        self->g_V_col_indices[this_ws->V_col_idx];

                    }

                    if (check_len > this_ws->v_edges.max_len){
                        this_ws->v_edges.max_len = check_len + 100;
                        self->g_V_ws_H1[ws_counter].v_edges.o_ab =\
                                                (EDGE_ID*)realloc(self->g_V_ws_H1[ws_counter].v_edges.o_ab\
                                                                , this_ws->v_edges.max_len*sizeof(EDGE_ID));

                    }

                    this_ws->v_edges.o_ab[this_ws->v_edges.last++] = this_ws->reduce_w_bndry;

                    coboundary_H1 ttemp;

                    ttemp.o_ab = this_ws->reduce_w_bndry;

                    //if (this_ws->edge == self->g_debug_edge){
                    //    printf("\n%d: Appending to v edge in parallel %d", this_ws->edge, ttemp.o_ab);
                    //    getchar();
                    //}

                    find_H1_cohom_greater(self, &(ttemp), &(this_ws->pivot));

                    insert_in_implicit_v(self, ws_counter, &(ttemp), 1);


                    // IF the V was recorded, add the bndries
                    if (this_ws->V_col_idx){

                        // We have to cycle through the col in V and add all the other boundary columns for reduction

                        EDGE_ID start = self->g_V_col_indices[this_ws->V_col_idx];
                        EDGE_ID end = self->g_V_col_indices[this_ws->V_col_idx+1];


                            for (EDGE_ID mm = start; mm < end; mm++){
                                  
                                this_ws->v_edges.o_ab[this_ws->v_edges.last++] = self->g_V_sparse_H1[mm];

                                ttemp.o_ab = self->g_V_sparse_H1[mm];

                                //if (this_ws->edge == self->g_debug_edge){
                                //    printf(", %d", ttemp.o_ab);
                                //}
                                // Find the first low greater than or equal pivot
                                find_H1_cohom_greater(self, &(ttemp), &(this_ws->pivot));

                                insert_in_implicit_v(self, ws_counter,  &(ttemp), 1);

                            }

                    }

                    //if (this_ws->edge == self->g_debug_edge){
                    //    getchar();
                    //}

                    reduce_hash_table_coH1(self, ws_counter);

                    //if (this_ws->edge == self->g_debug_edge){
                    //      printf("\nPivot after reduction in parallel is (%d, %d)", this_ws->pivot.key1\
                    //                                                  , this_ws->pivot.key2);
                    //}

                    if (!this_ws->flag_non_empty){
                        break;
                    }

                    // Check with trivial pair
                    if ((self->g_coH1_all_lows[this_ws->pivot.key1].low.key1 == this_ws->pivot.key1)\
                        && (self->g_coH1_all_lows[this_ws->pivot.key1].low.key2 == this_ws->pivot.key2)){


                            this_ws->reduce_w_bndry = this_ws->pivot.key1;
                            this_ws->V_col_idx = 0;
                            continue;

                    }

                    // If this low is not a pivot
                    if (!self->g_H1_cohom_pivots_len[this_ws->pivot.key1]){
                        break;
                    }


                    EDGE_ID idx = search_H1_cohom_pivots(self->g_H1_cohom_pivots[this_ws->pivot.key1]\
                                        , 0 \
                                        , self->g_H1_cohom_pivots_len[this_ws->pivot.key1] - 1\
                                        , this_ws->pivot.key2 \
                                        , self->g_n_valid_edges);




                    if (idx == self->g_n_valid_edges){
                      break;
                    }

                    this_ws->reduce_w_bndry = self->g_H1_cohom_pivots[this_ws->pivot.key1][idx].bndry;
                    this_ws->V_col_idx = self->g_H1_cohom_pivots[this_ws->pivot.key1][idx].col_idx;

              }

          }

      }
        
}


void reduce_with_self_coH1(filtration* self){


      // Now we have to reduce

      for (int ws_counter = 0; ws_counter < self->g_ws_counter; ws_counter++){

            coboundary_H1_ws* this_ws = self->g_V_ws_H1 + ws_counter;

            // If empty, then continue and don't append to complex
            if (!this_ws->flag_non_empty){
                  //this_ws->flag_append_to_complex = 0;
                  continue;
            }



            int m = 0;

            // Keep reducing if reduce with complex flag is 0 and reduce with trivial flag is 0
            while((m < ws_counter)\
                && (!this_ws->flag_red_w_complex)\
                && (!this_ws->flag_red_w_trivial)){

                  

                
                coboundary_H1_ws* m_ws = self->g_V_ws_H1 + m;

                // If m is empty, continue
                if (!m_ws->flag_non_empty){
                      m++;
                      continue;
                }


                int compare;

                if (m_ws->pivot.key1 > this_ws->pivot.key1) compare = 1;
                else if (m_ws->pivot.key1 < this_ws->pivot.key1) compare = 0;
                else{
                  if (m_ws->pivot.key2 > this_ws->pivot.key2) compare = 1;
                  else if (m_ws->pivot.key2 < this_ws->pivot.key2) compare = 0;
                  else compare = -1;
                }

                // If pivot of m is higher than pivot of ws_counter
                // then we don't care
                if (compare == 1){
                      m++;
                      continue;
                }

                // If pivot of m is lower than pivot of ws_counter
                // then if m has to be reduced, we have to hold ws_counter
                if (compare == 0){

                      if (m_ws->flag_red_w_complex || m_ws->flag_red_w_trivial){
                            
                            this_ws->flag_append_to_complex = 0;
                            break;
                      }
                      m++;
                      continue;
                }

                // At this point they have same low
                if (m_ws->flag_red_w_complex || m_ws->flag_red_w_trivial){
                      this_ws->flag_append_to_complex = 0;
                      //m++;
                      break;
                      //continue;
                }

                //printf("\nin serial reducing %d with %d", this_ws->edge, m_ws->edge);
                //getchar();

                // Merge m and this_ws
                //




                
                // Merge v_edges
                if (this_ws->v_edges.last + m_ws->v_edges.last > this_ws->v_edges.max_len - 1){
                      
                      this_ws->v_edges.max_len += m_ws->v_edges.last + 100;

                      self->g_V_ws_H1[ws_counter].v_edges.o_ab = (EDGE_ID*)realloc(\
                                                 self->g_V_ws_H1[ws_counter].v_edges.o_ab\
                                                 , this_ws->v_edges.max_len*sizeof(EDGE_ID));

                }

                // Add the original edge
                
                  
                //if (this_ws->edge == self->g_debug_edge){
                //    printf("\nAdding to v edge in serial %d", m_ws->edge);
                //}

                this_ws->v_edges.o_ab[this_ws->v_edges.last++] = m_ws->edge;
                
                for (EDGE_ID bb = 0; bb < m_ws->v_edges.last; bb++){

                      //if (this_ws->edge == self->g_debug_edge){
                      //    printf(", %d", m_ws->v_edges.o_ab[bb]);
                      //}

                      this_ws->v_edges.o_ab[this_ws->v_edges.last++] =\
                                                    m_ws->v_edges.o_ab[bb];
                      
                }

                // Merge hash tables
                coboundary_H1 ttemp;

                for (EDGE_ID bb = 0; bb < m_ws->last; bb++){

                      EDGE_ID m_start = 0;

                      if (bb == m_ws->k1_ptr){
                          m_start = m_ws->k2_ptr;
                      }
                      
                      ttemp.low.key1 = m_ws->keys1[bb].k1;

                      for (EDGE_ID mm = m_start; mm < m_ws->keys1[bb].last; mm++){
                            

                            ttemp.low.key2 = m_ws->keys1[bb].keys2[mm].k2;
                            ttemp.o_ab = m_ws->keys1[bb].keys2[mm].o_ab;
                            ttemp.a_ptr = m_ws->keys1[bb].keys2[mm].a_ptr;
                            ttemp.b_ptr = m_ws->keys1[bb].keys2[mm].b_ptr;

                            insert_in_implicit_v(self, ws_counter, &(ttemp)\
                                                  , m_ws->keys1[bb].keys2[mm].flag_next);
                            
                            
                      }
                      
                }



                // Now reduce


                reduce_hash_table_coH1(self, ws_counter);

                //if (this_ws->edge == self->g_debug_edge){
                //      printf("\nPivot after reduction in serial is (%d, %d)", this_ws->pivot.key1\
                //                                                  , this_ws->pivot.key2);
                //}

                if (!this_ws->flag_non_empty){
                  break;
                }

                // Check with trivial pair
                if ((self->g_coH1_all_lows[this_ws->pivot.key1].low.key1 == this_ws->pivot.key1)\
                    && (self->g_coH1_all_lows[this_ws->pivot.key1].low.key2 == this_ws->pivot.key2)){

                        this_ws->flag_red_w_trivial = 1;
                        this_ws->flag_red_w_complex = 0;
                        this_ws->flag_append_to_complex = 0;

                        this_ws->reduce_w_bndry = this_ws->pivot.key1;
                        this_ws->V_col_idx = 0;
                        break;

                }


                 // If this low is not a pivot
                 if (self->g_H1_cohom_pivots_len[this_ws->pivot.key1]){

                    EDGE_ID idx = search_H1_cohom_pivots(self->g_H1_cohom_pivots[this_ws->pivot.key1]\
                                        , 0 \
                                        , self->g_H1_cohom_pivots_len[this_ws->pivot.key1] - 1\
                                        , this_ws->pivot.key2 \
                                        , self->g_n_valid_edges);

                    if (idx != self->g_n_valid_edges){

                        this_ws->flag_red_w_trivial = 0;
                        this_ws->flag_red_w_complex = 1;
                        this_ws->flag_append_to_complex = 0;
                          
                        this_ws->reduce_w_bndry = self->g_H1_cohom_pivots[this_ws->pivot.key1][idx].bndry;
                        this_ws->V_col_idx = self->g_H1_cohom_pivots[this_ws->pivot.key1][idx].col_idx;

                        break;

                    }

                      
                 }


                 // Reset m after single reduction
                 m = 0;

                
            }

            
      }
      
      
    
}



void reduce_hash_table_coH1(filtration* self, int ws_counter){
      
      
      // Now we have to reduce
      int coeff = 1;
      
      coboundary_H1_ws* this_ws = self->g_V_ws_H1 + ws_counter;

      coboundary_H1 ttemp;
      
      EDGE_ID* k1_ptr = &(this_ws->k1_ptr);
      EDGE_ID* k2_ptr = &(this_ws->k2_ptr);
      
      while (1){
            
            
            if (this_ws->keys1[*k1_ptr].last == 1){
                this_ws->pivot.key1 = this_ws->keys1[*k1_ptr].k1;
                this_ws->pivot.key2 = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].k2;
                break;
            }
      
      
      
            if (this_ws->keys1[*k1_ptr].keys2[*k2_ptr].k2 ==\
                          this_ws->keys1[*k1_ptr].keys2[*k2_ptr+1].k2){
      
                  coeff = 1 - coeff;
      
                  if (this_ws->keys1[*k1_ptr].keys2[*k2_ptr].o_ab ==\
                            this_ws->keys1[*k1_ptr].keys2[*k2_ptr+1].o_ab){
      
                        if (this_ws->keys1[*k1_ptr].keys2[*k2_ptr].flag_next==\
                                                this_ws->keys1[*k1_ptr].keys2[*k2_ptr+1].flag_next){
      
                            this_ws->keys1[*k1_ptr].keys2[*k2_ptr].flag_next = 0;
                            this_ws->keys1[*k1_ptr].keys2[*k2_ptr+1].flag_next = 0;
      
                        }
                        
      
                  }
                  
              
            }
            else{
                  
                 if (coeff){
                    
                      this_ws->pivot.key1 = this_ws->keys1[*k1_ptr].k1;
                      this_ws->pivot.key2 = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].k2;
                      break;
      
                 }
                 else{
      
                      coeff = 1;
                 }
            }
      
      
      
            if (this_ws->keys1[*k1_ptr].keys2[*k2_ptr].flag_next){
      
                ttemp.o_ab = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].o_ab;
                ttemp.a_ptr = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].a_ptr;
                ttemp.b_ptr = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].b_ptr;
      
                ttemp.low.key1 = this_ws->keys1[*k1_ptr].k1;
                ttemp.low.key2 = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].k2;
      
                //if (this_ws->edge == self->g_debug_edge){
                //  printf("\nFinding next of %d:(%d, %d)", ttemp.o_ab\
                //                                        , ttemp.low.key1\
                //                                        , ttemp.low.key2\
                //                                        );
                //}

                find_H1_cohom_next(self, &(ttemp));

                this_ws->keys1[*k1_ptr].keys2[*k2_ptr].flag_next = 0;
      
                insert_in_implicit_v(self, ws_counter, &(ttemp), 1);
      
      
                // It is possible that last key1 and last key2 changed. Make sure last is consistent
                
                
            }
      
      
            *k2_ptr = *k2_ptr + 1;
      
            if (*k2_ptr == this_ws->keys1[*k1_ptr].last-1){

                if (coeff){
                          this_ws->pivot.key1 = this_ws->keys1[*k1_ptr].k1;
                          this_ws->pivot.key2 = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].k2;
                          break;
                }
      
                if (this_ws->keys1[*k1_ptr].keys2[*k2_ptr].flag_next){
                      
                      ttemp.o_ab  = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].o_ab;
                      ttemp.a_ptr = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].a_ptr;
                      ttemp.b_ptr = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].b_ptr;
      
                      ttemp.low.key1 = this_ws->keys1[*k1_ptr].k1;
                      ttemp.low.key2 = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].k2;
      
                      //if (this_ws->edge == self->g_debug_edge){
                      //  printf("\nFinding next of %d:(%d, %d)", ttemp.o_ab\
                      //                                        , ttemp.low.key1\
                      //                                        , ttemp.low.key2\
                      //                                        );
                      //}

                      find_H1_cohom_next(self, &(ttemp));
      
                      this_ws->keys1[*k1_ptr].keys2[*k2_ptr].flag_next = 0;

                      insert_in_implicit_v(self, ws_counter, &(ttemp), 1);
      
                      
                }
      
      
                if (*k2_ptr == this_ws->keys1[*k1_ptr].last-2){
      
                          *k2_ptr = *k2_ptr + 1;
                          this_ws->pivot.key1 = this_ws->keys1[*k1_ptr].k1;
                          this_ws->pivot.key2 = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].k2;
                          break;
                }
                else{
      
                          // Mark this key1 as empty
                          this_ws->keys1[*k1_ptr].flag_empty = 1;

                          // Reallocate to prune space
                          if (this_ws->keys1[*k1_ptr].max_len > 5){
                               this_ws->keys1[*k1_ptr].max_len = 5; 
                               self->g_V_ws_H1[ws_counter].keys1[*k1_ptr].keys2 = \
                                                        (implicit_keys2*)realloc\
                                                             (self->g_V_ws_H1[ws_counter].keys1[*k1_ptr].keys2\
                                                            , self->g_V_ws_H1[ws_counter].keys1[*k1_ptr].max_len\
                                                                  *sizeof(implicit_keys2));
                          }
      
                          EDGE_ID current_ptr = 0;
                          EDGE_ID minn = self->g_n_valid_edges;
      
                          for (EDGE_ID mm = 0; mm < this_ws->last; mm++){
      
                              if (this_ws->keys1[mm].flag_empty){
                                  continue;
                              }
      
                              implicit_keys1 ttemp = this_ws->keys1[current_ptr];
      
                              this_ws->keys1[current_ptr] = this_ws->keys1[mm];
                              
                              this_ws->keys1[mm] = ttemp;
      
                              if (this_ws->keys1[current_ptr].k1 < minn){
                                   minn = this_ws->keys1[current_ptr].k1;
                                   *k1_ptr = current_ptr;
                              }
      
                              current_ptr++;
                                
                                
                          }
      
                          this_ws->last = current_ptr;
      
                          if (minn == self->g_n_valid_edges){
                              this_ws->flag_non_empty = 0;
                              break;
                          }
      
                          coeff = 1;
                          *k2_ptr = 0;
      
                          sorter7_tim_sort(this_ws->keys1[*k1_ptr].keys2\
                                            , this_ws->keys1[*k1_ptr].last);          
      
                      
                }
      
                  
            }
      
      
      }
      
}




void reduce_ws_coH1(filtration* self){



      //printf("\nBefore parallel");
      //for (EDGE_ID mm = 0; mm < self->g_cohom_ws_size; mm++){
      //    printf("\n%d:(%d, %d) flag %d", self->g_V_ws_H1[mm].edge\
      //                                  , self->g_V_ws_H1[mm].pivot.key1\
      //                                  , self->g_V_ws_H1[mm].pivot.key2\
      //                                  , self->g_V_ws_H1[mm].flag_append_to_complex);
      //    printf(" v: ");
      //    for (EDGE_ID bb = 0; bb < self->g_V_ws_H1[mm].v_edges.last; bb++){
      //          printf("%d, ", self->g_V_ws_H1[mm].v_edges.o_ab[bb]);
      //    }
      //}

      //

      // PARALLEL
      self->g_processed_threads = 0;

      pthread_cond_broadcast(&(self->g_start_workers));

      while (self->g_processed_threads != self->g_cpu_count){
            
            pthread_cond_wait(&(self->g_start_boss) \
                            ,&(self->g_thread_lock));
      }


      //printf("\nAfter parallel");
      //for (EDGE_ID mm = 0; mm < self->g_cohom_ws_size; mm++){
      //    printf("\n%d:(%d, %d) flag %d", self->g_V_ws_H1[mm].edge\
      //                                  , self->g_V_ws_H1[mm].pivot.key1\
      //                                  , self->g_V_ws_H1[mm].pivot.key2\
      //                                  , self->g_V_ws_H1[mm].flag_append_to_complex);
      //    printf(" v: ");
      //    for (EDGE_ID bb = 0; bb < self->g_V_ws_H1[mm].v_edges.last; bb++){
      //          printf("%d, ", self->g_V_ws_H1[mm].v_edges.o_ab[bb]);
      //    }
      //}

      // SERIAL
      reduce_with_self_coH1(self);

      //printf("\nAfter serial");
      //for (EDGE_ID mm = 0; mm < self->g_cohom_ws_size; mm++){
      //    printf("\n%d:(%d, %d) flag %d", self->g_V_ws_H1[mm].edge\
      //                                  , self->g_V_ws_H1[mm].pivot.key1\
      //                                  , self->g_V_ws_H1[mm].pivot.key2\
      //                                  , self->g_V_ws_H1[mm].flag_append_to_complex);
      //    printf(" v: ");
      //    for (EDGE_ID bb = 0; bb < self->g_V_ws_H1[mm].v_edges.last; bb++){
      //          printf("%d, ", self->g_V_ws_H1[mm].v_edges.o_ab[bb]);
      //    }
      //}

      //getchar();


      // CLEARANCE
      int count_valid = 0;

      for (int ws_counter = 0; ws_counter < self->g_ws_counter; ws_counter++){

            
            if (!self->g_V_ws_H1[ws_counter].flag_non_empty){
                // Add the undead H1
                
                if (self->g_H1_pers_pairs_len+2 == self->g_H1_pers_pairs_max_len){
                      self->g_H1_pers_pairs_max_len += 1000;
                      self->g_H1_pers_pairs = (PAR*)realloc(self->g_H1_pers_pairs\
                                              , self->g_H1_pers_pairs_max_len*sizeof(PAR));
                
                }

                self->g_H1_pers_pairs[self->g_H1_pers_pairs_len++] = \
                                                     self->g_edge_parameter[self->g_V_ws_H1[ws_counter].edge];

                self->g_H1_pers_pairs[self->g_H1_pers_pairs_len++] = -1;

//#ifdef HOM_CYCLES
                if (self->g_compute_cycles){
                    self->g_H1_undead[self->g_H1_undead_ptr++] = self->g_V_ws_H1[ws_counter].edge;
                    if (self->g_H1_undead_ptr == self->g_H1_undead_max){
                        self->g_H1_undead_max += 100;
                        self->g_H1_undead = (EDGE_ID*)realloc(self->g_H1_undead\
                                                              , self->g_H1_undead_max*sizeof(EDGE_ID));
                    }
                }
//#endif
              
                continue;
            }

            if (self->g_V_ws_H1[ws_counter].flag_append_to_complex){
                  update_V_coH1(self, ws_counter);
                  continue;
            }


            // Swap V
            coboundary_H1_ws temp = self->g_V_ws_H1[count_valid];
            self->g_V_ws_H1[count_valid] = self->g_V_ws_H1[ws_counter];
            self->g_V_ws_H1[ws_counter] = temp;

            // At this point, this has to be a non-zero column
            self->g_V_ws_H1[count_valid].flag_non_empty = 1;

            // Run through parallel at least once
            self->g_V_ws_H1[count_valid].flag_append_to_complex = 0;
            
            count_valid++;

      }

      self->g_ws_counter = count_valid;
      
      
}


void reduce_with_self_coH2(filtration* self){


      // Now we have to reduce

      for (int ws_counter = 0; ws_counter < self->g_ws_counter; ws_counter++){

            coboundary_H2_ws* this_ws = self->g_V_ws_H2 + ws_counter;

            // If empty, then continue and don't append to complex
            if (!this_ws->flag_non_empty){
                  //this_ws->flag_append_to_complex = 0;
                  continue;
            }

            int m = 0;

            // Keep reducing if reduce with complex flag is 0 and reduce with trivial flag is 0
            while((m < ws_counter)\
                && (!this_ws->flag_red_w_complex)\
                && (!this_ws->flag_red_w_trivial)){
                
                coboundary_H2_ws* m_ws = self->g_V_ws_H2 + m;

                // If m is empty, continue
                if (!m_ws->flag_non_empty){
                      m++;
                      continue;
                }

                int compare;

                if (m_ws->pivot.key1 > this_ws->pivot.key1) compare = 1;
                else if (m_ws->pivot.key1 < this_ws->pivot.key1) compare = 0;
                else{
                  if (m_ws->pivot.key2 > this_ws->pivot.key2) compare = 1;
                  else if (m_ws->pivot.key2 < this_ws->pivot.key2) compare = 0;
                  else compare = -1;
                }

                // If pivot of m is higher than pivot of ws_counter
                // then we don't care
                if (compare == 1){
                      m++;
                      continue;
                }

                // If pivot of m is lower than pivot of ws_counter
                // then if m has to be reduced, we have to hold ws_counter
                if (compare == 0){

                      if (m_ws->flag_red_w_complex || m_ws->flag_red_w_trivial){
                            
                            this_ws->flag_append_to_complex = 0;
                            break;
                      }
                      m++;
                      continue;
                }

                // At this point they have same low
                if (m_ws->flag_red_w_complex || m_ws->flag_red_w_trivial){
                      this_ws->flag_append_to_complex = 0;
                      //m++;
                      break;
                      //continue;
                }

                //printf("\nin serial reducing %d with %d", this_ws->edge, m_ws->edge);
                //getchar();

                // Merge m and this_ws
                //

                
                // Merge v_triangles
                if (this_ws->v_triangles.last + m_ws->v_triangles.last > this_ws->v_triangles.max_len - 1){
                      
                      this_ws->v_triangles.max_len += m_ws->v_triangles.last + 100;

                      self->g_V_ws_H2[ws_counter].v_triangles.o_abc = (simplex*)realloc(\
                                                 self->g_V_ws_H2[ws_counter].v_triangles.o_abc\
                                                 , this_ws->v_triangles.max_len*sizeof(simplex));

                }

                // Add the original edge
                

                this_ws->v_triangles.o_abc[this_ws->v_triangles.last++] = m_ws->triangle;
                
                for (EDGE_ID bb = 0; bb < m_ws->v_triangles.last; bb++){

                      this_ws->v_triangles.o_abc[this_ws->v_triangles.last++] =\
                                                    m_ws->v_triangles.o_abc[bb];
                      
                }

                // Merge hash tables
                coboundary_H2 ttemp;

                for (EDGE_ID bb = 0; bb < m_ws->last; bb++){

                      EDGE_ID m_start = 0;

                      if (bb == m_ws->k1_ptr){
                          m_start = m_ws->k2_ptr;
                      }
                      
                      ttemp.low.key1 = m_ws->keys1[bb].k1;

                      for (EDGE_ID mm = m_start; mm < m_ws->keys1[bb].last; mm++){
                            

                            ttemp.low.key2 = m_ws->keys1[bb].keys2[mm].k2;
                            ttemp.triangle = m_ws->keys1[bb].keys2[mm].o_abc;
                            ttemp.a_ptr = m_ws->keys1[bb].keys2[mm].a_ptr;
                            ttemp.b_ptr = m_ws->keys1[bb].keys2[mm].b_ptr;
                            ttemp.c_ptr = m_ws->keys1[bb].keys2[mm].c_ptr;
                            ttemp.vertex = m_ws->keys1[bb].keys2[mm].vertex;

                            coH2_insert_in_implicit_v(self, ws_counter, &(ttemp)\
                                                  , m_ws->keys1[bb].keys2[mm].flag_next);
                            
                            
                      }
                      
                }



                // Now reduce

                reduce_hash_table_coH2(self, ws_counter);

                if (!this_ws->flag_non_empty){
                  break;
                }

                coboundary_H2 temptemp;

                // CHECK FOR TRIVIAL PAIR
                // Get low for maximum triangle <ab, d> in this_pivot
                temptemp.triangle.key1 = this_ws->pivot.key1;
                temptemp.triangle.key2 = self->g_edges_list[2*this_ws->pivot.key2+1];

                find_H2_cohom_low(self, &temptemp);


                // Check if the low of this triangle is same as self->g_this_pivot
                if ((temptemp.low.key1 == this_ws->pivot.key1)\
                    && (temptemp.low.key2 == this_ws->pivot.key2)){

                        this_ws->flag_red_w_trivial = 1;
                        this_ws->flag_red_w_complex = 0;
                        this_ws->flag_append_to_complex = 0;

                        this_ws->reduce_w_bndry = temptemp.triangle;
                        this_ws->V_col_idx = 0;
                        break;

                }


                 // If this low is not a pivot
                 if (self->g_H2_cohom_pivots_len[this_ws->pivot.key1]){

                    EDGE_ID idx = search_H2_cohom_pivots(self->g_H2_cohom_pivots[this_ws->pivot.key1]\
                                        , 0 \
                                        , self->g_H2_cohom_pivots_len[this_ws->pivot.key1] - 1\
                                        , this_ws->pivot.key2 \
                                        , self->g_n_valid_edges);



                    if (idx != self->g_n_valid_edges){

                        this_ws->flag_red_w_trivial = 0;
                        this_ws->flag_red_w_complex = 1;
                        this_ws->flag_append_to_complex = 0;
                          
                        this_ws->reduce_w_bndry = self->g_H2_cohom_pivots[this_ws->pivot.key1][idx].bndry;
                        this_ws->V_col_idx = self->g_H2_cohom_pivots[this_ws->pivot.key1][idx].col_idx;

                        break;

                    }

                      
                 }


                 // Reset m after single reduction
                 m = 0;
                
            }
            
      }
      
}


void* reduce_with_complex_coH2(void* arg){
        
      filtration* self = arg;


      pthread_mutex_lock(&(self->g_thread_lock));

      int tid = ++self->g_thread_id;

      coboundary_H2_ws* this_ws;

      coboundary_H2 temp;

      EDGE_ID idx, check_len, start, end;


      pthread_mutex_unlock(&(self->g_thread_lock));

      for (;;){

          pthread_mutex_lock(&(self->g_thread_lock));

          self->g_sleeping_threads++;
          self->g_processed_threads++;

          if (self->g_sleeping_threads == self->g_cpu_count){

              pthread_cond_signal(&(self->g_start_boss));
          }


          pthread_cond_wait(&(self->g_start_workers), &(self->g_thread_lock));

          if (self->g_delete_threads){
            pthread_mutex_unlock(&(self->g_thread_lock));
            pthread_exit(NULL);
          }

          self->g_sleeping_threads--;

          pthread_mutex_unlock(&(self->g_thread_lock));


          for (int ws_counter = self->g_jobs[tid - 1]; ws_counter < self->g_jobs[tid]; ws_counter++){


              this_ws = self->g_V_ws_H2 + ws_counter;

              //coboundary_H2 temp;

              //printf("\nProcessing (%d, %d)", this_ws->triangle.key1\
                                            , this_ws->triangle.key2);


              if (!this_ws->flag_non_empty){
                  // We are sure that we will exit only if there is no reduction
                  // required with existing complex or with trivial pair
                  //printf("\nEmpty. Skipping.");
                  this_ws->flag_red_w_complex = 0;
                  this_ws->flag_red_w_trivial = 0;
                  this_ws->flag_append_to_complex = 0;
                  continue;
              }

              if (this_ws->flag_append_to_complex){
                  //printf("\nAppend to complex. Nothing to do.");
                  continue;
              }

              if (this_ws->flag_first){

                  this_ws->flag_first = 0;

                  // CHECK WITH TRIVIAL
                  temp.triangle.key1 = this_ws->pivot.key1;
                  temp.triangle.key2 = self->g_edges_list[2*this_ws->pivot.key2+1];

                  find_H2_cohom_low(self, &temp);

                  if ((temp.low.key1 == this_ws->pivot.key1)\
                      && (temp.low.key2 == this_ws->pivot.key2)){

                        
                      this_ws->flag_red_w_trivial = 1;

                      this_ws->reduce_w_bndry = temp.triangle;
                      this_ws->V_col_idx = 0;


                  }
                  else{


                      if (!self->g_H2_cohom_pivots_len[this_ws->pivot.key1]){
                          
                          //printf("\npivot not in complex c1. append");
                          this_ws->flag_red_w_complex = 0;
                          this_ws->flag_red_w_trivial = 0;
                          this_ws->flag_append_to_complex = 1;
                          continue;
                      }
                      else{

                          idx = search_H2_cohom_pivots(self->g_H2_cohom_pivots[this_ws->pivot.key1]\
                                        , 0 \
                                        , self->g_H2_cohom_pivots_len[this_ws->pivot.key1] - 1\
                                        , this_ws->pivot.key2 \
                                        , self->g_n_valid_edges);
                          
                          // If this low is not a pivot
                          if (idx == self->g_n_valid_edges){

                                //printf("\npivot not in complex c2. append");
                                this_ws->flag_red_w_complex = 0;
                                this_ws->flag_red_w_trivial = 0;
                                this_ws->flag_append_to_complex = 1;
                                continue;


                          }
                          else{

                                this_ws->flag_red_w_complex = 1;
                                this_ws->reduce_w_bndry = self->g_H2_cohom_pivots[this_ws->pivot.key1][idx].bndry;
                                this_ws->V_col_idx = self->g_H2_cohom_pivots[this_ws->pivot.key1][idx].col_idx;

                          }

                      }

                  }

                    
              }



              if ((!this_ws->flag_red_w_trivial) && (!this_ws->flag_red_w_complex)){
                  //printf("\nno red with trivial and no red with complex");
                  this_ws->flag_append_to_complex = 1;
                  continue;
              }

              //printf("\nReducing...");

              // Presume that this will be flagged to be added to complex
              this_ws->flag_red_w_complex = 0;
              this_ws->flag_red_w_trivial = 0;
              this_ws->flag_append_to_complex = 1;
                    
              int flag = 0;

              while(1){

                    check_len = this_ws->v_triangles.last + 1;

                    if (this_ws->V_col_idx){

                          check_len += self->g_V_col_indices[this_ws->V_col_idx+1] -\
                                                        self->g_V_col_indices[this_ws->V_col_idx];

                    }

                    if (check_len > this_ws->v_triangles.max_len){
                        this_ws->v_triangles.max_len = check_len + 100;
                        self->g_V_ws_H2[ws_counter].v_triangles.o_abc = \
                                        (simplex*)realloc(self->g_V_ws_H2[ws_counter].v_triangles.o_abc\
                                                                , this_ws->v_triangles.max_len*sizeof(simplex));

                    }

                    this_ws->v_triangles.o_abc[this_ws->v_triangles.last++] = this_ws->reduce_w_bndry;

                    temp.triangle = this_ws->reduce_w_bndry;


                    find_H2_cohom_greater(self, &(temp), &(this_ws->pivot));

                    //if (flag){
                    //    printf("\ninserting (%d, %d):(%d, %d)"\
                    //                                , temp.triangle.key1\
                    //                                , temp.triangle.key2\
                    //                                , temp.low.key1\
                    //                                , temp.low.key2\
                    //                                );
                    //}

                    //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
                    //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

                    //    printf("\nBefore inserting (%d, %d):(%d, %d)"\
                    //                                , temp.triangle.key1\
                    //                                , temp.triangle.key2\
                    //                                , temp.low.key1\
                    //                                , temp.low.key2\
                    //                                );

                    //    coH2_print_v_implicit(self, ws_counter);

                    //}
                     
                    coH2_insert_in_implicit_v(self, ws_counter, &(temp), 1);

                    //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
                    //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

                    //    printf("\nAfter inserting");

                    //    coH2_print_v_implicit(self, ws_counter);

                    //}

                    // IF the V was recorded, add the bndries
                    if (this_ws->V_col_idx){

                        // We have to cycle through the col in V and add all the other boundary columns for reduction

                        start = self->g_V_col_indices[this_ws->V_col_idx];
                        end = self->g_V_col_indices[this_ws->V_col_idx+1];


                            for (EDGE_ID mm = start; mm < end; mm++){
                                  
                                this_ws->v_triangles.o_abc[this_ws->v_triangles.last++] = self->g_V_sparse_H2[mm];

                                temp.triangle = self->g_V_sparse_H2[mm];

                                // Find the first low greater than or equal pivot
                                find_H2_cohom_greater(self, &(temp), &(this_ws->pivot));

                                //if (flag){
                                //    printf("\ninserting (%d, %d):(%d, %d)"\
                                //                                , temp.triangle.key1\
                                //                                , temp.triangle.key2\
                                //                                , temp.low.key1\
                                //                                , temp.low.key2\
                                //                                );
                                //}

                                //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
                                //              (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

                                //    printf("\nBefore inserting (%d, %d):(%d, %d)"\
                                //                                , temp.triangle.key1\
                                //                                , temp.triangle.key2\
                                //                                , temp.low.key1\
                                //                                , temp.low.key2\
                                //                                );

                                //    coH2_print_v_implicit(self, ws_counter);

                                //}


                                coH2_insert_in_implicit_v(self, ws_counter,  &(temp), 1);

                                //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
                                //              (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

                                //    printf("\nAfter inserting");

                                //    coH2_print_v_implicit(self, ws_counter);

                                //}

                            }

                    }


                    //simplex test_low = this_ws->pivot;


                    //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
                    //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

                    //    printf("\nBefore reduction");
                    //    coH2_print_v_implicit(self, ws_counter);
                    //    printf("\nPivot is (%d, %d)"\
                    //                                      , this_ws->pivot.key1\
                    //                                      , this_ws->pivot.key2\
                    //                                      );
                    //                         
                    //}


                    reduce_hash_table_coH2(self, ws_counter);

                    //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
                    //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

                    //    printf("\nAfter reduction");

                    //    coH2_print_v_implicit(self, ws_counter);

                    //    printf("\nPivot is (%d, %d)"\
                    //                                      , this_ws->pivot.key1\
                    //                                      , this_ws->pivot.key2\
                    //                                      );
                    //}


                    //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
                    //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){
                    //    printf("\nNew pivot (%d, %d)", this_ws->pivot.key1\
                    //                                 , this_ws->pivot.key2);
                    //    getchar();
                    //}
                    

                    //if ((test_low.key1 == this_ws->pivot.key1)
                    //    &&(test_low.key2 == this_ws->pivot.key2)){

                    //    printf("\nprinting..");

                    //    flag = 1;
                    //    coH2_print_v_implicit(self, ws_counter);

                    //    getchar();

                    //}

                    //printf("\nNew low (%d, %d,)", this_ws->pivot.key1\
                    //                            , this_ws->pivot.key2);

                    if (!this_ws->flag_non_empty){
                        break;
                    }


                    // Check with trivial pair
                    temp.triangle.key1 = this_ws->pivot.key1;
                    temp.triangle.key2 = self->g_edges_list[2*this_ws->pivot.key2+1];

                    find_H2_cohom_low(self, &temp);

                    if ((temp.low.key1 == this_ws->pivot.key1)\
                        && (temp.low.key2 == this_ws->pivot.key2)){


                            //if (flag){
                            //    printf("\nReducing with trivial (%d, %d)"\
                            //                                , temp.triangle.key1\
                            //                                , temp.triangle.key2\
                            //                                );

                            //}

                            //printf("\nreduce with trivial");

                            this_ws->reduce_w_bndry = temp.triangle;
                            this_ws->V_col_idx = 0;
                            continue;
                          
                    }

                    // If this low is not a pivot
                    if (!self->g_H2_cohom_pivots_len[this_ws->pivot.key1]){
                        break;
                    }


                    idx = search_H2_cohom_pivots(self->g_H2_cohom_pivots[this_ws->pivot.key1]\
                                        , 0 \
                                        , self->g_H2_cohom_pivots_len[this_ws->pivot.key1] - 1\
                                        , this_ws->pivot.key2 \
                                        , self->g_n_valid_edges\
                                        );


                    if (idx == self->g_n_valid_edges){
                      break;
                    }

                    //if (flag){

                    //    for (EDGE_ID mm = 0; mm < self->g_H2_cohom_pivots_len[this_ws->pivot.key1]; mm++){
                    //          printf("\n(%d, %d, %d), "\
                    //                              , self->g_H2_cohom_pivots[this_ws->pivot.key1][mm].key2\
                    //                              , self->g_H2_cohom_pivots[this_ws->pivot.key1][mm].bndry.key1\
                    //                              , self->g_H2_cohom_pivots[this_ws->pivot.key1][mm].bndry.key2\
                    //                              );
                    //    }
                    //    printf("\nReducing with complex (%d, %d) at idx %d"\
                    //                                , self->g_H2_cohom_pivots[this_ws->pivot.key1][idx].bndry.key1\
                    //                                , self->g_H2_cohom_pivots[this_ws->pivot.key1][idx].bndry.key2\
                    //                                , idx\
                    //                                );

                    //}

                    //printf("\nreduce with complex");

                    this_ws->reduce_w_bndry = self->g_H2_cohom_pivots[this_ws->pivot.key1][idx].bndry;
                    this_ws->V_col_idx = self->g_H2_cohom_pivots[this_ws->pivot.key1][idx].col_idx;

              }

          }

      }
        
}


void update_V_coH2(filtration* self, int ws_counter){


    EDGE_ID red_col = 0;

    coboundary_H2_ws* this_ws = self->g_V_ws_H2 + ws_counter;

    //if ((this_ws->triangle.key1 == 227282)\
    //    &&(this_ws->triangle.key2 == 1807632)){
    //    

    //    self->g_new_debug = 1;
    //    

    //}
    //else{
    //    self->g_new_debug = 0;
    //}


    //if (self->g_new_debug){
    // 
    //    printf("\n ADDDDDING %d, %d, %d, %d", this_ws->triangle.key1\
    //                     , this_ws->triangle.key2\
    //                     , this_ws->pivot.key1\
    //                     , this_ws->pivot.key2\
    //                     );
    //    getchar();

    //}



    //printf("\nENTERING UPDATE_V_coH2");


    self->g_V_sparse_beg_ptr = self->g_V_sparse_ptr;


    if (this_ws->v_triangles.last){

        //printf("\n%d, %d, %d, %d", this_ws->triangle.key1\
        //                 , this_ws->triangle.key2\
        //                 , this_ws->pivot.key1\
        //                 , this_ws->pivot.key2\
        //                 );


        if ((this_ws->v_triangles.last + self->g_V_sparse_ptr) + 1 > self->g_V_sparse_max){

                    self->g_V_sparse_max = self->g_V_sparse_ptr + this_ws->v_triangles.last + 100000;
                    self->g_V_sparse_H2 = (simplex*)realloc(self->g_V_sparse_H2\
                                                        , self->g_V_sparse_max*sizeof(simplex));

        }


        // //TRYING EDIT
        if (this_ws->v_triangles.last > 1){

#ifdef VREDUCE2

            sorter4_tim_sort(this_ws->v_triangles.o_abc, this_ws->v_triangles.last);          
          
            int coeff = 1;

            for (EDGE_ID vv = 0; vv < this_ws->v_triangles.last-1; vv++){
                  
                if ((this_ws->v_triangles.o_abc[vv].key1 == this_ws->v_triangles.o_abc[vv+1].key1) &&
                    (this_ws->v_triangles.o_abc[vv].key2 == this_ws->v_triangles.o_abc[vv+1].key2))
                {
                    coeff = 1 - coeff;
                }
                else{
                    if (coeff){
                        self->g_V_sparse_H2[self->g_V_sparse_ptr++] = this_ws->v_triangles.o_abc[vv];
                    }

                    coeff = 1;
                }
                  
                  
            }

            if (coeff){
                 self->g_V_sparse_H2[self->g_V_sparse_ptr++] = this_ws->v_triangles.o_abc[this_ws->v_triangles.last-1];
            }

#else

            for (EDGE_ID vv = 0; vv < this_ws->v_triangles.last; vv++){

                 self->g_V_sparse_H2[self->g_V_sparse_ptr++] = this_ws->v_triangles.o_abc[vv];

            }

#endif


        }
        else if (this_ws->v_triangles.last == 1){

            self->g_V_sparse_H2[self->g_V_sparse_ptr++] = this_ws->v_triangles.o_abc[0];
            
        }



        // All have been added
        this_ws->v_triangles.last = 0;


        if ((self->g_V_sparse_ptr - self->g_V_sparse_beg_ptr) > 0){

              red_col = self->g_V_col_indices_ptr;

              if (self->g_V_col_indices_ptr+1 == self->g_V_col_indices_max){
                    
                    self->g_V_col_indices_max += 10000;
                    self->g_V_col_indices = (EDGE_ID*)realloc(self->g_V_col_indices\
                                                          , self->g_V_col_indices_max*sizeof(EDGE_ID));
              
              }

              self->g_V_col_indices[self->g_V_col_indices_ptr] = self->g_V_sparse_beg_ptr;
              self->g_V_col_indices[self->g_V_col_indices_ptr+1] = self->g_V_sparse_ptr;

              self->g_V_col_indices_ptr++;

        }

    }

    //printf("\n(%d, %d):(%d, %d)"\
    //                   , this_ws->triangle.key1\
    //                   , this_ws->triangle.key2\
    //                   , this_ws->pivot.key1\
    //                   , this_ws->pivot.key2\
    //                   );

    //printf("\nAdding v: ");
    //
    //for (EDGE_ID mm = self->g_V_sparse_beg_ptr; mm < self->g_V_sparse_ptr; mm++){
    //      printf("(%d, %d),", self->g_V_sparse_H2[mm].key1\
    //                        , self->g_V_sparse_H2[mm].key2);
    //}

    //getchar();

    // ADDING THE LOW
    //

    add_coH2_pivot(self, this_ws->triangle, this_ws->pivot, red_col);


}



void add_coH2_pivot (filtration* self, simplex triangle, simplex pivot, EDGE_ID red_col){
      

    if (!self->g_H2_cohom_pivots_max_len[pivot.key1]){
        
       self->g_H2_cohom_pivots_max_len[pivot.key1] = 2;
       self->g_H2_cohom_pivots[pivot.key1] = \
                              (H2_cohom_pivots*)malloc(self->g_H2_cohom_pivots_max_len[pivot.key1]*sizeof(H2_cohom_pivots));
    }
      
    if (self->g_H2_cohom_pivots_len[pivot.key1]\
                            == self->g_H2_cohom_pivots_max_len[pivot.key1]){
          
          self->g_H2_cohom_pivots_max_len[pivot.key1] += 5;
          self->g_H2_cohom_pivots[pivot.key1] = (H2_cohom_pivots*)realloc( \
                          self->g_H2_cohom_pivots[pivot.key1] \
                          , self->g_H2_cohom_pivots_max_len[pivot.key1]*sizeof(H2_cohom_pivots));

    }




    EDGE_ID old_ptr = self->g_H2_cohom_pivots_len[pivot.key1];
    EDGE_ID new_ptr = self->g_H2_cohom_pivots_len[pivot.key1];

    while (old_ptr){
          
          old_ptr--;
          
          if (self->g_H2_cohom_pivots[pivot.key1][old_ptr].key2 > pivot.key2){

                self->g_H2_cohom_pivots[pivot.key1][new_ptr--] =\
                                                       self->g_H2_cohom_pivots[pivot.key1][old_ptr];
                continue;

          }
          break;

    }


    self->g_H2_cohom_pivots[pivot.key1][new_ptr].key2 = pivot.key2;
    self->g_H2_cohom_pivots[pivot.key1][new_ptr].col_idx = red_col;
                            
    self->g_H2_cohom_pivots[pivot.key1][new_ptr].bndry = triangle;

    self->g_H2_cohom_pivots_len[pivot.key1]++;

    // PERS PAIRS
    // Add non-zero barcodes
        
    PAR birth = self->g_edge_parameter[triangle.key1];
    PAR death = self->g_edge_parameter[pivot.key1];
    if (birth != death){

           //printf("\nNon trivial pers pair (%f, %f)", birth, death);
           

           if (birth > death){
               printf("\nBirth, death (%lf, %lf)", birth, death);
               printf("\nError (%d, %d) at pair (%d, %d)", triangle.key1\
                                                    , triangle.key2\
                                                    , pivot.key1\
                                                    , pivot.key2);
               getchar();
             
           }


          if (self->g_H2_pers_pairs_len+2 == self->g_H2_pers_pairs_max_len){
                self->g_H2_pers_pairs_max_len += 1000;
                self->g_H2_pers_pairs = (PAR*)realloc(self->g_H2_pers_pairs\
                                              , self->g_H2_pers_pairs_max_len*sizeof(PAR));
          
          }
          self->g_H2_pers_pairs[self->g_H2_pers_pairs_len++] = birth;
          self->g_H2_pers_pairs[self->g_H2_pers_pairs_len++] = death;
          
    }
      
      
      
      
}


void reduce_ws_coH2(filtration* self){



      //printf("\nBefore parallel");
      //for (EDGE_ID mm = 0; mm < self->g_cohom_ws_size; mm++){
      //    printf("\n%d:(%d, %d) flag %d", self->g_V_ws_H1[mm].edge\
      //                                  , self->g_V_ws_H1[mm].pivot.key1\
      //                                  , self->g_V_ws_H1[mm].pivot.key2\
      //                                  , self->g_V_ws_H1[mm].flag_append_to_complex);
      //    printf(" v: ");
      //    for (EDGE_ID bb = 0; bb < self->g_V_ws_H1[mm].v_edges.last; bb++){
      //          printf("%d, ", self->g_V_ws_H1[mm].v_edges.o_ab[bb]);
      //    }
      //}

      //

      // PARALLEL
      self->g_processed_threads = 0;

      pthread_cond_broadcast(&(self->g_start_workers));

      while (self->g_processed_threads != self->g_cpu_count){
            
            pthread_cond_wait(&(self->g_start_boss) \
                            ,&(self->g_thread_lock));
      }


      //printf("\nAfter parallel");
      //for (EDGE_ID mm = 0; mm < self->g_cohom_ws_size; mm++){
      //    printf("\n%d:(%d, %d) flag %d", self->g_V_ws_H1[mm].edge\
      //                                  , self->g_V_ws_H1[mm].pivot.key1\
      //                                  , self->g_V_ws_H1[mm].pivot.key2\
      //                                  , self->g_V_ws_H1[mm].flag_append_to_complex);
      //    printf(" v: ");
      //    for (EDGE_ID bb = 0; bb < self->g_V_ws_H1[mm].v_edges.last; bb++){
      //          printf("%d, ", self->g_V_ws_H1[mm].v_edges.o_ab[bb]);
      //    }
      //}

      // SERIAL
      reduce_with_self_coH2(self);

      //printf("\nAfter serial");
      //for (EDGE_ID mm = 0; mm < self->g_cohom_ws_size; mm++){
      //    printf("\n%d:(%d, %d) flag %d", self->g_V_ws_H1[mm].edge\
      //                                  , self->g_V_ws_H1[mm].pivot.key1\
      //                                  , self->g_V_ws_H1[mm].pivot.key2\
      //                                  , self->g_V_ws_H1[mm].flag_append_to_complex);
      //    printf(" v: ");
      //    for (EDGE_ID bb = 0; bb < self->g_V_ws_H1[mm].v_edges.last; bb++){
      //          printf("%d, ", self->g_V_ws_H1[mm].v_edges.o_ab[bb]);
      //    }
      //}

      //getchar();


      // CLEARANCE
      int count_valid = 0;

      for (int ws_counter = 0; ws_counter < self->g_ws_counter; ws_counter++){

            
            if (!self->g_V_ws_H2[ws_counter].flag_non_empty){
                // Add the undead H2

                //printf("\nAdding undead for H2");
                
                if (self->g_H2_pers_pairs_len+2 == self->g_H2_pers_pairs_max_len){
                      self->g_H2_pers_pairs_max_len += 1000;
                      self->g_H2_pers_pairs = (PAR*)realloc(self->g_H2_pers_pairs\
                                              , self->g_H2_pers_pairs_max_len*sizeof(PAR));
                
                }

                self->g_H2_pers_pairs[self->g_H2_pers_pairs_len++] = \
                                                     self->g_edge_parameter[self->g_V_ws_H2[ws_counter].triangle.key1];

                self->g_H2_pers_pairs[self->g_H2_pers_pairs_len++] = -1;

//#ifdef HOM_CYCLES
                if (self->g_compute_cycles){
                    self->g_H2_undead[self->g_H2_undead_ptr++] = self->g_V_ws_H2[ws_counter].triangle;
                    if (self->g_H2_undead_ptr == self->g_H2_undead_max){
                        self->g_H2_undead_max += 100;
                        self->g_H2_undead = (simplex*)realloc(self->g_H2_undead\
                                                              , self->g_H2_undead_max*sizeof(simplex));
                    }
                }
//#endif
                //printf("\nEmpty. Not adding to complex.");
              
                continue;
            }

            if (self->g_V_ws_H2[ws_counter].flag_append_to_complex){
                  //printf("\nAdding to complex.");
                  update_V_coH2(self, ws_counter);
                  continue;
            }

            //printf("\nSwapping...");

            // Swap V
            coboundary_H2_ws temp = self->g_V_ws_H2[count_valid];
            self->g_V_ws_H2[count_valid] = self->g_V_ws_H2[ws_counter];
            self->g_V_ws_H2[ws_counter] = temp;

            // At this point, this has to be a non-zero column
            self->g_V_ws_H2[count_valid].flag_non_empty = 1;

            // Run through parallel at least once
            self->g_V_ws_H2[count_valid].flag_append_to_complex = 0;
            
            count_valid++;

      }

      self->g_ws_counter = count_valid;
      
      
}


void reduce_hash_table_coH2(filtration* self, int ws_counter){
      
      
      // Now we have to reduce
      int coeff = 1;
      
      coboundary_H2_ws* this_ws = self->g_V_ws_H2 + ws_counter;

      coboundary_H2 ttemp;
      
      EDGE_ID* k1_ptr = &(this_ws->k1_ptr);
      EDGE_ID* k2_ptr = &(this_ws->k2_ptr);
      
      while (1){
            
            
            if (this_ws->keys1[*k1_ptr].last == 1){
                this_ws->pivot.key1 = this_ws->keys1[*k1_ptr].k1;
                this_ws->pivot.key2 = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].k2;
                break;
            }
      
      
      
            if (this_ws->keys1[*k1_ptr].keys2[*k2_ptr].k2 ==\
                          this_ws->keys1[*k1_ptr].keys2[*k2_ptr+1].k2){
      
                  coeff = 1 - coeff;
      
                  if ((this_ws->keys1[*k1_ptr].keys2[*k2_ptr].o_abc.key1 ==\
                            this_ws->keys1[*k1_ptr].keys2[*k2_ptr+1].o_abc.key1) &&
                          (this_ws->keys1[*k1_ptr].keys2[*k2_ptr].o_abc.key2 ==\
                            this_ws->keys1[*k1_ptr].keys2[*k2_ptr+1].o_abc.key2))
                  {
      
                        if (this_ws->keys1[*k1_ptr].keys2[*k2_ptr].flag_next ==\
                                                this_ws->keys1[*k1_ptr].keys2[*k2_ptr+1].flag_next){
      
                            this_ws->keys1[*k1_ptr].keys2[*k2_ptr].flag_next = 0;
                            this_ws->keys1[*k1_ptr].keys2[*k2_ptr+1].flag_next = 0;
      
                        }
                        
      
                  }
                  
              
            }
            else{
                  
                 if (coeff){
                    
                      this_ws->pivot.key1 = this_ws->keys1[*k1_ptr].k1;
                      this_ws->pivot.key2 = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].k2;
                      break;
      
                 }
                 else{
      
                      coeff = 1;
                 }
            }
      
      
            if (this_ws->keys1[*k1_ptr].keys2[*k2_ptr].flag_next){
      
                ttemp.triangle = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].o_abc;
                ttemp.a_ptr = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].a_ptr;
                ttemp.b_ptr = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].b_ptr;
                ttemp.c_ptr = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].c_ptr;
                ttemp.vertex = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].vertex;
      
                ttemp.low.key1 = this_ws->keys1[*k1_ptr].k1;
                ttemp.low.key2 = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].k2;
      
                //if (this_ws->edge == self->g_debug_edge){
                //  printf("\nFinding next of %d:(%d, %d)", ttemp.o_ab\
                //                                        , ttemp.low.key1\
                //                                        , ttemp.low.key2\
                //                                        );
                //}

                //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
                //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

                //    printf("\nFinding next of (%d, %d):(%d, %d)"\
                //                                      , ttemp.triangle.key1\
                //                                      , ttemp.triangle.key2\
                //                                      , ttemp.low.key1\
                //                                      , ttemp.low.key2\
                //                                      );

                //}

                find_H2_cohom_next(self, &(ttemp));

                //printf("\nInserting (%d, %d):(%d, %d)"\
                //                                  , ttemp.triangle.key1\
                //                                  , ttemp.triangle.key2\
                //                                  , ttemp.low.key1\
                //                                  , ttemp.low.key2\
                //                                  );
      
                this_ws->keys1[*k1_ptr].keys2[*k2_ptr].flag_next = 0;

                coH2_insert_in_implicit_v(self, ws_counter, &(ttemp), 1);

      
                
                //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
                //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

                //    printf("\nAfter inserting");
                //    coH2_print_v_implicit(self, ws_counter);

                //}
      
                // It is possible that last key1 and last key2 changed. Make sure last is consistent
                
                
            }
      
      
            *k2_ptr = *k2_ptr + 1;
      
            if (*k2_ptr == this_ws->keys1[*k1_ptr].last-1){
      
                if (coeff){
                          this_ws->pivot.key1 = this_ws->keys1[*k1_ptr].k1;
                          this_ws->pivot.key2 = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].k2;
                          break;
                }

                if (this_ws->keys1[*k1_ptr].keys2[*k2_ptr].flag_next){
                      
                      ttemp.triangle  = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].o_abc;
                      ttemp.a_ptr = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].a_ptr;
                      ttemp.b_ptr = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].b_ptr;
                      ttemp.c_ptr = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].c_ptr;
                      ttemp.vertex = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].vertex;
      
                      ttemp.low.key1 = this_ws->keys1[*k1_ptr].k1;
                      ttemp.low.key2 = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].k2;
      

                      find_H2_cohom_next(self, &(ttemp));
      
                      this_ws->keys1[*k1_ptr].keys2[*k2_ptr].flag_next = 0;

                      coH2_insert_in_implicit_v(self, ws_counter, &(ttemp), 1);
      
                      
                }
      
      
                if (*k2_ptr == this_ws->keys1[*k1_ptr].last-2){
                      
                          *k2_ptr = *k2_ptr + 1;
                          this_ws->pivot.key1 = this_ws->keys1[*k1_ptr].k1;
                          this_ws->pivot.key2 = this_ws->keys1[*k1_ptr].keys2[*k2_ptr].k2;
                          break;
                }
      
                else{
      
                          // Mark this key1 as empty
                          this_ws->keys1[*k1_ptr].flag_empty = 1;

                          // Reallocate to prune space
                          if (this_ws->keys1[*k1_ptr].max_len > 5){
                               this_ws->keys1[*k1_ptr].max_len = 5; 
                               self->g_V_ws_H2[ws_counter].keys1[*k1_ptr].keys2 = \
                                                        (coH2_implicit_keys2*)realloc\
                                                             (self->g_V_ws_H2[ws_counter].keys1[*k1_ptr].keys2\
                                                            , self->g_V_ws_H2[ws_counter].keys1[*k1_ptr].max_len\
                                                                  *sizeof(coH2_implicit_keys2));
                          }


      
                          EDGE_ID current_ptr = 0;
                          EDGE_ID minn = self->g_n_valid_edges;
      
                          for (EDGE_ID mm = 0; mm < this_ws->last; mm++){
      
                              if (this_ws->keys1[mm].flag_empty){
                                  continue;
                              }
      
                              coH2_implicit_keys1 ttemp = this_ws->keys1[current_ptr];
      
                              this_ws->keys1[current_ptr] = this_ws->keys1[mm];
                              
                              this_ws->keys1[mm] = ttemp;
      
                              if (this_ws->keys1[current_ptr].k1 < minn){
                                   minn = this_ws->keys1[current_ptr].k1;
                                   *k1_ptr = current_ptr;
                              }
      
                              current_ptr++;
                                
                                
                          }
      
                          this_ws->last = current_ptr;
      
                          if (minn == self->g_n_valid_edges){
                              this_ws->flag_non_empty = 0;
                              break;
                          }
      
                          // Otherwise reset coefficient and begin reduction
                          coeff = 1;
                          *k2_ptr = 0;
      
                          sorter9_tim_sort(this_ws->keys1[*k1_ptr].keys2\
                                            , this_ws->keys1[*k1_ptr].last);          
      
                }
      
                  
            }
      
      
      }
      
}



void coH2_insert_in_implicit_v(filtration* self, int ws_counter, coboundary_H2* phi, int flag_next){


      if (phi->low.key1 == self->g_n_valid_edges){
          return;
      }

      coboundary_H2_ws* this_ws = self->g_V_ws_H2 + ws_counter;

      //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
      //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

      //    printf("\nINSERTING (%d, %d):(%d, %d)"\
      //                                , phi->triangle.key1\
      //                                , phi->triangle.key2\
      //                                , phi->low.key1\
      //                                , phi->low.key2\
      //                                );
      //}


      if (phi->low.key1 == this_ws->keys1[this_ws->k1_ptr].k1){


            if (this_ws->keys1[this_ws->k1_ptr].last ==\
                                            this_ws->keys1[this_ws->k1_ptr].max_len){

                  self->g_V_ws_H2[ws_counter].keys1[self->g_V_ws_H2[ws_counter].k1_ptr].max_len += 10;
                  self->g_V_ws_H2[ws_counter].keys1[self->g_V_ws_H2[ws_counter].k1_ptr].keys2 = \
                                                        (coH2_implicit_keys2*)realloc\
                                                             (self->g_V_ws_H2[ws_counter].keys1[self->g_V_ws_H2[ws_counter].k1_ptr].keys2\
                                                            , self->g_V_ws_H2[ws_counter].keys1[self->g_V_ws_H2[ws_counter].k1_ptr].max_len\
                                                                  *sizeof(coH2_implicit_keys2));
                  
            }


            EDGE_ID mm = this_ws->keys1[this_ws->k1_ptr].last;

            coH2_implicit_keys2* this_key2 = &(this_ws->keys1[this_ws->k1_ptr].keys2[mm-1]);

            int compare;

            while (1){


                //int compare = coH2_compare_implicit(this_ws->keys1[this_ws->k1_ptr].keys2[mm-1], *phi);
                  

                //if (this_ws->keys1[this_ws->k1_ptr].keys2[mm-1].k2 < phi->low.key2) compare = 0;
                //else if (this_ws->keys1[this_ws->k1_ptr].keys2[mm-1].k2 > phi->low.key2) compare = 1;
                if (this_key2->k2 < phi->low.key2) compare = 0;
                else if (this_key2->k2 > phi->low.key2) compare = 1;
                else{

                      if (this_key2->o_abc.key1 < phi->triangle.key1) compare = 0;
                      else if (this_key2->o_abc.key1 > phi->triangle.key1) compare = 1;
                      else{

                          if (this_key2->o_abc.key2 < phi->triangle.key2) compare = 0;
                          else compare = 1;

                      }

                }



                if (compare){

                      this_ws->keys1[this_ws->k1_ptr].keys2[mm] =\
                                                        this_ws->keys1[this_ws->k1_ptr].keys2[mm-1];
                }
                else{
                    
                      //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
                      //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

                      //    printf("\nin c1");
                      //    coH2_print_v_implicit(self, ws_counter);

                      //}
                  
                      coH2_implicit_keys1* this_key1 = &(this_ws->keys1[this_ws->k1_ptr]);
                      coH2_implicit_keys2* this_key2 = &(this_key1->keys2[mm]);

                      this_key2->k2 = phi->low.key2; 
                      this_key2->o_abc = phi->triangle; 
                      this_key2->a_ptr = phi->a_ptr; 
                      this_key2->b_ptr = phi->b_ptr; 
                      this_key2->c_ptr = phi->c_ptr; 
                      this_key2->vertex = phi->vertex; 
                      this_key2->flag_next = flag_next; 

                      this_ws->keys1[this_ws->k1_ptr].last++;

                      //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
                      //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

                      //    printf("\nin c2");
                      //    coH2_print_v_implicit(self, ws_counter);

                      //}

                      return;

                }

                this_key2--;
                mm--;

                //// ERROR CHECKING, REMOVE LATER
                //if (!mm){
                //    printf("\nk2_ptr %d", v_implicit->k2_ptr);
                //    printf("\nADDING %d:(%d, %d) to ", phi->o_ab, phi->low.key1, phi->low.key2);
                //    print_v_implicit(self);
                //    exit(0);
                //    
                //}

                if (mm == this_ws->k2_ptr){

                      //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
                      //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

                      //    printf("\nin c3, idx %d, last is %d"\
                      //                            , self->g_V_ws_H2[ws_counter].k1_ptr\
                      //                            , self->g_V_ws_H2[ws_counter].keys1[self->g_V_ws_H2[ws_counter].k1_ptr].last);
                      //    coH2_print_v_implicit(self, ws_counter);

                      //}

                      coH2_implicit_keys1* this_key1 = &(this_ws->keys1[this_ws->k1_ptr]);
                      coH2_implicit_keys2* this_key2 = &(this_key1->keys2[mm]);

                      this_key2->k2 = phi->low.key2; 
                      this_key2->o_abc = phi->triangle; 
                      this_key2->a_ptr = phi->a_ptr; 
                      this_key2->b_ptr = phi->b_ptr; 
                      this_key2->c_ptr = phi->c_ptr; 
                      this_key2->vertex = phi->vertex; 
                      this_key2->flag_next = flag_next; 

                      this_ws->keys1[this_ws->k1_ptr].last++;

                      //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
                      //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

                      //    printf("\nin c4, idx %d, last is %d"\
                      //                            , self->g_V_ws_H2[ws_counter].k1_ptr\
                      //                            , self->g_V_ws_H2[ws_counter].keys1[self->g_V_ws_H2[ws_counter].k1_ptr].last);
                      //    coH2_print_v_implicit(self, ws_counter);

                      //}

                      return;

                }

            }
            
      }

      for (EDGE_ID mm = 0; mm < self->g_V_ws_H2[ws_counter].last; mm++){

            
            if (self->g_V_ws_H2[ws_counter].keys1[mm].k1 == phi->low.key1){
                
                  
                  //check_space_implicit_keys2(&(v_implicit->keys1[mm]));

                  if (self->g_V_ws_H2[ws_counter].keys1[mm].last ==\
                                                  self->g_V_ws_H2[ws_counter].keys1[mm].max_len){

                        self->g_V_ws_H2[ws_counter].keys1[mm].max_len += 10;
                        self->g_V_ws_H2[ws_counter].keys1[mm].keys2 = (coH2_implicit_keys2*)realloc\
                                                                 (self->g_V_ws_H2[ws_counter].keys1[mm].keys2\
                                                                , self->g_V_ws_H2[ws_counter].keys1[mm].max_len*sizeof(coH2_implicit_keys2));
                        
                  }

                  //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
                  //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

                  //    printf("\nin c5, idx %d, last is %d"\
                  //                            , self->g_V_ws_H2[ws_counter].k1_ptr\
                  //                            , self->g_V_ws_H2[ws_counter].keys1[self->g_V_ws_H2[ws_counter].k1_ptr].last);
                  //    coH2_print_v_implicit(self, ws_counter);
                  //}

                  coH2_implicit_keys1* this_key1 = &(this_ws->keys1[mm]);
                  coH2_implicit_keys2* this_key2 = &(this_key1->keys2[this_key1->last]);

                  this_ws->keys1[mm].flag_empty = 0;
                  

                  this_key2->k2 = phi->low.key2;
                  this_key2->o_abc = phi->triangle;
                  this_key2->a_ptr = phi->a_ptr;
                  this_key2->b_ptr = phi->b_ptr;
                  this_key2->c_ptr = phi->c_ptr;
                  this_key2->vertex = phi->vertex;
                  this_key2->flag_next = flag_next;
                
                  this_key1->last++;

                  //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
                  //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

                  //    printf("\nin c6, idx %d, last is %d"\
                  //                            , self->g_V_ws_H2[ws_counter].k1_ptr\
                  //                            , self->g_V_ws_H2[ws_counter].keys1[self->g_V_ws_H2[ws_counter].k1_ptr].last);
                  //    coH2_print_v_implicit(self, ws_counter);
                  //}
                  return;

                  
            }
            
            
      }


      if (self->g_V_ws_H2[ws_counter].last == self->g_V_ws_H2[ws_counter].max_len){
            
            EDGE_ID mm = self->g_V_ws_H2[ws_counter].max_len;

            self->g_V_ws_H2[ws_counter].max_len += 10;
            self->g_V_ws_H2[ws_counter].keys1 = (coH2_implicit_keys1*)realloc(self->g_V_ws_H2[ws_counter].keys1\
                                                                  , self->g_V_ws_H2[ws_counter].max_len*sizeof(coH2_implicit_keys1));

            while (mm < self->g_V_ws_H2[ws_counter].max_len){
                  
                  this_ws->keys1[mm].flag_empty = 1;
                  this_ws->keys1[mm].max_len = 10;
                  this_ws->keys1[mm].last = 0;
                  self->g_V_ws_H2[ws_counter].keys1[mm].keys2 = (coH2_implicit_keys2*)malloc(10*sizeof(coH2_implicit_keys2));

                  mm++;
                
                  
            }


      }


      //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
      //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

      //    printf("\nin c7, idx %d, last is %d"\
      //                            , self->g_V_ws_H2[ws_counter].k1_ptr\
      //                            , self->g_V_ws_H2[ws_counter].keys1[self->g_V_ws_H2[ws_counter].k1_ptr].last);
      //    coH2_print_v_implicit(self, ws_counter);
      //}

      coH2_implicit_keys1* this_key1 = &(this_ws->keys1[this_ws->last]);
      coH2_implicit_keys2* this_key2 = &(this_key1->keys2[0]);

      this_key1->flag_empty = 0;
      this_key1->k1 = phi->low.key1;
      this_key2->k2 = phi->low.key2;

      this_key2->o_abc = phi->triangle;
      this_key2->a_ptr = phi->a_ptr;
      this_key2->b_ptr = phi->b_ptr;
      this_key2->c_ptr = phi->c_ptr;
      this_key2->vertex = phi->vertex;
      this_key2->flag_next = flag_next;

      this_key1->last = 1;
      this_ws->last++;


      //if ((this_ws->triangle.key1 == self->g_debug_triangle.key1) &&\
      //            (this_ws->triangle.key2 == self->g_debug_triangle.key2)){

      //    printf("\nin c7, idx %d, last is %d"\
      //                            , self->g_V_ws_H2[ws_counter].k1_ptr\
      //                            , self->g_V_ws_H2[ws_counter].keys1[self->g_V_ws_H2[ws_counter].k1_ptr].last);
      //    coH2_print_v_implicit(self, ws_counter);
      //}

      return;


        
}





BIGINT compute_num_simplices(filtration* self){

    
      self->g_n_all_simp  = (BIGINT)(self->g_n_vert) + (BIGINT)(self->g_n_valid_edges);
      printf("\n");

      for (EDGE_ID o_ab = 0; o_ab < self->g_n_valid_edges; o_ab++){


        printf("\redge%d", o_ab);
            
         VERT_ID a = self->g_edges_list[2*o_ab];
         VERT_ID b = self->g_edges_list[2*o_ab+1];


         VERT_ID a_ptr = 0;
         VERT_ID b_ptr = 0;



         while ((a_ptr < self->g_Neigh_len[a])\
              && (b_ptr < self->g_Neigh_len[b])){

              
              if (self->g_Neighbors[a][a_ptr].neighbor < self->g_Neighbors[b][b_ptr].neighbor){

                    a_ptr++;

              }
              else if (self->g_Neighbors[a][a_ptr].neighbor > self->g_Neighbors[b][b_ptr].neighbor){

                    b_ptr++;
              }
              else{
                    
                    
                    VERT_ID c = self->g_Neighbors[a][a_ptr].neighbor;
                    
                    EDGE_ID o_ac = self->g_Neighbors[a][a_ptr].order;

                    if (o_ac > o_ab){
                        a_ptr++;
                        b_ptr++;
                        continue;
                    }

                    EDGE_ID o_bc = self->g_Neighbors[b][b_ptr].order;

                    if (o_bc > o_ab){
                        a_ptr++;
                        b_ptr++;
                        continue;
                    }

                    // This is a valid triangle
                    self->g_n_all_simp++;
                    


                    for (VERT_ID mm = 0; mm < self->g_Neigh_len[c]; mm++){

                          if (self->g_Neighbors[c][mm].neighbor < c){
                            continue;
                          }

                          VERT_ID d = self->g_Neighbors[c][mm].neighbor;
                          
                          VERT_ID idx = search_Neighbors(self, a, d, 0, self->g_Neigh_len[a]-1);
                          if (idx == self->g_n_vert) continue;
                          EDGE_ID o_ad = self->g_Neighbors[a][idx].order;
                          if (o_ad > o_ab) continue;

                          idx = search_Neighbors(self, b, d, 0, self->g_Neigh_len[b]-1);
                          if (idx == self->g_n_vert) continue;
                          EDGE_ID o_bd = self->g_Neighbors[b][idx].order;
                          if (o_bd > o_ab) continue;

                          idx = search_Neighbors(self, c, d, 0, self->g_Neigh_len[c]-1);
                          if (idx == self->g_n_vert) continue;
                          EDGE_ID o_cd = self->g_Neighbors[c][idx].order;
                          if (o_cd > o_ab) continue;

                          self->g_n_all_simp++;

                          //printf("\n %d, %d, %d, %d", a, b, c, d);
                          
                          
                    }
                    a_ptr++;
                    b_ptr++;


              }
            
         }
             
      }

      printf("\nNumber of simplices %llu\n", self->g_n_all_simp);

      
}



void compute_H1_homology_cycles(filtration* self){

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
//            STEP H1.1: Find homology now for the triangles
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


     if (!self->g_suppress_output){
        printf("\n\n---------------");
        printf("\nComputing H1...");
        printf("\n---------------\n");
     }



     struct timespec start_wall_clock, finish_wall_clock;

     clock_gettime(CLOCK_MONOTONIC, &start_wall_clock);



     self->g_R_max_len_H1 = 100;
     self->g_R_H1 = (EDGE_ID*)malloc(self->g_R_max_len_H1*sizeof(EDGE_ID));
     self->g_R_len_H1 = 0;


     self->g_R_col_idx_max_len_H1 = 100;
     self->g_R_col_idx_H1 = (EDGE_ID*)malloc(self->g_R_col_idx_max_len_H1*sizeof(EDGE_ID));

     self->g_R_col_idx_H1_ptr = 0;


     self->g_workspace_size = 1000;

     self->g_ws_pre_alloc = 1000;
     // Initialize ws counter
     self->g_ws_counter = 0;

     // H1 workspace structures
     self->g_workspace_H1 = (EDGE_ID**)malloc(self->g_workspace_size*sizeof(EDGE_ID*));

     // H1 workspace info
     self->g_workspace_H1_info = (boundary_H1_ws*)malloc(self->g_workspace_size*sizeof(boundary_H1_ws));



     for (int i = 0; i < self->g_workspace_size; i++){

         self->g_workspace_H1_info[i].max_len = self->g_ws_pre_alloc;

         self->g_workspace_H1[i] = (EDGE_ID*)malloc(2*self->g_workspace_H1_info[i].max_len*sizeof(EDGE_ID));

         self->g_workspace_H1_info[i].trivial_boundary = (EDGE_ID*)malloc(3*sizeof(EDGE_ID));


     }

     // Pivots
     self->g_pivots_H1 = (EDGE_ID*)calloc(self->g_n_valid_edges, sizeof(EDGE_ID));

     // Convenient info for pers pairs
     self->g_homH1_pers_len = 0;
     self->g_homH1_pers_max_len = 100;
     self->g_homH1_pers = (homH1_pers*)malloc(self->g_homH1_pers_max_len*sizeof(homH1_pers));

     // Temporary space for birth cycles
     self->g_temp_V_primary.max_len = 10;
     self->g_temp_V_primary.VV = (EDGE_ID*)malloc(self->g_temp_V_primary.max_len*sizeof(EDGE_ID));
     self->g_temp_V_primary.len = 0;

     // Temporary space for birth cycles
     self->g_temp_R_birth_cycles.max_len = 100;
     self->g_temp_R_birth_cycles.RR = (EDGE_ID*)malloc(2*self->g_temp_R_birth_cycles.max_len*sizeof(EDGE_ID));
     self->g_temp_R_birth_cycles.original = 0;
     self->g_temp_R_birth_cycles.len = 0;


//#ifdef HOM_CYCLES
     // Need this info for birth voids
     if (self->g_compute_cycles){
          self->g_H1_pivot_of = (V_H1*)malloc(self->g_n_valid_edges*sizeof(V_H1));
     }

//#endif

#ifdef ADAPTIVE_V_STORAGE
     // Create pointers to store V

     if (self->g_compute_cycles){
        for (EDGE_ID mm = 0; mm < self->g_n_vert; mm++){

           self->g_H0_pivot_of[mm].V_usage = 0;
           self->g_H0_pivot_of[mm].V_stored = 0;
           self->g_H0_pivot_of[mm].V_len = 0;
           self->g_H0_pivot_of[mm].VV = NULL;

        }
     }


     // Create pointers to store V per extraction call
     self->g_store_V_for_len = 0;
     self->g_store_V_for_max_len = 10;
     self->g_store_V_for = (EDGE_ID*)malloc(self->g_store_V_for_max_len*sizeof(EDGE_ID));


#endif

//#ifdef MINIMIZE_BIRTH_CYCLES
     //if (self->g_reduce_cyc_lengths){

        self->g_all_V_stored_num = 0;
        self->g_all_V_stored_max_num = 10;
        //self->g_all_V_stored_len = (EDGE_ID*)calloc(self->g_all_V_stored_max_num, sizeof(EDGE_ID));

        //self->g_all_V_H0_stored = (EDGE_ID**)malloc(self->g_all_V_stored_max_num*sizeof(EDGE_ID*));
        self->g_all_V_H0_stored = (cyc_info*)malloc(self->g_all_V_stored_max_num*sizeof(cyc_info));
         
        self->g_edges_in_cycles = (EDGE_ID**)malloc(self->g_n_valid_edges*sizeof(EDGE_ID*));
        self->g_edges_in_cycles_len = (EDGE_ID*)calloc(self->g_n_valid_edges, sizeof(EDGE_ID));

     //}

//#endif


#ifdef MINIMIZE_HOM_CYCLES
     self->g_all_R_hom_stored_num = 0;
     self->g_all_R_hom_stored_max_num = 10;

     //self->g_all_R_hom_stored_len = (EDGE_ID*)calloc(self->g_all_V_hom_stored_max_num, sizeof(EDGE_ID));
     //self->g_all_V_hom_H1_stored = (EDGE_ID**)malloc(self->g_all_V_hom_stored_max_num*sizeof(EDGE_ID*));

     self->g_all_R_hom_H1_stored = (cyc_info*)malloc(self->g_all_R_hom_stored_max_num*sizeof(cyc_info));

#endif

     ////////////////////////////////////////////////////////////////
     //
     //   Allocate jobs for parallel H1
     // 
     ////////////////////////////////////////////////////////////////

     self->g_jobs = (int*)malloc((self->g_cpu_count + 1)*sizeof(int));

     allocate_jobs(self, self->g_workspace_size);

     self->g_threads = (pthread_t *)malloc(self->g_cpu_count*sizeof(pthread_t));

     int rtn;

     if ((rtn = pthread_mutex_init(&(self->g_thread_lock), NULL)) !=0)
        fprintf(stderr, "pthread_mutex_init %s", strerror(rtn)), exit(-1);

     if ((rtn = pthread_cond_init(&(self->g_start_boss), NULL)) !=0)
        fprintf(stderr, "pthread_cond_init %s", strerror(rtn)), exit(-1);

     if ((rtn = pthread_cond_init(&(self->g_start_workers), NULL)) !=0)
        fprintf(stderr, "pthread_cond_init %s", strerror(rtn)), exit(-1);


     // Initialize thread creation
     self->g_thread_id = 0;
     self->g_sleeping_threads = 0;
     self->g_delete_threads = 0;

     for (int i = 0; i < self->g_cpu_count; i++){

        if ((rtn = pthread_create( \
                                &(self->g_threads[i]) \
                                , NULL \
                                , reduce_with_complex_H1 \
                                , (void*)self)!= 0))
          fprintf(stderr, "pthread_create %d", rtn), exit(-1);
      
     }

     // Wait for threads to be initialized
     pthread_mutex_lock(&(self->g_thread_lock));

     while(self->g_sleeping_threads != self->g_cpu_count){
        
          pthread_cond_wait(&(self->g_start_boss) \
                          , &(self->g_thread_lock));

     }

     ////////////////////////////////

     // STEP 1: Compute R (Note: Do not include trivial pairs)


     for (EDGE_ID o_ab = 0; o_ab < self->g_n_valid_edges; o_ab++){
          
        if (!self->g_H1_cohom_pivots_len[o_ab]){
            continue;
        }

        VERT_ID a = self->g_edges_list[2*o_ab];
        VERT_ID b = self->g_edges_list[2*o_ab+1];
        
        for (VERT_ID mm = 0; mm < self->g_H1_cohom_pivots_len[o_ab]; mm++){
              
            
              //This triangle is in a non-trivial persistence pair

              // Workspace attributes
              boundary_H1_ws* this_ws = self->g_workspace_H1_info + self->g_ws_counter;

              // Initially, the original is at 0
              this_ws->original = 0;

              this_ws->flag_first = 1;

              // Parallel control flags
              this_ws->flag_empty = 0;
              this_ws->flag_red_w_complex = 0;
              this_ws->flag_append_to_complex = 1;

              this_ws->triangle.key1 = o_ab;
              this_ws->triangle.key2 = self->g_H1_cohom_pivots[o_ab][mm].key2;

              // Initial length of boundary
              this_ws->len = 3;

              compute_boundary_triangle(self, this_ws->triangle, self->g_workspace_H1[self->g_ws_counter]);

              this_ws->pivot = o_ab;

              //printf("\nOutside the boundary is (%d, %d, %d)"\
              //                                      , self->g_workspace_H1[self->g_ws_counter][0]\
              //                                      , self->g_workspace_H1[self->g_ws_counter][1]\
              //                                      , self->g_workspace_H1[self->g_ws_counter][2]\
              //                                      );


              self->g_ws_counter++;


              if (self->g_ws_counter == self->g_workspace_size){
                    
                    reduce_ws_H1(self);
              }
              
        }

          
     }




     //printf("\n press key for the last batch");
     //self->g_new_debug2 = 1;
     //getchar();

     // Reduction of final batch
     while (self->g_ws_counter){
          
          allocate_jobs(self, self->g_ws_counter);
          reduce_ws_H1(self);

     }


     //printf("\nComputed H1.");
     //getchar();


     /////////////////////////
     // Cancel the threads
     /////////////////////////

     self->g_delete_threads = 1;

     pthread_cond_broadcast(&(self->g_start_workers));

     pthread_mutex_unlock(&(self->g_thread_lock));

     for (int i = 0; i < self->g_cpu_count; i++){

        pthread_join(self->g_threads[i], NULL);
      
     }



     free(self->g_jobs);
     free(self->g_threads);



     // CANNOT FREE THIS IS H2 CYCLES ARE NEEDED
     //free(self->g_R_H1);
     //free(self->g_R_col_idx_H1);
     //free(self->g_pivots_H1);


     for (int i = 0; i < self->g_workspace_size; i++){

         free(self->g_workspace_H1_info[i].trivial_boundary);
         free(self->g_workspace_H1[i]);

     }


     free(self->g_workspace_H1);
     free(self->g_workspace_H1_info);



     if (!self->g_suppress_output){
        clock_gettime(CLOCK_MONOTONIC, &finish_wall_clock);
        self->g_timer_computeH1 = (finish_wall_clock.tv_sec - start_wall_clock.tv_sec);
        self->g_timer_computeH1 += (finish_wall_clock.tv_nsec - start_wall_clock.tv_nsec) / 1000000000.0;

        printf("\nComputed H1.");

     ////////////////////////
     // HOMOLOGY CYCLES
     ////////////////////////
     
        clock_gettime(CLOCK_MONOTONIC, &start_wall_clock);

     }

     FILE* fp2 = fopen(self->g_homH1_cycles_file, "w");

     PAR birth, death;

     self->g_n_H1_birth_cycles = 0;
     self->g_n_H0_stored_V = 0;

     int add_flag = 0;
     // Go over the pers pairs of features that died and compute the cycles
     for (EDGE_ID mm = 0; mm < self->g_homH1_pers_len; mm++){
          
          //if (n_cycles%1000 == 0) printf("\r Computing cycle num %llu", n_cycles);

          self->g_n_H1_birth_cycles++;

          if (self->g_filetype == 1){
              birth = sqrt(self->g_edge_parameter[self->g_homH1_pers[mm].birth_edge]);
              death = sqrt(self->g_edge_parameter[self->g_homH1_pers[mm].death_triangle_key1]);

          }
          else{
              birth = self->g_edge_parameter[self->g_homH1_pers[mm].birth_edge];
              death = self->g_edge_parameter[self->g_homH1_pers[mm].death_triangle_key1];
          }
        
          fprintf(fp2, "%lf, %lf", birth , death);

          fprintf(fp2, "\nhomology cycle");

#ifdef MINIMIZE_HOM_CYCLES
          self->g_all_V_hom_stored_len[self->g_all_V_hom_stored_num] =\
                                                self->g_R_col_idx_H1[self->g_homH1_pers[mm].R_col_idx + 1]\
                                                - self->g_R_col_idx_H1[self->g_homH1_pers[mm].R_col_idx];

          self->g_all_V_hom_H1_stored[self->g_all_V_hom_stored_num] =\
                                             (EDGE_ID*)\
                                             malloc(self->g_all_V_hom_stored_len[self->g_all_V_hom_stored_num]\
                                                                      *sizeof(EDGE_ID));
          EDGE_ID bbb = 0;
#endif

          for (EDGE_ID bb = self->g_R_col_idx_H1[self->g_homH1_pers[mm].R_col_idx]\
              ; bb < self->g_R_col_idx_H1[self->g_homH1_pers[mm].R_col_idx + 1]\
              ; bb++){
                
                fprintf(fp2, ", %d, %d", self->g_edges_list[2*self->g_R_H1[bb]]\
                                       , self->g_edges_list[2*self->g_R_H1[bb]+1]);

#ifdef MINIMIZE_HOM_CYCLES
                  
                self->g_all_V_hom_H1_stored[self->g_all_V_hom_stored_num][bbb++] = self->g_R_H1[bb];
#endif

          }

            
#ifdef MINIMIZE_HOM_CYCLES

          self->g_all_V_hom_stored_num++;
          if (self->g_all_V_hom_stored_num == self->g_all_V_hom_stored_max_num){
                self->g_all_V_hom_stored_max_num += 100;
                self->g_all_V_hom_H1_stored = (EDGE_ID**)realloc(self->g_all_V_hom_H1_stored\
                                                        , self->g_all_V_hom_stored_max_num*sizeof(EDGE_ID*));

                self->g_all_V_hom_stored_len = (EDGE_ID*)realloc(self->g_all_V_hom_stored_len\
                                                        , self->g_all_V_hom_stored_max_num*sizeof(EDGE_ID));
          }

#endif
          
          // Always store the birth cycles for now
          fprintf(fp2, "\nbirth cycle");

          //printf("\nGetting birth cycle %d out of %d", mm, self->g_homH1_pers_len);
          //getchar();


#ifdef ADAPTIVE_V_STORAGE
          self->g_store_V_for_len = 0;
#endif

          get_birth_cycle(self, self->g_homH1_pers[mm].birth_edge);


//#ifdef MINIMIZE_BIRTH_CYCLES
          //if (self->g_reduce_cyc_lengths){

              ///self->g_all_V_stored_len[self->g_all_V_stored_num] = self->g_temp_V_primary.len;
              ///self->g_all_V_H0_stored[self->g_all_V_stored_num] =\
              ///                                   (EDGE_ID*)malloc(self->g_temp_V_primary.len*sizeof(EDGE_ID));


              //add_flag = 1;
              //if (death - birth > 4){
              //    add_flag = 1;
              //}

              //if (add_flag){
                  self->g_all_V_H0_stored[self->g_all_V_stored_num].boundary =\
                                                     (EDGE_ID*)malloc(self->g_temp_V_primary.len*sizeof(EDGE_ID));

                  self->g_all_V_H0_stored[self->g_all_V_stored_num].len = self->g_temp_V_primary.len;


                  //self->g_all_V_H0_stored[self->g_all_V_stored_num].ops_len = 1;
                  //self->g_all_V_H0_stored[self->g_all_V_stored_num].ops = (EDGE_ID*)malloc(sizeof(EDGE_ID));
                  //self->g_all_V_H0_stored[self->g_all_V_stored_num].ops[0] = self->g_all_V_stored_num;


                  self->g_all_V_H0_stored[self->g_all_V_stored_num].perspair[0] = birth;
                  self->g_all_V_H0_stored[self->g_all_V_stored_num].perspair[1] = death;

                  self->g_all_V_H0_stored[self->g_all_V_stored_num].updated_birth = birth;


                  //printf("\n%d idx is pers pair (%lf, %lf) has len %d"\
                  //                                          , self->g_all_V_stored_num\
                  //                                          , birth\
                  //                                          , death\
                  //                                          , self->g_temp_V_primary.len);

              //}


          //}

//#endif



          for (EDGE_ID nn = 0; nn < self->g_temp_V_primary.len; nn++){
               
                 fprintf(fp2, ", %d, %d", self->g_edges_list[2*self->g_temp_V_primary.VV[nn]]\
                                        , self->g_edges_list[2*self->g_temp_V_primary.VV[nn]+1]);

//#ifdef MINIMIZE_BIRTH_CYCLES
                 //if (self->g_reduce_cyc_lengths){
                      self->g_all_V_H0_stored[self->g_all_V_stored_num].boundary[nn] = self->g_temp_V_primary.VV[nn];
                 //}
//#endif
               
          }


          fprintf(fp2, "\n");


//#ifdef MINIMIZE_BIRTH_CYCLES

          //if (add_flag){

          //if (self->g_reduce_cyc_lengths){

              self->g_all_V_stored_num++;
              if (self->g_all_V_stored_num == self->g_all_V_stored_max_num){
                    self->g_all_V_stored_max_num += 100;
                    self->g_all_V_H0_stored = (cyc_info*)realloc(self->g_all_V_H0_stored\
                                                            , self->g_all_V_stored_max_num*sizeof(cyc_info));

                    //self->g_all_V_stored_len = (EDGE_ID*)realloc(self->g_all_V_stored_len\
                    //                                        , self->g_all_V_stored_max_num*sizeof(EDGE_ID));
              }

          //}

//#endif





#ifdef ADAPTIVE_V_STORAGE

          store_V_H0(self);

#endif

          
     }



     // Go over the pers pairs of undead features and compute the birth cycles
     for (EDGE_ID mm = 0; mm < self->g_H1_undead_ptr; mm++){

          self->g_n_H1_birth_cycles++;

          if (self->g_filetype == 1){
              birth = sqrt(self->g_edge_parameter[self->g_H1_undead[mm]]);

          }
          else{
              birth = self->g_edge_parameter[self->g_H1_undead[mm]];
          }
        

          fprintf(fp2, "%lf, -1", birth);
          
#ifdef ADAPTIVE_V_STORAGE
          self->g_store_V_for_len = 0;
#endif

          get_birth_cycle(self, self->g_H1_undead[mm]);

//#ifdef MINIMIZE_BIRTH_CYCLES
          //if (self->g_reduce_cyc_lengths){

              //self->g_all_V_stored_len[self->g_all_V_stored_num] = self->g_temp_V_primary.len;
              //self->g_all_V_H0_stored[self->g_all_V_stored_num] =\
              //                                   (EDGE_ID*)malloc(self->g_temp_V_primary.len*sizeof(EDGE_ID));
              self->g_all_V_H0_stored[self->g_all_V_stored_num].boundary =\
                                                 (EDGE_ID*)malloc(self->g_temp_V_primary.len*sizeof(EDGE_ID));

              self->g_all_V_H0_stored[self->g_all_V_stored_num].len = self->g_temp_V_primary.len;

              //self->g_all_V_H0_stored[self->g_all_V_stored_num].ops_len = 1;
              //self->g_all_V_H0_stored[self->g_all_V_stored_num].ops = (EDGE_ID*)malloc(sizeof(EDGE_ID));
              //self->g_all_V_H0_stored[self->g_all_V_stored_num].ops[0] = self->g_all_V_stored_num;

              self->g_all_V_H0_stored[self->g_all_V_stored_num].perspair[0] = birth;
              self->g_all_V_H0_stored[self->g_all_V_stored_num].perspair[1] = -1;

              self->g_all_V_H0_stored[self->g_all_V_stored_num].updated_birth = birth;

          //}
//#endif
          
          fprintf(fp2, "\nbirth cycle");

          for (EDGE_ID nn = 0; nn < self->g_temp_V_primary.len; nn++){
               
                 fprintf(fp2, ", %d, %d", self->g_edges_list[2*self->g_temp_V_primary.VV[nn]]\
                                        , self->g_edges_list[2*self->g_temp_V_primary.VV[nn]+1]);

//#ifdef MINIMIZE_BIRTH_CYCLES
                 //if (self->g_reduce_cyc_lengths){
                      self->g_all_V_H0_stored[self->g_all_V_stored_num].boundary[nn] = self->g_temp_V_primary.VV[nn];
                 //}
//#endif
               
          }


          fprintf(fp2, "\n");


//#ifdef MINIMIZE_BIRTH_CYCLES

          //if (self->g_reduce_cyc_lengths){

              self->g_all_V_stored_num++;
              if (self->g_all_V_stored_num == self->g_all_V_stored_max_num){
                    self->g_all_V_stored_max_num += 100;
                    self->g_all_V_H0_stored = (cyc_info*)realloc(self->g_all_V_H0_stored\
                                                            , self->g_all_V_stored_max_num*sizeof(cyc_info));
                    //self->g_all_V_stored_len = (EDGE_ID*)realloc(self->g_all_V_stored_len\
                    //                                        , self->g_all_V_stored_max_num*sizeof(EDGE_ID));
              }

          //}

//#endif



#ifdef ADAPTIVE_V_STORAGE

          store_V_H0(self);

#endif


     }




     fclose(fp2);


     free(self->g_homH1_pers);
     free(self->g_temp_V_primary.VV);
     free(self->g_temp_R_birth_cycles.RR);

 

#ifdef RECORD_V_USAGE
     FILE* fp3 = fopen(self->g_V_H0_usage_file, "w");
#endif


     for (EDGE_ID mm = 0; mm < self->g_n_vert; mm++){
          
        if (self->g_H0_pivot_of[mm].V_len){
#ifdef RECORD_V_USAGE
            fprintf(fp3, "%d, %d\n"\
                                  , self->g_H0_pivot_of[mm].V_usage\
                                  , self->g_H0_pivot_of[mm].V_depth);
#endif
            //printf("\n%d is used %d times with depth %d"\
            //                        , mm\
            //                        , self->g_H0_pivot_of[mm].V_usage\
            //                        , self->g_H0_pivot_of[mm].V_depth\
            //                        );
            free(self->g_H0_pivot_of[mm].VV);
        }
          
     }

#ifdef RECORD_V_USAGE
     fclose(fp3);
#endif

     free(self->g_H0_pivot_of);
     //printf("\nPress key to continue...");
     //getchar();

#ifdef ADAPTIVE_V_STORAGE
     free(self->g_store_V_for);
#endif

     if (!self->g_suppress_output){
        clock_gettime(CLOCK_MONOTONIC, &finish_wall_clock);
        self->g_timer_H1cycles = (finish_wall_clock.tv_sec - start_wall_clock.tv_sec);
        self->g_timer_H1cycles += (finish_wall_clock.tv_nsec - start_wall_clock.tv_nsec) / 1000000000.0;
     }



     //if (self->g_reduce_cyc_lengths){

     struct timespec start_wall_clock2, finish_wall_clock2;

     if (!self->g_suppress_output){
        printf("\nMinimizing birth cycles...");

        clock_gettime(CLOCK_MONOTONIC, &start_wall_clock2);
     }


     minimize_birth_cycles_H0_v3(self\
                            , self->g_all_V_H0_stored\
                            , self->g_all_V_stored_num\
                            , self->g_minimal_V_H0_file\
                            , self->g_V_H0_birthcyc_lens_file\
                            , self->g_minimal_V_H0_birthcyc_lens_file\
                            , self->g_birth_subset_points_file_H0\
                            );

                            //, self->g_minimal_V_H0_in_cycles_file\

     if (!self->g_suppress_output){
        clock_gettime(CLOCK_MONOTONIC, &finish_wall_clock2);
        self->g_timer_minimize_H1cycles = (finish_wall_clock2.tv_sec - start_wall_clock2.tv_sec);
        self->g_timer_minimize_H1cycles += (finish_wall_clock2.tv_nsec - start_wall_clock2.tv_nsec) / 1000000000.0;
     }

     //}
     //else{

     //   for (EDGE_ID ci = 0; ci < self->g_all_V_stored_num; ci++){
     //      free(self->g_all_V_H0_stored[ci].boundary);
     //   }

     //}



#ifdef MINIMIZE_HOM_CYCLES

     if (!self->g_suppress_output){
        printf("\nMinimizing hom cycles...");
        clock_gettime(CLOCK_MONOTONIC, &start_wall_clock2);
     }

     minimize_birth_cycles_H0_v3(self, self->g_all_V_hom_H1_stored\
                            , self->g_all_V_hom_stored_len\
                            , self->g_all_V_hom_stored_num\
                            , self->g_minimal_V_hom_H1_file\
                            , self->g_birth_subset_points_file_H0\
                            );

     if (!self->g_suppress_output){
        clock_gettime(CLOCK_MONOTONIC, &finish_wall_clock2);
        self->g_timer_minimize_H1_homcycles = (finish_wall_clock2.tv_sec - start_wall_clock2.tv_sec);
        self->g_timer_minimize_H1_homcycles += (finish_wall_clock2.tv_nsec - start_wall_clock2.tv_nsec) / 1000000000.0;
     }

#endif





     if (!self->g_suppress_output){
        printf("\nComputed homology and birth cycles for H1.");
     }
     //getchar();
     

}





void reduce_ws_H1(filtration* self){


      if (self->g_new_debug2){
            for (int kk = 0; kk < self->g_ws_counter; kk++){
                  
                  printf("\n%d has triangle (%d, %d) with pivot %d", kk\
                                                                  , self->g_workspace_H1_info[kk].triangle.key1\
                                                                  , self->g_workspace_H1_info[kk].triangle.key2\
                                                                  , self->g_workspace_H1_info[kk].pivot\
                                                                  );
            }
            printf("\nbefore parallel. press key to start parallel");
            //getchar();
      }
      self->g_processed_threads = 0;

      //printf("\npress key to reduce with complex");
      //getchar();

      pthread_cond_broadcast(&(self->g_start_workers));

      while (self->g_processed_threads != self->g_cpu_count){
            
            pthread_cond_wait(&(self->g_start_boss) \
                            ,&(self->g_thread_lock));
      }

      //printf("\npress key to reduce with self");
      //getchar();


      if (self->g_new_debug2){
            for (int kk = 0; kk < self->g_ws_counter; kk++){
                  
                  printf("\n%d has triangle (%d, %d) with pivot %d", kk\
                                                                  , self->g_workspace_H1_info[kk].triangle.key1\
                                                                  , self->g_workspace_H1_info[kk].triangle.key2\
                                                                  , self->g_workspace_H1_info[kk].pivot\
                                                                  );
            }
            printf("\nafter parallel. press key to start serial");
            //getchar();
      }

      reduce_with_self_H1( \
                            self \
                            );

      if (self->g_new_debug2){
            for (int kk = 0; kk < self->g_ws_counter; kk++){
                  
                  printf("\n%d has triangle (%d, %d) with pivot %d", kk\
                                                                  , self->g_workspace_H1_info[kk].triangle.key1\
                                                                  , self->g_workspace_H1_info[kk].triangle.key2\
                                                                  , self->g_workspace_H1_info[kk].pivot\
                                                                  );
            }
            printf("\nafter serial. press key to update ");
            //getchar();
      }
      int count_valid = 0;

      for (int ws_counter=0; ws_counter < self->g_ws_counter; ws_counter++){

            if (self->g_workspace_H1_info[ws_counter].flag_append_to_complex){
                  

                  update_R_H1(self \
                                , ws_counter\
                                );
                  continue;
                    
            }
            
            //if (!self->g_workspace_H1_info[ws_counter].len){continue;}
            if (self->g_workspace_H1_info[ws_counter].flag_empty){

                  continue;
            }


            // Swap R
            EDGE_ID* temp = self->g_workspace_H1[count_valid];
            self->g_workspace_H1[count_valid] = self->g_workspace_H1[ws_counter];
            self->g_workspace_H1[ws_counter] = temp;

            // Swap R info
            boundary_H1_ws temp2 = self->g_workspace_H1_info[count_valid];
            self->g_workspace_H1_info[count_valid] = self->g_workspace_H1_info[ws_counter];
            self->g_workspace_H1_info[ws_counter] = temp2;


            // At this point, this has to be a non-zero column
            self->g_workspace_H1_info[count_valid].flag_empty = 0;

            count_valid += 1;

      }

      self->g_ws_counter = count_valid;

      if (self->g_new_debug2){
            for (int kk = 0; kk < self->g_ws_counter; kk++){
                  
                  printf("\n%d has triangle (%d, %d) with pivot %d", kk\
                                                                  , self->g_workspace_H1_info[kk].triangle.key1\
                                                                  , self->g_workspace_H1_info[kk].triangle.key2\
                                                                  , self->g_workspace_H1_info[kk].pivot\
                                                                  );
            }
            printf("\nafter update. press key to continue ");
            //getchar();
      }
      //if (dim)
      //  self->g_H0_MAX = self->g_n_reduced_simplex[dim];

}




void reduce_with_self_H1( \
                      filtration* self \
                      ){

    int compare;

    int i, m;
    int idx;

    EDGE_ID count, j, k;
    count = 0;


    for (i=0; i < self->g_ws_counter; i++){

        boundary_H1_ws* this_ws = self->g_workspace_H1_info + i;


        this_ws->flag_reduce = 0;

        // If the simplex has already been reduced to 0
        // then continue
        if (this_ws->flag_empty){ 

          this_ws->flag_append_to_complex = 0;
          continue;

        }

        
        EDGE_ID* orig = self->g_workspace_H1[i] + this_ws->original*this_ws->max_len;

        m = 0;
        while (m < i){

            boundary_H1_ws* m_ws = self->g_workspace_H1_info + m;

            if (m_ws->flag_empty){
                m++;
                continue;
            }


            EDGE_ID* original_m = self->g_workspace_H1[m] + m_ws->original*m_ws->max_len;

            orig = self->g_workspace_H1[i] + this_ws->original*this_ws->max_len;


            if (m_ws->pivot > this_ws->pivot){
                  
                  if (m_ws->flag_red_w_complex){
                        
                        this_ws->flag_append_to_complex = 0;
                        break;
                  }
                  m++;
                  continue;
                  
            }


            if (m_ws->pivot < this_ws->pivot){
                    m++;
                    continue;
            }

            if (m_ws->flag_red_w_complex){
                    this_ws->flag_append_to_complex = 0;
                    //m++;
                    //continue;
                    break;
            }


            if (this_ws->len + m_ws->len > this_ws->max_len ){

                //printf("\nReallocating inside self");
                //simp_max_len_i = len_i + len_m + 1000;


                if (this_ws->original){
                      
                      for (EDGE_ID mm = 0; mm < this_ws->len; mm++){
                            self->g_workspace_H1[i][mm] = self->g_workspace_H1[i][mm + this_ws->max_len];
                      }

                      this_ws->original = 0;

                }

                this_ws->max_len = this_ws->len + m_ws->len + 1000;

                self->g_workspace_H1[i] = (EDGE_ID*)realloc(self->g_workspace_H1[i]\
                                                      , 2*this_ws->max_len*sizeof(EDGE_ID));

                orig  = self->g_workspace_H1[i];


            }

            EDGE_ID* scratch = self->g_workspace_H1[i] + (1-this_ws->original)*this_ws->max_len;

            // Store the result in scratch

            count = 0;

            j = 0;
            k = 0;

            while ((j < this_ws->len) && (k < m_ws->len)){

                if (orig[j] < original_m[k]){
                    scratch[count++] = orig[j++];
                }
                else if (orig[j] > original_m[k]){
                    scratch[count++] = original_m[k++];
                }
                else{
                    j++;
                    k++;
                }

            }

            while (j < this_ws->len){
                    scratch[count++] = orig[j++];
            }

            while (k < m_ws->len){
                    scratch[count++] = original_m[k++];
            }

            this_ws->len = count;

            this_ws->original = 1 - this_ws->original;

            if (!count){
                  
                  this_ws->flag_append_to_complex = 0;
                  this_ws->flag_empty = 1;
                  break;
                  
            }

            this_ws->pivot = scratch[this_ws->len-1];
            //printf("\npivot is %d", this_ws->pivot);

            // Check if pivot is trivial or if it is pivot in H1
            if (self->g_coH1_all_lows[this_ws->pivot].low.key1 == this_ws->pivot){


                compute_boundary_triangle(self, self->g_coH1_all_lows[this_ws->pivot].low, this_ws->trivial_boundary);

                this_ws->R_col_idx = 0;
                this_ws->reduce_with_len = 3;
                this_ws->flag_reduce = 1;

                this_ws->flag_red_w_complex = 1;
                this_ws->flag_append_to_complex = 0;
                break;

                
            }
            else{

                idx = self->g_pivots_H1[this_ws->pivot];

                if (idx){
                    this_ws->R_col_idx = idx;
                    this_ws->reduce_with_len = self->g_R_col_idx_H1[idx+1] - self->g_R_col_idx_H1[idx];
                    this_ws->flag_reduce = 1;

                    this_ws->flag_red_w_complex = 1;
                    this_ws->flag_append_to_complex = 0;
                    break;
                }
                    


            }



            //idx = self->g_pivots_H1[this_ws->pivot];
            //// If the pivot is in red complex, then this has to be reduced w/ complex
            ////if (idx != self->g_n_reduced_simplex[self->g_dim_now]){
            //if (idx){
            //      
            //      this_ws->flag_red_w_complex = 1;
            //      this_ws->flag_append_to_complex = 0;
            //      break;

            //}


            //}

            m = 0;

        }//End of m loop

    }

}//End of red_ws_w_self_single



void* reduce_with_complex_H1(void* arg){
      
      filtration* self = arg;

      pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, 0);

      pthread_mutex_lock(&(self->g_thread_lock));

      int tid = ++self->g_thread_id;


      EDGE_ID *simp, *original_simp, *scratch_simp, *red_start;
      EDGE_ID j, k ,count;

      //EDGE_ID reduced_col;
      int i ,reduced_col ,idx;


      for (;;){

          self->g_sleeping_threads++;
          
          if (self->g_sleeping_threads == self->g_cpu_count)
              pthread_cond_signal(&(self->g_start_boss));

          pthread_cond_wait(&(self->g_start_workers), &(self->g_thread_lock));

          if (self->g_delete_threads){
            //printf("\nexiting from thread %d", tid);
            pthread_mutex_unlock(&(self->g_thread_lock));
            pthread_exit(NULL);
          }

          self->g_sleeping_threads--;

          pthread_mutex_unlock(&(self->g_thread_lock));

          for (i = self->g_jobs[tid - 1]; i < self->g_jobs[tid]; i++){

              boundary_H1_ws* this_ws = self->g_workspace_H1_info + i;

              if (this_ws->flag_first){

                  this_ws->flag_first = 0;

                  // Check if this is part of trivial pair

                  if (self->g_coH1_all_lows[this_ws->pivot].low.key1 == this_ws->pivot){
                        

                        compute_boundary_triangle(self\
                                            , self->g_coH1_all_lows[this_ws->pivot].low\
                                            , this_ws->trivial_boundary);

                        this_ws->R_col_idx = 0;
                        this_ws->reduce_with_len = 3;

                        this_ws->flag_reduce = 1;

                  }
                  else{

                        idx = self->g_pivots_H1[this_ws->pivot];

                        if (idx){

                            this_ws->R_col_idx = idx;
                            this_ws->reduce_with_len = self->g_R_col_idx_H1[idx+1] - self->g_R_col_idx_H1[idx];
                            this_ws->flag_reduce = 1;

                        }

                  }

              }


              if (this_ws->flag_empty){
                  // We are sure that we will exit only if there is no reduction
                  // required with existing complex or with trivial pair
                  this_ws->flag_red_w_complex = 0;
                  this_ws->flag_append_to_complex = 0;
                  continue;
              }

              this_ws->flag_red_w_complex = 0;
              this_ws->flag_append_to_complex = 1;


              EDGE_ID* orig = self->g_workspace_H1[i] \
                                  + this_ws->original*this_ws->max_len;

              EDGE_ID* scratch = self->g_workspace_H1[i] \
                                  + (1-this_ws->original)*this_ws->max_len;
                                    

              //printf("\nreducing with %d", idx);

              while(this_ws->flag_reduce){

                    //printf("\nreducing with %d", idx);
                    
                    //red_start = self->g_R_H1 + self->g_R_col_idx_H1[idx];

                    //len_a2 = self->g_R_col_idx_H1[idx+1] - self->g_R_col_idx_H1[idx];


                    if (this_ws->len + this_ws->reduce_with_len > this_ws->max_len){


                        if (this_ws->original){
                              
                              
                              for (EDGE_ID mm = 0; mm < this_ws->len; mm++){
                                    
                                    self->g_workspace_H1[i][mm] = self->g_workspace_H1[i][mm + this_ws->max_len];
                              }

                              this_ws->original = 0;
                              
                        }

                        this_ws->max_len = this_ws->len + this_ws->reduce_with_len + 1000;

                        self->g_workspace_H1[i] = (EDGE_ID*)realloc(self->g_workspace_H1[i]\
                                                            , 2*this_ws->max_len*sizeof(EDGE_ID));


                        orig = self->g_workspace_H1[i];

                        
                    }

                    if (!this_ws->R_col_idx){
                        red_start = this_ws->trivial_boundary;
                    }
                    else{

                        red_start = self->g_R_H1 + self->g_R_col_idx_H1[this_ws->R_col_idx];

                    }

                    
                    
                    scratch = self->g_workspace_H1[i] \
                                    + (1-this_ws->original)*this_ws->max_len;

                    count = 0;
                    j = 0;
                    k = 0;


                    while ((j < this_ws->len) && (k < this_ws->reduce_with_len)){

                        if (orig[j] < red_start[k]){
                            scratch[count++] = orig[j++];
                        }
                        else if (orig[j] > red_start[k]){
                            scratch[count++] = red_start[k++];
                        }
                        else{
                            j++;
                            k++;
                        }

                    }

                    while (j < this_ws->len){

                        scratch[count++] = orig[j++];
                    }

                    while (k < this_ws->reduce_with_len){

                        scratch[count++] = red_start[k++];
                    }

                    this_ws->original = 1 - this_ws->original;

                    this_ws->len = count;

                    if (!this_ws->len){

                        //idx = self->g_n_reduced_simplex[self->g_dim_now];
                        //idx = -1;
                        this_ws->flag_empty = 1;
                        break;

                    }

                    orig = self->g_workspace_H1[i] + this_ws->original*this_ws->max_len;

                    this_ws->pivot = orig[this_ws->len-1];

                    //printf("\npivot is %d", this_ws->pivot);

                    this_ws->flag_reduce = 0;

                    // First check if this is trivial pair
                    if (self->g_coH1_all_lows[this_ws->pivot].low.key1 == this_ws->pivot){

                        compute_boundary_triangle(self\
                                            , self->g_coH1_all_lows[this_ws->pivot].low\
                                            , this_ws->trivial_boundary);

                        this_ws->R_col_idx = 0;

                        this_ws->reduce_with_len = 3;
                        this_ws->flag_reduce = 1;
                        
                    }
                    else{

                        idx = self->g_pivots_H1[this_ws->pivot];

                        if (idx){
                            this_ws->R_col_idx = idx;
                            this_ws->reduce_with_len = self->g_R_col_idx_H1[idx+1] - self->g_R_col_idx_H1[idx];
                            this_ws->flag_reduce = 1;
                        }


                    }
                    

                    
                    //printf("\nidx is %d", idx);

                    

              }


          }

          pthread_mutex_lock(&(self->g_thread_lock));

          self->g_processed_threads++;

            
      }


}





void update_R_H1 (filtration* self, int ws_counter){
      
      //printf("\nupdating R H1");

      boundary_H1_ws* this_ws = self->g_workspace_H1_info + ws_counter;
      EDGE_ID* orig = self->g_workspace_H1[ws_counter] \
                                    + this_ws->original*this_ws->max_len;
      
      // Update R
      if ((self->g_R_len_H1 + this_ws->len) > self->g_R_max_len_H1){
            
            self->g_R_max_len_H1 += 1000 + this_ws->len;

            self->g_R_H1 = (EDGE_ID*)realloc(self->g_R_H1, self->g_R_max_len_H1*sizeof(EDGE_ID));
            
      }


      // Update R col idx
      self->g_R_col_idx_H1_ptr++;

      if (self->g_R_col_idx_H1_ptr == self->g_R_col_idx_max_len_H1 - 1){
              
            self->g_R_col_idx_max_len_H1 += 1000;

            self->g_R_col_idx_H1 = (EDGE_ID*)realloc(self->g_R_col_idx_H1\
                                          , self->g_R_col_idx_max_len_H1*sizeof(EDGE_ID));


      }

      //if (this_ws->pivot == 12631){
      //      printf("\nR is stored at R_col_idx %d: ", self->g_R_col_idx_H1_ptr);
      //      for (EDGE_ID mm = 0; mm < this_ws->len; mm++){
      //          printf("%d, ", orig[mm]);
      //      }
      //      getchar();

      //}


      //printf("\nAdding pivot %d at %d", this_ws->pivot, self->g_R_col_idx_H1_ptr);
      self->g_pivots_H1[this_ws->pivot] = self->g_R_col_idx_H1_ptr;

      self->g_R_col_idx_H1[self->g_R_col_idx_H1_ptr] = self->g_R_len_H1;


      for (EDGE_ID mm = 0; mm < this_ws->len; mm++){
            
            self->g_R_H1[self->g_R_len_H1++] = orig[mm];
            //printf("\nadded %d", self->g_R_H1[self->g_R_len_H1 - 1]);

      }

      
      self->g_R_col_idx_H1[self->g_R_col_idx_H1_ptr+1] = self->g_R_len_H1;

      PAR birth = self->g_edge_parameter[this_ws->pivot];
      PAR death = self->g_edge_parameter[this_ws->triangle.key1];

//#ifdef HOM_CYCLES
      if (self->g_compute_cycles){
            self->g_H1_pivot_of[this_ws->pivot].coface.key1 = this_ws->triangle.key1;
            self->g_H1_pivot_of[this_ws->pivot].coface.key2 = this_ws->triangle.key2;
            //if (this_ws->pivot == 2879)
            //printf("\nAdding pivot %d with simplex (%d, %d)", this_ws->pivot\
            //                                                , this_ws->triangle.key1\
            //                                                , this_ws->triangle.key2\
            //                                                );
      }
//#endif

      //printf("\n%lf, %lf", birth, death);

      if (death != birth){
          //printf("\n(%d, %d) has pivot %d at (%lf, %lf)", this_ws->triangle.key1\
          //                                , this_ws->triangle.key2\
          //                                , this_ws->pivot\
          //                                , birth\
          //                                , death\
          //                                );

          
          self->g_homH1_pers[self->g_homH1_pers_len].birth_edge = this_ws->pivot;
          self->g_homH1_pers[self->g_homH1_pers_len].death_triangle_key1 = this_ws->triangle.key1;
          self->g_homH1_pers[self->g_homH1_pers_len++].R_col_idx = self->g_R_col_idx_H1_ptr;

          if (self->g_homH1_pers_len == self->g_homH1_pers_max_len){
                self->g_homH1_pers_max_len += 100;
                self->g_homH1_pers = (homH1_pers*)realloc(self->g_homH1_pers\
                                                    , self->g_homH1_pers_max_len*sizeof(homH1_pers));
                 
          }


          //if (self->g_filetype == 1){
          //    fprintf(self->g_homH1_pers_file, "%0.12lf, %0.12lf\n", sqrt(birth), sqrt(death));
          //}
          //else{
          //    fprintf(self->g_homH1_pers_file, "%0.12lf, %0.12lf\n", birth, death);
          //}

          //for (EDGE_ID mm = 0; mm < this_ws->len; mm++){
          //      
          //      fprintf(self->g_homH1_cycles_file, "%d, %d,", self->g_edges_list[orig[mm]][0]\
          //                                                  , self->g_edges_list[orig[mm]][1]);
          //      //printf("\nadded %d", self->g_R_H1[self->g_R_len_H1 - 1]);

          //}
          //
          //fprintf(self->g_homH1_cycles_file, "\n");

      }

      

}



void get_birth_cycle(filtration* self, EDGE_ID bo_idx){

    
    self->g_temp_V_primary.len = 1;
    self->g_temp_V_primary.VV[0] = bo_idx;

    self->g_temp_R_birth_cycles.original = 0;

    self->g_temp_R_birth_cycles.RR[0] = self->g_edges_list[2*bo_idx];
    self->g_temp_R_birth_cycles.RR[1] = self->g_edges_list[2*bo_idx+1];
    self->g_temp_R_birth_cycles.len = 2;

    EDGE_ID* original_result;

    EDGE_ID* scratch_result;

    EDGE_ID j, k, count, possible_len, ro, red_simp_len, pivot;

    EDGE_ID* red_start;


    self->g_depth = 0;

#ifdef ADAPTIVE_V_STORAGE
    self->g_store_V_for_len = 0;
#endif
    

    while (self->g_temp_R_birth_cycles.len){
        

          original_result = self->g_temp_R_birth_cycles.RR \
                       + (self->g_temp_R_birth_cycles.original)*self->g_temp_R_birth_cycles.max_len;

          pivot = original_result[self->g_temp_R_birth_cycles.len-1];

          bo_idx = self->g_H0_pivot_of[pivot].coface;

          // FIND THE V RECURSIVELY
          find_V_recursively_edges(self, bo_idx, pivot);

          ro = (EDGE_ID)self->g_pivots_H0[pivot];

          red_simp_len = self->g_R_col_indices_H0[ro+1] - \
                           self->g_R_col_indices_H0[ro];


          red_start = self->g_R_sparse_H0 + self->g_R_col_indices_H0[ro];
        

          // Check for overflow
          possible_len = self->g_temp_R_birth_cycles.len + red_simp_len;
          if (possible_len > self->g_temp_R_birth_cycles.max_len - 1){

                if (self->g_temp_R_birth_cycles.original){

                      for (EDGE_ID k = 0; k < self->g_temp_R_birth_cycles.len; k++){
                          self->g_temp_R_birth_cycles.RR[k] =\
                                        self->g_temp_R_birth_cycles.RR[k + self->g_temp_R_birth_cycles.max_len];
                      }
                      
                }

                self->g_temp_R_birth_cycles.max_len = possible_len + 1000;
                self->g_temp_R_birth_cycles.RR = (EDGE_ID*)realloc(self->g_temp_R_birth_cycles.RR\
                                        , (2*self->g_temp_R_birth_cycles.max_len)*sizeof(EDGE_ID));

                original_result = self->g_temp_R_birth_cycles.RR;
                self->g_temp_R_birth_cycles.original = 0;
                
          }


          scratch_result = self->g_temp_R_birth_cycles.RR \
                       + (1-self->g_temp_R_birth_cycles.original)*self->g_temp_R_birth_cycles.max_len;

          // Reduce
          j = 0;
          k = 0;
          count = 0;

          while ((j < self->g_temp_R_birth_cycles.len) && (k < red_simp_len)){

              if (original_result[j] < red_start[k]){
                   
                  scratch_result[count] = original_result[j];
                  count = count + 1;
                  j = j + 1;
                   
              }
              else if (original_result[j] > red_start[k]){

                  scratch_result[count] = red_start[k];
                  count = count + 1;
                  k = k + 1;

              }
              else{

                  j = j + 1;
                  k = k + 1;

              }

          }

          while (j < self->g_temp_R_birth_cycles.len){

               scratch_result[count++] = original_result[j++];

          }

          while (k < red_simp_len){

               scratch_result[count++] = red_start[k++];

          }

          self->g_temp_R_birth_cycles.len = count;

          self->g_temp_R_birth_cycles.original = 1 - self->g_temp_R_birth_cycles.original;
          
          
    }

    //// Reduce V
    reduce_temp_V_H0(self);

}



void find_V_recursively_edges(filtration* self, EDGE_ID bo_idx, EDGE_ID bo_pivot){

#ifdef RECORD_V_USAGE
      self->g_H0_pivot_of[bo_pivot].V_usage++;
#endif

#ifdef ADAPTIVE_V_STORAGE

      if (self->g_H0_pivot_of[bo_pivot].V_len){

          //printf("\nUsing stored cycle for %d", bo_idx);
          
          if ((self->g_temp_V_primary.len + self->g_H0_pivot_of[bo_pivot].V_len) > self->g_temp_V_primary.max_len - 1){
                
                self->g_temp_V_primary.max_len = self->g_temp_V_primary.len + self->g_H0_pivot_of[bo_pivot].V_len + 1000;
                self->g_temp_V_primary.VV = (EDGE_ID*)realloc(self->g_temp_V_primary.VV\
                                                            , self->g_temp_V_primary.max_len*sizeof(EDGE_ID));

          }

          for (EDGE_ID mm = 0; mm < self->g_H0_pivot_of[bo_pivot].V_len; mm++){
                
                self->g_temp_V_primary.VV[self->g_temp_V_primary.len++] =\
                                           self->g_H0_pivot_of[bo_pivot].VV[mm];
          }
          return;
          
      }
#endif

      EDGE_ID res_max_len = 100;
      EDGE_ID* result = (EDGE_ID*)malloc((2*res_max_len)*sizeof(EDGE_ID));
      int res_original = 0;
      EDGE_ID* original_result;
      EDGE_ID* scratch_result;

      EDGE_ID res_len = 2;

      result[0] = self->g_edges_list[2*bo_idx];
      result[1] = self->g_edges_list[2*bo_idx+1];

      EDGE_ID possible_len;

      EDGE_ID* red_start;
      EDGE_ID red_simp_len;
      EDGE_ID ro;
      EDGE_ID j, k, count;
      EDGE_ID pivot;

      // RECORD THIS IN REDUCTION OPERATION for e_o
      self->g_temp_V_primary.VV[self->g_temp_V_primary.len++] = bo_idx;
      
      // Check for overflow
      if (self->g_temp_V_primary.len == self->g_temp_V_primary.max_len){

           //printf("\nReallocating");
           //getchar();
           self->g_temp_V_primary.max_len += 100;
           self->g_temp_V_primary.VV = (EDGE_ID*)realloc(self->g_temp_V_primary.VV\
                                                   , self->g_temp_V_primary.max_len*sizeof(EDGE_ID));

      }

      while(res_len != 0){

#ifdef ADAPTIVE_V_STORAGE
          self->g_depth++;
#endif

          original_result = result + (res_original*res_max_len);
          scratch_result = result + ((1-res_original)*res_max_len);

          pivot = original_result[res_len-1];

          // The new pivot is pivot of R(bo_idx)
          bo_idx = self->g_H0_pivot_of[pivot].coface;

          ro = (EDGE_ID)self->g_pivots_H0[pivot];


          red_simp_len = self->g_R_col_indices_H0[ro+1] - \
                           self->g_R_col_indices_H0[ro];


          red_start = self->g_R_sparse_H0 + self->g_R_col_indices_H0[ro];

          // Check for overflow
          possible_len = res_len + red_simp_len;
          if (possible_len > res_max_len - 1){

                if (res_original){

                      for (k = 0; k < res_len; k++){
                          result[k] = result[k + res_max_len];
                      }
                      
                }

                res_max_len = possible_len + 1000;
                result = (EDGE_ID*)realloc(result, (2*res_max_len)*sizeof(EDGE_ID));

                original_result = result;
                scratch_result = result + res_max_len;
                res_original = 0;
                
          }


          // Reduce
          j = 0;
          k = 0;
          count = 0;

          while ((j < res_len) && (k < red_simp_len)){

              if (original_result[j] < red_start[k]){
                   
                  scratch_result[count] = original_result[j];
                  count = count + 1;
                  j = j + 1;
                   
              }
              else if (original_result[j] > red_start[k]){

                  scratch_result[count] = red_start[k];
                  count = count + 1;
                  k = k + 1;

              }
              else{

                  j = j + 1;
                  k = k + 1;

              }

          }

          while (j < res_len){

               scratch_result[count++] = original_result[j++];

          }

          while (k < red_simp_len){

               scratch_result[count++] = red_start[k++];

          }

          res_len = count;

          res_original = 1 - res_original;

          if (res_len != 0){
                
                    find_V_recursively_edges(self, bo_idx, pivot);

          }

      }


#ifndef RECORD_V_USAGE
      self->g_H0_pivot_of[bo_pivot].V_usage++;
#endif

#ifdef ADAPTIVE_V_STORAGE


      if ((self->g_H0_pivot_of[bo_pivot].V_usage > self->g_cycle_usage_thresh)\
           &&(!self->g_H0_pivot_of[bo_pivot].V_stored)){
          self->g_store_V_for[self->g_store_V_for_len++] = bo_pivot;
          if (self->g_store_V_for_len == self->g_store_V_for_max_len){
                self->g_store_V_for_max_len += 100;
                self->g_store_V_for = (EDGE_ID*)realloc(self->g_store_V_for\
                                                    , self->g_store_V_for_max_len*sizeof(EDGE_ID));
          }
      }


#endif

      free(result);
  
}



void compute_H2_homology_cycles(filtration* self){

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
//            STEP H2.1: Find homology now for the tetrahedrons (H2)
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


     if (!self->g_suppress_output){
          printf("\n\n---------------");
          printf("\nComputing H2...");
          printf("\n---------------\n");
     }


     struct timespec start_wall_clock, finish_wall_clock;

     clock_gettime(CLOCK_MONOTONIC, &(self->g_start_wall_clock));

     self->g_R_max_len_H2 = 100;
     self->g_R_H2 = (simplex*)malloc(self->g_R_max_len_H2*sizeof(simplex));
     self->g_R_len_H2 = 0;


     self->g_R_col_idx_max_len_H2 = 100;
     self->g_R_col_idx_H2 = (EDGE_ID*)malloc(self->g_R_col_idx_max_len_H2*sizeof(EDGE_ID));

     self->g_R_col_idx_H2_ptr = 0;

     self->g_workspace_size = 1000;

     self->g_ws_pre_alloc = 1000;
     // Initialize ws counter
     self->g_ws_counter = 0;

     // H1 workspace structures
     self->g_workspace_H2 = (simplex**)malloc(self->g_workspace_size*sizeof(simplex*));

     // H1 workspace info
     self->g_workspace_H2_info = (boundary_H2_ws*)malloc(self->g_workspace_size*sizeof(boundary_H2_ws));


     for (int i = 0; i < self->g_workspace_size; i++){

         self->g_workspace_H2_info[i].max_len = self->g_ws_pre_alloc;

         self->g_workspace_H2[i] = (simplex*)malloc(2*self->g_workspace_H2_info[i].max_len*sizeof(simplex));

         self->g_workspace_H2_info[i].trivial_boundary = (simplex*)malloc(4*sizeof(simplex));

     }

     // Pivots
     self->g_H2_pivots = (H2_pivots**)malloc(self->g_n_valid_edges*sizeof(H2_pivots*));
     self->g_H2_pivots_len = (EDGE_ID*)calloc(self->g_n_valid_edges, sizeof(EDGE_ID));
     self->g_H2_pivots_max_len = (EDGE_ID*)calloc(self->g_n_valid_edges, sizeof(EDGE_ID));

     for (EDGE_ID mm = 0; mm < self->g_n_valid_edges; mm++){
          self->g_H2_pivots_max_len[mm] = 5;
          self->g_H2_pivots[mm] = (H2_pivots*)malloc(self->g_H2_pivots_max_len[mm]*sizeof(H2_pivots));
     }

     // Convenient info for pers pairs
     self->g_homH2_pers_len = 0;
     self->g_homH2_pers_max_len = 100;
     self->g_homH2_pers = (homH2_pers*)malloc(self->g_homH2_pers_max_len*sizeof(homH2_pers));

     // Temporary space for birth cycles
     self->g_temp_V_H2_primary.max_len = 10;
     self->g_temp_V_H2_primary.VV = (simplex*)malloc(self->g_temp_V_H2_primary.max_len*sizeof(simplex));
     self->g_temp_V_H2_primary.len = 0;

     // Temporary space for birth cycles
     self->g_temp_R_H2_birth_cycles.max_len = 100;
     self->g_temp_R_H2_birth_cycles.RR = (EDGE_ID*)malloc(2*self->g_temp_R_H2_birth_cycles.max_len*sizeof(EDGE_ID));
     self->g_temp_R_H2_birth_cycles.original = 0;
     self->g_temp_R_H2_birth_cycles.len = 0;

#ifdef ADAPTIVE_V_STORAGE
     // Create pointers to store V

     for (EDGE_ID mm = 0; mm < self->g_n_valid_edges; mm++){

        self->g_H1_pivot_of[mm].V_usage = 0;
        self->g_H1_pivot_of[mm].V_stored = 0;
        self->g_H1_pivot_of[mm].V_len = 0;
        self->g_H1_pivot_of[mm].VV = NULL;

     }

     // Create pointers to store V per extraction call
     self->g_store_V_for_len = 0;
     self->g_store_V_for_max_len = 10;
     self->g_store_V_for = (EDGE_ID*)malloc(self->g_store_V_for_max_len*sizeof(EDGE_ID));


#endif

//#ifdef MINIMIZE_BIRTH_CYCLES
     if (self->g_reduce_cyc_lengths){

        self->g_all_V_stored_num = 0;
        self->g_all_V_stored_max_num = 10;

        self->g_all_V_H1_stored = (cyc_info_H2*)malloc(self->g_all_V_stored_max_num*sizeof(cyc_info_H2));

     }

//#endif

     ////////////////////////////////////////////////////////////////
     //
     //   Allocate jobs for parallel H1
     // 
     ////////////////////////////////////////////////////////////////

     self->g_jobs = (int*)malloc((self->g_cpu_count + 1)*sizeof(int));

     allocate_jobs(self, self->g_workspace_size);

     self->g_threads = (pthread_t *)malloc(self->g_cpu_count*sizeof(pthread_t));

     int rtn;

     if ((rtn = pthread_mutex_init(&(self->g_thread_lock), NULL)) !=0)
        fprintf(stderr, "pthread_mutex_init %s", strerror(rtn)), exit(-1);

     if ((rtn = pthread_cond_init(&(self->g_start_boss), NULL)) !=0)
        fprintf(stderr, "pthread_cond_init %s", strerror(rtn)), exit(-1);

     if ((rtn = pthread_cond_init(&(self->g_start_workers), NULL)) !=0)
        fprintf(stderr, "pthread_cond_init %s", strerror(rtn)), exit(-1);


     // Initialize thread creation
     self->g_thread_id = 0;
     self->g_sleeping_threads = 0;
     self->g_delete_threads = 0;

     for (int i = 0; i < self->g_cpu_count; i++){

        if ((rtn = pthread_create( \
                                &(self->g_threads[i]) \
                                , NULL \
                                , reduce_with_complex_H2 \
                                , (void*)self)!= 0))
          fprintf(stderr, "pthread_create %d", rtn), exit(-1);
      
     }

     // Wait for threads to be initialized
     pthread_mutex_lock(&(self->g_thread_lock));

     while(self->g_sleeping_threads != self->g_cpu_count){
        
          pthread_cond_wait(&(self->g_start_boss) \
                          , &(self->g_thread_lock));

     }

     ////////////////////////////////


     H2_preprocess* temp_tetra = (H2_preprocess*)malloc(self->g_n_valid_edges*sizeof(H2_preprocess));
     EDGE_ID temp_tetra_len = 0;

     int ccounter = 0;
     self->g_new_debug2 = 0;


     EDGE_ID o_cd;
    
     self->g_ws_counter = 0;

     for (EDGE_ID o_ab = 0; o_ab < self->g_n_valid_edges; o_ab++){

        if (!self->g_H2_cohom_pivots_len[o_ab]){
            continue;
        }


        for (VERT_ID mm = 0; mm < self->g_H2_cohom_pivots_len[o_ab]; mm++){

            o_cd = self->g_H2_cohom_pivots[o_ab][mm].key2;

            // Workspace attributes
            boundary_H2_ws* this_ws = self->g_workspace_H2_info + self->g_ws_counter;

            // Initially, the original is at 0
            this_ws->original = 0;

            this_ws->flag_first = 1;

            // Parallel control flags
            this_ws->flag_empty = 0;
            this_ws->flag_red_w_complex = 0;
            this_ws->flag_append_to_complex = 1;


            this_ws->tetrahedron.key1 = o_ab;
            this_ws->tetrahedron.key2 = o_cd;

            // Initial length of boundary
            this_ws->len = 4;

            compute_boundary_tetra(self, this_ws->tetrahedron, self->g_workspace_H2[self->g_ws_counter]);

            this_ws->pivot = self->g_workspace_H2[self->g_ws_counter][3];

            self->g_ws_counter++;


            if (self->g_ws_counter == self->g_workspace_size){
                  
                  reduce_ws_H2(self);
            }


        }


     }


     // Reduction of final batch
     while (self->g_ws_counter){
          
          allocate_jobs(self, self->g_ws_counter);
          reduce_ws_H2(self);

     }


     clock_gettime(CLOCK_MONOTONIC, &finish_wall_clock);
     self->g_timer_computeH2 = (finish_wall_clock.tv_sec - start_wall_clock.tv_sec);
     self->g_timer_computeH2 += (finish_wall_clock.tv_nsec - start_wall_clock.tv_nsec) / 1000000000.0;

     if (!self->g_suppress_output){
        printf("\nComputed H2.");
     }

     ////////////////////////
     // HOMOLOGY CYCLES
     ////////////////////////

     clock_gettime(CLOCK_MONOTONIC, &start_wall_clock);
     
     FILE* fp2 = fopen(self->g_homH2_cycles_file, "w");

     PAR birth, death;

     self->g_n_H1_stored_V = 0;

     // Go over the pers pairs of features that died and compute the cycles
     for (EDGE_ID mm = 0; mm < self->g_homH2_pers_len; mm++){

          self->g_n_H2_birth_cycles++;

          //printf("\nFinding cycle for (%d, %d)", self->g_homH2_pers[mm].birth_simplex.key1\
          //                                     , self->g_homH2_pers[mm].birth_simplex.key2);
          //getchar();
          
          if (self->g_filetype == 1){
              birth = sqrt(self->g_edge_parameter[self->g_homH2_pers[mm].birth_simplex.key1]);
              death = sqrt(self->g_edge_parameter[self->g_homH2_pers[mm].death_edge]);

          }
          else{
              birth = self->g_edge_parameter[self->g_homH2_pers[mm].birth_simplex.key1];
              death = self->g_edge_parameter[self->g_homH2_pers[mm].death_edge];
          }
        
          fprintf(fp2, "%lf, %lf", birth, death);

          fprintf(fp2, "\nhomology cycle");

          for (EDGE_ID bb = self->g_R_col_idx_H2[self->g_homH2_pers[mm].R_col_idx]\
              ; bb < self->g_R_col_idx_H2[self->g_homH2_pers[mm].R_col_idx + 1]\
              ; bb++){
                
                fprintf(fp2, ", %d, %d, %d", self->g_edges_list[2*self->g_R_H2[bb].key1]\
                                           , self->g_edges_list[2*self->g_R_H2[bb].key1+1]\
                                           , self->g_R_H2[bb].key2\
                                           );
          }

            
          
          // Always write the birth cycles to file for now
          fprintf(fp2, "\nbirth cycle");

          //printf("\nGetting birth cycle %d out of %d", mm, self->g_homH1_pers_len);
          //getchar();

#ifdef ADAPTIVE_V_STORAGE
          self->g_store_V_for_len = 0;
#endif

          //printf("\nGetting void");
          get_birth_void(self, self->g_homH2_pers[mm].birth_simplex);


//#ifdef MINIMIZE_BIRTH_CYCLES
          if (self->g_reduce_cyc_lengths){

              //self->g_all_V_stored_len[self->g_all_V_stored_num] = self->g_temp_V_H2_primary.len;
              //self->g_all_V_H1_stored[self->g_all_V_stored_num] =\
              //                                   (simplex*)malloc(self->g_temp_V_H2_primary.len*sizeof(simplex));
              
              self->g_all_V_H1_stored[self->g_all_V_stored_num].boundary =\
                                                     (simplex*)malloc(self->g_temp_V_H2_primary.len*sizeof(simplex));

              self->g_all_V_H1_stored[self->g_all_V_stored_num].len = self->g_temp_V_H2_primary.len;


              self->g_all_V_H1_stored[self->g_all_V_stored_num].perspair[0] = birth;
              self->g_all_V_H1_stored[self->g_all_V_stored_num].perspair[1] = death;

              self->g_all_V_H1_stored[self->g_all_V_stored_num].updated_birth = birth;
          }

//#endif


          for (EDGE_ID nn = 0; nn < self->g_temp_V_H2_primary.len; nn++){
               
                 fprintf(fp2, ", %d, %d, %d", self->g_edges_list[2*self->g_temp_V_H2_primary.VV[nn].key1]\
                                            , self->g_edges_list[2*self->g_temp_V_H2_primary.VV[nn].key1+1]\
                                            , self->g_temp_V_H2_primary.VV[nn].key2\
                                            );

//#ifdef MINIMIZE_BIRTH_CYCLES
                 if (self->g_reduce_cyc_lengths){
                      self->g_all_V_H1_stored[self->g_all_V_stored_num].boundary[nn] = self->g_temp_V_H2_primary.VV[nn];
                 }
//#endif

               
          }

          fprintf(fp2, "\n");


//#ifdef MINIMIZE_BIRTH_CYCLES

          if (self->g_reduce_cyc_lengths){

              self->g_all_V_stored_num++;
              if (self->g_all_V_stored_num == self->g_all_V_stored_max_num){
                    self->g_all_V_stored_max_num += 100;
                    self->g_all_V_H1_stored = (cyc_info_H2*)realloc(self->g_all_V_H1_stored\
                                                            , self->g_all_V_stored_max_num*sizeof(cyc_info_H2));

                    //self->g_all_V_stored_len = (EDGE_ID*)realloc(self->g_all_V_stored_len\
                    //                                        , self->g_all_V_stored_max_num*sizeof(EDGE_ID));
              }

          }

//#endif





#ifdef ADAPTIVE_V_STORAGE

          //printf("\nPress key to store V at line 10412");
          //getchar();
          store_V_H1(self);
          //printf("\nstored V at line 10412. Pres key to continue.");
          //getchar();

#endif
          
     }
     
     // Go over the pers pairs of undead features and compute the birth cycles
     for (EDGE_ID mm = 0; mm < self->g_H2_undead_ptr; mm++){

          self->g_n_H2_birth_cycles++;

          if (self->g_filetype == 1){
              birth = sqrt(self->g_edge_parameter[self->g_H2_undead[mm].key1]);

          }
          else{
              birth = self->g_edge_parameter[self->g_H2_undead[mm].key1];
          }

          fprintf(fp2, "%lf, -1", birth);
          
#ifdef ADAPTIVE_V_STORAGE
          self->g_store_V_for_len = 0;
#endif

          get_birth_void(self, self->g_H2_undead[mm]);

//#ifdef MINIMIZE_BIRTH_CYCLES
          if (self->g_reduce_cyc_lengths){

              //self->g_all_V_stored_len[self->g_all_V_stored_num] = self->g_temp_V_H2_primary.len;
              self->g_all_V_H1_stored[self->g_all_V_stored_num].boundary =\
                                                 (simplex*)malloc(self->g_temp_V_H2_primary.len*sizeof(simplex));

              self->g_all_V_H1_stored[self->g_all_V_stored_num].len = self->g_temp_V_H2_primary.len;

              self->g_all_V_H1_stored[self->g_all_V_stored_num].perspair[0] = birth;
              self->g_all_V_H1_stored[self->g_all_V_stored_num].perspair[1] = -1;

              self->g_all_V_H1_stored[self->g_all_V_stored_num].updated_birth = birth;

          }

//#endif
          
          fprintf(fp2, "\nbirth cycle");

          for (EDGE_ID nn = 0; nn < self->g_temp_V_H2_primary.len; nn++){
               
                 fprintf(fp2, ", %d, %d, %d", self->g_edges_list[2*self->g_temp_V_H2_primary.VV[nn].key1]\
                                            , self->g_edges_list[2*self->g_temp_V_H2_primary.VV[nn].key1+1]\
                                            , self->g_temp_V_H2_primary.VV[nn].key2\
                                            );

//#ifdef MINIMIZE_BIRTH_CYCLES
                 if (self->g_reduce_cyc_lengths){

                    self->g_all_V_H1_stored[self->g_all_V_stored_num].boundary[nn] = self->g_temp_V_H2_primary.VV[nn];

                 }
//#endif
               
          }

          fprintf(fp2, "\n");


//#ifdef MINIMIZE_BIRTH_CYCLES

          if (self->g_reduce_cyc_lengths){

              self->g_all_V_stored_num++;
              if (self->g_all_V_stored_num == self->g_all_V_stored_max_num){
                    self->g_all_V_stored_max_num += 100;
                    //self->g_all_V_H1_stored = (simplex**)realloc(self->g_all_V_H1_stored\
                    //                                        , self->g_all_V_stored_max_num*sizeof(simplex*));
                    self->g_all_V_H1_stored = (cyc_info_H2*)realloc(self->g_all_V_H1_stored\
                                                            , self->g_all_V_stored_max_num*sizeof(cyc_info_H2));
                    //self->g_all_V_stored_len = (EDGE_ID*)realloc(self->g_all_V_stored_len\
                    //                                        , self->g_all_V_stored_max_num*sizeof(EDGE_ID));
              }

          }

//#endif




#ifdef ADAPTIVE_V_STORAGE

          store_V_H1(self);

#endif

     }


     if (!self->g_suppress_output){
        printf("\nComputed birth cycles.");
     }
     //getchar();



     /////////////////////////
     // Cancel the threads used in getting next during reduction
     /////////////////////////

     self->g_delete_threads = 1;

     pthread_cond_broadcast(&(self->g_start_workers));

     pthread_mutex_unlock(&(self->g_thread_lock));

     for (int i = 0; i < self->g_cpu_count; i++){

        pthread_join(self->g_threads[i], NULL);
      
     }

     free(self->g_jobs);
     free(self->g_threads);


     for (int i = 0; i < self->g_workspace_size; i++){

         free(self->g_workspace_H2_info[i].trivial_boundary);
         free(self->g_workspace_H2[i]);

     }


     free(self->g_workspace_H2);
     free(self->g_workspace_H2_info);


     free(self->g_R_H2);
     free(self->g_R_col_idx_H2);

     free(self->g_H2_pivots);
     free(self->g_homH2_pers);
     free(self->g_temp_V_H2_primary.VV);
     free(self->g_temp_R_H2_birth_cycles.RR);

     free(self->g_R_H1);
     free(self->g_R_col_idx_H1);
     free(self->g_pivots_H1);


#ifdef ADAPTIVE_V_STORAGE
     free(self->g_store_V_for);

#ifdef RECORD_V_USAGE
     FILE* fp3 = fopen(self->g_V_H1_usage_file, "w");
#endif

     for (EDGE_ID mm = 0; mm < self->g_n_valid_edges; mm++){
          if (self->g_H1_pivot_of[mm].V_len){

#ifdef RECORD_V_USAGE
              fprintf(fp3, "%d, %d\n"\
                                  , self->g_H1_pivot_of[mm].V_usage\
                                  , self->g_H1_pivot_of[mm].V_depth);
#endif
              free(self->g_H1_pivot_of[mm].VV);
          }
     }

#ifdef RECORD_V_USAGE
     fclose(fp3);
#endif


#endif


 
//#ifdef HOM_CYCLES
     free(self->g_H1_pivot_of);

     clock_gettime(CLOCK_MONOTONIC, &finish_wall_clock);
     self->g_timer_H2cycles = (finish_wall_clock.tv_sec - start_wall_clock.tv_sec);
     self->g_timer_H2cycles += (finish_wall_clock.tv_nsec - start_wall_clock.tv_nsec) / 1000000000.0;

//#endif

     if (!self->g_suppress_output){
        printf("\nQUITTING H2 cycle computation");
     }
     //getchar();

//#ifdef MINIMIZE_BIRTH_CYCLES

     if (self->g_reduce_cyc_lengths){

          if (!self->g_suppress_output){
              printf("\nMinimizing birth cycles...");
          }
          //getchar();

          clock_gettime(CLOCK_MONOTONIC, &start_wall_clock);

          minimize_birth_cycles_H1_v2(self\
                                 , self->g_all_V_H1_stored\
                                 , self->g_all_V_stored_num\
                                 , self->g_minimal_V_H1_file\
                                 , self->g_V_H1_birthcyc_lens_file\
                                 , self->g_minimal_V_H1_birthcyc_lens_file\
                                 );

          clock_gettime(CLOCK_MONOTONIC, &finish_wall_clock);
          self->g_timer_minimize_H2cycles = (finish_wall_clock.tv_sec - start_wall_clock.tv_sec);
          self->g_timer_minimize_H2cycles += (finish_wall_clock.tv_nsec - start_wall_clock.tv_nsec) / 1000000000.0;

     }

//#endif



}

void reduce_ws_H2(filtration* self){


      //if (self->g_new_debug2){
            //printf("\nBEFORE parallel.");
            //for (int kk = 0; kk < self->g_ws_counter; kk++){
            //      
            //      printf("\n%d has triangle (%d, %d) with pivot (%d, %d) and reduce_with %p", kk\
            //                                                      , self->g_workspace_H2_info[kk].tetrahedron.key1\
            //                                                      , self->g_workspace_H2_info[kk].tetrahedron.key2\
            //                                                      , self->g_workspace_H2_info[kk].pivot.key1\
            //                                                      , self->g_workspace_H2_info[kk].pivot.key2\
            //                                                      , self->g_workspace_H2_info[kk].reduce_with\
            //                                                      );
            //}
            //getchar();
      //}
      
      self->g_processed_threads = 0;

      //printf("\npress key to reduce with complex");
      //getchar();

      pthread_cond_broadcast(&(self->g_start_workers));

      while (self->g_processed_threads != self->g_cpu_count){
            
            pthread_cond_wait(&(self->g_start_boss) \
                            ,&(self->g_thread_lock));
      }

      //printf("\npress key to reduce with self");
      //getchar();


      //if (self->g_new_debug2){
      //      for (int kk = 0; kk < self->g_ws_counter; kk++){
      //            
      //            printf("\n%d has triangle (%d, %d) with pivot %d", kk\
      //                                                            , self->g_workspace_H1_info[kk].triangle.key1\
      //                                                            , self->g_workspace_H1_info[kk].triangle.key2\
      //                                                            , self->g_workspace_H1_info[kk].pivot\
      //                                                            );
      //      }
      //      printf("\nafter parallel. press key to start serial");
      //      //getchar();
      //}

      reduce_with_self_H2( \
                            self \
                            );

            //printf("\nAFTER SERIAL parallel.");
            //for (int kk = 0; kk < self->g_ws_counter; kk++){
            //      
            //      printf("\n%d has triangle (%d, %d) with pivot (%d, %d) and reduce_with %p", kk\
            //                                                      , self->g_workspace_H2_info[kk].tetrahedron.key1\
            //                                                      , self->g_workspace_H2_info[kk].tetrahedron.key2\
            //                                                      , self->g_workspace_H2_info[kk].pivot.key1\
            //                                                      , self->g_workspace_H2_info[kk].pivot.key2\
            //                                                      , self->g_workspace_H2_info[kk].reduce_with\
            //                                                      );
            //}
      //if (self->g_new_debug2){
      //      for (int kk = 0; kk < self->g_ws_counter; kk++){
      //            
      //            printf("\n%d has triangle (%d, %d) with pivot %d", kk\
      //                                                            , self->g_workspace_H1_info[kk].triangle.key1\
      //                                                            , self->g_workspace_H1_info[kk].triangle.key2\
      //                                                            , self->g_workspace_H1_info[kk].pivot\
      //                                                            );
      //      }
      //      printf("\nafter serial. press key to update ");
      //      //getchar();
      //}
      
      int count_valid = 0;

      for (int ws_counter=0; ws_counter < self->g_ws_counter; ws_counter++){

            if (self->g_workspace_H2_info[ws_counter].flag_append_to_complex){

                  update_R_H2(self \
                                , ws_counter\
                                );
                  continue;
                    
            }
            
            //if (!self->g_workspace_H1_info[ws_counter].len){continue;}
            if (self->g_workspace_H2_info[ws_counter].flag_empty){

                  continue;
            }


            // Swap R
            simplex* temp = self->g_workspace_H2[count_valid];
            self->g_workspace_H2[count_valid] = self->g_workspace_H2[ws_counter];
            self->g_workspace_H2[ws_counter] = temp;

            // Swap R info
            boundary_H2_ws temp2 = self->g_workspace_H2_info[count_valid];
            self->g_workspace_H2_info[count_valid] = self->g_workspace_H2_info[ws_counter];
            self->g_workspace_H2_info[ws_counter] = temp2;


            // At this point, this has to be a non-zero column
            self->g_workspace_H2_info[count_valid].flag_empty = 0;

            count_valid += 1;

      }

      self->g_ws_counter = count_valid;

            //printf("\nAFTER UPDATE .");
            //for (int kk = 0; kk < self->g_ws_counter; kk++){
            //      
            //      printf("\n%d has triangle (%d, %d) with pivot (%d, %d)", kk\
            //                                                      , self->g_workspace_H2_info[kk].tetrahedron.key1\
            //                                                      , self->g_workspace_H2_info[kk].tetrahedron.key2\
            //                                                      , self->g_workspace_H2_info[kk].pivot.key1\
            //                                                      , self->g_workspace_H2_info[kk].pivot.key2\
            //                                                      );
            //}
      //if (self->g_new_debug2){
      //      for (int kk = 0; kk < self->g_ws_counter; kk++){
      //            
      //            printf("\n%d has triangle (%d, %d) with pivot %d", kk\
      //                                                            , self->g_workspace_H1_info[kk].triangle.key1\
      //                                                            , self->g_workspace_H1_info[kk].triangle.key2\
      //                                                            , self->g_workspace_H1_info[kk].pivot\
      //                                                            );
      //      }
      //      printf("\nafter update. press key to continue ");
      //      //getchar();
      //}
      
      //if (dim)
      //  self->g_H0_MAX = self->g_n_reduced_simplex[dim];

}



void reduce_with_self_H2( \
                      filtration* self \
                      ){

    int compare;

    int i, m;
    int idx;

    EDGE_ID count, j, k;


    for (i=0; i < self->g_ws_counter; i++){

        boundary_H2_ws* this_ws = self->g_workspace_H2_info + i;

        this_ws->flag_reduce = 0;

        // If the simplex has already been reduced to 0
        // then continue
        if (this_ws->flag_empty){ 

          this_ws->flag_append_to_complex = 0;
          continue;

        }

        
        simplex* orig = self->g_workspace_H2[i] + this_ws->original*this_ws->max_len;

        m = 0;
        while (m < i){

            boundary_H2_ws* m_ws = self->g_workspace_H2_info + m;

            if (m_ws->flag_empty){
                m++;
                continue;
            }


            simplex* original_m = self->g_workspace_H2[m] + m_ws->original*m_ws->max_len;

            orig = self->g_workspace_H2[i] + this_ws->original*this_ws->max_len;

            int compare;

            if (m_ws->pivot.key1 < this_ws->pivot.key1) compare = 1;
            else if (m_ws->pivot.key1 > this_ws->pivot.key1) compare = 0;
            else{

                if (m_ws->pivot.key2 < this_ws->pivot.key2) compare = 1;
                else if (m_ws->pivot.key2 > this_ws->pivot.key2) compare = 0;
                else compare = -1;
            }

            //if (m_ws->pivot > this_ws->pivot){
            if (compare == 0){
                  
                  if (m_ws->flag_red_w_complex){
                        
                        this_ws->flag_append_to_complex = 0;
                        break;
                  }
                  m++;
                  continue;
                  
            }


            //if (m_ws->pivot < this_ws->pivot){
            if (compare == 1){
                    m++;
                    continue;
            }

            if (m_ws->flag_red_w_complex){
                    this_ws->flag_append_to_complex = 0;
                    //m++;
                    //continue;
                    break;
            }


            if (this_ws->len + m_ws->len > this_ws->max_len ){

                //printf("\nReallocating inside self");
                //simp_max_len_i = len_i + len_m + 1000;


                if (this_ws->original){
                      
                      for (EDGE_ID mm = 0; mm < this_ws->len; mm++){
                            self->g_workspace_H2[i][mm] = self->g_workspace_H2[i][mm + this_ws->max_len];
                      }

                      this_ws->original = 0;

                }

                this_ws->max_len = this_ws->len + m_ws->len + 1000;

                self->g_workspace_H2[i] = (simplex*)realloc(self->g_workspace_H2[i]\
                                                      , 2*this_ws->max_len*sizeof(simplex));

                orig  = self->g_workspace_H2[i];


            }

            simplex* scratch = self->g_workspace_H2[i] + (1-this_ws->original)*this_ws->max_len;

            // Store the result in scratch

            count = 0;

            j = 0;
            k = 0;

            while ((j < this_ws->len) && (k < m_ws->len)){

                
                if (orig[j].key1 < original_m[k].key1) compare = 1;
                else if (orig[j].key1 > original_m[k].key1) compare = 0;
                else{

                    if (orig[j].key2 < original_m[k].key2) compare = 1;
                    else if (orig[j].key2 > original_m[k].key2) compare = 0;
                    else compare = -1;
                }
                

                //if (orig[j] < original_m[k]){
                if (compare == 1){
                    scratch[count++] = orig[j++];
                }
                //else if (orig[j] > original_m[k]){
                else if (compare == 0){
                    scratch[count++] = original_m[k++];
                }
                else{
                    j++;
                    k++;
                }

            }

            while (j < this_ws->len){
                    scratch[count++] = orig[j++];
            }

            while (k < m_ws->len){
                    scratch[count++] = original_m[k++];
            }

            this_ws->len = count;

            this_ws->original = 1 - this_ws->original;

            if (!count){
                  
                  this_ws->flag_append_to_complex = 0;
                  this_ws->flag_empty = 1;
                  break;
                  
            }

            this_ws->pivot = scratch[this_ws->len-1];


            coboundary_H2 temp;

            temp.triangle.key1 = this_ws->pivot.key1;
            temp.triangle.key2 = this_ws->pivot.key2;

            find_H2_cohom_low(self, &temp);


            // Check if pivot is trivial or if it is pivot in H2
            // To check trivial, check if maximum triangle in low of cob of triangle is triangle
            if ((temp.low.key1 == temp.triangle.key1)\
              &&(self->g_edges_list[2*temp.low.key2+1] == temp.triangle.key2)){

                // In which case, reduce with boundary of temp.low
                
                
                compute_boundary_tetra(self, temp.low, this_ws->trivial_boundary);

                this_ws->R_col_idx = 0;

                this_ws->reduce_with_len = 4;
                this_ws->flag_reduce = 1;


                this_ws->flag_red_w_complex = 1;
                this_ws->flag_append_to_complex = 0;
                break;



            }
            else{

                if (self->g_H2_pivots_len[this_ws->pivot.key1]){

                    EDGE_ID idx = search_H2_pivots(self->g_H2_pivots[this_ws->pivot.key1]\
                                              , 0\
                                              , self->g_H2_pivots_len[this_ws->pivot.key1] - 1\
                                              , this_ws->pivot.key2 \
                                              , self->g_n_valid_edges\
                                              );

                    if (idx != self->g_n_valid_edges){
                            //this_ws->flag_red_w_trivial = 0;
                            
                            this_ws->R_col_idx = self->g_H2_pivots[this_ws->pivot.key1][idx].col_idx;

                            this_ws->reduce_with_len = self->g_R_col_idx_H2[this_ws->R_col_idx+1] -\
                                                        self->g_R_col_idx_H2[this_ws->R_col_idx];

                            this_ws->flag_reduce = 1;

                            
                            this_ws->flag_red_w_complex = 1;
                            this_ws->flag_append_to_complex = 0;
                            break;
                    }

                }


            }
            



            //}

            m = 0;

        }//End of m loop

    }

}//End of red_ws_w_self_single



void* reduce_with_complex_H2(void* arg){
      
      filtration* self = arg;

      pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, 0);

      pthread_mutex_lock(&(self->g_thread_lock));

      int tid = ++self->g_thread_id;

      simplex *red_start;
      EDGE_ID j, k ,count;


      for (;;){

          self->g_sleeping_threads++;
          
          if (self->g_sleeping_threads == self->g_cpu_count)
              pthread_cond_signal(&(self->g_start_boss));

          pthread_cond_wait(&(self->g_start_workers), &(self->g_thread_lock));

          if (self->g_delete_threads){
            //printf("\nexiting from thread %d", tid);
            pthread_mutex_unlock(&(self->g_thread_lock));
            pthread_exit(NULL);
          }

          self->g_sleeping_threads--;

          pthread_mutex_unlock(&(self->g_thread_lock));

          for (int i = self->g_jobs[tid - 1]; i < self->g_jobs[tid]; i++){

              boundary_H2_ws* this_ws = self->g_workspace_H2_info + i;


              if (this_ws->flag_empty){
                  // We are sure that we will exit only if there is no reduction
                  // required with existing complex or with trivial pair
                  this_ws->flag_red_w_complex = 0;
                  this_ws->flag_append_to_complex = 0;
                  continue;
              }



              if (this_ws->flag_first){

                  //if ((this_ws->tetrahedron.key1 == self->g_debug_tetra.key1)\
                  //  &&(this_ws->tetrahedron.key2 == self->g_debug_tetra.key2))
                  //{
                  //    printf("\nProcessing tetra: (%d, %d) first pivot is (%d, %d)"\
                  //                                   ,self->g_debug_tetra.key1\
                  //                                   ,self->g_debug_tetra.key2\
                  //                                   ,this_ws->pivot.key1\
                  //                                   ,this_ws->pivot.key2\
                  //                                   );
                  //    
                  //}

                  this_ws->flag_first = 0;

                  coboundary_H2 temp;

                  temp.triangle.key1 = this_ws->pivot.key1;
                  temp.triangle.key2 = this_ws->pivot.key2;

                  find_H2_cohom_low(self, &temp);

                  this_ws->flag_reduce = 0;


                  // Check if pivot is trivial or if it is pivot in H2
                  // To check trivial, check if maximum triangle in low of cob of triangle is triangle
                  if ((temp.low.key1 == temp.triangle.key1)\
                    &&(self->g_edges_list[2*temp.low.key2+1] == temp.triangle.key2)){

                      // In which case, reduce with boundary of temp.low
                      
                      
                      compute_boundary_tetra(self, temp.low, this_ws->trivial_boundary);

                      this_ws->R_col_idx = 0;
                      this_ws->reduce_with_len = 4;
                      this_ws->flag_reduce = 1;


                  }
                  else{
                        
                      
                      if (self->g_H2_pivots_len[this_ws->pivot.key1]){

                          EDGE_ID idx = search_H2_pivots(self->g_H2_pivots[this_ws->pivot.key1]\
                                                    , 0\
                                                    , self->g_H2_pivots_len[this_ws->pivot.key1] - 1\
                                                    , this_ws->pivot.key2 \
                                                    , self->g_n_valid_edges\
                                                    );

                          if (idx != self->g_n_valid_edges){
                                  //this_ws->flag_red_w_trivial = 0;
                                  
                                  this_ws->R_col_idx = self->g_H2_pivots[this_ws->pivot.key1][idx].col_idx;

                                  //this_ws->reduce_with = self->g_R_H2 + self->g_R_col_idx_H2[R_col_idx];

                                  this_ws->reduce_with_len = self->g_R_col_idx_H2[this_ws->R_col_idx+1] -\
                                                              self->g_R_col_idx_H2[this_ws->R_col_idx];
                                  this_ws->flag_reduce = 1;

                          }

                      }
                      
                  }


              }



              //if ((this_ws->tetrahedron.key1 == 368174)\
              //  &&(this_ws->tetrahedron.key2 == 333271)){

              //    printf("\npivot before reduction is (%d, %d)"\
              //                        , this_ws->pivot.key1\
              //                        , this_ws->pivot.key2\
              //                        );
              //}

              this_ws->flag_red_w_complex = 0;
              this_ws->flag_append_to_complex = 1;


              //simplex* orig = self->g_workspace_H2[i] \
              //                    + this_ws->original*this_ws->max_len;

              //simplex* scratch = self->g_workspace_H2[i] \
              //                    + (1-this_ws->original)*this_ws->max_len;
                                    



              //printf("\nreducing with %d", idx);

              while(this_ws->flag_reduce){


                    simplex* orig = self->g_workspace_H2[i] + this_ws->original*this_ws->max_len;

                    if ( this_ws->len + this_ws->reduce_with_len > this_ws->max_len){

                        //printf("\nREALLOCATING");

                        if (this_ws->original){
                              
                              //printf("\nCOPYING");
                              
                              for (EDGE_ID mm = 0; mm < this_ws->len; mm++){
                                    
                                    self->g_workspace_H2[i][mm] = self->g_workspace_H2[i][mm + this_ws->max_len];
                                    
                              }

                              this_ws->original = 0;
                              
                        }

                        this_ws->max_len = this_ws->len + this_ws->reduce_with_len + 1000;

                        self->g_workspace_H2[i] = (simplex*)realloc(self->g_workspace_H2[i]\
                                                            , 2*this_ws->max_len*sizeof(simplex));


                        orig = self->g_workspace_H2[i];

                        
                    }
                    
                    
                    simplex* scratch = self->g_workspace_H2[i] \
                                    + (1-this_ws->original)*this_ws->max_len;

                    if (!this_ws->R_col_idx){
                          //if ((this_ws->tetrahedron.key1 == 368174)\
                          //  &&(this_ws->tetrahedron.key2 == 333271)){

                          //    printf("\nreducing with trivial boundary "\
                          //                        );
                          //}

                          red_start = this_ws->trivial_boundary;
                    }
                    else{
                          //if ((this_ws->tetrahedron.key1 == 368174)\
                          //  &&(this_ws->tetrahedron.key2 == 333271)){

                          //    printf("\nreducing with red complex "\
                          //                        );
                          //}

                          red_start = self->g_R_H2 + self->g_R_col_idx_H2[this_ws->R_col_idx];
                    }

                    //if ((this_ws->tetrahedron.key1 == 368174)\
                    //  &&(this_ws->tetrahedron.key2 == 333271)){

                    //    printf("\nReducing ");
                    //    for (EDGE_ID mm = 0; mm < this_ws->len; mm++){
                    //        printf("(%d, %d), "\
                    //                        ,orig[mm].key1\
                    //                        ,orig[mm].key2\
                    //                        );
                    //    }

                    //    printf("\nwith ");
                    //    for (EDGE_ID mm = 0; mm < this_ws->reduce_with_len; mm++){
                    //        printf("(%d, %d), "\
                    //                        ,red_start[mm].key1\
                    //                        ,red_start[mm].key2\
                    //                        );
                    //    }
                    //}

                    count = 0;
                    j = 0;
                    k = 0;

                    int compare;

                    while ((j < this_ws->len) && (k < this_ws->reduce_with_len)){


                        if (orig[j].key1 < red_start[k].key1) compare = 1;
                        else if (orig[j].key1 > red_start[k].key1) compare = 0;
                        else{

                            if (orig[j].key2 < red_start[k].key2) compare = 1;
                            else if (orig[j].key2 > red_start[k].key2) compare = 0;
                            else compare = -1;
                        }



                        if (compare == 1){
                            scratch[count++] = orig[j++];
                        }
                        else if (compare == 0){
                            scratch[count++] = red_start[k++];
                        }
                        else{
                            j++;
                            k++;
                        }

                    }

                    while (j < this_ws->len){

                        scratch[count++] = orig[j++];
                    }

                    while (k < this_ws->reduce_with_len){

                        scratch[count++] = red_start[k++];
                    }

                    this_ws->original = 1 - this_ws->original;

                    this_ws->len = count;

                    //if ((this_ws->tetrahedron.key1 == 368174)\
                    //  &&(this_ws->tetrahedron.key2 == 333271)){

                    //    printf("\nAfter reduction ");
                    //    for (EDGE_ID mm = 0; mm < this_ws->len; mm++){
                    //        printf("(%d, %d), "\
                    //                        ,scratch[mm].key1\
                    //                        ,scratch[mm].key2\
                    //                        );
                    //    }

                    //}

                    if (!this_ws->len){

                        //idx = self->g_n_reduced_simplex[self->g_dim_now];
                        //idx = -1;
                        this_ws->flag_empty = 1;
                        break;

                    }
                    //else{
                    

                    this_ws->pivot = scratch[this_ws->len-1];

                    //printf("\nNew pivot is (%d, %d)"\
                    //                              , this_ws->pivot.key1\
                    //                              , this_ws->pivot.key2\
                    //                              );

                    // Check trivial pers pair


                    //if ((this_ws->tetrahedron.key1 == 368174)\
                    //  &&(this_ws->tetrahedron.key2 == 333271)){

                    //    printf("\npivot after reduction is (%d, %d)"\
                    //                  , this_ws->pivot.key1\
                    //                  , this_ws->pivot.key2\
                    //                  );
                    //    getchar();

                    //}


                    coboundary_H2 temp;

                    temp.triangle.key1 = this_ws->pivot.key1;
                    temp.triangle.key2 = this_ws->pivot.key2;

                    find_H2_cohom_low(self, &temp);

                    this_ws->flag_reduce = 0;

                    // Check if pivot is trivial or if it is pivot in H2
                    // To check trivial, check if maximum triangle in low of cob of triangle is triangle
                    if ((temp.low.key1 == temp.triangle.key1)\
                      &&(self->g_edges_list[2*temp.low.key2+1] == temp.triangle.key2)){

                        // In which case, reduce with boundary of temp.low
                        
                        //if ((this_ws->tetrahedron.key1 == 368174)\
                        //  &&(this_ws->tetrahedron.key2 == 333271)){
                        //    printf("\nreducing with trivial");
                        //}
                        
                        compute_boundary_tetra(self, temp.low, this_ws->trivial_boundary);


                        //this_ws->reduce_with = this_ws->trivial_boundary;
                        
                        this_ws->R_col_idx = 0;
                        this_ws->reduce_with_len = 4;
                        this_ws->flag_reduce = 1;

                    }
                    else{

                        if (self->g_H2_pivots_len[this_ws->pivot.key1]){

                            EDGE_ID idx = search_H2_pivots(self->g_H2_pivots[this_ws->pivot.key1]\
                                                      , 0\
                                                      , self->g_H2_pivots_len[this_ws->pivot.key1] - 1\
                                                      , this_ws->pivot.key2 \
                                                      , self->g_n_valid_edges\
                                                      );

                            if (idx != self->g_n_valid_edges){
                                    //this_ws->flag_red_w_trivial = 0;
                                    
                                    this_ws->R_col_idx = self->g_H2_pivots[this_ws->pivot.key1][idx].col_idx;

                                    this_ws->reduce_with_len = self->g_R_col_idx_H2[this_ws->R_col_idx+1] -\
                                                                self->g_R_col_idx_H2[this_ws->R_col_idx];
                                    this_ws->flag_reduce = 1;

                                    //if ((this_ws->tetrahedron.key1 == 368174)\
                                    //  &&(this_ws->tetrahedron.key2 == 333271)){
                                    //    printf("\nreducing with complex");
                                    //}


                            }

                        }

                    }


                    //printf("\nidx is %d, pivot key2 %d: in red at col idx %d, in R at %d to %d"\
                    //                                      , idx\
                    //                                      , self->g_H2_pivots[this_ws->pivot.key1][idx].key2\
                    //                                      , self->g_H2_pivots[this_ws->pivot.key1][idx].col_idx\
                    //                                      , self->g_R_col_idx_H2[self->g_H2_pivots[this_ws->pivot.key1][idx].col_idx]\
                    //                                      , self->g_R_col_idx_H2[self->g_H2_pivots[this_ws->pivot.key1][idx].col_idx+1]\
                    //                                      );

                    //}
                    
              }


              //if ((this_ws->tetrahedron.key1 == 368174)\
              //  &&(this_ws->tetrahedron.key2 == 333271)){

              //    printf("\nQuitting parallel");
              //}


          }

          pthread_mutex_lock(&(self->g_thread_lock));

          self->g_processed_threads++;
            
      }


}



void update_R_H2 (filtration* self, int ws_counter){
      

      boundary_H2_ws* this_ws = self->g_workspace_H2_info + ws_counter;
      simplex* orig = self->g_workspace_H2[ws_counter] \
                                    + this_ws->original*this_ws->max_len;
      
      // Update R
      if ((self->g_R_len_H2 + this_ws->len) > self->g_R_max_len_H2){
            
            self->g_R_max_len_H2 += 1000 + this_ws->len;

            self->g_R_H2 = (simplex*)realloc(self->g_R_H2, self->g_R_max_len_H2*sizeof(simplex));
            
      }


      // Update R col idx
      self->g_R_col_idx_H2_ptr++;

      if (self->g_R_col_idx_H2_ptr == self->g_R_col_idx_max_len_H2 - 1){
              
            self->g_R_col_idx_max_len_H2 += 1000;

            self->g_R_col_idx_H2 = (EDGE_ID*)realloc(self->g_R_col_idx_H2\
                                          , self->g_R_col_idx_max_len_H2*sizeof(EDGE_ID));


      }


      self->g_R_col_idx_H2[self->g_R_col_idx_H2_ptr] = self->g_R_len_H2;


      //printf("\nAdding pivot %d at %d", this_ws->pivot, self->g_R_col_idx_H2_ptr);
      // ////////////////////
      // ADD PIVOT
      // ////////////////////
      
      add_H2_pivot(self, this_ws->tetrahedron, this_ws->pivot, self->g_R_col_idx_H2_ptr);

      //printf("\nAdding to R with pivot (%d, %d): "\
      //                                    , this_ws->pivot.key1\
      //                                    , this_ws->pivot.key2\
      //                                    );
      
      //////////////////////



      for (EDGE_ID mm = 0; mm < this_ws->len; mm++){
            
            self->g_R_H2[self->g_R_len_H2++] = orig[mm];
            //printf("(%d, %d), ", self->g_R_H2[self->g_R_len_H2 - 1].key1\
            //           , self->g_R_H2[self->g_R_len_H2 - 1].key2);

      }


      
      self->g_R_col_idx_H2[self->g_R_col_idx_H2_ptr+1] = self->g_R_len_H2;

      //if (this_ws->pivot.key1 == 24782){

            //printf("\nR_col_idx at %d, in R from from %d to %d"\
            //                                  , self->g_R_col_idx_H2_ptr\
            //                                  , self->g_R_col_idx_H2[self->g_R_col_idx_H2_ptr]\
            //                                  , self->g_R_col_idx_H2[self->g_R_col_idx_H2_ptr+1]);
            //getchar();
      //}



      //if (this_ws->pivot.key1 == 7173){
      //    printf("\nThe R is from %d to %d"\
      //                    , self->g_R_col_idx_H2[self->g_R_col_idx_H2_ptr]\
      //                    , self->g_R_col_idx_H2[self->g_R_col_idx_H2_ptr+1]\
      //                    );
      //    getchar();
      //}

      //PAR birth = self->g_edge_parameter[this_ws->pivot.key1];
      //PAR death = self->g_edge_parameter[this_ws->tetrahedron.key1];

//#ifdef BIRTH_HOM_CYCLES
//      self->g_H1_pivot_of[this_ws->pivot].key1 = this_ws->triangle.key1;
//      self->g_H1_pivot_of[this_ws->pivot].key2 = this_ws->triangle.key2;
//#endif

      //printf("\n%lf, %lf", birth, death);

}



void add_H2_pivot (filtration* self, simplex tetrahedron, simplex pivot, EDGE_ID red_col){
      
      
    if (self->g_H2_pivots_len[pivot.key1]\
                            == self->g_H2_pivots_max_len[pivot.key1]){
          
          self->g_H2_pivots_max_len[pivot.key1] += 5;
          self->g_H2_pivots[pivot.key1] = (H2_pivots*)realloc( \
                          self->g_H2_pivots[pivot.key1] \
                          , self->g_H2_pivots_max_len[pivot.key1]*sizeof(H2_pivots));

    }


    EDGE_ID old_ptr = self->g_H2_pivots_len[pivot.key1];
    EDGE_ID new_ptr = self->g_H2_pivots_len[pivot.key1];

    while (old_ptr){
          
          old_ptr--;
          
          if (self->g_H2_pivots[pivot.key1][old_ptr].key2 > pivot.key2){

                self->g_H2_pivots[pivot.key1][new_ptr--] =\
                                                       self->g_H2_pivots[pivot.key1][old_ptr];
                continue;

          }
          break;

    }


    self->g_H2_pivots[pivot.key1][new_ptr].key2 = pivot.key2;
    self->g_H2_pivots[pivot.key1][new_ptr].col_idx = red_col;

                      
    self->g_H2_pivots[pivot.key1][new_ptr].tetrahedron = tetrahedron;

    self->g_H2_pivots_len[pivot.key1]++;

    //if (pivot.key1 == 7173){
    //  printf("\nAdding (%d, %d) at R col idx %d"\
    //                        , pivot.key1\
    //                        , pivot.key2\
    //                        , self->g_H2_pivots[pivot.key1][new_ptr].col_idx\
    //                        );
    //  getchar();
    //}

    // PERS PAIRS
    // Add non-zero barcodes
        
    PAR birth = self->g_edge_parameter[pivot.key1];
    PAR death = self->g_edge_parameter[tetrahedron.key1];

    if (birth != death){

           //printf("\nNon trivial pers pair (%f, %f)", birth, death);
           

           //if (birth > death){
           //    printf("\nBirth, death (%lf, %lf)", birth, death);
           //    printf("\nError (%d, %d) at pair (%d, %d)", triangle.key1\
           //                                         , triangle.key2\
           //                                         , pivot.key1\
           //                                         , pivot.key2);
           //    getchar();
           //  
           //}


           self->g_homH2_pers[self->g_homH2_pers_len].birth_simplex = pivot;
           self->g_homH2_pers[self->g_homH2_pers_len].death_edge = tetrahedron.key1;
           self->g_homH2_pers[self->g_homH2_pers_len++].R_col_idx = red_col;

           //printf("\nAdding (%d, %d) at %d in homH2pers for tetra (%d, %d)"\
           //                           , pivot.key1\
           //                           , pivot.key2\
           //                           , self->g_homH2_pers_len-1\
           //                           , tetrahedron.key1\
           //                           , tetrahedron.key2\
           //                           );


           if (self->g_homH2_pers_len == self->g_homH2_pers_max_len){
                 self->g_homH2_pers_max_len += 100;
                 self->g_homH2_pers = (homH2_pers*)realloc(self->g_homH2_pers\
                                                     , self->g_homH2_pers_max_len*sizeof(homH2_pers));
                  
           }



    }
      
}



EDGE_ID search_H2_pivots(H2_pivots* arr, EDGE_ID l, EDGE_ID r, EDGE_ID key2, EDGE_ID max) 
{ 
    if (r >= l) { 
        EDGE_ID mid = l + (r - l) / 2; 

        if (arr[mid].key2 == key2) 
            return mid; 
  
        // If element is smaller than mid, then 
        // it can only be present in left subarray 
        if (arr[mid].key2 > key2) 
        {

          /// PRECAUTIONARY: CAN REMOVE LATER
            if (!mid){
                return max; 
                //printf("\nMID 0 WILL GIVE ERROR FOR UNSIGNED NEXT");
                //getchar();
            }
          ///////////////////
            return search_H2_pivots(arr, l, mid - 1, key2, max); 
        }
  
        // Else the element can only be present 
        // in right subarray 
        return search_H2_pivots(arr, mid + 1, r, key2, max); 
    } 
  
    // We reach here when element is not 
    // present in array 
    //printf("\nNOT FOUND");
    return max; 
} 




void get_birth_void(filtration* self, simplex bo_idx){

    
    self->g_temp_V_H2_primary.len = 1;
    self->g_temp_V_H2_primary.VV[0] = bo_idx;

    // Initiate R with the boundary of bo_idx
    compute_boundary_triangle(self, bo_idx, self->g_temp_R_H2_birth_cycles.RR);

    //printf("\nComputed cob at line 11582");
    //getchar();


    self->g_temp_R_H2_birth_cycles.original = 0;
    self->g_temp_R_H2_birth_cycles.len = 3;

    EDGE_ID* original_result;

    EDGE_ID* scratch_result;

    EDGE_ID j, k, count, possible_len, ro, red_simp_len;

    EDGE_ID* red_start;

    EDGE_ID* trivial_boundary = (EDGE_ID*)malloc(3*sizeof(EDGE_ID));

    EDGE_ID bo_pivot;
    
#ifdef ADAPTIVE_V_STORAGE
    self->g_store_V_for_len = 0;
#endif

    self->g_depth = 0;

    while (self->g_temp_R_H2_birth_cycles.len){

          //printf("\nV len is %d", self->g_temp_V_H2_primary.len);
        

          original_result = self->g_temp_R_H2_birth_cycles.RR \
                       + (self->g_temp_R_H2_birth_cycles.original)*self->g_temp_R_H2_birth_cycles.max_len;


          bo_pivot = original_result[self->g_temp_R_H2_birth_cycles.len-1];

          // Check if the pivot-edge is in a trivial pers pair
          if (self->g_coH1_all_lows[bo_pivot].low.key1 == bo_pivot){

              bo_idx = self->g_coH1_all_lows[bo_pivot].low;

              // Get the boundary of R
              compute_boundary_triangle(self, bo_idx, trivial_boundary);
                
              red_simp_len = 3;
              red_start = trivial_boundary;

              // The V-operation for this is exactly bo_idx 
              self->g_temp_V_H2_primary.VV[self->g_temp_V_H2_primary.len++] = bo_idx;
              // Check for overflow
              if (self->g_temp_V_H2_primary.len == self->g_temp_V_H2_primary.max_len){

                   //getchar();
                   self->g_temp_V_H2_primary.max_len += 100;
                   //printf("\nReallocating %d", self->g_temp_V_H2_primary.max_len);
                   self->g_temp_V_H2_primary.VV = (simplex*)realloc(self->g_temp_V_H2_primary.VV\
                                                           , self->g_temp_V_H2_primary.max_len*sizeof(simplex));

              }
                
          }
          else{

              bo_idx = self->g_H1_pivot_of[bo_pivot].coface;

              ro = self->g_pivots_H1[bo_pivot];

              red_simp_len = self->g_R_col_idx_H1[ro+1] - \
                               self->g_R_col_idx_H1[ro];


              red_start = self->g_R_H1 + self->g_R_col_idx_H1[ro];

              // For this one we have to find V recursively
              // FIND THE V RECURSIVELY
              find_V_recursively_triangles(self, bo_idx, bo_pivot);
               
          }


          // Check for overflow of R_H2
          possible_len = self->g_temp_R_H2_birth_cycles.len + red_simp_len;
          if (possible_len > self->g_temp_R_H2_birth_cycles.max_len - 1){

                if (self->g_temp_R_H2_birth_cycles.original){

                      for (EDGE_ID k = 0; k < self->g_temp_R_H2_birth_cycles.len; k++){
                          self->g_temp_R_H2_birth_cycles.RR[k] =\
                                        self->g_temp_R_H2_birth_cycles.RR[k + self->g_temp_R_H2_birth_cycles.max_len];
                      }
                      
                }

                self->g_temp_R_H2_birth_cycles.max_len = possible_len + 1000;
                self->g_temp_R_H2_birth_cycles.RR = (EDGE_ID*)realloc(self->g_temp_R_H2_birth_cycles.RR\
                                        , (2*self->g_temp_R_H2_birth_cycles.max_len)*sizeof(EDGE_ID));

                original_result = self->g_temp_R_H2_birth_cycles.RR;
                self->g_temp_R_H2_birth_cycles.original = 0;
                
          }

          scratch_result = self->g_temp_R_H2_birth_cycles.RR \
                       + (1-self->g_temp_R_H2_birth_cycles.original)*self->g_temp_R_H2_birth_cycles.max_len;

          // Reduce
          j = 0;
          k = 0;
          count = 0;

          while ((j < self->g_temp_R_H2_birth_cycles.len) && (k < red_simp_len)){

              if (original_result[j] < red_start[k]){
                   
                  scratch_result[count] = original_result[j];
                  count = count + 1;
                  j = j + 1;
                   
              }
              else if (original_result[j] > red_start[k]){

                  scratch_result[count] = red_start[k];
                  count = count + 1;
                  k = k + 1;

              }
              else{

                  j = j + 1;
                  k = k + 1;

              }

          }

          while (j < self->g_temp_R_H2_birth_cycles.len){

               scratch_result[count++] = original_result[j++];

          }

          while (k < red_simp_len){

               scratch_result[count++] = red_start[k++];

          }

          self->g_temp_R_H2_birth_cycles.len = count;

          self->g_temp_R_H2_birth_cycles.original = 1 - self->g_temp_R_H2_birth_cycles.original;
          
          
    }


    //printf("\nreducing temp V at line 11735");
    //getchar();
    reduce_temp_V_H1(self);
    //printf("\nreduced V at line 11735. Press key to continue");
    //getchar();

    free(trivial_boundary);

}



void find_V_recursively_triangles(filtration* self, simplex bo_idx, EDGE_ID bo_pivot){


#ifdef RECORD_V_USAGE
      self->g_H1_pivot_of[bo_pivot].V_usage++;
#endif


#ifdef ADAPTIVE_V_STORAGE
      if (self->g_H1_pivot_of[bo_pivot].V_stored == 1){

          //printf("\nUsing stored cycle for %d of length %d, temp_V_H2 length is %d, max %d"\
          //                                                            , bo_idx\
          //                                                            , self->g_H1_pivot_of[bo_pivot].V_len\
          //                                                            , self->g_temp_V_H2_primary.len\
          //                                                            , self->g_temp_V_H2_primary.max_len\
          //                                                            );
          
          if ((self->g_temp_V_H2_primary.len + self->g_H1_pivot_of[bo_pivot].V_len) > self->g_temp_V_H2_primary.max_len - 1){
                
                self->g_temp_V_H2_primary.max_len =\
                                            self->g_temp_V_H2_primary.len + self->g_H1_pivot_of[bo_pivot].V_len + 1000;
                //printf("\nReallocating %d", self->g_temp_V_H2_primary.max_len);
                self->g_temp_V_H2_primary.VV = (simplex*)realloc(self->g_temp_V_H2_primary.VV\
                                                            , self->g_temp_V_H2_primary.max_len*sizeof(simplex));

          }

          for (EDGE_ID mm = 0; mm < self->g_H1_pivot_of[bo_pivot].V_len; mm++){
                
                self->g_temp_V_H2_primary.VV[self->g_temp_V_H2_primary.len++] =\
                                           self->g_H1_pivot_of[bo_pivot].VV[mm];

          }
          return;
          
      }
      else if (self->g_H1_pivot_of[bo_pivot].V_stored != -1){

#ifndef RECORD_V_USAGE
          self->g_H1_pivot_of[bo_pivot].V_usage++;
#endif

          if ((self->g_H1_pivot_of[bo_pivot].V_usage > self->g_cycle_usage_thresh)\
               &&(!self->g_H1_pivot_of[bo_pivot].V_stored)){

              //printf("\nstore simplex (%d, %d) that has pivot %d", bo_idx.key1\
              //                                                    , bo_idx.key2\
              //                                                    , bo_pivot);

              self->g_store_V_for[self->g_store_V_for_len++] = bo_pivot;
              if (self->g_store_V_for_len == self->g_store_V_for_max_len){
                    self->g_store_V_for_max_len += 100;
                    self->g_store_V_for = (EDGE_ID*)realloc(self->g_store_V_for\
                                                        , self->g_store_V_for_max_len*sizeof(EDGE_ID));
              }

          }

      }

#endif


      //printf("\nComputing V for %d, temp_V_H2 length is %d, max %d" , bo_idx\
      //                                                              , self->g_temp_V_H2_primary.len\
      //                                                              , self->g_temp_V_H2_primary.max_len\
                                                                    );


      EDGE_ID res_max_len = 100;
      EDGE_ID* result = (EDGE_ID*)malloc((2*res_max_len)*sizeof(EDGE_ID));
      int res_original = 0;
      EDGE_ID* original_result;
      EDGE_ID* scratch_result;

      //printf("\nstarting red for depth %d", self->g_depth);

      // Initiate R with the boundary of bo_idx
      compute_boundary_triangle(self, bo_idx, result);

      EDGE_ID res_len = 3;

      EDGE_ID possible_len;

      EDGE_ID* red_start;
      EDGE_ID red_simp_len;
      EDGE_ID ro;
      EDGE_ID j, k, count;


      EDGE_ID* trivial_boundary  = (EDGE_ID*)malloc(3*sizeof(EDGE_ID));

      self->g_temp_V_H2_primary.VV[self->g_temp_V_H2_primary.len++] = bo_idx;
      
      // Check for overflow
      if (self->g_temp_V_H2_primary.len == self->g_temp_V_H2_primary.max_len){

           //printf("\nReallocating");
           //getchar();
           self->g_temp_V_H2_primary.max_len += 100;
           //printf("\nReallocating %d", self->g_temp_V_H2_primary.max_len);
           self->g_temp_V_H2_primary.VV = (simplex*)realloc(self->g_temp_V_H2_primary.VV\
                                                   , self->g_temp_V_H2_primary.max_len*sizeof(simplex));

      }


      //int flag_recursive = 0;

      while(res_len != 0){

#ifdef ADAPTIVE_V_STORAGE
          self->g_depth++;
#endif

          original_result = result + (res_original*res_max_len);

          bo_pivot = original_result[res_len-1];

          // Check for trivial pers pair
          if (self->g_coH1_all_lows[bo_pivot].low.key1 == bo_pivot){

              //printf("\nreducing with trivial");
              bo_idx = self->g_coH1_all_lows[bo_pivot].low;

              // Get the boundary of R
              compute_boundary_triangle(self, bo_idx, trivial_boundary);
                
              red_simp_len = 3;
              red_start = trivial_boundary;


          }
          else{

              //printf("\nreducing with complex");
              bo_idx = self->g_H1_pivot_of[bo_pivot].coface;

              ro = self->g_pivots_H1[bo_pivot];

              //printf("\n r_col_idx is %d", ro);

              red_simp_len = self->g_R_col_idx_H1[ro+1] - \
                               self->g_R_col_idx_H1[ro];


              red_start = self->g_R_H1 + self->g_R_col_idx_H1[ro];


          }


          // Check for overflow
          possible_len = res_len + red_simp_len;
          if (possible_len > res_max_len - 1){

                if (res_original){

                      for (k = 0; k < res_len; k++){
                          result[k] = result[k + res_max_len];
                      }
                      
                }

                res_max_len = possible_len + 1000;
                result = (EDGE_ID*)realloc(result, (2*res_max_len)*sizeof(EDGE_ID));

                original_result = result;
                res_original = 0;
                
          }

          scratch_result = result + ((1-res_original)*res_max_len);

          //printf("\nReducing ");
          //for (EDGE_ID mm = 0; mm < res_len; mm++){
          //    printf("%d, ", original_result[mm]);
          //}
          //printf("\nwith ");
          //for (EDGE_ID mm = 0; mm < red_simp_len; mm++){
          //    printf("%d, ", red_start[mm]);
          //}

          // Reduce
          j = 0;
          k = 0;
          count = 0;


          while ((j < res_len) && (k < red_simp_len)){

              if (original_result[j] < red_start[k]){
                   
                  scratch_result[count] = original_result[j];
                  count = count + 1;
                  j = j + 1;
                   
              }
              else if (original_result[j] > red_start[k]){

                  scratch_result[count] = red_start[k];
                  count = count + 1;
                  k = k + 1;

              }
              else{

                  j = j + 1;
                  k = k + 1;

              }

          }

          while (j < res_len){

               scratch_result[count++] = original_result[j++];

          }

          while (k < red_simp_len){

               scratch_result[count++] = red_start[k++];

          }


          res_len = count;

          res_original = 1 - res_original;

          if (res_len != 0){

                find_V_recursively_triangles(self, bo_idx, bo_pivot);

          }

      }



      free(result);
      free(trivial_boundary);
  
}



void compute_boundary_tetra(filtration* self, simplex tetra, simplex* simp){

        

        EDGE_ID o_ab = tetra.key1;
        
        VERT_ID a = self->g_edges_list[2*o_ab];
        VERT_ID b = self->g_edges_list[2*o_ab+1];
    
        EDGE_ID o_cd = tetra.key2;

        VERT_ID c = self->g_edges_list[2*o_cd];
        VERT_ID d = self->g_edges_list[2*o_cd+1];


        // Since d > c, abd is the maximum triangle
        simp[3].key1 = o_ab;
        simp[3].key2 = d;

        // abc
        simp[2].key1 = o_ab;
        simp[2].key2 = c;

        // Remaining triangles are bcd and acd
        // Get all edges first
        VERT_ID idx = search_Neighbors(self, b, d, 0, self->g_Neigh_len[b]-1);
        EDGE_ID o_bd = self->g_Neighbors[b][idx].order;

        idx = search_Neighbors(self, a, d, 0, self->g_Neigh_len[a]-1);
        EDGE_ID o_ad = self->g_Neighbors[a][idx].order;


        idx = search_Neighbors(self, a, c, 0, self->g_Neigh_len[a]-1);
        EDGE_ID o_ac = self->g_Neighbors[a][idx].order;

        idx = search_Neighbors(self, b, c, 0, self->g_Neigh_len[b]-1);
        EDGE_ID o_bc = self->g_Neighbors[b][idx].order;



        // Triangle bcd
        EDGE_ID o_max = o_bc;
        VERT_ID v3 = d;

        if (o_cd > o_max){
            o_max = o_cd;
            v3 = b;
        }

        if (o_bd > o_max){
            o_max = o_bd;
            v3 = c;
        }

        // Triangle acd
        EDGE_ID o_max_2 = o_ac;
        VERT_ID v3_2 = d;

        if (o_cd > o_max_2){
            o_max_2 = o_cd;
            v3_2 = a;
        }

        if (o_ad > o_max_2){
            o_max_2 = o_ad;
            v3_2 = c;
        }


        int compare = -1;
        if (o_max > o_max_2) compare = 1;
        else if (o_max < o_max_2) compare = 0;
        else{
            if (v3 > v3_2) compare = 1;
            else if (v3 < v3_2) compare = 0;
            else compare = -1;
        }

        if (compare == 1){
              
            simp[1].key1 = o_max;
            simp[1].key2 = v3;

            simp[0].key1 = o_max_2;
            simp[0].key2 = v3_2;
              
        }
        else if (compare == 0){

            simp[1].key1 = o_max_2;
            simp[1].key2 = v3_2;

            simp[0].key1 = o_max;
            simp[0].key2 = v3;

        }
        else{
              printf("\nERROR? 10119");
              getchar();
        }
      

}



void compute_boundary_triangle(filtration* self, simplex triangle, EDGE_ID* simp){

      
        VERT_ID a = self->g_edges_list[2*triangle.key1];
        VERT_ID b = self->g_edges_list[2*triangle.key1+1];
        VERT_ID c = triangle.key2;
        
        
        VERT_ID idx = search_Neighbors(self, a, c, 0, self->g_Neigh_len[a]-1);
        EDGE_ID o_ac = self->g_Neighbors[a][idx].order;
        
        idx = search_Neighbors(self, b, c, 0, self->g_Neigh_len[b]-1);
        EDGE_ID o_bc = self->g_Neighbors[b][idx].order;
        
        if (o_ac > o_bc){
              simp[0] = o_bc;
              simp[1] = o_ac;
        }
        else{
              simp[0] = o_ac;
              simp[1] = o_bc;
        }
        
        simp[2] = triangle.key1;

        //printf("\nComputed boundary of triangle (%d, %d, %d)"\
        //                                    , simp[0]\
        //                                    , simp[1]\
        //                                    , simp[2]\
        //                                    );

}






void store_V_H0(filtration* self){
    
      for (EDGE_ID nn = 0; nn < self->g_store_V_for_len; nn++){
          
            EDGE_ID bo_pivot = self->g_store_V_for[nn];
      
            if (self->g_H0_pivot_of[bo_pivot].V_len){
              continue;
            }
            //printf("\nstoring %d", bo_pivot);
      
            self->g_temp_V_primary.len = 0;
      
            EDGE_ID bo_idx = self->g_H0_pivot_of[bo_pivot].coface;
      
            self->g_depth = 0;
      
            find_V_recursively_edges(self, bo_idx, bo_pivot);
      
#ifdef RECORD_V_USAGE
            self->g_H0_pivot_of[bo_pivot].V_depth = self->g_depth;
#endif
            if (self->g_depth < self->g_cycle_depth_thresh){
                self->g_H0_pivot_of[bo_pivot].V_stored = -1;
                continue;
            }


            self->g_n_H0_stored_V++;
      
            //printf("\nstoring %d with depth %d", bo_pivot, self->g_depth);

            self->g_H0_pivot_of[bo_pivot].V_stored = 1;

            // Reduce V
            reduce_temp_V_H0(self);
      
            //printf("\nlen is %d", self->g_temp_V_primary.len);
      
            self->g_H0_pivot_of[bo_pivot].V_len = self->g_temp_V_primary.len;
            self->g_H0_pivot_of[bo_pivot].VV = \
                                  (EDGE_ID*)malloc(self->g_temp_V_primary.len*sizeof(EDGE_ID));
      
            for (EDGE_ID oo = 0; oo < self->g_temp_V_primary.len; oo++){
                  self->g_H0_pivot_of[bo_pivot].VV[oo] = self->g_temp_V_primary.VV[oo];
            }
            
            
      }

}

void reduce_temp_V_H0(filtration* self){

      //// Reduce V
      sorter8_tim_sort(self->g_temp_V_primary.VV, self->g_temp_V_primary.len);          

      int coeff = 1;
      EDGE_ID idx = 0;

      for (EDGE_ID vv = 0; vv < self->g_temp_V_primary.len-1; vv++){
            
          if (self->g_temp_V_primary.VV[vv] == self->g_temp_V_primary.VV[vv+1])
          {
              coeff = 1 - coeff;
          }
          else{
              if (coeff){
                  //self->g_V_sparse_H1[self->g_V_sparse_ptr++] = this_ws->v_edges.o_ab[vv];
                  self->g_temp_V_primary.VV[idx++] = self->g_temp_V_primary.VV[vv];
              }

              coeff = 1;
          }
            
            
      }

      if (coeff){
           self->g_temp_V_primary.VV[idx++] = self->g_temp_V_primary.VV[self->g_temp_V_primary.len-1];
      }

      self->g_temp_V_primary.len = idx;

}



void store_V_H1(filtration* self){
    
      for (EDGE_ID nn = 0; nn < self->g_store_V_for_len; nn++){
          
            EDGE_ID bo_pivot = self->g_store_V_for[nn];

            //printf("\nPress key to store V for pivot %d", bo_pivot);
            //getchar();
      
            if (self->g_H1_pivot_of[bo_pivot].V_len){
              continue;
            }

            // Check for trivial pers pair
            if (self->g_coH1_all_lows[bo_pivot].low.key1 == bo_pivot){
                self->g_H1_pivot_of[bo_pivot].V_stored = -1;
                //printf("\nNot storing because trivial, %d", bo_pivot);
                //getchar();
                continue;
            }

            self->g_n_H1_stored_V++;

            //printf("\nchecking depth for pivot %d", bo_pivot);
            //getchar();
      
            self->g_temp_V_H2_primary.len = 0;
      
            simplex bo_idx = self->g_H1_pivot_of[bo_pivot].coface;

            //printf(", simplex is (%d, %d)", bo_idx.key1, bo_idx.key2);
            //getchar();
      
            self->g_depth = 0;
      
            find_V_recursively_triangles(self, bo_idx, bo_pivot);

#ifdef RECORD_V_USAGE
            self->g_H1_pivot_of[bo_pivot].V_depth = self->g_depth;
#endif
      
            if (self->g_depth < self->g_cycle_depth_thresh){
                self->g_H1_pivot_of[bo_pivot].V_stored = -1;
                continue;
            }

            //printf("\nstoring %d with depth %d", bo_pivot, self->g_depth);
      
            self->g_H1_pivot_of[bo_pivot].V_stored = 1;
            //self->g_H1_pivot_of[bo_pivot].V_depth = self->g_depth;

            // Reduce V
            reduce_temp_V_H1(self);
      
      
            self->g_H1_pivot_of[bo_pivot].V_len = self->g_temp_V_H2_primary.len;
            self->g_H1_pivot_of[bo_pivot].VV = \
                                  (simplex*)malloc(self->g_temp_V_H2_primary.len*sizeof(simplex));
      
            for (EDGE_ID oo = 0; oo < self->g_temp_V_H2_primary.len; oo++){
                  self->g_H1_pivot_of[bo_pivot].VV[oo] = self->g_temp_V_H2_primary.VV[oo];
            }
            
            
      }

}


void reduce_temp_V_H1(filtration* self){
    

    //// Reduce V
    sorter4_tim_sort(self->g_temp_V_H2_primary.VV, self->g_temp_V_H2_primary.len);          

    int coeff = 1;
    EDGE_ID idx = 0;

    for (EDGE_ID vv = 0; vv < self->g_temp_V_H2_primary.len-1; vv++){
          
        if ((self->g_temp_V_H2_primary.VV[vv].key1 == self->g_temp_V_H2_primary.VV[vv+1].key1)\
          &&(self->g_temp_V_H2_primary.VV[vv].key2 == self->g_temp_V_H2_primary.VV[vv+1].key2))
        {
            coeff = 1 - coeff;
        }
        else{
            if (coeff){
                //self->g_V_sparse_H1[self->g_V_sparse_ptr++] = this_ws->v_edges.o_ab[vv];
                self->g_temp_V_H2_primary.VV[idx++] = self->g_temp_V_H2_primary.VV[vv];
            }

            coeff = 1;
        }
          
          
    }

    if (coeff){
         self->g_temp_V_H2_primary.VV[idx++] = self->g_temp_V_H2_primary.VV[self->g_temp_V_H2_primary.len-1];
    }

    self->g_temp_V_H2_primary.len = idx;


}






EDGE_ID bin_search_max_less_V(EDGE_ID* arr, EDGE_ID l, EDGE_ID r, EDGE_ID x, EDGE_ID MAX){

    //printf("\nl %d, r %d", l, r);

    if (arr[r] >= x){
      //printf("\nreturing max");
      return MAX;
    }

    if (arr[l] < x){
      //printf("\nreturing l %d", l);
      return l;
    }

    EDGE_ID mid = l + (r-l)/2;
    //printf("\nmid is %d", mid);

    if (arr[mid] >= x){
        
        l = mid+1;
        bin_search_max_less_V(arr, l, r, x, MAX);

    }
    else{
        r = mid-1;
        bin_search_max_less_V(arr, l , r, x, MAX);

    }
     
}




EDGE_ID bin_search_max_less_updated_V(min_update_V* arr, EDGE_ID l, EDGE_ID r, EDGE_ID x, EDGE_ID MAX){

    //printf("\nl %d, r %d", l, r);

    if (arr[r].cycid >= x){
      //printf("\nreturing max");
      return MAX;
    }

    if (arr[l].cycid < x){
      //printf("\nreturing l %d", l);
      return l;
    }

    EDGE_ID mid = l + (r-l)/2;
    //printf("\nmid is %d", mid);

    if (arr[mid].cycid >= x){
        
        l = mid+1;
        bin_search_max_less_updated_V(arr, l, r, x, MAX);

    }
    else{
        r = mid-1;
        bin_search_max_less_updated_V(arr, l , r, x, MAX);

    }
     
}


EDGE_ID bin_search_min_greater_updated_V_byLidx(EDGE_ID* arr, EDGE_ID l, EDGE_ID r, EDGE_ID x, EDGE_ID MAX){

    //printf("\nl %d, r %d", l, r);

    if (arr[r] <= x){
      //printf("\nreturing max");
      return MAX;
    }

    if (arr[l] > x){
      //printf("\nreturing l %d", l);
      return l;
    }

    EDGE_ID mid = l + (r-l)/2;
    //printf("\nmid is %d", mid);

    if (arr[mid] <= x){
        
        l = mid+1;
        bin_search_min_greater_updated_V_byLidx(arr, l, r, x, MAX);

    }
    else{
        r = mid;
        if (arr[r-1] <= x) return r;
        bin_search_min_greater_updated_V_byLidx(arr, l , r, x, MAX);

    }
     
}



void minimize_birth_cycles_H0_v2(filtration* self\
                              , EDGE_ID** stored_boundaries\
                              , EDGE_ID* len_boundaries\
                              , EDGE_ID stored_num\
                              , char* filename\
                              ){
      

     EDGE_ID** all_diff = (EDGE_ID**)malloc(self->g_all_V_stored_num*sizeof(EDGE_ID*));

     for (EDGE_ID mm = 0; mm < self->g_all_V_stored_num; mm++){
          all_diff[mm] = (EDGE_ID*)malloc(2*sizeof(EDGE_ID));
     }


     EDGE_ID* updated = (EDGE_ID*)calloc(stored_num, sizeof(EDGE_ID));
     EDGE_ID* original_id = (EDGE_ID*)malloc(stored_num*sizeof(EDGE_ID));
     for (EDGE_ID mm = 0; mm < stored_num; mm++){
         original_id[mm] = mm;
     }

     mergeSort_V_H0(len_boundaries, stored_boundaries, updated, original_id \
                , 0, stored_num-1);

     //printf("\nInitial sums");

     //omp_set_num_threads(2*self->g_cpu_count-1);

     #pragma omp parallel for schedule(guided) shared(self, stored_boundaries, len_boundaries, stored_num, all_diff)
     for (EDGE_ID mm = 0; mm < stored_num; mm++){

          if (!self->g_suppress_output){
              if (mm%1000 == 0)
              printf("\nDoing %d", mm);
          }
          
          all_diff[mm][0] = 0;
          all_diff[mm][1] = 0;


          for (EDGE_ID nn = mm+1; nn < stored_num; nn++){

                if (len_boundaries[nn] < all_diff[mm][1]){
                    break;
                }


                EDGE_ID j, k, count;
                j = 0;
                k = 0;
                count = 0;

                int quit_flag = 0;

                while ((j < len_boundaries[mm]) && (k < len_boundaries[nn])){
                      
                      if (stored_boundaries[mm][j] < stored_boundaries[nn][k]){
                          count++;
                          j++;
                      }
                      else if (stored_boundaries[mm][j] > stored_boundaries[nn][k]){
                          count++;
                          k++;
                      }
                      else{
                          j++;
                          k++;
                      }

                      if (count > len_boundaries[mm] - all_diff[mm][1]){
                          quit_flag = 1;
                          break;
                      }
                      
                }

                if (quit_flag) continue;

                if (j < len_boundaries[mm]){
                    count += len_boundaries[mm] - j;
                }

                if (count > len_boundaries[mm] - all_diff[mm][1]){
                    continue;
                }

                if (k < len_boundaries[nn]){
                    count += len_boundaries[nn] - k;
                }

                if (count > len_boundaries[mm] - all_diff[mm][1]){
                    continue;
                }


                if (count < len_boundaries[mm]){

                    if ((len_boundaries[mm] - count) > all_diff[mm][1]){

                          all_diff[mm][0] = nn;
                          all_diff[mm][1] = len_boundaries[mm] - count;
                        
                    }
                    
                }

          }
          
     }


    
      
}



//////////////////////////////////////////////////////////
// MERGE SORT ALGORITHMS FOR update_V
//////////////////////////////////////////////////////////

// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge_update_V(min_update_V* arr, EDGE_ID l, EDGE_ID m, EDGE_ID r) 
{ 
    EDGE_ID i, j, k; 
    EDGE_ID n1 = m - l + 1; 
    EDGE_ID n2 =  r - m; 
    //printf("\nn1, n2: %u, %u", n1, n2);
  
    /* create temp arrays */
    min_update_V *L, *R;
    L = (min_update_V*)malloc(n1*sizeof(min_update_V));
    R = (min_update_V*)malloc(n2*sizeof(min_update_V));

    
    //int L[n1], R[n2]; 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){ 
        L[i] = arr[l + i]; 
    }

    for (j = 0; j < n2; j++) {
        R[j] = arr[m + 1+ j]; 
    }
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 

    while (i < n1 && j < n2) 
    { 
          
          if (L[i].cycid > R[j].cycid) 
          { 
              arr[k] = L[i]; 
              i++; 
          } 
          else
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

    free(L);
    free(R);
} 
  
/* l is for left index and r is right index of the 
   sub-array of arr to be sorted */
void mergeSort_update_V(min_update_V* arr, EDGE_ID l, EDGE_ID r) 
{ 
    if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        EDGE_ID m = l+(r-l)/2; 
  
        // Sort first and second halves 
        mergeSort_update_V(arr, l, m); 
        mergeSort_update_V(arr, m+1, r); 
  
        merge_update_V(arr, l, m, r); 
    } 

} 


// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge_update_V_byLidx(min_update_V* arr, EDGE_ID l, EDGE_ID m, EDGE_ID r) 
{ 
    EDGE_ID i, j, k; 
    EDGE_ID n1 = m - l + 1; 
    EDGE_ID n2 =  r - m; 
    //printf("\nn1, n2: %u, %u", n1, n2);
  
    /* create temp arrays */
    min_update_V *L, *R;
    L = (min_update_V*)malloc(n1*sizeof(min_update_V));
    R = (min_update_V*)malloc(n2*sizeof(min_update_V));

    
    //int L[n1], R[n2]; 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){ 
        L[i] = arr[l + i]; 
    }

    for (j = 0; j < n2; j++) {
        R[j] = arr[m + 1+ j]; 
    }
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 

    while (i < n1 && j < n2) 
    { 
          
          if (L[i].Lidx < R[j].Lidx) 
          { 
              arr[k] = L[i]; 
              i++; 
          } 
          else
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

    free(L);
    free(R);
} 
  
/* l is for left index and r is right index of the 
   sub-array of arr to be sorted */
void mergeSort_update_V_byLidx(min_update_V* arr, EDGE_ID l, EDGE_ID r) 
{ 
    if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        EDGE_ID m = l+(r-l)/2; 
  
        // Sort first and second halves 
        mergeSort_update_V_byLidx(arr, l, m); 
        mergeSort_update_V_byLidx(arr, m+1, r); 
  
        merge_update_V_byLidx(arr, l, m, r); 
    } 

} 








EDGE_ID find_first_diff_H0(EDGE_ID* len_boundaries\
                          , EDGE_ID stored_num\
                          , EDGE_ID** stored_boundaries){
      
    EDGE_ID difff = 0;
      
    for (EDGE_ID mm = 0; mm < stored_num; mm++){
    
          if (len_boundaries[mm] < difff){
              break;
          }
    
    
          for (EDGE_ID nn = mm+1; nn < stored_num; nn++){
    
                if (len_boundaries[nn] < difff){
                    break;
                }
    
                EDGE_ID j, k, count;
                j = 0;
                k = 0;
                count = 0;
    
                int quit_flag = 0;
    
                while ((j < len_boundaries[mm]) && (k < len_boundaries[nn])){
                      
                      if (stored_boundaries[mm][j] < stored_boundaries[nn][k]){
                          count++;
                          j++;
                      }
                      else if (stored_boundaries[mm][j] > stored_boundaries[nn][k]){
                          count++;
                          k++;
                      }
                      else{
                          j++;
                          k++;
                      }
    
                      if (count > len_boundaries[mm] - difff){
                          quit_flag = 1;
                          break;
                      }
                      
                }
    
                if (quit_flag) continue;
    
                if (j < len_boundaries[mm]){
                    count += len_boundaries[mm] - j;
                }
                
                if (count > len_boundaries[mm] - difff){
                    continue;
                }
    
                if (k < len_boundaries[nn]){
                    count += len_boundaries[nn] - k;
                }
    
                if (count > len_boundaries[mm] - difff){
                    continue;
                }
    
                //printf("\n%d, %d, %d",len_boundaries[mm], count, difff);
    
                if (count < len_boundaries[mm]){
    
                    if ((len_boundaries[mm] - count) > difff){
    
                          difff = len_boundaries[mm] - count;
    
                    }

                }
    
          }
    
    }

    return difff;

}





void minimize_birth_cycles_H0_v3(filtration* self\
                              , cyc_info* CC\
                              , EDGE_ID stored_num\
                              , char* filename\
                              , char* lens_filename\
                              , char* minimal_lens_filename\
                              , char* subset_points_file\
                              ){

                              //, char* filename2\

      if (!self->g_suppress_output){
          printf("\nNumber of cycles %d", stored_num);
      }
      //getchar();
      

#ifdef STORE_LENGTHS_CYCLES
      FILE* fp0 = fopen(lens_filename, "w");
      for (EDGE_ID ci = 0; ci < stored_num; ci++){
           fprintf(fp0, "%d, ", CC[ci].len);
           if (!self->g_reduce_cyc_lengths){
                free(CC[ci].boundary);
           }
      }
      if (!self->g_reduce_cyc_lengths){
        free(CC);
      }
      fclose(fp0);
#endif

      if (!self->g_reduce_cyc_lengths){
            return;
      }



      omp_set_num_threads(2*self->g_cpu_count - 1);

      EDGE_ID* Lcycid = (EDGE_ID*)malloc(stored_num*sizeof(EDGE_ID));
      EDGE_ID* Llen = (EDGE_ID*)malloc(stored_num*sizeof(EDGE_ID));
      EDGE_ID* Lupdated = (EDGE_ID*)calloc(stored_num, sizeof(EDGE_ID));


      update_in_cyc** update_in_cycle = (update_in_cyc**)malloc(self->g_n_valid_edges*sizeof(update_in_cyc*));
      EDGE_ID* update_in_cycle_len = (EDGE_ID*)calloc(self->g_n_valid_edges, sizeof(EDGE_ID));
      EDGE_ID* update_in_cycle_max_len = (EDGE_ID*)calloc(self->g_n_valid_edges, sizeof(EDGE_ID));

      for (EDGE_ID bb = 0; bb < self->g_n_valid_edges; bb++){
            update_in_cycle_max_len[bb] = 1;
            update_in_cycle[bb] = (update_in_cyc*)malloc(sizeof(update_in_cyc));
      }


      EDGE_ID* update_edges = (EDGE_ID*)malloc(self->g_n_valid_edges*sizeof(EDGE_ID));
      EDGE_ID update_edges_num = 0;


      EDGE_ID* case2a2bb = (EDGE_ID*)malloc(stored_num*sizeof(EDGE_ID));
      EDGE_ID case2a2bb_num = 0;

      EDGE_ID* case2ba = (EDGE_ID*)malloc(stored_num*sizeof(EDGE_ID));
      EDGE_ID case2ba_num = 0;

      //printf("\nIntializing L, C.Lidx...");

      // Step 1. Initialize L and C.Lidx
      for (EDGE_ID i = 0; i < stored_num; i++){
          Lcycid[i] = i;
          Llen[i] = CC[i].len;
          Lupdated[i] = 0;
      }

      
      //printf("\nSorting Llen...");

      // Step 2(a): Sort Llen, Lcycid, Lupdated by Llen
      mergeSort_Llen(Llen, Lcycid, Lupdated, 0, stored_num - 1);

      //printf("\nInitializing C.Lidx...");
      // Step 2(b): Initialize C.Lidx
      for (EDGE_ID li = 0; li < stored_num; li++){
            CC[Lcycid[li]].Lidx = li;
      }


      //// Cycle-intersects-cycle information
      //cyc_in_cyc** cyc_inter_cyc = (cyc_in_cyc**)malloc(stored_num*sizeof(cyc_in_cyc*));
      //EDGE_ID* cyc_inter_cyc_len = (EDGE_ID*)calloc(stored_num, sizeof(EDGE_ID));
      //EDGE_ID* cyc_inter_cyc_max_len = (EDGE_ID*)calloc(stored_num, sizeof(EDGE_ID));


      // BUILD VERT to cycle association
      for (EDGE_ID ci = 0; ci < stored_num; ci++){

            EDGE_ID li = CC[ci].Lidx;
      
            for (EDGE_ID idx = 0; idx < CC[ci].len; idx++){
                
                  EDGE_ID edge = CC[ci].boundary[idx];

                  //printf("\nedge is %d with %d cycles ", edge, self->g_edges_in_cycles_len[edge]);

                  //for (EDGE_ID nn = 0; nn < self->g_edges_in_cycles_len[edge]; nn++){
                  //    printf("%d, ", self->g_edges_in_cycles[edge][nn]);
                  //}
                  
                  if (!self->g_edges_in_cycles_len[edge]){
                      self->g_edges_in_cycles[edge] = (EDGE_ID*)malloc(sizeof(EDGE_ID));
                  }
                  else{
                      self->g_edges_in_cycles[edge] = (EDGE_ID*)realloc(\
                                                  self->g_edges_in_cycles[edge]\
                                                  ,(self->g_edges_in_cycles_len[edge]+1)*sizeof(EDGE_ID));
                  }

                  self->g_edges_in_cycles[edge][self->g_edges_in_cycles_len[edge]++] = ci;

            }
            
      }

      

      // NOTE: At this point g_edges_in_cycles is sorted by cyc_id

      if (!self->g_suppress_output){
          printf("\nInitializing diff...");
      }

      #pragma omp parallel for schedule(static, 1000) shared(Lcycid, CC, Llen)
      for (EDGE_ID li = 0; li < stored_num; li++){
            
            EDGE_ID ci = Lcycid[li];
            CC[ci].diff = 0;

            if (!self->g_suppress_output){
                if (li %1000 == 0){
                    printf("\n%d", li);
                }
            }

            for (EDGE_ID lj = li + 1; lj < stored_num; lj++){
            //for (EDGE_ID lj = 0; lj < stored_num; lj++){

                  if (Llen[lj] < CC[ci].diff){
                      break;
                  }

                  //if (Llen[lj] > Llen[li]){
                  //    continue;
                  //}

                  EDGE_ID cj = Lcycid[lj];

                  //if (CC[cj].perspair[0] > CC[ci].perspair[0]){
                  //    continue;
                  //}

                  //if (CC[cj].updated_birth > self->g_cycle_min_birth_thresh){
                  //    continue;
                  //}

                  EDGE_ID j, k, count;
                  j = 0;
                  k = 0;
                  count = 0;
                  
                  int quit_flag = 0;

                  while ((j < Llen[li]) && (k < Llen[lj])){
                        
                      if (CC[ci].boundary[j] < CC[cj].boundary[k]){
                            j++;
                            count++;
                      }
                      else if (CC[ci].boundary[j] > CC[cj].boundary[k]){
                            k++;
                            count++;
                      }
                      else{
                            j++;
                            k++;
                      }

                      if (count >= Llen[li]){
                          quit_flag = 1;
                          break;
                      }

                      if ((Llen[li] - CC[ci].diff) <= count){
                          quit_flag = 1;
                          break;
                      }
                        
                  }

                  if (quit_flag){
                      continue;
                  }

                  if (j < Llen[li]){
                      count += Llen[li] - j;
                  }
                  if (k < Llen[lj]){
                      count += Llen[lj] - k;
                  }

                  if (count >= Llen[li]){
                      continue;
                  }

                  if ((Llen[li] - CC[ci].diff) <= count){
                      continue;
                  }


                  // At this point len - count < diff?
                  CC[ci].diff = Llen[li] - count;
                  CC[ci].redw = cj;

            }


      }

      //for (EDGE_ID ci = 0; ci < stored_num; ci++){
      //    if (CC[ci].diff)
      //    printf("\nmax diff for %d is %d", ci, CC[ci].diff);
      //}
      //getchar();


      // Define V that will store summations to be done
      EDGE_ID V_len = 0;
      EDGE_ID V_max_len = 10;
      min_update_V* update_V = (min_update_V*)malloc(V_max_len*sizeof(min_update_V));
      EDGE_ID* update_v_indices = (EDGE_ID*)malloc(stored_num*sizeof(EDGE_ID));


      int* update_edges_flag = (int*)calloc(self->g_n_valid_edges, sizeof(int));
      update_edges_num = 0;

      EDGE_ID it_counter = 0;

      // Step 4: Loop for minimization
      while (1){

            struct timespec ss0, ss1, ss2, ss3, ss4, ss5, ss6;
            struct timespec ff0, ff1, ff2, ff3, ff4, ff5, ff6;
            double cc0 = 0;
            double cc1 = 0;
            double cc2 = 0;
            double cc3 = 0;
            double cc4 = 0;
            double cc5 = 0;
            double cc6 = 0;

            if (!self->g_suppress_output){
                printf("\n\nIteration %d", it_counter++);
                clock_gettime(CLOCK_MONOTONIC, &ss0);
                clock_gettime(CLOCK_MONOTONIC, &ss1);

            }

            //printf("\nFINDING MAX DIFF");
            // Step 4(a): Find max diff
            EDGE_ID dm = 1;
            V_len = 0;
            for (EDGE_ID ci = 0; ci < stored_num; ci++){
                  if (CC[ci].diff > dm){
                        dm = CC[ci].diff;
                        update_V[0].cycid = ci;
                        V_len = 1;
                  }
                  else if (CC[ci].diff == dm){

                        update_V[V_len++].cycid = ci;

                        if (V_len == V_max_len){
                              V_max_len += 100;
                              update_V = (min_update_V*)realloc(update_V\
                                                      , V_max_len*sizeof(min_update_V));
                        }
                  }
            }

            if (!self->g_suppress_output){
                clock_gettime(CLOCK_MONOTONIC, &ff1);
                cc1 += (ff1.tv_sec - ss1.tv_sec);
                cc1 += (ff1.tv_nsec - ss1.tv_nsec) / 1000000000.0;
            }



            if (!V_len){
                //printf("\nDiff 0. EXITING.");
                break;
            }


            //printf("\nMax diff %d in pairs %d", dm, V_len);
            
            //printf("\ndo 15453");
            //getchar();
            
            if (!self->g_suppress_output){
                clock_gettime(CLOCK_MONOTONIC, &ss2);
            }

            //#pragma omp parallel for schedule(static, 50) shared(update_V, CC)
            for (EDGE_ID vi = 0; vi < V_len; vi++){

                //if (vi % 1000 == 0){
                //    printf("\nreducing %d", vi);
                //}
                  
                EDGE_ID ci = update_V[vi].cycid;
                EDGE_ID cj = CC[ci].redw;

                //printf("\nReducing %d with %d", ci, cj);

                EDGE_ID* scratch = (EDGE_ID*)malloc((CC[ci].len + CC[cj].len)*sizeof(EDGE_ID));


                EDGE_ID edge;

                EDGE_ID j = 0;
                EDGE_ID k = 0;
                EDGE_ID count = 0;

                while ((j < CC[ci].len) && (k < CC[cj].len)){
                      
                      if (CC[ci].boundary[j] < CC[cj].boundary[k]){

                            scratch[count++] = CC[ci].boundary[j++];

                      }
                      else if (CC[ci].boundary[j] > CC[cj].boundary[k]){

                            scratch[count] = CC[cj].boundary[k];

                            // This edge is new in ci
                            // So, add this edge to ci-info
                            edge = CC[cj].boundary[k];

                            if (!update_edges_flag[edge]){
                                update_edges_flag[edge] = 1;
                                update_edges[update_edges_num++] = edge;
                            }

                            update_in_cycle[edge][update_in_cycle_len[edge]].cyc = ci;
                            // Add ci in edge
                            update_in_cycle[edge][update_in_cycle_len[edge]].flag = 1;

                            update_in_cycle_len[edge]++;
                            if (update_in_cycle_len[edge] == update_in_cycle_max_len[edge]){
                                  update_in_cycle_max_len[edge] += 100;
                                  update_in_cycle[edge] = (update_in_cyc*)realloc(update_in_cycle[edge]\
                                                                      , update_in_cycle_max_len[edge]*sizeof(update_in_cyc));
                            }

                            count++;
                            k++;
                      }
                      else{

                            // This edge is not there anymore in ci
                            // so, remove ci from this edge's info
                            edge = CC[ci].boundary[j];
                            
                            if (!update_edges_flag[edge]){
                                update_edges_flag[edge] = 1;
                                update_edges[update_edges_num++] = edge;
                            }

                            update_in_cycle[edge][update_in_cycle_len[edge]].cyc = ci;
                            // Remove ci from edge
                            update_in_cycle[edge][update_in_cycle_len[edge]].flag = 0;

                            update_in_cycle_len[edge]++;
                            if (update_in_cycle_len[edge] == update_in_cycle_max_len[edge]){
                                  update_in_cycle_max_len[edge] += 100;
                                  update_in_cycle[edge] = (update_in_cyc*)realloc(update_in_cycle[edge]\
                                                                      , update_in_cycle_max_len[edge]*sizeof(update_in_cyc));
                            }

                            j++;
                            k++;

                      }

                }

                // update_in_cyc will have unique

                while (j < CC[ci].len){

                    // No change
                    scratch[count++] = CC[ci].boundary[j++];

                }

                while (k < CC[cj].len){

                    // These are all new edges
                    // So, add ci to every edge-info
                    scratch[count] = CC[cj].boundary[k];

                    edge = CC[cj].boundary[k];

                    if (!update_edges_flag[edge]){
                        update_edges_flag[edge] = 1;
                        update_edges[update_edges_num++] = edge;
                    }

                    update_in_cycle[edge][update_in_cycle_len[edge]].cyc = ci;
                    // Add ci in edge
                    update_in_cycle[edge][update_in_cycle_len[edge]].flag = 1;

                    update_in_cycle_len[edge]++;
                    if (update_in_cycle_len[edge] == update_in_cycle_max_len[edge]){
                          update_in_cycle_max_len[edge] += 100;
                          update_in_cycle[edge] = (update_in_cyc*)realloc(update_in_cycle[edge]\
                                                              , update_in_cycle_max_len[edge]*sizeof(update_in_cyc));
                    }

                    count++;
                    k++;

                }


                scratch = (EDGE_ID*)realloc(scratch, count*sizeof(EDGE_ID));


                update_V[vi].VV = scratch;
                update_V[vi].V_len = count;


            }

            if (!self->g_suppress_output){

                clock_gettime(CLOCK_MONOTONIC, &ff2);
                cc2 += (ff2.tv_sec - ss2.tv_sec);
                cc2 += (ff2.tv_nsec - ss2.tv_nsec) / 1000000000.0;

                clock_gettime(CLOCK_MONOTONIC, &ss6);

                clock_gettime(CLOCK_MONOTONIC, &ff6);
                cc6 += (ff6.tv_sec - ss6.tv_sec);
                cc6 += (ff6.tv_nsec - ss6.tv_nsec) / 1000000000.0;

                clock_gettime(CLOCK_MONOTONIC, &ss3);

                printf("\nNumber of edges to update %d", update_edges_num);

            }

            #pragma omp parallel for schedule(static, 50) shared(self, stored_num\
                                                            , CC\
                                                            , Lcycid, Lupdated\
                                                            , update_V, V_len\
                                                            , update_edges\
                                                            , update_in_cycle\
                                                            , update_in_cycle_len\
                                                            )
            for (EDGE_ID iddx = 0; iddx < update_edges_num; iddx++){
                  
                  EDGE_ID ei = update_edges[iddx];

                  update_edges_flag[ei] = 0;
                
                  EDGE_ID* scratch = (EDGE_ID*)malloc(\
                                            (self->g_edges_in_cycles_len[ei]+update_in_cycle_len[ei])*sizeof(EDGE_ID));
                  
                  EDGE_ID o_ptr = 0;
                  EDGE_ID u_ptr = 0;

                  EDGE_ID s_ptr = 0;

                  while ((o_ptr < self->g_edges_in_cycles_len[ei]) && (u_ptr < update_in_cycle_len[ei])){
                        
                        if (self->g_edges_in_cycles[ei][o_ptr] < update_in_cycle[ei][u_ptr].cyc){
                              
                              scratch[s_ptr++] = self->g_edges_in_cycles[ei][o_ptr++];
                              
                        }
                        else if (self->g_edges_in_cycles[ei][o_ptr] > update_in_cycle[ei][u_ptr].cyc){

                              // Check if update_in_cycle is flagged for addition
                              if (update_in_cycle[ei][u_ptr].flag){
                                  
                                  scratch[s_ptr++] = update_in_cycle[ei][u_ptr++].cyc;

                              }
                              else{
                                  // Sanity check
                                  printf("\nERROR 15604");
                                  getchar();
                              }
                              
                        }
                        else{
                              
                              // Check if update_in_cycle is flagged for addition
                              if (update_in_cycle[ei][u_ptr].flag){

                                  // Flagged for addition: means, it is
                                  // already there in original. Just copy and increment all pointers
                                  
                                  scratch[s_ptr++] = self->g_edges_in_cycles[ei][o_ptr++];
                                  u_ptr++;

                              }
                              else{
                                  // Flagged for removal: means, it is
                                  // to be skipped 
                                  o_ptr++;
                                  u_ptr++;
                              }
                              
                        }

                  }

                  // If o_ptr did not reach end: means copy all that are remaining as nothing to
                  // update
                  while (o_ptr < self->g_edges_in_cycles_len[ei]){

                      scratch[s_ptr++] = self->g_edges_in_cycles[ei][o_ptr++];

                  }

                  // If u_ptr did not reach end: means all of these should have been flagged for addition
                  while (u_ptr < update_in_cycle_len[ei]){

                      scratch[s_ptr++] = update_in_cycle[ei][u_ptr++].cyc;

                  }


                  free(self->g_edges_in_cycles[ei]);

                  self->g_edges_in_cycles[ei] = scratch;
                  self->g_edges_in_cycles_len[ei] = s_ptr;

                  update_in_cycle_len[ei] = 0;
                  update_in_cycle_max_len[ei] = 1;
                  update_in_cycle[ei] = (update_in_cyc*)realloc(\
                                                  update_in_cycle[ei]\
                                                 , sizeof(update_in_cyc));

            }

            update_edges_num = 0;


            if (!self->g_suppress_output){

                clock_gettime(CLOCK_MONOTONIC, &ff3);
                cc3 += (ff3.tv_sec - ss3.tv_sec);
                cc3 += (ff3.tv_nsec - ss3.tv_nsec) / 1000000000.0;

                clock_gettime(CLOCK_MONOTONIC, &ss4);

            }

            // Update CC with new boundaries, update Llen and Lupdated
            for (EDGE_ID vi = 0; vi < V_len; vi++){
                  
                  EDGE_ID ci = update_V[vi].cycid;
                  EDGE_ID li = CC[ci].Lidx;

                  free(CC[ci].boundary);

                  CC[ci].boundary = update_V[vi].VV;
                  CC[ci].len = update_V[vi].V_len;

                  Llen[li] = update_V[vi].V_len;
                  Lupdated[li] = 1;

            }



            // Sort Llen, Lupdated, Lcycid
            mergeSort_Llen(Llen, Lcycid, Lupdated, 0, stored_num - 1);

            // Update C.Lidx
            for (EDGE_ID li = 0; li < stored_num; li++){
                  CC[Lcycid[li]].Lidx = li;
            }

            // Update V.Lidx
            for (EDGE_ID vi = 0; vi < V_len; vi++){
                  update_V[vi].Lidx = CC[update_V[vi].cycid].Lidx;
            }


            // Sort update_V by Lidx
            if (V_len > 1){
                  mergeSort_update_V_byLidx(update_V, 0, V_len-1);
            }

            for (EDGE_ID vi = 0; vi < V_len; vi++){
                  update_v_indices[vi] = update_V[vi].Lidx; 
            }

            //// Sort each g_edges_in_cycles by len of cycles
            //for (EDGE_ID ei = 0; ei < self->g_n_valid_edges; ei++){
            //      
            //      if (self->g_edges_in_cycles_len[ei] < 2) continue;

            //      mergeSort_edges_in_cycles(self->g_edges_in_cycles[ei], CC, 0, self->g_edges_in_cycles_len[ei]-1);

            //      
            //}

            if (!self->g_suppress_output){

                clock_gettime(CLOCK_MONOTONIC, &ff4);
                cc4 += (ff4.tv_sec - ss4.tv_sec);
                cc4 += (ff4.tv_nsec - ss4.tv_nsec) / 1000000000.0;

                
                clock_gettime(CLOCK_MONOTONIC, &ss5);
            }

            // Get the cases
            // Case 1: x is updated -> check with all with diff = 0 (this is n update_V)
            // CASE 2a: x is not updated + diff is 0 -> only check with updated
            // CASE 2ba: x is not updated + diff is not 0 + y is updated -> check with all after new diff
            // Case 2bb: x is not updated + diff is not 0 + y is not updated -> only check with updated
            
            case2a2bb_num = 0;
            case2ba_num = 0;

            for (EDGE_ID li = 0; li < stored_num; li++){
                  
                  EDGE_ID ci = Lcycid[li];

                  if (Lupdated[li]) continue;

                  // x is not updated
                  if (!CC[ci].diff){

                      // diff is 0 -> only check with updated -> case 2
                      case2a2bb[case2a2bb_num++] = li;

                  }
                  else{

                        EDGE_ID cj = CC[ci].redw;
                        EDGE_ID lj = CC[cj].Lidx;
                        if (Lupdated[lj]){

                              // diff is not 0 and y is updated -> check with all, but diff is not set to 0 -> case 1
                              case2ba[case2ba_num++] = li;
                              
                        }
                        else{
                              // diff is not 0 and y is not updated -> only check with updated -> case 2
                              case2a2bb[case2a2bb_num++] = li;

                        }

                  }

            }

            if (!self->g_suppress_output){
                clock_gettime(CLOCK_MONOTONIC, &ff5);
                cc5 += (ff5.tv_sec - ss5.tv_sec);
                cc5 += (ff5.tv_nsec - ss5.tv_nsec) / 1000000000.0;
            }


            //printf("\nStarting 15790");
            //getchar();

            //////////////////////////////////////////////////////////
            //printf("\nUpdating diffs...");
            // Update diffs
            //////////////////////////////////////////////////////////
            
            // NOTE: updated_V is sorted by increasing Lidx at this point


            struct timespec start1, start2, start3;
            struct timespec finish1, finish2, finish3;
            double c1 = 0;
            double c2 = 0;
            double c3 = 0;

            // NEW NEW 
            // CASE 1

            //printf("\nStarting 15810");
            //getchar();

            if (!self->g_suppress_output){
                clock_gettime(CLOCK_MONOTONIC, &start1);
            }

            #pragma omp parallel for schedule(static, 50) shared(stored_num\
                                                            , CC\
                                                            , Lcycid, Lupdated\
                                                            , update_V, V_len\
                                                            )
            for (EDGE_ID vi = 0; vi < V_len; vi++){

                  EDGE_ID li = update_V[vi].Lidx;
                  EDGE_ID ci = Lcycid[li];

            //for (EDGE_ID ci = 0; ci < stored_num; ci++){

                  CC[ci].diff = 0;
                  //EDGE_ID li = CC[ci].Lidx;
                  //printf("\ncase 1 diff before for cycle %d is %d", ci, CC[ci].diff);
                  //// flag_case = 0 to check with all
                  update_diff(self, li, Lupdated, 0, Lcycid, CC, stored_num);
                  //printf("\ncase 1 diff after for cycle %d is %d", ci, CC[ci].diff);

            }

            if (!self->g_suppress_output){
                clock_gettime(CLOCK_MONOTONIC, &finish1);
                c1 += (finish1.tv_sec - start1.tv_sec);
                c1 += (finish1.tv_nsec - start1.tv_nsec) / 1000000000.0;
                clock_gettime(CLOCK_MONOTONIC, &start2);
            }
            // CASE 2a2bb: Only check with updated cycles
            
            #pragma omp parallel for schedule(static, 50) shared(stored_num\
                                                            , CC\
                                                            , Lcycid, Lupdated\
                                                            , case2a2bb_num, case2a2bb)
            for (EDGE_ID idx = 0; idx < case2a2bb_num; idx++){
                  
                  EDGE_ID li = case2a2bb[idx];
                  //EDGE_ID ci = Lcycid[li];
                  //printf("\ncase 2 diff before for cycle %d is %d", ci, CC[ci].diff);
                  
                  // flag_case = 1 to check with only updated
                  
                  //update_diff(self, li, Lupdated, 1, Lcycid, CC, stored_num);


                  minimal_CASE2(self, li, CC, Lcycid, Llen, update_v_indices, V_len);

                  //printf("\ncase 2 diff after for cycle %d is %d", ci, CC[ci].diff);
                  
            }


            if (!self->g_suppress_output){
                clock_gettime(CLOCK_MONOTONIC, &finish2);
                c2 += (finish2.tv_sec - start2.tv_sec);
                c2 += (finish2.tv_nsec - start2.tv_nsec) / 1000000000.0;
                clock_gettime(CLOCK_MONOTONIC, &start3);
            }

            // CASE 2ba
            
            #pragma omp parallel for schedule(static, 50) shared(stored_num\
                                                            , CC\
                                                            , Lcycid, Lupdated\
                                                            , case2ba, case2ba_num)
            for(EDGE_ID idx = 0; idx < case2ba_num; idx++){
                  
                  EDGE_ID li = case2ba[idx];
                  EDGE_ID ci = Lcycid[li];
                  EDGE_ID cj = CC[ci].redw;

                  //printf("\ncase 3 diff before for cycle %d is %d", ci, CC[ci].diff);

                  //////////////
                  // CASE 2ba: x is NOT updated and diff is NOT 0 and y is updated
                  //////////////
                  //printf(" Case 2ba");

                  EDGE_ID j = 0;
                  EDGE_ID k = 0;
                  EDGE_ID count = 0;
                  while ((j < CC[ci].len) && (k < CC[cj].len)){
                      if (CC[ci].boundary[j] < CC[cj].boundary[k]){
                            j++;
                            count++;
                      }
                      else if (CC[ci].boundary[j] > CC[cj].boundary[k]){
                            k++;
                            count++;
                      }
                      else{
                            j++;
                            k++;
                      }

                      if (count >= CC[ci].len){
                          break;
                      }
                  }

                  if (j < CC[ci].len){
                      count += CC[ci].len - j;
                  }
                  if (k < CC[cj].len){
                      count += CC[cj].len - k;
                  }

                  if (count < CC[ci].len){
                        CC[ci].diff = CC[ci].len - count;
                  }
                  else{
                        CC[ci].diff = 0;
                  }

                  // flag_case = 0 to check with all
                  update_diff(self, li, Lupdated, 0, Lcycid, CC, stored_num);
                                    
                  //printf("\ncase 3 diff after for cycle %d is %d", ci, CC[ci].diff);

            }

            if (!self->g_suppress_output){
                clock_gettime(CLOCK_MONOTONIC, &finish3);
                c3 += (finish3.tv_sec - start3.tv_sec);
                c3 += (finish3.tv_nsec - start3.tv_nsec) / 1000000000.0;

                clock_gettime(CLOCK_MONOTONIC, &ff0);
                cc0 += (ff0.tv_sec - ss0.tv_sec);
                cc0 += (ff0.tv_nsec - ss0.tv_nsec) / 1000000000.0;

                printf("\nmax diff %d, num pairs %d, case 1 %lf, case 2a2bb %lf, case 2ba %lf, (%lf,%lf, %lf, %lf, %lf, %lf), %lf"\
                                                                    , dm\
                                                                    , V_len\
                                                                    , c1\
                                                                    , c2\
                                                                    , c3\
                                                                    , cc1\
                                                                    , cc2\
                                                                    , cc3\
                                                                    , cc4\
                                                                    , cc5\
                                                                    , cc6\
                                                                    , cc0\
                                                                    );
            }

            // Reset Lupdated
            for (EDGE_ID li = 0; li < stored_num; li++){
                  Lupdated[li] = 0;
            }

            //getchar();
            
      }



//////////////////////////////////////////////////////////////
      // TESTING
//////////////////////////////////////////////////////////////
//

      //printf("\nPress key to test");
      //getchar();
      

      //printf("\nTESTING...");

      //#pragma omp parallel for schedule(static) shared(stored_num\
      //                                                  , CC)
      //for (EDGE_ID ci = 0; ci < stored_num; ci++){

      //      EDGE_ID diff = 0;
      //      
      //      for (EDGE_ID cj = 0; cj < stored_num; cj++){
      //            
      //            if (cj == ci) continue;

      //            EDGE_ID j = 0;
      //            EDGE_ID k = 0;
      //            EDGE_ID count = 0;

      //            int quit_flag = 0;

      //            while ((j < CC[ci].len) && (k < CC[cj].len)){
      //                  if (CC[ci].boundary[j] < CC[cj].boundary[k]){
      //                        j++;
      //                        count++;
      //                  }
      //                  else if (CC[ci].boundary[j] > CC[cj].boundary[k]){
      //                        k++;
      //                        count++;
      //                  }
      //                  else{
      //                        j++;
      //                        k++;
      //                  }


      //            }

      //            if (j < CC[ci].len){
      //                  count += CC[ci].len - j;
      //            }

      //            if (k < CC[cj].len){
      //                  count += CC[cj].len - k;
      //            }

      //            
      //            if (count < CC[ci].len){
      //                  if (CC[ci].len - count > diff){
      //                      printf("\nImprovement possible by reducing ci, li %d, %d with cj, lj %d, %d!"\
      //                                                                              , ci, CC[ci].Lidx\
      //                                                                              , cj, CC[cj].Lidx);
      //                      printf("\nTESTING FAILED!!!");
      //                      getchar();
      //                  }
      //            }
      //            
      //      }
      //      
      //}

      //printf("\nTESTED OK.");


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


      // Sort Llen, Lupdated, Lcycid
      mergeSort_Llen(Llen, Lcycid, Lupdated, 0, stored_num - 1);
      // Update C.Lidx
      for (EDGE_ID li = 0; li < stored_num; li++){
            CC[Lcycid[li]].Lidx = li;
      }


#ifdef STORE_LENGTHS_CYCLES
      FILE* fp1 = fopen(minimal_lens_filename, "w");
      for (EDGE_ID li = 0; li < stored_num; li++){
           fprintf(fp1, "%d, ", Llen[li]);
      }

      //for (EDGE_ID ci = 0; ci < stored_num; ci++){
      //     fprintf(fp1, "%d, ", CC[ci].len);
      //}
      fclose(fp1);
#endif


      FILE* fp2 = fopen(filename, "w");

      if (!self->g_suppress_output){
          printf("\n");
      }


      for (EDGE_ID li = 0; li < stored_num; li++){

           if (Llen[li] < 5){
                break;
           }

           EDGE_ID ci = Lcycid[li];

           //if (li < 15)
           //printf("\nlen %d", Llen[li]);


           for (EDGE_ID nn = 0; nn < CC[ci].len; nn++){
                
                 fprintf(fp2, "%d, %d, ", self->g_edges_list[2*CC[ci].boundary[nn]]\
                                        , self->g_edges_list[2*CC[ci].boundary[nn]+1]);

           }
           
           fprintf(fp2, "\n");

           free(CC[ci].boundary);


      }

      fclose(fp2);


      free(CC);

      free(update_V);

      free(Lcycid);
      free(Llen);
      free(Lupdated);
      //free(update_v_indices);
      //
      free(case2ba);
      free(case2a2bb);

      free(update_edges);
      free(update_edges_flag);

      free(update_in_cycle_len);
      free(update_in_cycle_max_len);

      for (EDGE_ID ii = 0; ii < self->g_n_valid_edges; ii++){
          free(update_in_cycle[ii]);
      }


      if (!self->g_suppress_output){
          printf("\nDone. Press key to quit H1 cycle shortening");
      }

}

void minimize_birth_cycles_H0_v4(filtration* self\
                              , cyc_info* CC\
                              , EDGE_ID stored_num\
                              , char* filename\
                              , char* filename2\
                              ){
      
      
      //printf("\nNumber of cycles %d", stored_num);
      //getchar();

      omp_set_num_threads(2*self->g_cpu_count - 1);

      EDGE_ID* Lcycid = (EDGE_ID*)malloc(stored_num*sizeof(EDGE_ID));
      EDGE_ID* Llen = (EDGE_ID*)malloc(stored_num*sizeof(EDGE_ID));
      EDGE_ID* Lupdated = (EDGE_ID*)calloc(stored_num, sizeof(EDGE_ID));

      printf("\nIntializing L, C.Lidx...");

      // Step 1. Initialize L and C.Lidx
      for (EDGE_ID i = 0; i < stored_num; i++){
          Lcycid[i] = i;
          Llen[i] = CC[i].len;
          Lupdated[i] = 0;
      }

      
      printf("\nSorting Llen...");

      // Step 2(a): Sort Llen, Lcycid, Lupdated by Llen
      mergeSort_Llen(Llen, Lcycid, Lupdated, 0, stored_num - 1);

      printf("\nInitializing C.Lidx...");
      // Step 2(b): Initialize C.Lidx
      for (EDGE_ID li = 0; li < stored_num; li++){
            CC[Lcycid[li]].Lidx = li;
      }


      PAR* sorted_par = (PAR*)malloc(stored_num*sizeof(PAR));
      EDGE_ID* Pcycid = (EDGE_ID*)malloc(stored_num*sizeof(EDGE_ID));

      for (EDGE_ID ci = 0; ci < stored_num; ci++){
        
          if (CC[ci].perspair[1] != -1){
              sorted_par[ci] = CC[ci].perspair[1] - CC[ci].perspair[0];
          }
          else{
              sorted_par[ci] = self->g_thresh - CC[ci].perspair[0];
          }
          Pcycid[ci] = ci;
      }

      // SORT temp_par by pers barcode as follows: 
      // Sorted in decreasing order of parameter where:
      // For dead features: parameter = death - birth
      // For undead features: parameter = thresh - birth
      mergeSort_temp_par(sorted_par, Pcycid, 0, stored_num-1);


      for (EDGE_ID pi = 0; pi < stored_num; pi++){
            
            EDGE_ID ci = Pcycid[pi];

            while(1){
               
                  EDGE_ID li = CC[ci].Lidx;
                  CC[ci].diff = 0;

                  for (EDGE_ID lj = li + 1; lj < stored_num; lj++){
                        
                        if (Llen[lj] < CC[ci].diff) break;

                        EDGE_ID cj = Lcycid[lj];

                        EDGE_ID j = 0;
                        EDGE_ID k = 0;
                        EDGE_ID count = 0;


                        while ((j < CC[ci].len) && (k < CC[cj].len)){
                              
                              if (CC[ci].boundary[j] < CC[cj].boundary[k]){
                                    j++;
                                    count++;
                              }
                              else if (CC[ci].boundary[j] > CC[cj].boundary[k]){
                                    k++;
                                    count++;
                              }
                              else{
                                    j++;
                                    k++;
                              }

                        }

                        if (j < CC[ci].len){
                            count += CC[ci].len - j;
                        }

                        if (k < CC[cj].len){
                            count += CC[cj].len - k;
                        }

                        if (count >= CC[ci].len - CC[ci].diff){
                            continue;
                        }

                        CC[ci].redw = cj;
                        CC[ci].diff = CC[ci].len - count;

                        
                  }

                  if (!CC[ci].diff){
                        printf("\nFinal new len of (%lf, %lf) is %d"\
                                                            , CC[ci].perspair[0]\
                                                            , CC[ci].perspair[1]\
                                                            , CC[ci].len);
                        getchar();
                        break;
                  }


                  EDGE_ID j = 0;
                  EDGE_ID k = 0;
                  EDGE_ID count = 0;
                  EDGE_ID cj = CC[ci].redw;

                  printf("\nreducing %d (len %d) with %d (len %d, (%lf, %lf))"\
                                                                 , ci\
                                                                 , CC[ci].len\
                                                                 , cj\
                                                                 , CC[cj].len\
                                                                 , CC[cj].perspair[0]\
                                                                 , CC[cj].perspair[1]\
                                                                 );

                  EDGE_ID* scratch = (EDGE_ID*)malloc((CC[ci].len+CC[cj].len)*sizeof(EDGE_ID));

                  while ((j < CC[ci].len) && (k < CC[cj].len)){
                        
                        if (CC[ci].boundary[j] < CC[cj].boundary[k]){
                              scratch[count++] = CC[ci].boundary[j++];
                        }
                        else if (CC[ci].boundary[j] > CC[cj].boundary[k]){
                              scratch[count++] = CC[cj].boundary[k++];
                        }
                        else{
                              j++;
                              k++;
                        }

                  }

                  while (j < CC[ci].len){
                      scratch[count++] = CC[ci].boundary[j++];
                  }

                  while (k < CC[cj].len){
                      scratch[count++] = CC[cj].boundary[k++];
                  }


                  free(CC[ci].boundary);

                  printf("\nfor pers pair (%lf, %lf), Old len is %d and new len is %d"\
                                                                          , CC[ci].perspair[0]\
                                                                          , CC[ci].perspair[1]\
                                                                          , CC[ci].len\
                                                                          , count);

                  CC[ci].boundary = scratch;
                  CC[ci].len = count;

                  Llen[ci] = count;

                  printf("\nUpdate lengths");
                  // Sort Llen, Lupdated, Lcycid
                  mergeSort_Llen(Llen, Lcycid, Lupdated, 0, stored_num - 1);
                  // Update C.Lidx
                  for (EDGE_ID li = 0; li < stored_num; li++){
                        CC[Lcycid[li]].Lidx = li;
                  }


            }

      }



      
}

void minimal_CASE1(EDGE_ID li, cyc_info* CC, EDGE_ID* Lcycid, EDGE_ID* Llen, EDGE_ID stored_num){
      

      EDGE_ID ci = Lcycid[li];

      //CC[ci].diff = 0;

      for (EDGE_ID lj = li + 1; lj < stored_num; lj++){
      //for (EDGE_ID lj = 0; lj < stored_num; lj++){
            
            if (Llen[lj] < CC[ci].diff){
                  break;
            }

            //if (lj == li) continue;

            EDGE_ID cj = Lcycid[lj];

            EDGE_ID j = 0;
            EDGE_ID k = 0;
            EDGE_ID count = 0;

            int quit_flag = 0;

            while ((j < CC[ci].len) && (k < CC[cj].len)){
                  
                  if (CC[ci].boundary[j] < CC[cj].boundary[k]){
                        j++;
                        count++;
                  }
                  else if (CC[ci].boundary[j] > CC[cj].boundary[k]){
                        k++;
                        count++;
                  }
                  else{
                        j++;
                        k++;
                  }

                  if (count >= CC[ci].len){
                      quit_flag = 1;
                      break;
                  }

                  if ((CC[ci].len - count) <= CC[ci].diff){
                      quit_flag = 1;
                      break;
                  }
                  
            }

            if (quit_flag){
                continue;
            }

            if (j < CC[ci].len){
                count += CC[ci].len - j;
            }

            if (k < CC[cj].len){
                count += CC[cj].len - k;
            }

            // NEED TO CHECK LOGIC HERE!!!!!!!!!!!!!!!!!!!!!!
            if (count >= CC[ci].len){
                continue;
            }

            // NEED TO CHECK THIS <= OR < !!!!!!!!!!!!!!!!
            if ((CC[ci].len - count) <= CC[ci].diff){
                  continue;
            }

            CC[ci].diff = CC[ci].len - count;
            CC[ci].redw = cj;
            //printf("\ncase1 updating diff to %d", CC[ci].diff);
            //getchar();
            
      }


}



void minimal_CASE2(filtration* self, EDGE_ID li, cyc_info* CC, EDGE_ID* Lcycid, EDGE_ID* Llen\
                  , EDGE_ID* update_v_indices, EDGE_ID V_len){

      if (update_v_indices[V_len-1] <= li) return;
      
      EDGE_ID idx = bin_search_min_greater_updated_V_byLidx(update_v_indices\
                                            , 0, V_len-1\
                                            , li\
                                            , V_len);
      
      EDGE_ID ci = Lcycid[li];
      
      for (EDGE_ID vj = idx; vj < V_len; vj++){
            
            EDGE_ID lj = update_v_indices[vj];
            if (lj <= li) continue;
            if (Llen[lj] < CC[ci].diff){
                break;
            }
      
            EDGE_ID cj = Lcycid[lj];
             
            //if (CC[cj].perspair[0] > CC[ci].perspair[0]){
            //    continue;
            //}

            // CHECK BIRTH THRESH
            //if (CC[cj].updated_birth > self->g_cycle_min_birth_thresh){
            //    continue;
            //}

      
            EDGE_ID j = 0;
            EDGE_ID k = 0;
            EDGE_ID count = 0;
            EDGE_ID common = 0;
      
            int quit_flag = 0;
      
            while ((j < CC[ci].len) && (k < CC[cj].len)){
                  
                  if (CC[ci].boundary[j] < CC[cj].boundary[k]){
                        j++;
                        count++;
                  }
                  else if (CC[ci].boundary[j] > CC[cj].boundary[k]){
                        k++;
                        count++;
                  }
                  else{
                        j++;
                        k++;
                  }
      
                  if (count >= CC[ci].len){
                      quit_flag = 1;
                      break;
                  }
      
                  if ((CC[ci].len - count) <= CC[ci].diff){
                      quit_flag = 1;
                      break;
                  }
      
                   
            }
      
            if (quit_flag){
                continue;
            }
      
            if (j < CC[ci].len){
                count += CC[ci].len - j;
            }
      
            if (k < CC[cj].len){
                count += CC[cj].len - k;
            }
      
            // NEED TO CHECK LOGIC HERE!!!!!!!!!!!!!!!!!!!!!!
            if (count >= CC[ci].len){
                continue;
            }
            
            // NEED TO CHECK LOGIC HERE!!!!!!!!!!!!!!!!!!!!!!
            if ((CC[ci].len - count) <= CC[ci].diff){
                continue;
            }
      
            CC[ci].diff = CC[ci].len - count;
            CC[ci].redw = cj;
            
            
      }

      
}

void shuffle_cyc(cyc_info* CC, EDGE_ID num){
      

    int n = (int)num;
    
    srand((unsigned)time(NULL));
    for (int i = 0; i < n - 1; i++) {
        size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
        cyc_info t = CC[j];
        CC[j] = CC[i];
        CC[i] = t;
    } 
      
      
}


void update_diff(filtration* self, EDGE_ID li, EDGE_ID* Lupdated, int flag_case\
                , EDGE_ID* Lcycid, cyc_info* CC, EDGE_ID stored_num){


    EDGE_ID ci = Lcycid[li];
    
    EDGE_ID* cj_diff = (EDGE_ID*)calloc(stored_num, sizeof(EDGE_ID));

    EDGE_ID max_diff = 0;
    
    for (EDGE_ID idx = 0; idx < CC[ci].len; idx++){
    
          EDGE_ID edge = CC[ci].boundary[idx];
    
          for (EDGE_ID idx2 = 0; idx2 < self->g_edges_in_cycles_len[edge]; idx2++){
                
                EDGE_ID cj = self->g_edges_in_cycles[edge][idx2];

                //if (CC[cj].perspair[0] > CC[ci].perspair[0]){
                //    continue;
                //}

                // CHECK BIRTH THRESH
                //if (CC[cj].updated_birth > self->g_cycle_min_birth_thresh){
                //    continue;
                //}

                EDGE_ID lj = CC[cj].Lidx;

                //if (CC[cj].len < max_diff) break;

                if (lj <= li) continue;

                if (CC[cj].len < CC[ci].diff) continue;
    
                if (flag_case){
                    if (!Lupdated[lj]) continue;
                }

                cj_diff[cj] += 2;



                if (CC[cj].len < cj_diff[cj]){
                      
                      if ((cj_diff[cj] - CC[cj].len) > CC[ci].diff){

                          CC[ci].diff = cj_diff[cj] - CC[cj].len;
                          CC[ci].redw = cj;

                          max_diff = CC[ci].len + CC[cj].len - cj_diff[cj];

                      }

                }

    
          }
    
    }


    free(cj_diff);
    
      
}

void find_first_diff(filtration* self, EDGE_ID li, EDGE_ID* Lupdated\
                , EDGE_ID* Lcycid, cyc_info* CC, EDGE_ID stored_num){


    EDGE_ID ci = Lcycid[li];
    
    EDGE_ID* cj_diff = (EDGE_ID*)calloc(stored_num, sizeof(EDGE_ID));

    for (EDGE_ID idx = 0; idx < CC[ci].len; idx++){
    
          EDGE_ID edge = CC[ci].boundary[idx];
    
          for (EDGE_ID idx2 = 0; idx2 < self->g_edges_in_cycles_len[edge]; idx2++){
                
                EDGE_ID cj = self->g_edges_in_cycles[edge][idx2];
                EDGE_ID lj = CC[cj].Lidx;

                if (CC[cj].len < CC[ci].diff) break;

                if (lj <= li) continue;
    
                cj_diff[cj] += 2;

                if (CC[cj].len < cj_diff[cj]){
                      
                      if ((cj_diff[cj] - CC[cj].len) > CC[ci].diff){

                          CC[ci].diff = cj_diff[cj] - CC[cj].len;
                          CC[ci].redw = cj;

                      }

                }

    
          }
    
    }


    free(cj_diff);
    
      
}







void minimize_birth_cycles_H1_v2(filtration* self\
                              , cyc_info_H2* CC\
                              , EDGE_ID stored_num\
                              , char* filename\
                              , char* lens_filename\
                              , char* minimal_lens_filename\
                              ){

                              //, char* filename2\

      if (!self->g_suppress_output){
          printf("\nNumber of cycles %d", stored_num);
      }

#ifdef STORE_LENGTHS_CYCLES
      FILE* fp0 = fopen(lens_filename, "w");
      for (EDGE_ID ci = 0; ci < stored_num; ci++){
           fprintf(fp0, "%d, ", CC[ci].len);
           if (!self->g_reduce_cyc_lengths){
                free(CC[ci].boundary);
           }
      }
      fclose(fp0);
#endif

      if (!self->g_reduce_cyc_lengths){
            free(CC);
            return;
      }


      omp_set_num_threads(2*self->g_cpu_count - 1);


      EDGE_ID* Lcycid = (EDGE_ID*)malloc(stored_num*sizeof(EDGE_ID));
      EDGE_ID* Llen = (EDGE_ID*)malloc(stored_num*sizeof(EDGE_ID));
      EDGE_ID* Lupdated = (EDGE_ID*)calloc(stored_num, sizeof(EDGE_ID));


      // Step 1. Initialize L and C.Lidx
      for (EDGE_ID i = 0; i < stored_num; i++){
          Lcycid[i] = i;
          Llen[i] = CC[i].len;
          Lupdated[i] = 0;
      }

      
      //printf("\nSorting Llen...");

      // Step 2(a): Sort Llen, Lcycid, Lupdated by Llen
      mergeSort_Llen(Llen, Lcycid, Lupdated, 0, stored_num - 1);

      //printf("\nInitializing C.Lidx...");
      // Step 2(b): Initialize C.Lidx
      for (EDGE_ID li = 0; li < stored_num; li++){
            CC[Lcycid[li]].Lidx = li;
      }


      // Define V that will store summations to be done
      EDGE_ID V_len = 0;
      EDGE_ID V_max_len = 10;
      min_update_V_H2* update_V = (min_update_V_H2*)malloc(V_max_len*sizeof(min_update_V_H2));

      EDGE_ID it_counter = 0;


      // Step 4: Loop for minimization
      while (1){

            if (!self->g_suppress_output){
                printf("\n\nIteration %d", it_counter++);
            }
            
            // Step 4(a): Find max diff
            EDGE_ID dm = 1;
            V_len = 0;
            for (EDGE_ID li = 0; li < stored_num; li++){

                  //printf("\nli %d", li);

                  if (Llen[li] < dm) break;

                  int add_flag = 0;

                  EDGE_ID ci = Lcycid[li];

                  for (EDGE_ID lj = li + 1; lj < stored_num; lj++){
                        
                        if (Llen[lj] < dm) break;
                        //printf("\nlj %d", lj);

                        EDGE_ID cj = Lcycid[lj];

                        EDGE_ID i = 0;
                        EDGE_ID j = 0;
                        EDGE_ID count = 0;

                        int quit_flag = 0;

                        while (i < CC[ci].len && j < CC[cj].len){

                              int compare;

                              if (CC[ci].boundary[i].key1 < CC[cj].boundary[j].key1){
                                  compare = 1;
                              }
                              else if (CC[ci].boundary[i].key1 > CC[cj].boundary[j].key1){
                                  compare = 0;
                              }
                              else{
                                  if (CC[ci].boundary[i].key2 < CC[cj].boundary[j].key2){
                                      compare = 1;
                                  }
                                  else if (CC[ci].boundary[i].key2 > CC[cj].boundary[j].key2){
                                      compare = 0;
                                  }
                                  else{
                                      compare = -1;
                                  }
                                  
                              }


                              
                              if (compare == 1){
                                    i++;
                                    count++;
                              }
                              else if (!compare){
                                    j++;
                                    count++;
                              }
                              else{
                                    i++;
                                    j++;
                              }

                              if (count > CC[ci].len - dm){
                                    quit_flag = 1;
                                    break;
                              }

                              
                        }

                        if (quit_flag) continue;

                        if (i < CC[ci].len){
                            count += CC[ci].len - i;
                        }
                        if (j < CC[cj].len){
                            count += CC[cj].len - j;
                        }
                        
                        if (count > CC[ci].len - dm){
                              continue;
                        }

                        //printf("\nREACHED HERE");
                        //getchar();


                        if (CC[ci].len - count > dm){
                              dm = CC[ci].len - count;
                              update_V[0].cycid = ci;
                              CC[ci].redw = cj;
                              V_len = 1;
                              add_flag = 1;
                        }
                        else if ((CC[ci].len - count == dm) && (!add_flag)){
                            
                              add_flag = 1;
                              update_V[V_len++].cycid = ci;
                              CC[ci].redw = cj;

                              if (V_len == V_max_len){
                                    V_max_len += 100;
                                    update_V = (min_update_V_H2*)realloc(update_V\
                                                            , V_max_len*sizeof(min_update_V_H2));
                              }

                        }
                        
                  }

            }


            if (!V_len){
                //printf("\nDiff 0. EXITING.");
                break;
            }

            if (!self->g_suppress_output){
                printf("\nmaxdiff %d, num pairs %d", dm, V_len);
            }


            #pragma omp parallel for schedule(static, 50) shared(update_V, CC)
            for (EDGE_ID vi = 0; vi < V_len; vi++){

                //if (vi % 1000 == 0){
                    //printf("\nreducing %d", vi);
                //}
                  
                EDGE_ID ci = update_V[vi].cycid;
                EDGE_ID cj = CC[ci].redw;

                //printf("\nReducing %d with %d", ci, cj);

                simplex* scratch = (simplex*)malloc((CC[ci].len + CC[cj].len)*sizeof(simplex));


                EDGE_ID edge;

                EDGE_ID i = 0;
                EDGE_ID j = 0;
                EDGE_ID count = 0;

                while ((i < CC[ci].len) && (j < CC[cj].len)){
                      
                      int compare;

                      if (CC[ci].boundary[i].key1 < CC[cj].boundary[j].key1){
                          compare = 1;
                      }
                      else if (CC[ci].boundary[i].key1 > CC[cj].boundary[j].key1){
                          compare = 0;
                      }
                      else{
                          if (CC[ci].boundary[i].key2 < CC[cj].boundary[j].key2){
                              compare = 1;
                          }
                          else if (CC[ci].boundary[i].key2 > CC[cj].boundary[j].key2){
                              compare = 0;
                          }
                          else{
                              compare = -1;
                          }
                          
                      }


                      if (compare==1){

                            scratch[count++] = CC[ci].boundary[i++];

                      }
                      else if (!compare){

                            scratch[count++] = CC[cj].boundary[j++];

                      }
                      else{

                            i++;
                            j++;

                      }

                }

                // update_in_cyc will have unique

                while (i < CC[ci].len){

                    // No change
                    scratch[count++] = CC[ci].boundary[i++];

                }

                while (j < CC[cj].len){

                    scratch[count++] = CC[cj].boundary[j++];

                }

                scratch = (simplex*)realloc(scratch, count*sizeof(simplex));


                update_V[vi].VV = scratch;
                update_V[vi].V_len = count;


            }


            // Update CC with new boundaries, update Llen and Lupdated
            for (EDGE_ID vi = 0; vi < V_len; vi++){
                  
                  EDGE_ID ci = update_V[vi].cycid;
                  EDGE_ID li = CC[ci].Lidx;

                  free(CC[ci].boundary);

                  CC[ci].boundary = update_V[vi].VV;
                  CC[ci].len = update_V[vi].V_len;

                  Llen[li] = update_V[vi].V_len;
                  Lupdated[li] = 1;

            }



            // Sort Llen, Lupdated, Lcycid
            mergeSort_Llen(Llen, Lcycid, Lupdated, 0, stored_num - 1);

            // Update C.Lidx
            for (EDGE_ID li = 0; li < stored_num; li++){
                  CC[Lcycid[li]].Lidx = li;
            }


      }



//////////////////////////////////////////////////////////////
      // TESTING
//////////////////////////////////////////////////////////////
//

      //printf("\nPress key to test");
      //getchar();
      

      //printf("\nTESTING...");

      //#pragma omp parallel for schedule(static) shared(stored_num\
      //                                                  , CC)
      //for (EDGE_ID ci = 0; ci < stored_num; ci++){

      //      EDGE_ID diff = 0;
      //      
      //      for (EDGE_ID cj = 0; cj < stored_num; cj++){
      //            
      //            if (cj == ci) continue;

      //            EDGE_ID j = 0;
      //            EDGE_ID k = 0;
      //            EDGE_ID count = 0;

      //            int quit_flag = 0;

      //            while ((j < CC[ci].len) && (k < CC[cj].len)){
      //                  if (CC[ci].boundary[j] < CC[cj].boundary[k]){
      //                        j++;
      //                        count++;
      //                  }
      //                  else if (CC[ci].boundary[j] > CC[cj].boundary[k]){
      //                        k++;
      //                        count++;
      //                  }
      //                  else{
      //                        j++;
      //                        k++;
      //                  }


      //            }

      //            if (j < CC[ci].len){
      //                  count += CC[ci].len - j;
      //            }

      //            if (k < CC[cj].len){
      //                  count += CC[cj].len - k;
      //            }

      //            
      //            if (count < CC[ci].len){
      //                  if (CC[ci].len - count > diff){
      //                      printf("\nImprovement possible by reducing ci, li %d, %d with cj, lj %d, %d!"\
      //                                                                              , ci, CC[ci].Lidx\
      //                                                                              , cj, CC[cj].Lidx);
      //                      printf("\nTESTING FAILED!!!");
      //                      getchar();
      //                  }
      //            }
      //            
      //      }
      //      
      //}

      //printf("\nTESTED OK.");


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


#ifdef STORE_LENGTHS_CYCLES
      FILE* fp1 = fopen(minimal_lens_filename, "w");
      for (EDGE_ID ci = 0; ci < stored_num; ci++){
           fprintf(fp1, "%d, ", CC[ci].len);
      }
      fclose(fp1);
#endif



      FILE* fp2 = fopen(filename, "w");

      for (EDGE_ID li = 0; li < stored_num; li++){

           EDGE_ID ci = Lcycid[li];
           //if (li < 15)
           //printf("\nlen %d", Llen[li]);

           for (EDGE_ID nn = 0; nn < CC[ci].len; nn++){
                
                 fprintf(fp2, "%d, %d, %d, ", self->g_edges_list[2*CC[ci].boundary[nn].key1]\
                                            , self->g_edges_list[2*CC[ci].boundary[nn].key1+1]\
                                            , CC[ci].boundary[nn].key2\
                                            );
                
           }
           fprintf(fp2, "\n");

           free(CC[ci].boundary);
      }

      fclose(fp2);

      //}





      free(CC);

      free(update_V);

      free(Lcycid);
      free(Llen);
      free(Lupdated);

}




//static PyMethodDef DoryMethods[] = {
//      
//      {"compute_PH", compute_PH, METH_VARARGS, "Compute PH"},
//      {NULL, NULL, 0, NULL}
//
//};
//
//static struct PyModuleDef dorymodule = {
//      PyModuleDef_HEAD_INIT,
//      "pydory", /* name of module*/
//      NULL, /* documentation */
//      -1, /* ??? */
//      DoryMethods
//};
//
//
//PyMODINIT_FUNC PyInit_pydory(void){
//      return PyModule_Create(&dorymodule);
//}
























