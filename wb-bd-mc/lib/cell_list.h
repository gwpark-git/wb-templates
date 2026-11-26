
#ifndef CELL_LIST_H
#define CELL_LIST_H

#define TRUE 1
#define FALSE 0
#define BIT 64

#include <iostream>
#include <math.h>
#include "matrix.h" 
#include "trajectory.h" // the converting interface will be implemented
/* #include "read_file_condition.h" */
#include "file_IO.h"

class CLIST
{
 private:
 public:

  bool CELL_LIST_BOOST;
  long N_dimension; // spatial dimension
  double *box_dimension; // box dimension for one axis
  bool INITIALIZATION = false;
  /*
    Based on the given critical radius, rc, the real length for cell is computed based on proper divisor of whole box. For instance, if the given box length is 10 and the critical length scale is 0.9, then the length for cell becomes 1.0, which is the slight larger than the critical length scale.
  */
  double cut_off_radius; // cut_off_radius
  double *cell_length; // length for cell
  long *N_div; // number of cells per axis
  long N_cells; // number of cells
  long N_neighbor_cells; // including itself
  long MAX_IN_CELL; // it state how many particles are allowed for one cell. Set as all the number of particles as default.
  
  long **CELL;
  long *cell_index;
  long **NEIGHBOR_CELLS;


  
  long ***BEYOND_BOX;
  long *TOKEN;
  long Np;

  // following are related with mechanical perturbations
  bool SIMPLE_SHEAR;
  bool STEP_SHEAR;
  long shear_axis;
  long shear_grad_axis;
  double map_to_central_box_image;
  
  // when simple shear is implemented
  // long N_offset;
  // long ***NEIGHBOR_CELLS_OFFSET; // the first dimensionality suggest that up and bottom check (right and bottom in general description)
  

  
  // mapping function
  long
    index_vec2sca(const long* index_vec, long& index_sca);
  long
    index_sca2vec(const long& index_sca, long* index_vec);
  long
    get_neighbor_cell_list(const long& index_sca,
			   long* index_neighbor_cells,
			   long* self_index_vec_boost,
			   long* sf_vec_boost);
  long
    allocate_index_neighbor_cell_list();
  long
    identify_cell_from_given_position(TRAJECTORY& TRAJ, long index_t_now,
				      long index_particle, long *index_vec_boost);
  long
    allocate_cells_from_positions(TRAJECTORY& TRAJ, long index_t_now,
				  long *index_vec_boost);
  long
    identify_cell_from_given_position(TRAJECTORY_HDF5& TRAJ, long index_t_now,
				      long index_particle, long *index_vec_boost);
  long
    allocate_cells_from_positions(TRAJECTORY_HDF5& TRAJ, long index_t_now,
				  long *index_vec_boost);

    
  // operator overloading
  long&
    operator()
    (long i, long j)
    {
      return CELL[i][j];
    }
  
  // constructor
  CLIST(){}
  CLIST(COND& given_condition);

  // destructor
  virtual
    ~CLIST()
    {
      if(INITIALIZATION)
        {
          delete[] TOKEN;
	  delete[] box_dimension;
	  delete[] cell_length;
	  // delete[] N_cells;
	  delete[] N_div;
	  
          for(long i=0; i<N_cells; i++)
            {
              delete[] CELL[i];
              delete[] NEIGHBOR_CELLS[i];
              for(long j=0; j<N_neighbor_cells; j++)
                delete[] BEYOND_BOX[i][j];
              delete[] BEYOND_BOX[i];
            }
          delete[] CELL;
          delete[] NEIGHBOR_CELLS;
          delete[] BEYOND_BOX;
          delete[] cell_index;

          // if(SIMPLE_SHEAR)
          //   {
          //     for(long i=0; i<2; i++)
          //       {
          //         for(long j=0; j<N_cells + 1; j++)
          //           {
          //             for(long k=0; k<N_neighbor_cells + N_offset; k++)
          //               {
          //                 delete[] NEIGHBOR_CELLS_OFFSET[i][j][k];
          //               }
          //             delete[] NEIGHBOR_CELLS_OFFSET[i][j];
          //           }
          //         delete[] NEIGHBOR_CELLS_OFFSET[i];
          //       }
          //   }
        }
    }
};

namespace UTILITY
{
  long
  index_vec2sca(const long* index_vec, long& index_sca,
                const long N_dimension, const long *N_div);
  long
  index_sca2vec(const long& index_sca, long* index_vec,
                const long N_dimension, const long *N_div);
}

#endif
