
#include "geometry.h"

// H-D07/H-D08 (docs/DEVELOPMENT_HISTORY.org): a physical bridge survives leaving
// the repulsive cell neighborhood. Refresh its geometry from bond adjacency,
// even when the spatial candidate pass did not visit that pair. The independent
// non-neighbor-bond force fixture in test_legacy_physics.cpp guards this seam.
long GEOMETRY::compute_bonded_RDIST_particle(RDIST& distances, CONNECTIVITY& bonds,
                                            long particle, TRAJECTORY& trajectory, long time)
{
  for(long h=1; h<bonds.TOKEN[particle]; ++h) {
    const long target=static_cast<long>(bonds.HASH[particle](h));
    MATRIX& vector=distances.Rvec[particle][target];
    distances.Rsca[particle](target)=(distances.SIMPLE_SHEAR || distances.STEP_SHEAR)
      ? get_minimum_distance_simple_shear_fixed_axis(trajectory,time,particle,target,
                                                    vector,distances.map_to_central_box_image)
      : get_minimum_distance(trajectory,time,particle,target,vector);
  }
  return 0;
}

long
GEOMETRY::
compute_RDIST_particle(RDIST& R_boost, const long index_particle,
		       TRAJECTORY& TRAJ, long index_t)
{
  long cell_index_particle = R_boost.cell_index[index_particle];
  for(long k=0; k<R_boost.N_neighbor_cells; k++)
    {
      long cell_index_neighbor = R_boost.NEIGHBOR_CELLS[cell_index_particle][k];
      for(long p=0; p<R_boost.TOKEN[cell_index_neighbor]; p++)
        {
          long index_target = R_boost(cell_index_neighbor, p);
          double distance = R_boost.measure_minimum_distance(R_boost, TRAJ, index_t, index_particle, index_target, R_boost.BEYOND_BOX[cell_index_particle][k]);
          R_boost.Rsca[index_particle](index_target) = distance;
        }
    }
  return 0;
}


RDIST::
RDIST(COND& given_condition)
  : CLIST(given_condition)
{
  Rvec = new MATRIX* [Np];
  Rsca = new MATRIX [Np];
    
  for(long i=0; i<Np; i++)
    {
      // since the space complexity is not the matter for our simulation (at this moment),
      // the Rvec have Np*Np array that is much larger when we used cell-list approaches
      Rvec[i] = new MATRIX [Np];
      for(long j=0; j<Np; j++)
        {
          Rvec[i][j].initial(N_dimension, 1, 0.);
        }
      Rsca[i].initial(Np, 1, 0.);
    }
  // if(SIMPLE_SHEAR) // it already declare in CLIST class
  if(SIMPLE_SHEAR || STEP_SHEAR) 
    {
      // note that both of the cases are applicable based on the following functions
      // hence the name of functions are subject to change from SIMPLE_SHEAR to SHEAR
      if(CELL_LIST_BOOST)
        {
          printf("ERR: CELL LIST with SIMPLE SHEAR FLOW is not implemented. \n");
        }
      if (!(shear_axis == 0 && shear_grad_axis == 1) && !CELL_LIST_BOOST)
        {
          printf("ERR: SIMPLE SHEAR FLOW with different direction without CELL LIST is not implemented\n");
        }
      // note again, the cell list implementation is not work at this moment

      measure_minimum_distance = GEOMETRY::measure_minimum_distance_simple_shear_fixed_axis;
      // check
      printf("DISTANCE is using the function GEOMETRY::measure_minimum_distance_simple_shear_fixed_axis\n");
      
    }
  else
    {
      if(CELL_LIST_BOOST)
        measure_minimum_distance = GEOMETRY::measure_minimum_distance_cell_list;
      else
        measure_minimum_distance = GEOMETRY::measure_minimum_distance_default;
    }
}

double
GEOMETRY::
get_minimum_distance_for_particle(TRAJECTORY& TRAJ, long index_t,
				  long index_particle,
				  MATRIX& R_minimum_boost_particle,
				  MATRIX** R_minimum_vec_boost)
{
  // this function generating minimum relative vector into R_minimum_vec_boost,
  // and distance of it into R_minimum_boost.
  for(long i=0; i<TRAJ.Np; i++)
    {
      R_minimum_boost_particle(i) = GEOMETRY::get_minimum_distance(TRAJ, index_t, index_particle, i, R_minimum_vec_boost[index_particle][i]);
    }
  return 0;
}

double
GEOMETRY::
minimum_image_convention(TRAJECTORY& TRAJ, long target_t)
{
  return minimum_image_convention_loop(TRAJ,target_t);
}

double
GEOMETRY::
minimum_image_convention_loop(TRAJECTORY& TRAJ, long target_t)
{
  double time_st = omp_get_wtime();
  for(long i=0; i<TRAJ.Np; i++)
    {
      for(long k=0; k<TRAJ.N_dimension; k++)
        {
          double coord = TRAJ(target_t, i, k);
          TRAJ(target_t,i,k)=wrapped_periodic_coordinate(coord,TRAJ.box_dimension[k]);
        }
    }
  return omp_get_wtime() - time_st;

}

double
GEOMETRY::
apply_shear_boundary_condition(TRAJECTORY& TRAJ, long target_t,
                               const long shear_axis, const long shear_grad_axis, const double shift_factor)
{
  double time_st = omp_get_wtime();
  if(shear_axis<0 || shear_grad_axis<0 || shear_axis>=TRAJ.N_dimension ||
     shear_grad_axis>=TRAJ.N_dimension || shear_axis==shear_grad_axis)
    throw std::invalid_argument("invalid legacy shear axes");
  if(!std::isfinite(shift_factor))
    throw std::domain_error("shear boundary shift must be finite");
  for (long i=0; i<TRAJ.Np; i++)
    {
      double coord = TRAJ(target_t, i, shear_grad_axis);
      const double quotient=checked_periodic_quotient(coord,TRAJ.box_dimension[shear_grad_axis]);
      checked_periodic_quotient(TRAJ(target_t,i,shear_axis),TRAJ.box_dimension[shear_axis]);
      if(coord < 0 || coord >= TRAJ.box_dimension[shear_grad_axis])
        {
          // Count all signed gradient crossings before ordinary wrapping. A
          // single +/- box correction loses the corresponding flow displacement.
          // Use the destination frame's shear shift, supplied by the driver.
          const double winding=std::floor(quotient);
          const double shifted=TRAJ(target_t,i,shear_axis)-winding*shift_factor;
          checked_periodic_quotient(shifted,TRAJ.box_dimension[shear_axis]);
          TRAJ(target_t,i,shear_axis)=shifted;
          
        }
    }
  return omp_get_wtime() - time_st;
}

double
GEOMETRY::
apply_step_shear(TRAJECTORY& TRAJ, long target_t,
                       const long shear_axis, const long shear_grad_axis,
                       const double gamma_0, const double box_dimension)
{
  double time_st = omp_get_wtime();
  printf("APPLY STEP DEFORMATION (STEP_SHEAR is on)\n");
  for(long i=0; i<TRAJ.Np; i++)
    {
      // apply initial step deformation
      // note that shift in x-direction is given by: y*gamma_0 <= 1
      TRAJ(target_t, i, shear_axis) += gamma_0*TRAJ(target_t, i, shear_grad_axis);
    }
  // // the following is not necessary at this moment since the initial box deformation only apply for shear direction (shear gradient direction still in minimum image convention)
  // GEOMETRY::
  //   apply_shear_boundary_condition(TRAJ, target_t, shear_axis, shear_grad_axis, gamma_0*box_dimension);

  // following will apply minimum image convention (for shear direction)
  GEOMETRY::
    minimum_image_convention_loop(TRAJ, target_t);
  return omp_get_wtime() - time_st;
}

double
GEOMETRY::
set_box_shift_factor(double &shear_PBC_shift,
                     double const maximum_displacement,
                     double const box_dimension_shear_axis)
{
  double time_st = omp_get_wtime();
  shear_PBC_shift = fmod(maximum_displacement, box_dimension_shear_axis);
  return omp_get_wtime() - time_st;
}

double
GEOMETRY::
set_box_shift_factor(double &shear_PBC_shift,
                     RDIST& R_boost,
                     double const maximum_displacement,
                     double const box_dimension_shear_axis)
{
  /* // Original Form
     VAR.shear_PBC_shift = fmod(VAR.gamma_0*TRAJ.box_dimension[VAR.shear_grad_axis], TRAJ.box_dimension[VAR.shear_axis]);
     // VAR.shear_PBC_shift = VAR.gamma_0*TRAJ.box_dimension[VAR.shear_grad_axis];
     R_boost.map_to_central_box_image = fmod(VAR.shear_PBC_shift, TRAJ.box_dimension[VAR.shear_axis]);
     long central_standard = (long)(2*R_boost.map_to_central_box_image/TRAJ.box_dimension[VAR.shear_axis]);
     R_boost.map_to_central_box_image -= TRAJ.box_dimension[VAR.shear_axis]*(double)central_standard;
  */
  double time_st = omp_get_wtime();
  shear_PBC_shift = fmod(maximum_displacement, box_dimension_shear_axis);
  R_boost.map_to_central_box_image = fmod(shear_PBC_shift, box_dimension_shear_axis);
  long central_standard = (long)(2*R_boost.map_to_central_box_image/box_dimension_shear_axis);
  R_boost.map_to_central_box_image -= box_dimension_shear_axis*(double)central_standard;
  return omp_get_wtime() - time_st;
}
