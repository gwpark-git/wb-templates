
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include "trajectory.h"
/* #include "read_file_condition.h" */
#include "file_IO.h"
#include "cell_list.h"
#include "numerics.h"
#include "connectivity.h"
#include <limits>

class RDIST : public CLIST
{
public:

  MATRIX** Rvec; // Rvec[i][j] will be the relative vector between i- and j-th particles
  MATRIX* Rsca; // Rsca[i](j) will be the relative distance between i- and j-th particles. The form is differ from the original.
  
  
  RDIST()
  {
    std::cout << "There is no empty constructor for RDIST class\n" << std::endl;
  }
  RDIST(COND& given_condition);
  virtual
    ~RDIST()
  {
    if(INITIALIZATION)
      {
        delete[] Rsca;
        for(long i=0; i<Np; i++)
          {
            delete[] Rvec[i];
          }
        delete[] Rvec;
      }
  }

  // member function case dependency
  // function point to measure minimum distance between pair of micelle

  // in order to use pointer-to-function rather than pointer-to-member_function, the measuring_minimum_distance have the type argument with itself: RDIST&. Note that the allocation is related with inline wrapping function defined in GEOMETRY namespace
  double
  (*measure_minimum_distance)(RDIST&, TRAJECTORY&, const long, const long, const long, long*);


  // the following are related with pointer-to-member_function.
  // in consequence of poingter-to-member_function, the fellowing have low readability which is not satisfactory by myself. Hence, the fellowing methods are deleted and left as recorded history to explain. (in future, these will be deleted completely) 

  // by some reason, the pointer-to-member_function is slightly differ from pointer-to-function.
  // it is not sure what exactly differ each other, but if this allocation step comapred with member function for POTENTIAL_SET, it uses typical allocation for pointer and uses it as name (as reference variable)
  // while the case on here is just works as pointer.
  // In consequence, allocation uses reference type rather than typical pointer type.
  // in future, the explanation will be revisited

  
};

namespace UTIL_ARR
{
  double
    get_minimum_image_k_from_x(double x, double k, double dimension);

  long
    get_index_minimum_abs(double *k, long N);
}


namespace GEOMETRY
{
  // Bound work and floating-point image-count precision explicitly. No images
  // are silently truncated when a requested domain exceeds these limits.
  // H-D08: these are supported arithmetic/work domains, not a license to
  // approximate an extreme lattice with fewer images. Extreme-aspect-ratio
  // counterexamples in test_legacy_physics.cpp must fail before cast or scan.
  constexpr long max_shear_gradient_images = 4096;
  constexpr double max_periodic_winding = 1125899906842624.; // 2^50

  inline void validate_shear_box(double lx,double ly,double lz)
  {
    if(!std::isfinite(lx) || !std::isfinite(ly) || !std::isfinite(lz) ||
       lx<=0. || ly<=0. || lz<=0.)
      throw std::domain_error("shear box lengths must be positive and finite");
    const double search_bound=std::hypot(lx/ly,1.)+2.;
    if(!std::isfinite(search_bound) || search_bound>max_shear_gradient_images)
      throw std::domain_error("legacy shear lattice exceeds 4096 gradient images; unsupported box aspect ratio");
  }

  inline double checked_periodic_quotient(double coordinate,double length)
  {
    if(!std::isfinite(coordinate) || !std::isfinite(length) || length<=0.)
      throw std::domain_error("periodic coordinates and positive box lengths must be finite");
    const double quotient=coordinate/length;
    if(!std::isfinite(quotient) || std::abs(quotient)>max_periodic_winding)
      throw std::domain_error("periodic winding exceeds the supported 2^50 image-count domain");
    return quotient;
  }

  inline double wrapped_periodic_coordinate(double coordinate,double length)
  {
    const double quotient=checked_periodic_quotient(coordinate,length);
    if(coordinate>=0. && coordinate<length) return coordinate;
    const double wrapped=coordinate-length*std::floor(quotient);
    if(!std::isfinite(wrapped) || wrapped<0. || wrapped>length)
      throw std::domain_error("periodic remapping is outside finite representable coordinate domain");
    // The two endpoints identify the same periodic point; retain [0,L).
    return wrapped==length ? 0. : wrapped;
  }

  long
    compute_RDIST_particle(RDIST& R_boost, const long index_particle,
			   TRAJECTORY& TRAJ, long index_t);
  long compute_bonded_RDIST_particle(RDIST& distances, CONNECTIVITY& bonds,
                                    long particle, TRAJECTORY& trajectory, long time);
  // the original set
  double
    get_minimum_distance(TRAJECTORY& TRAJ, long index_t,
			 long index_i, long index_j,
			 MATRIX& given_vec);
  double
    get_minimum_distance_for_particle(TRAJECTORY& TRAJ, long index_t,
				      long index_particle,
				      MATRIX& R_minimum_boost_particle,
				      MATRIX** R_minimum_vec_boost); 
  double
    return_minimum_distance(TRAJECTORY& TRAJ, long index_t,
			    long index_i, long index_j);
  double
    minimum_image_convention(TRAJECTORY& TRAJ, long target_t);
  double
    minimum_image_convention_loop(TRAJECTORY& TRAJ, long target_t);
  double
  apply_shear_boundary_condition(TRAJECTORY& TRAJ, long target_t,
				   const long shear_axis, const long shear_grad_axis, const double shift_factor);

  double
  set_box_shift_factor(double &shear_PBC_shift,
                       double const maximum_displacement,
                       double const box_dimension_shear_axis);
  
  double
  set_box_shift_factor(double &shear_PBC_shift,
                       RDIST& R_boost,
                       double const maximum_displacement,
                       double const box_dimension_shear_axis);
  
  double
  apply_step_shear(TRAJECTORY& TRAJ, long target_t,
                   const long shear_axis, const long shear_grad_axis,
                   const double gamma_0, const double box_dimension);
  
  long
    get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, long index_t,
				    long given_index, long target_index,
				    MATRIX& given_vec);
  long
    get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, long index_t,
				    long given_index, long target_index,
				    MATRIX& given_vec);
  long
    get_minimum_distance_rel_vector(TRAJECTORY& TRAJ, long index_t,
				    long given_index, long target_index,
				    MATRIX& given_vec);
  double
    get_simple_distance(TRAJECTORY& TRAJ, long index_t,
			long index_i, long i);

  double
    get_simple_distance(TRAJECTORY& TRAJ, long index_t,
			long index_i, long i);


  // inlined functions
  double
    measure_minimum_distance_default(RDIST& R_boost,
				     TRAJECTORY& TRAJ, const long index_t,
				     const long index_particle, const long index_target,
				     long* beyond_box_check);
  double
    measure_minimum_distance_cell_list(RDIST& R_boost,
				       TRAJECTORY& TRAJ, const long index_t,
				       const long index_particle, const long index_target,
				       long* beyond_box_check);
  double
    measure_minimum_distance_simple_shear_fixed_axis(RDIST& R_boost,
						     TRAJECTORY& TRAJ, const long index_t,
						     const long index_particle, const long index_target,
						     long* beyond_box_check);

  
  long
    get_minimum_distance_rel_vector_cell_list(TRAJECTORY& TRAJ, long index_t,
					      long given_index, long target_index,
					      MATRIX& given_vec,
					      long* beyond_box_check);
  double
    get_minimum_distance_cell_list(TRAJECTORY& TRAJ, long index_t,
				   long index_i, long index_j,
				   MATRIX& given_vec,
				   long* beyond_box_check);
  double
    get_minimum_distance_rel_vector_simple_shear_fixed_axis(TRAJECTORY& TRAJ, const long index_t,
							    const long index_i, const long index_j,
							    MATRIX& given_vec,
							    const double map_to_central_box_image);
  double
    get_minimum_distance_simple_shear_fixed_axis(TRAJECTORY& TRAJ, long index_t,
						 long index_i, long index_j,
						 MATRIX& given_vec,
						 const double map_to_central_box_image);
  long
    minimum_image_convention_particle(TRAJECTORY& TRAJ, long target_t,
				      long index_particle);
  double
    get_minimum_distance_rel_vector_equilibrium_test(TRAJECTORY& TRAJ, const long index_t,
						     const long index_i, const long index_j,
						     MATRIX& given_vec);
  double
    get_minimum_distance_equilibrium_test(TRAJECTORY& TRAJ, long index_t,
					  long index_i, long index_j,
					  MATRIX& given_vec);
  long
    minimum_image_convention_particle(TRAJECTORY& TRAJ, long target_t,
				      long index_particle);
  
} // namespace GEOMETRY

// the following are related with inlined function
// it is extracted from above for readability
inline double
UTIL_ARR::
get_minimum_image_k_from_x(double x, double k, double dimension)
{
  double kd[3] = {k-dimension - x, k - x, k+dimension - x};
  double re= kd[get_index_minimum_abs(kd, 3)] + x;
  return re;
}

inline long
UTIL_ARR::
get_index_minimum_abs(double *k, long N)
{
  long re = 0;
  for(long i=1; i<N; i++)
    {
      if (fabs(k[i]) < fabs(k[re]))
        re = i;
    }
  return re;
}

inline long
GEOMETRY::
get_minimum_distance_rel_vector_cell_list(TRAJECTORY& TRAJ, long index_t,
					  long given_index, long target_index,
					  MATRIX& given_vec,
					  long* beyond_box_check)
{
  // direction convention:
  // +: direction to the given bead
  // -: direction to the target bead

  for(long k=0; k<TRAJ.N_dimension; k++)
    {
      // the scheme is differ from the default minimum distance definition
      // since there are beyond_box_check will identify a pair of particles have different beloning for PBC box
      // the cell list already have information related with a pair of micelles wheather the micelles are belong to different PBC box or not
      // in equilibrium simulation, beyond_box_check is constant, so there are no meaning to minimize distance for given all possible images of the pair of micelles.
      double PBC_coord_target = TRAJ(index_t, target_index, k) + (double)beyond_box_check[k]*TRAJ.box_dimension[k];
      given_vec(k) = PBC_coord_target - TRAJ(index_t, given_index, k);
    }
  return 0;
}
  
inline double
GEOMETRY::
get_minimum_distance_cell_list(TRAJECTORY& TRAJ, long index_t,
			       long index_i, long index_j,
			       MATRIX& given_vec,
			       long* beyond_box_check)
{
  GEOMETRY::get_minimum_distance_rel_vector_cell_list(TRAJ, index_t, index_i, index_j, given_vec, beyond_box_check);

  return given_vec.norm();
  /* return cblas_dnrm2(given_vec.size, given_vec.data, 1); // 1 denote increment for x */
  /* the cblas_dnrm2 is slower than inlined norm */
}


inline double
GEOMETRY::
get_minimum_distance_rel_vector_equilibrium_test(TRAJECTORY& TRAJ, const long index_t,
						 const long index_i, const long index_j,
						 MATRIX& given_vec)
{
  double min_norm = std::numeric_limits<double>::infinity();
  long min_sf_0 = 0, min_sf_1 = 0, min_sf_2 = 0;
  for(long sf_0=-1; sf_0<=+1; sf_0++)
    {
      for(long sf_1=-1; sf_1<=+1; sf_1++)
        {
          for(long sf_2=-1; sf_2<=+1; sf_2++)
            {
              double rpos_ij_0 = TRAJ(index_t, index_j, 0) + sf_0*TRAJ.box_dimension[0] - TRAJ(index_t, index_i, 0);
              double rpos_ij_1 = TRAJ(index_t, index_j, 1) + sf_1*TRAJ.box_dimension[1] - TRAJ(index_t, index_i, 1);
              double rpos_ij_2 = TRAJ(index_t, index_j, 2) + sf_2*TRAJ.box_dimension[2] - TRAJ(index_t, index_i, 2);
              double norm_pos = sqrt(rpos_ij_0*rpos_ij_0 + rpos_ij_1*rpos_ij_1 + rpos_ij_2*rpos_ij_2);
              if(norm_pos < min_norm)
                {
                  min_sf_0 = sf_0;
                  min_sf_1 = sf_1;
                  min_sf_2 = sf_2;
                  min_norm = norm_pos;
                }
            } // it is already assumed for 3-dimensional case
        }
    }
  given_vec(0) = TRAJ(index_t, index_j, 0) + min_sf_0*TRAJ.box_dimension[0] - TRAJ(index_t, index_i, 0);
  given_vec(1) = TRAJ(index_t, index_j, 1) + min_sf_1*TRAJ.box_dimension[1] - TRAJ(index_t, index_i, 1);
  given_vec(2) = TRAJ(index_t, index_j, 2) + min_sf_2*TRAJ.box_dimension[2] - TRAJ(index_t, index_i, 2);
  return min_norm;
}



inline double
GEOMETRY::
get_minimum_distance_equilibrium_test(TRAJECTORY& TRAJ, const long index_t,
				      const long index_i, const long index_j,
				      MATRIX& given_vec)
{
  return GEOMETRY::get_minimum_distance_rel_vector_equilibrium_test(TRAJ, index_t, index_i, index_j, given_vec);
}

inline double
GEOMETRY::
get_minimum_distance_rel_vector_simple_shear_fixed_axis(TRAJECTORY& TRAJ, const long index_t,
							const long index_i, const long index_j,
							MATRIX& given_vec,
							const double map_to_central_box_image)
{
  if (TRAJ.N_dimension != 3)
    throw std::invalid_argument("legacy shear distance requires three dimensions");
  const double lx=TRAJ.box_dimension[0], ly=TRAJ.box_dimension[1], lz=TRAJ.box_dimension[2];
  validate_shear_box(lx,ly,lz);
  if(!std::isfinite(map_to_central_box_image))
    throw std::domain_error("shear image displacement must be finite");
  const double dx=TRAJ(index_t,index_j,0)-TRAJ(index_t,index_i,0);
  const double dy=TRAJ(index_t,index_j,1)-TRAJ(index_t,index_i,1);
  const double dz=TRAJ(index_t,index_j,2)-TRAJ(index_t,index_i,2);
  checked_periodic_quotient(dx,lx);
  checked_periodic_quotient(dy,ly);
  checked_periodic_quotient(dz,lz);
  const double z=dz-lz*std::floor(dz/lz+.5);
  const double central=-std::floor(dy/ly+.5);
  auto x_image = [&](double winding) {
    const double x=dx+winding*map_to_central_box_image;
    checked_periodic_quotient(x,lx);
    return x-lx*std::floor(x/lx+.5);
  };
  double best_x=x_image(central), best_y=dy+central*ly;
  double best=best_x*best_x+best_y*best_y;
  // Preserve ordinary-case arithmetic, but reject unsafe or excessive work
  // before floating-to-integer conversion or the bounded lattice loop.
  // A candidate with a larger |y| than this existing 2D norm cannot improve it.
  // This bounds a complete oblique-lattice search; nearest y alone, or only
  // {-1,0,1} images, can miss the closest pair in an elongated/sheared box.
  const double bound=std::sqrt(best);
  const double lower_value=std::ceil((-bound-dy)/ly);
  const double upper_value=std::floor((bound-dy)/ly);
  const double long_lower=static_cast<double>(std::numeric_limits<long>::min());
  const double long_upper_exclusive=-long_lower;
  if(!std::isfinite(best) || !std::isfinite(lower_value) || !std::isfinite(upper_value) ||
     lower_value<long_lower || upper_value>=long_upper_exclusive ||
     upper_value-lower_value+1.>max_shear_gradient_images || lower_value>upper_value)
    throw std::domain_error("legacy shear lattice exceeds finite representable 4096-image search domain");
  const long lower=static_cast<long>(lower_value), upper=static_cast<long>(upper_value);
  for(long winding=lower;; ++winding) {
    const double x=x_image(static_cast<double>(winding)), y=dy+winding*ly;
    const double norm=x*x+y*y;
    if(!std::isfinite(norm)) throw std::domain_error("nonfinite shear lattice candidate distance");
    if (norm < best) { best=norm; best_x=x; best_y=y; }
    if(winding==upper) break; // do not increment a terminal LONG_MAX value
  }
  given_vec(0)=best_x; given_vec(1)=best_y; given_vec(2)=z;
  const double distance=std::sqrt(best+z*z);
  if(!std::isfinite(distance)) throw std::domain_error("nonfinite shear lattice distance");
  return distance;
}

inline double
GEOMETRY::
get_minimum_distance_simple_shear_fixed_axis(TRAJECTORY& TRAJ, long index_t,
					     long index_i, long index_j,
					     MATRIX& given_vec,
					     const double map_to_central_box_image)
{
  double distance = GEOMETRY::get_minimum_distance_rel_vector_simple_shear_fixed_axis(TRAJ, index_t, index_i, index_j, given_vec, map_to_central_box_image);
  return distance;
}

  
inline long
GEOMETRY::
minimum_image_convention_particle(TRAJECTORY& TRAJ, long target_t,
				  long index_particle)
{
  for(long k=0; k<TRAJ.N_dimension; k++)
    {
      double& coordinate=TRAJ(target_t,index_particle,k);
      coordinate=wrapped_periodic_coordinate(coordinate,TRAJ.box_dimension[k]);
    }
  return 0;
}

inline double
GEOMETRY::
get_minimum_distance(TRAJECTORY& TRAJ, long index_t,
		     long index_i, long index_j,
		     MATRIX& given_vec)
{
  // this comentation is related with the assumes that
  // the given dimensionality is properly working
  // therefore, it is assumed that there is no mistake
  // and no internal check for the dimensionality.
  GEOMETRY::get_minimum_distance_rel_vector(TRAJ, index_t, index_i, index_j, given_vec);
  return given_vec.norm();
}

inline double
GEOMETRY::
return_minimum_distance(TRAJECTORY& TRAJ, long index_t,
			long index_i, long index_j)
{
  MATRIX rel_vec(TRAJ.N_dimension, 1, 0.);
  GEOMETRY::get_minimum_distance(TRAJ, index_t, index_i, index_j, rel_vec);
  return rel_vec.norm();
}

inline long
GEOMETRY::
get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, long index_t,
				long given_index, long target_index,
				MATRIX& given_vec)
{
  for(long k=0; k<TRAJ.N_dimension; k++)
    {
      given_vec(k) = UTIL_ARR::get_minimum_image_k_from_x(TRAJ(index_t, given_index, k), TRAJ(index_t, target_index, k), TRAJ.box_dimension[k]);
    }
  return 0;
}

inline long
GEOMETRY::
get_minimum_distance_rel_vector(TRAJECTORY& TRAJ, long index_t,
				long given_index, long target_index,
				MATRIX& given_vec)
{
  GEOMETRY::get_minimum_distance_pos_vector(TRAJ, index_t, given_index, target_index, given_vec);
  for(long k=0; k<TRAJ.N_dimension; k++)
    {
      // direction convention:
      // +: direction to the given bead
      // -: direction to the target bead
      given_vec(k) -= TRAJ(index_t, given_index, k);
    }
  return 0;
}

inline double
GEOMETRY::
get_simple_distance(TRAJECTORY& TRAJ, long index_t,
		    long index_i, long i)
{
  double distance = 0.;
  for(long k=0; k<TRAJ.N_dimension; k++)
    {
      distance += pow(TRAJ(index_t, index_i, k) - TRAJ(index_t, i, k), 2.0);
    }
  return sqrt(distance);

}

inline double
GEOMETRY::
measure_minimum_distance_default(RDIST& R_boost,
				 TRAJECTORY& TRAJ, const long index_t,
				 const long index_particle, const long index_target,
				 long* beyond_box_check)
{
  return GEOMETRY::get_minimum_distance(TRAJ, index_t, index_particle, index_target, R_boost.Rvec[index_particle][index_target]);
  // the following function is reference (and very slow).
  // to be on safe side, it still alive (just commented)
  // return GEOMETRY::get_minimum_distance_equilibrium_test(TRAJ, index_t, index_particle, index_target, R_boost.Rvec[index_particle][index_target]);
}

inline double
GEOMETRY::
measure_minimum_distance_cell_list(RDIST& R_boost,
				   TRAJECTORY& TRAJ, const long index_t,
				   const long index_particle, const long index_target,
				   long* beyond_box_check)
{
  return GEOMETRY::get_minimum_distance_cell_list(TRAJ, index_t, index_particle, index_target, R_boost.Rvec[index_particle][index_target], beyond_box_check);
}

inline double
GEOMETRY::
measure_minimum_distance_simple_shear_fixed_axis(RDIST& R_boost,
						 TRAJECTORY& TRAJ, const long index_t,
						 const long index_particle, const long index_target,
						 long* beyond_box_check)
{
  return GEOMETRY::get_minimum_distance_simple_shear_fixed_axis(TRAJ, index_t, index_particle, index_target, R_boost.Rvec[index_particle][index_target], R_boost.map_to_central_box_image);
}

#endif
