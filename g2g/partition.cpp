#include <iostream>
#include <limits>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include "common.h"
#include "init.h"
#include "matrix.h"
#include "partition.h"
#include "timer.h"
#include "mkl.h"
using namespace std;

namespace G2G {
Partition partition;

ostream& operator<<(ostream& io, const Timers& t) {
#ifdef TIMINGS
    ostringstream ss;
    ss << "memcpys: " << t.memcpy << "trmms: " << t.trmms << "density_calcs: " << t.density_calcs << "rmm: " << t.rmm << " density: " 
       << t.density << " pot: " << t.pot << " forces: " << t.forces << " resto: " << t.resto << " functions: " << t.functions;
    io << ss.str() << endl;
#endif
  return io;
}

/********************
 * PointGroup
 ********************/

template<class scalar_type>
void PointGroup<scalar_type>::get_rmm_input(HostMatrix<scalar_type>& rmm_input) const {
  rmm_input.zero();
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
          rmm_input(ii, jj) = (scalar_type)fortran_vars.rmm_input_ndens1.data[big_index];
          rmm_input(jj, ii) = rmm_input(ii, jj);
        }
      }
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output(const HostMatrix<scalar_type>& rmm_output) const {
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
          fortran_vars.rmm_output(big_index) += (double)rmm_output(ii, jj);
        }
      }
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::compute_nucleii_maps(void)
{
  if (total_functions_simple() != 0) {
    func2global_nuc.resize(total_functions_simple());
    for (uint i = 0; i < total_functions_simple(); i++) {
      func2global_nuc(i) = fortran_vars.nucleii(local2global_func[i]) - 1;
    }

    func2local_nuc.resize(total_functions());
    uint ii = 0;
    for (uint i = 0; i < total_functions_simple(); i++) {
      uint global_atom = func2global_nuc(i);
      uint local_atom = std::distance(local2global_nuc.begin(), std::find(local2global_nuc.begin(), local2global_nuc.end(), global_atom));
      uint inc = small_function_type(i);
      for (uint k = 0; k < inc; k++, ii++) func2local_nuc(ii) = local_atom;
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::add_point(const Point& p) {
  points.push_back(p);
  number_of_points++;
}

#define EXP_PREFACTOR 1.01057089636005 // (2 * pow(4, 1/3.0)) / M_PI

template<class scalar_type>
bool PointGroup<scalar_type>::is_significative(FunctionType type, double exponent, double coeff, double d2) {
  switch(type) {
    case FUNCTION_S:
      return (exponent * d2 < max_function_exponent-log(pow((2.*exponent/M_PI),3))/4);
    break;
    default:
    {
      double x = 1;
      double delta;
      double e = 0.1;
      double factor = pow((2.0*exponent/M_PI),3);
      factor = sqrt(factor*4.0*exponent) ;
      double norm = (type == FUNCTION_P ? sqrt(factor) : abs(factor)) ;
      do {
        double div = (type == FUNCTION_P ? log(x) : 2 * log(x));
        double x1 = sqrt((max_function_exponent - log(norm) + div) / exponent);
        delta = abs(x-x1);
        x = x1;
      } while (delta > e);
      return (sqrt(d2) < x);
    }
    break;
  }
}

template<class scalar_type>
long long PointGroup<scalar_type>::cost() const {
    long long np = number_of_points, gm = total_functions();
    static const long long MIN_COST = 0;
    // Primer termino: multiplicaciones de matrices.
    // Segundo termino: Calcular rmm
    // Tercer termino: overhead 
    return (10 * np * gm * gm) / 2 + np * 2 * gm * gm + MIN_COST;
}
template<class scalar_type>
bool PointGroup<scalar_type>::operator<(const PointGroup<scalar_type>& T) const{
    return cost() < T.cost();
}
template<class scalar_type>
int PointGroup<scalar_type>::pool_elements() const {
    int t = total_functions(), n = number_of_points;
    return t * n;
}
template<class scalar_type>
size_t PointGroup<scalar_type>::size_in_gpu() const
{
    uint total_cost=0;
    uint single_matrix_cost = COALESCED_DIMENSION(number_of_points) * total_functions();

    total_cost += single_matrix_cost;       //1 scalar_type functions
    if (fortran_vars.do_forces || fortran_vars.gga)
      total_cost += (single_matrix_cost*4); //4 vec_type gradient
    if (fortran_vars.gga)
      total_cost+= (single_matrix_cost*8);  //2*4 vec_type hessian
    return total_cost*sizeof(scalar_type);  // size in bytes according to precision
}
template<class scalar_type>
PointGroup<scalar_type>::~PointGroup<scalar_type>()
{

#if !CPU_KERNELS
    if(inGlobal)
    {
      globalMemoryPool::dealloc(size_in_gpu());
      function_values.deallocate();
      gradient_values.deallocate();
      hessian_values.deallocate();
    }

#endif
}

void Partition::solve(Timers& timers, bool compute_rmm,bool lda,bool compute_forces, 
                      bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr){
  double energy = 0.0;
  Timer total; total.start();

  omp_set_num_threads(cube_outer_threads);

  int matrices = cube_outer_threads; 
  if(sphere_outer_threads > matrices) 
    matrices = sphere_outer_threads;

  HostMatrix<double> fort_forces_ms[matrices];
  if (compute_forces) {
      for(int i = 0; i < matrices; i++) {
          fort_forces_ms[i].resize(fortran_vars.max_atoms, 3);
          fort_forces_ms[i].zero();
      }
  }

  HostMatrix<base_scalar_type> rmm_outputs[matrices];
  if (compute_rmm) {
      for(int i = 0; i < matrices; i++) {
          rmm_outputs[i].resize(fortran_vars.rmm_output.width, fortran_vars.rmm_output.height);
          rmm_outputs[i].zero();
      }
  }

  #pragma omp parallel for reduction(+:energy) 
  for(int i = 0; i< cube_work.size(); i++) {
      ThreadBufferPool<base_scalar_type> pool(10, cube_pool_sizes[i]);
      double local_energy = 0; Timers ts; Timer t;
      int id = omp_get_thread_num();

      omp_set_num_threads(cube_inner_threads);
      t.start();
      long long cost = 0;
      for(int j = 0; j < cube_work[i].size(); j++) {
         pool.reset();
         int ind = cube_work[i][j];

         cubes[ind].solve(ts, compute_rmm,lda,compute_forces, compute_energy, 
             local_energy, fort_forces_ms[i], pool, cube_inner_threads, rmm_outputs[i]);
         cost += cubes[ind].cost();
      }

      t.stop();
      printf("Cube workload %d took %ds %dms and it has %d elements (%lld nanounits) (%d)\n", i, 
        t.getSec(), t.getMicrosec(), cube_work[i].size(), cost, id);
      cout << ts;

      energy += local_energy;
  }

  omp_set_num_threads(sphere_outer_threads);

  #pragma omp parallel for reduction(+:energy) 
  for(int i = 0; i< sphere_work.size(); i++) {
      ThreadBufferPool<base_scalar_type> pool(10, sphere_pool_sizes[i]);
      double local_energy = 0; Timers ts; Timer t;
      int id = omp_get_thread_num();

      omp_set_num_threads(sphere_inner_threads);

      t.start();
      long long cost = 0;
      for(int j = 0; j < sphere_work[i].size(); j++) {
          pool.reset();

          int ind = sphere_work[i][j];
          spheres[ind].solve(ts, compute_rmm,lda,compute_forces, compute_energy, 
              local_energy, fort_forces_ms[i], pool, sphere_inner_threads, rmm_outputs[i]);
          cost += spheres[ind].cost();
      }

      t.stop();
      printf("Sphere workload %d took %ds %dms and it has %d elements (%lld nanounits) (%d)\n", i, 
        t.getSec(), t.getMicrosec(), sphere_work[i].size(), cost, id);
      cout << ts;

      energy += local_energy;
  }

  if (compute_forces) {
      FortranMatrix<double> fort_forces_out(fort_forces_ptr, fortran_vars.atoms, 3, fortran_vars.max_atoms);
      for(int k = 0; k < matrices; k++) {
          for(int i = 0; i < fortran_vars.atoms; i++) {
              for(int j = 0; j < 3; j++) {
                  fort_forces_out(i,j) += fort_forces_ms[k](i,j);
              }
          }
      }
  }

  if (compute_rmm) {
      for(int k = 0; k < matrices; k++) {
          for(int i = 0; i < rmm_outputs[k].width; i++) {
              for(int j = 0; j < rmm_outputs[k].height; j++) {
                  fortran_vars.rmm_output(i,j) += rmm_outputs[k](i,j);
              }
          }
      }
  }
    
  total.stop();
  cout << "iteracion total: " << total << endl;
  *fort_energy_ptr = energy;
  if(*fort_energy_ptr != *fort_energy_ptr) {
      std::cout << "I see dead peaple " << std::endl;
#ifndef CPU_KERNELS
    cudaDeviceReset();
#endif
     exit(1);
   }
}

/**********************
 * Sphere
 **********************/
Sphere::Sphere(void) : atom(0), radius(0) { }
Sphere::Sphere(uint _atom, double _radius) : atom(_atom), radius(_radius) { }

/**********************
 * Cube
 **********************/

template class PointGroup<double>;
template class PointGroup<float>;
}
