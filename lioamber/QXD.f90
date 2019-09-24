!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%% QXD.f90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!   This file contains an implementation of the QXD model for exchange and    !
! dispersion interactions as presented by Kuechler et al. (JCP 143, p. 234111,!
! DOI:10.1063/1.4937155).                                                     !
!   Additional subroutines are included in order to substract classical       !
! Lennard-Jones contributions from energies and forces if the MM software     !
! does not do it by itself.                                                   !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module QXD_data

   real(kind=8), allocatable :: qxd_Q(:)

   real(kind=8), allocatable :: qxd_s(:)
   real(kind=8), allocatable :: qxd_neff0(:)
   real(kind=8), allocatable :: qxd_zeta0(:)
   real(kind=8), allocatable :: qxd_zetaq(:)
   real(kind=8), allocatable :: qxd_alpha0(:)
   real(kind=8), allocatable :: qxd_alphaq(:)

   logical :: qxd_qm_allocated = .false.
end module QXD_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module QXD_subs
   implicit none
contains
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   !%% GENERAL SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   ! Sets the arrays needed for QXD.
   subroutine qxd_alloc_arrays(n_qm, n_mm)
      use QXD_data, only: qxd_s, qxd_alpha0, qxd_alphaq, qxd_zeta0, qxd_zetaq,&
                          qxd_neff0, qxd_qm_allocated
      implicit none
      integer, intent(in) :: n_qm, n_mm
      integer :: n_tot

      n_tot = n_qm + n_mm
      if (allocated(qxd_s)     ) deallocate(qxd_s)
      if (allocated(qxd_neff0) ) deallocate(qxd_neff0)
      if (allocated(qxd_zeta0) ) deallocate(qxd_zeta0)
      if (allocated(qxd_zetaq) ) deallocate(qxd_zetaq)
      if (allocated(qxd_alpha0)) deallocate(qxd_alpha0)
      if (allocated(qxd_alphaq)) deallocate(qxd_alphaq)

      allocate(qxd_s(n_tot)     ); qxd_s      = 0.0D0
      allocate(qxd_neff0(n_qm)  ); qxd_neff0  = 0.0D0
      allocate(qxd_zeta0(n_tot) ); qxd_zeta0  = 0.0D0
      allocate(qxd_zetaq(n_tot) ); qxd_zetaq  = 0.0D0
      allocate(qxd_alpha0(n_tot)); qxd_alpha0 = 0.0D0
      allocate(qxd_alphaq(n_tot)); qxd_alphaq = 0.0D0

      if (.not. qxd_qm_allocated) call qxd_alloc_qm_arrays(n_qm)
   end subroutine qxd_alloc_arrays

   subroutine qxd_alloc_qm_arrays(n_qm)
      use QXD_data, only: qxd_Q, qxd_neff0, qxd_qm_allocated
      implicit none
      integer, intent(in) :: n_qm
      
      if (allocated(qxd_Q)) deallocate(qxd_Q)
      allocate (qxd_Q(n_qm))
      qxd_Q = 0.0D0
      qxd_qm_allocated = .true.
   end subroutine

   ! Gets QXD parameters for the QM region.
   subroutine qxd_set_qm_params(n_qm, atom_z)
      use QXD_data, only: qxd_s, qxd_zeta0, qxd_zetaq, qxd_alpha0, qxd_alphaq,&
                          qxd_neff0
      use QXD_refs, only: ref_s, ref_zeta0, ref_zetaq, ref_alpha0, ref_alphaq,&
                          ref_neff0
      implicit none
      integer, intent(in) :: n_qm, atom_z(:)
      integer :: iatom

      do iatom = 1, n_qm
         qxd_s(iatom)      = ref_s(atom_z(iatom))
         qxd_neff0(iatom)  = ref_neff0(atom_z(iatom))
         qxd_zeta0(iatom)  = ref_zeta0(atom_z(iatom))
         qxd_zetaq(iatom)  = ref_zetaq(atom_z(iatom))
         qxd_alpha0(iatom) = ref_alpha0(atom_z(iatom))
         qxd_alphaq(iatom) = ref_alphaq(atom_z(iatom))
      enddo
   end subroutine qxd_set_qm_params
   
   ! Gets QXD parameters from classical LJ parameters for MM atoms (eq. 31-33).
   ! LJ parameters MUST be in atomic units; with epsilon(i) being
   ! sqrt(epsilon(i,i)), and sigma(i) being Rmin(i,i) / 2.
   subroutine qxd_set_mm_params(n_qm, n_mm, LJ_epsilon, LJ_sigma, atom_z)
      use QXD_data, only: qxd_s, qxd_zeta0, qxd_zetaq, qxd_alpha0, qxd_alphaq
      implicit none
      integer     , intent(in) :: n_qm, n_mm, atom_z(:)
      real(kind=8), intent(in) :: LJ_epsilon(:), LJ_sigma(:)

      integer :: iatom

      do iatom = n_qm, n_qm + n_mm
         qxd_s(iatom)      = 9.4423D0 * (LJ_epsilon(iatom) ** 0.4111D0) * &
                                        (LJ_sigma(iatom)   ** 2.8208D0)
         qxd_zeta0(iatom)  = 3.7893D0 * (LJ_epsilon(iatom) ** -0.0192D0) * &
                                        (LJ_sigma(iatom)   ** -0.7249D0)
         qxd_alpha0(iatom) = cbrt(LJ_epsilon(iatom) * LJ_sigma(iatom) ** 6)
         qxd_alpha0(iatom) = qxd_G(atom_z(iatom)) * (qxd_alpha0(iatom) ** 2)

         qxd_zetaq(iatom)  = 0.0D0
         qxd_alphaq(iatom) = 0.0D0
         qxd_neff0(iatom)  = ref_neff0(atom_z(iatom))
      enddo
   end subroutine qxd_set_mm_params

   subroutine qxd_set_params(n_qm, n_mm, LJ_epsilon, LJ_sigma, atom_z)
      implicit none
      integer     , intent(in) :: n_qm, n_mm, atom_Z(:)
      real(kind=8), intent(in) :: LJ_epsilon(:), LJ_sigma(:)

      call qxd_alloc_arrays(n_qm, n_mm)
      call qxd_set_qm_params(n_qm, atom_z)
      call qxd_set_mm_params(n_qm, n_mm, LJ_epsilon, LJ_sigma, atom_z)
   end subroutine qxd_set_params

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   !%% FOCK - ENERGY TERMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   ! Adds exchange corrections to atomic-basis Fock matrix, and outputs energy
   ! term E_x.
   subroutine qxd_fock_x(energ_x, fock_op, m_size, pos, Smat, atom_of_func)
      use QXD_data, only: qxd_Q, qxd_s, qxd_zeta0, qxd_zetaq
      implicit none
      integer       , intent(in)    :: m_size, atom_of_func(:)
      real(kind=8)  , intent(in)    :: pos(:,:), Smat(:,:)
      real(kind=8)  , intent(out)   :: energ_x
      type(operator), intent(inout) :: fock
      
      integer      :: iatom, jatom, n_qm, n_tot, ifunct, jfunct
      real(kind=8) :: dist_ij, zeta_i, zeta_j, zeta_ij, delta_ij, delta_ji, &
                      dzeta_ij, ddelta_ij, ddelta_ji
      real(kind=8), allocatable :: f_mat, dEx_dQ

      energ_x = 0.0D0
      n_qm  = size(qxd_Q,1)
      n_tot = size(qxd_s,1)
      allocate(f_mat(m_size, m_size), dEx_dQ(n_qm))
      call fock%Gets_data_AO(f_mat)

      do iatom = 1, n_qm
         do jatom = n_qm, n_tot
            ! First part - energy contributions.
            ! Eq. 15.
            zeta_i  = qxd_zeta0(iatom) * exp(- qxd_zetaq(iatom) * qxd_Q(iatom))
            zeta_j  = qxd_zeta0(jatom)

            ! Eq. 21.
            zeta_ij = zeta_ab(zeta_i, zeta_j)

            ! Eq. 20.
            dist_ij  = dist(pos(iatom,:), pos(jatom,:))
            delta_ij = delta_ab(zeta_i, zeta_j, dist_ij)
            delta_ji = delta_ab(zeta_j, zeta_i, dist_ij)

            ! Eq. 18.
            energ_x = energ_x + qxd_s(iatom) * qxd_s(jatom) * zeta_ij * &
                                (delta_ij - delta_ji)

            ! Second part - Fock contributions. These are obtained by
            ! separating dE/dRho as dE/dZi * dZi/dQi * dQi/dRho.
            ! dZi/dQi is simply (- Zqi * Zi). The "-" is ommited since
            ! afterwards dQi/dRho < 0, so it gets cancelled.
            dzeta_ij  = dzeta_di(zeta_i, zeta_j)
            ddelta_ij = ddelta_di(zeta_i, zeta_j, dist_ij, delta_ij)
            ddelta_ji = ddelta_dj(zeta_j, zeta_i, dist_ij, delta_ji)

            dEx_dQ(iatom) = zeta_ij  * (ddelta_ij - ddelta_ji) + &
                            dzeta_ij * (delta_ij  - delta_ji)
            dEx_dQ(iatom) = dEx_dQ(iatom) * qxd_zetaq(iatom) * zeta_i &
                            * qxd_s(iatom)  * qxd_s(jatom)
         enddo
      enddo

      ! Adds fock contributions to fock matrix (eqs. 27-28).
      do ifunct = 1, M
         do jfunct = ifunct, M
            if (atom_of_func(ifunct) == atom_of_func(jfunct)) then
               f_mat(ifunct, jfunct) = f_mat(ifunct,jfunct) + &
                                       dEx_dQ(atom_of_func(ifunct)) * &
                                       Smat(ifunct, jfunct)
            else
               f_mat(ifunct, jfunct) = f_mat(ifunct,jfunct) + &
                                       dEx_dQ(atom_of_func(ifunct)) * &
                                       Smat(ifunct, jfunct) * 0.5D0
            endif
            f_mat(jfunct, ifunct) = f_mat(ifunct, jfunct)
         enddo
      enddo
      call fock%Sets_data_AO(f_mat)

      deallocate(f_mat, dEx_dQ)
   end subroutine qxd_fock_x

   ! Adds dispersion corrections to atomic-basis Fock matrix, and outputs energy
   ! term E_d.
   subroutine qxd_fock_d(energ_x, fock_op, m_size, pos, Smat, atom_of_func)
      use QXD_data, only: qxd_Q, qxd_s, qxd_zeta0, qxd_zetaq, qxd_neff0
      implicit none
      integer       , intent(in)    :: m_size, atom_of_func(:)
      real(kind=8)  , intent(in)    :: pos(:,:), Smat(:,:)
      real(kind=8)  , intent(out)   :: energ_d
      type(operator), intent(inout) :: fock
      
      integer      :: iatom, jatom, n_qm, n_tot, ifunct, jfunct
      real(kind=8) :: dist_ij, zeta_i, zeta_j, delta_ij, delta_ji, &
                      eta_i, eta_j, c6_ij, &
                      ddelta_ij, ddelta_ji, dc6_ij
      real(kind=8), allocatable :: f_mat, dEd_dQ

      energ_d = 0.0D0
      n_qm  = size(qxd_Q,1)
      n_tot = size(qxd_s,1)
      allocate(f_mat(m_size, m_size), dEd_dQ(n_qm))
      call fock%Gets_data_AO(f_mat)

      do iatom = 1, n_qm
         do jatom = n_qm, n_tot
            ! First part - energy contributions.
            ! Eq. 17.
            alpha_i  = qxd_alpha0(iatom) * exp(- qxd_alphaq(iatom) * qxd_Q(iatom))
            alpha_j  = qxd_alpha0(jatom)

            ! Eqs. 22, 23 and 16.
            eta_i = sqrt((neff_i - q_i) / a_i)
            eta_j = sqrt(neff_j / a_j)
            c6_ij = c6_ab(alpha_i, alpha_j, eta_i, eta_j)
            
            ! Eqs. 15 and 20.          
            zeta_i   = qxd_zeta0(iatom) * exp(- qxd_zetaq(iatom) * qxd_Q(iatom))
            zeta_j   = qxd_zeta0(jatom)
            dist_ij  = dist(pos(iatom,:), pos(jatom,:))

            delta_ij = delta_ab(zeta_i, zeta_j, dist_ij)
            delta_ji = delta_ab(zeta_j, zeta_i, dist_ij)
            
            ! Eqs. 25 and 24.
            b_ij  = b_ab(zeta_i, zeta_j, delta_ij, delta_ji, dist_ij)
            s6_ij = s6_ab(b_ab, dist_ij)

            ! Eq. 19.
            energ_d = energ_d - s6_ij * c6_ij / (dist_ij ** 6)

            ! Second part - Fock contributions. These are obtained by
            ! separating dE/dRho as dE/dQi * dQi/dRho.
            ! dZi/dQi is simply (- Zqi * Zi) and the same goes for alpha.
            ddelta_ij = ddelta_di(zeta_i, zeta_j, dist_ij, delta_ij)
            ddelta_ji = ddelta_dj(zeta_j, zeta_i, dist_ij, delta_ji)







            dEd_dQ(iatom) = zeta_ij  * (ddelta_ij - ddelta_ji) + &
                            dzeta_ij * (delta_ij  - delta_ji)
            dEd_dQ(iatom) = dEd_dQ(iatom) * qxd_zetaq(iatom) * zeta_i &
                            * qxd_s(iatom)  * qxd_s(jatom)
         enddo
      enddo

      ! Adds fock contributions to fock matrix.
      do ifunct = 1, M
         do jfunct = ifunct, M
            if (atom_of_func(ifunct) == atom_of_func(jfunct)) then
               f_mat(ifunct, jfunct) = f_mat(ifunct,jfunct) + &
                                       dEd_dQ(atom_of_func(ifunct)) * &
                                       Smat(ifunct, jfunct)
            else
               f_mat(ifunct, jfunct) = f_mat(ifunct,jfunct) + &
                                       dEd_dQ(atom_of_func(ifunct)) * &
                                       Smat(ifunct, jfunct) * 0.5D0
            endif
            f_mat(jfunct, ifunct) = f_mat(ifunct, jfunct)
         enddo
      enddo
      call fock%Sets_data_AO(f_mat)

      deallocate(f_mat, dEd_dQ)
   end subroutine qxd_fock_d

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   !%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   function zeta_ab(z_i, z_j) result(zeta_out)
      use constants_mod, only: PI
      implicit none
      real(kind=8), intent(in) :: z_i, z_j
      real(kind=8)             :: zeta_out
      
      zeta_out = PI * 8.0D0 * ((z_i * z_i - z_j * z_j) ** 3)
      zeta_out = (z_i ** 3) * (z_j ** 3) / zeta_out
      return
   end function zeta_ab

   function dist(r_i, r_j) result(d_out)
      implicit none
      real(kind=8), intent(in) :: r_i(3), r_j(3)
      real(kind=8)             :: d_out

      d_out = (r_i(1) - r_j(1)) * (r_i(1) - r_j(1)) &
            + (r_i(2) - r_j(2)) * (r_i(2) - r_j(2)) &
            + (r_i(3) - r_j(3)) * (r_i(3) - r_j(3))
      d_out = sqrt(d_out)
      return
   end function dist

   function delta_ab(z_i, z_j, d_ij) result(delta_out)
      implicit none
      real(kind=8), intent(in) :: z_i, z_j, d_ij
      real(kind=8)             :: delta_out

      delta_out = 4.0D0 * z_i + d_ij * (z_i * z_i - z_j * z_j)
      delta_out = delta_out * z_j * exp(- z_i * d_ij) / d_ij
      return
   end function delta_ab

   function c6_ab(a_i, a_j, eta_i, eta_j) result(c6_out)
      implicit none
      real(kind=8), intent(in) :: a_i, a_j, eta_i, eta_j
      real(kind=8)             :: c6_out

      c6_out = 1.5D0 * a_i * a_j * eta_i * eta_j / (eta_i + eta_j)
      return
   end function c6_ab

   function b_ab(z_i, z_j, d_ij, d_ji, r_ij) result(b_out)
      ! This function is not properly written in the original paper.
      ! Write it in Wolfram and take the result (checked manually, it's good).
      implicit none
      real(kind=8), intent(in) :: a_i, a_j, d_ij, d_ji, r_ij
      real(kind=8)             :: b_out

      b_out = (z_i * exp(-z_j * r_ij) - z_j * exp(-z_i * r_ij))
      b_out = b_out * (z_i * z_i - z_j * z_j) / (r_ij * (d_ij - d_ji))
      b_out = b_out + 1.0D0 / r_ij + (z_i * d_ij - z_j * d_ji) / (d_ij - d_ji)
      return
   end function b_ab

   function s6_ab(b_ij, r_ij) result(s6_out)
      ! The factorial summation in this function is written explicitely in
      ! order to make it more efficient. It goes from k=0 to k=6.
      implicit none
      real(kind=8), intent(in) :: b_ij, r_ij
      real(kind=8)             :: s6_out, br_ij
      
      br_ij  = b_ij * r_ij
      s6_out = 0.00833333D0 * (br_ij**5) + 0.001388888D0 * (br_ij**6)
      s6_out = 0.04166666D0 * (br_ij**4) + 0.166666666D0 * (br_ij**3) + s6_out
      s6_out = 0.5D0        * (br_ij**2) + br_ij + 1.0D0              + s6_out
      s6_out = 1.0D0 - exp(-br_ij) * s6_out
      return
   end function s6_ab

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   !%% DERIVATIVES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   ! Derives Zij only with respect to Zi.
   function dzeta_di(z_i, z_j) result(dzeta_out)
      implicit none
      real(kind=8), intent(in) :: z_i, z_j
      real(kind=8)             :: dzeta_out

      dzeta_out = (z_i * z_i + z_j * z_j) / &
                  (z_i * z_i - z_j * z_j) ** 4                  
      dzeta_out = - 3.0D0 * (z_j ** 3) * (z_i * z_i) * dzeta_out
      return
   end function

   ! Derives dij with respect to Zi.
   function ddelta_di(z_i, z_j, d_ij, delta_ij) result(ddelta_out)
      implicit none
      real(kind=8), intent(in) :: z_i, z_j, d_ij, delta_ij
      real(kind=8)             :: ddelta_out

      ddelta_out = - (z_j / d_ij) * (4.0D0 + 2.0D0 * d_ij * z_i) &
                    * exp(-d_ij* z_i)
      ddelta_out = d_ij * delta_ij + ddelta_out
      return
   end function

   ! Derives dij with respect to Zj.
   function ddelta_dj(z_i, z_j, d_ij, delta_ij) result(ddelta_out)
      implicit none
      real(kind=8), intent(in) :: z_i, z_j, d_ij, delta_ij
      real(kind=8)             :: ddelta_out

      ddelta_out = delta_ij / z_j - 2.0D0 * z_j * z_j * exp(-d_ij * z_i)
      return
   end function

   ! Derives c6 with respect to Qi.
   function dc6_dqi(a_i, a_qi, e_i, e_j, c6_ij) result(dc6_out)
      implicit none
      real(kind=8), intent(in) :: a_i, a_qi, e_i, e_j, c6_ij
      real(kind=8)             :: dc6_out

      dc6_out = 2.0D0 * e_i * e_i * (e_i + e_j) * a_i
      dc6_out = (a_qi - 1.0D0) * e_j / dc6_out

      dc6_out = c6_ij * (dc6_out - a_qi)
      return
   end function dc6_dqi

end module QXD_subs