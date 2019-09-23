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

   ! Sets the arrays needed for QXD.
   subroutine qxd_alloc_arrays(n_qm, n_mm)
      use QXD_data, only: qxd_s, qxd_alpha0, qxd_alphaq, qxd_zeta0, qxd_zetaq,&
                          qxd_qm_allocated
      implicit none
      integer, intent(in) :: n_qm, n_mm
      integer :: n_tot

      n_tot = n_qm + n_mm
      if (allocated(qxd_s)     ) deallocate(qxd_s)
      if (allocated(qxd_zeta0) ) deallocate(qxd_zeta0)
      if (allocated(qxd_zetaq) ) deallocate(qxd_zetaq)
      if (allocated(qxd_alpha0)) deallocate(qxd_alpha0)
      if (allocated(qxd_alphaq)) deallocate(qxd_alphaq)

      allocate(qxd_s(n_tot)     ); qxd_s      = 0.0D0
      allocate(qxd_zeta0(n_tot) ); qxd_zeta0  = 0.0D0
      allocate(qxd_zetaq(n_tot) ); qxd_zetaq  = 0.0D0
      allocate(qxd_alpha0(n_tot)); qxd_alpha0 = 0.0D0
      allocate(qxd_alphaq(n_tot)); qxd_alphaq = 0.0D0

      if (.not. qxd_qm_allocated) call qxd_alloc_qm_arrays(n_qm)
   end subroutine qxd_alloc_arrays

   subroutine qxd_alloc_qm_arrays(n_qm)
      use QXD_data, only: qxd_Q, qxd_qm_allocated
      implicit none
      integer, intent(in) :: n_qm
      
      if (allocated(qxd_Q)) deallocate(qxd_Q)
      allocate (qxd_Q(n_qm))
      qxd_Q = 0.0D0
      qxd_qm_allocated = .true.
   end subroutine

   ! Gets QXD parameters for the QM region.
   subroutine qxd_set_qm_params(n_qm, atom_z)
      use QXD_data, only: qxd_s, qxd_zeta0, qxd_zetaq, qxd_alpha0, qxd_alphaq
      use QXD_refs, only: ref_s, ref_zeta0, ref_zetaq, ref_alpha0, ref_alphaq
      implicit none
      integer, intent(in) :: n_qm, atom_z(:)
      integer :: iatom

      do iatom = 1, n_qm
         qxd_s(iatom)      = ref_s(atom_z(iatom))
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

   ! Adds QXD to atomic-basis Fock matrix, and outputs energy term E_D.
   subroutine qxd_fock_x(energ_x, fock_op, m_size, pos)
      use QXD_data, only: qxd_Q, qxd_s, qxd_zeta0, qxd_zetaq
      implicit none
      real(kind=8)  , intent(in)    :: pos(:,:)
      real(kind=8)  , intent(out)   :: energ_x
      type(operator), intent(inout) :: fock
      
      integer      :: iatom, jatom, n_qm, n_tot
      real(kind=8) :: dist_ij, zeta_i, zeta_j, zeta_ij, delta_ij, delta_ji
      real(kind=8), allocatable :: f_mat

      energ_x = 0.0D0
      allocate(f_mat(m_size, m_size))
      call fock%Gets_data_AO(f_mat)

      n_qm  = size(qxd_Q,1)
      n_tot = size(qxd_s,1)
      do iatom = 1, n_qm
         do jatom = n_qm, n_tot
            ! First part - energy contributions.
            zeta_i  = qxd_zeta0(iatom) * exp(- qxd_zetaq(iatom) * qxd_Q(iatom))
            zeta_j  = qxd_zeta0(jatom)
            zeta_ij = zeta_ab(zeta_i, zeta_j)

            dist_ij  = dist(pos(iatom,:), pos(jatom,:))
            delta_ij = delta_ab(zeta_i, zeta_j, dist_ij)
            delta_ji = delta_ab(zeta_j, zeta_i, dist_ij)

            energ_x = energ_x + qxd_s(iatom) * qxd_s(iatom) * zeta_ij * &
                                (delta_ij - delta_ji)

            ! Second part - Fock contributions.
         enddo
      enddo


      deallocate(f_mat)
   end subroutine qxd_fock_x

   ! Helper functions
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
end module QXD_subs