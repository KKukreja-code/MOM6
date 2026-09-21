! This module contains subroutines to calculate how horizontal and vertical diffusion change tracer variance locally.
module MOM_tracer_parameterised_variance_production

use MOM_grid,          only : ocean_grid_type
use MOM_tracer_types,  only : tracer_type
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_domains,       only : pass_var

implicit none ; private

#include <MOM_memory.h>

public compute_cell_hordiff_variance_production, hor_interface_variance_production, T_cell_diabatic_variance_production

contains

!< Subroutine to calculate variance production at grid cell interfaces
subroutine hor_interface_variance_production(G, GV, Tr, Idt, Ihdxdy, h, Coef_x, Coef_y, k)

  type(ocean_grid_type),                        intent(in) :: G       !< Ocean grid structure
  type(verticalGrid_type),                      intent(in) :: GV      !< Ocean vertical grid structure
  type(tracer_type),                            intent(in) :: Tr      !< Pointer to the tracer regsitry
  real,                                         intent(in) :: Idt     !< Inverse time interval [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G)),             intent(in) :: Ihdxdy  !<  The inverse of the volume or mass of fluid
                                                                      !! in a layer in a grid cell
                                                                      !! [H-1 L-2 ~> m-3 or kg-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in) :: h       !< thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in) :: Coef_x  !< The coefficients relating zonal tracer
                                                                      !! differences to time-integrated fluxes, in
                                                                      !! [L2 ~> m2] for some schemes and
                                                                      !! [H L2 ~> m3 or kg] for others.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in) :: Coef_y  !< The coefficients relating meridional tracer
                                                                      !! differences to time-integrated fluxes, in
                                                                      !! [L2 ~> m2] for some schemes and
                                                                      !! [H L2 ~> m3 or kg] for others.
  integer,                                      intent(in) :: k       !< vertical index

  ! Local variables
  integer :: is, ie, js, je                        !< Grid cell centre and layer indexes
  integer :: i, j                                  !< Counters
  real :: dtr_left, dtr_right, dtr_top, dtr_bottom !< Changes in tracer content at neighbouring interfaces
                                                   !! due to diffusion, for variance production calculation [Conc]
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do i = is,ie ; do j = js,je
    ! Calculating the change in tracer caused by diffusive fluxes at the left, right, top, and bottom
    ! interfaces of a cell
    dtr_left = Ihdxdy(i,j) * (Coef_x(I-1,j,1) * (Tr%t(i-1,j,k) - Tr%t(i,j,k)))
    dtr_right = Ihdxdy(i,j) * (Coef_x(I,j,1) * (Tr%t(i,j,k) - Tr%t(i+1,j,k)))
    dtr_top = Ihdxdy(i,j) * (Coef_y(i,J,1) * (Tr%t(i,j,k) - Tr%t(i,j+1,k)))
    dtr_bottom = Ihdxdy(i,j) * (Coef_y(i,J-1,1) * (Tr%t(i,j-1,k) - Tr%t(i,j,k)))

    ! Assigning terms in expansion of thickness-weighted squared tracer concentration for a cell to each
    ! interface; see manuscript for derivation.
    Tr%leftint_hordiff_var_prod(i,j,k) = Tr%leftint_hordiff_var_prod(i,j,k) + &
    h(i,j,k) * Idt * G%areaT(i,j) * (2*Tr%t(i,j,k)*dtr_left - dtr_left*dtr_right - dtr_left*dtr_top + &
    dtr_left*dtr_bottom + dtr_left*dtr_left)
    Tr%rightint_hordiff_var_prod(i,j,k) = Tr%rightint_hordiff_var_prod(i,j,k) + &
    h(i,j,k) * Idt * G%areaT(i,j) * (-2*Tr%t(i,j,k)*dtr_right - dtr_right*dtr_left - dtr_right*dtr_bottom + &
    dtr_right*dtr_top + dtr_right*dtr_right)
    Tr%bottomint_hordiff_var_prod(i,j,k) = Tr%bottomint_hordiff_var_prod(i,j,k) + &
    h(i,j,k) * Idt * G%areaT(i,j) * (2*Tr%t(i,j,k)*dtr_bottom + dtr_bottom*dtr_left - dtr_bottom*dtr_right - &
    dtr_bottom*dtr_top + dtr_bottom*dtr_bottom)
    Tr%topint_hordiff_var_prod(i,j,k) = Tr%topint_hordiff_var_prod(i,j,k) + &
    h(i,j,k) * Idt * G%areaT(i,j) * (-2*Tr%t(i,j,k)*dtr_top - dtr_top*dtr_left + dtr_top*dtr_right - &
    dtr_top*dtr_bottom + dtr_top*dtr_top)
  enddo ; enddo

end subroutine hor_interface_variance_production

!< Subroutine to sum the contirbutions that contribute to changes in variance due to horizontal diffusion
subroutine compute_cell_hordiff_variance_production(G, GV, Tr)

  type(ocean_grid_type),                        intent(in) :: G         !< Ocean grid structure
  type(verticalGrid_type),                      intent(in) :: GV        !< Ocean vertical grid structure
  type(tracer_type),                            intent(in) :: Tr        !< Pointer to the tracer regsitry

  integer :: is, ie, js, je, nz  !< Grid cell centre and layer indexes
  integer :: i, j, k             !< Counters
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call pass_var(Tr%leftint_hordiff_var_prod, G%Domain, halo=1)
  call pass_var(Tr%rightint_hordiff_var_prod, G%Domain, halo=1)
  call pass_var(Tr%topint_hordiff_var_prod, G%Domain, halo=1)
  call pass_var(Tr%bottomint_hordiff_var_prod, G%Domain, halo=1)

  ! Computing a grid-cell variance production diagnostic by splitting variance change at each
  ! interface equally between its two neighbouring cells.
  do i=is,ie ; do j=js,je ; do k=1,nz
    Tr%cell_hordiff_var_prod(i,j,k) = &
    0.5 * (Tr%leftint_hordiff_var_prod(i,j,k)   + Tr%rightint_hordiff_var_prod(i,j,k) + &
           Tr%topint_hordiff_var_prod(i,j,k)    + Tr%bottomint_hordiff_var_prod(i,j,k) + &
           Tr%leftint_hordiff_var_prod(i+1,j,k) + Tr%rightint_hordiff_var_prod(i-1,j,k) + &
           Tr%topint_hordiff_var_prod(i,j-1,k)  + Tr%bottomint_hordiff_var_prod(i,j+1,k))
  enddo ; enddo ; enddo

end subroutine compute_cell_hordiff_variance_production

!< Subroutine to compute the variance production at cell centres due to diabatic diffusion.
subroutine T_cell_diabatic_variance_production(G, GV, Tdif_flx, temp_new, temp_diag, T_cell_var_prod)

  type(ocean_grid_type),                        intent(in) :: G         !< Ocean grid structure
  type(verticalGrid_type),                      intent(in) :: GV        !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(in) :: Tdif_flx  !< diffusive diapycnal heat flux across
                                                                        !! interfaces
                                                                        !! [C H T-1 ~> degC m s-1 or degC kg m-2 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in) :: temp_new  !! Updated temperatures [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in) :: temp_diag !< Diagnostic array of previous
                                                                        !! temperatures [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: T_cell_var_prod !< Averaged variance production in cell
                                                                              !! due to diabatic diffusion
                                                                              !! [CU2 H T-1 ~> conc2 m s-1]

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: T_int_var_prod !< Variance production due to diabatic diffusion
                                                                !! at an interface [CU2 H T-1 ~> conc2 m s-1]
  integer :: is, ie, js, je, nz  !< Grid cell centre and layer indexes
  integer :: i, j, k             !< Counters

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  T_int_var_prod(:,:,:) = 0.0

  ! See manuscript for deriving interface variance production by expanding the change in thickness-weighted squared
  ! tracer concentration in a cell due to vertical mixing, and assigning resulting terms to the relevant interface.
  ! This calculation involves substituting the update equation for grid cell tracer concentration due to vertical
  ! mixing to simplify the expression for interface variance production.

  do i=is,ie ; do j=js,je ; do K=2,nz
    T_int_var_prod(i,j,K) = Tdif_flx(i,j,K)*((temp_diag(i,j,K)- temp_diag(i,j,K-1)) + &
                                             (temp_new(i,j,K) - temp_new(i,j,K-1)))
  enddo ; enddo ; enddo

  ! Equally split variance production at an interface across its two neighbouring cells to get a grid cell measure.
  do i=is,ie ; do j=js,je ; do K=1,nz
    T_cell_var_prod(i,j,k) = 0.5*(T_int_var_prod(i,j,k)+T_int_var_prod(i,j,k+1))
  enddo ; enddo ; enddo

end subroutine T_cell_diabatic_variance_production

end module MOM_tracer_parameterised_variance_production
