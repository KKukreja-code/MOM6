! MOM6 module to calculate the variance production due to remapping.
module MOM_tracer_remap_variance_production

implicit none ; private

public remapping_variance_production

contains

subroutine remapping_variance_production(n0, h0, u0, n1, h1, u1, itgt_start, itgt_end, &
                                     isrc_start, isrc_end, h_sub, u_sub, col_var_production)

  integer,                intent(in) :: n0 !< Number of cells on source grid
  real, dimension(n0),    intent(in) :: h0 !< Cell widths on source grid [H]
  real, dimension(n0),    intent(in) :: u0 !< Cell averages on source grid [A]
  integer,                intent(in) :: n1 !< Number of cells on target grid
  real, dimension(n1),    intent(in) :: h1 !< Cell widths on target grid [H]
  real, dimension(n1),    intent(in) :: u1 !< Cell averages on target grid [A]
  integer,                intent(in) :: itgt_start(n1) !< Index of first sub-cell within each target cell
  integer,                intent(in) :: itgt_end(n1) !< Index of last sub-cell within each target cell
  integer, dimension(n0), intent(in) :: isrc_start ! Index of first sub-cell within each source cell
  integer, dimension(n0), intent(in) :: isrc_end ! Index of last sub-cell within each source cell
  real,                   intent(in) :: h_sub(n0+n1+1) !< Overlapping sub-cell thicknesses, h_sub [H]
  real,                   intent(in) :: u_sub(n0+n1+1) !< Sub-cell cell averages (size n1) [A]
  real, dimension(n1), intent(inout) :: col_var_production !< Remap variance production for a column

  ! Local variables
  integer :: i1, j2 !< Loop counters

  if (n0 >= n1) then
    do i1 = 1, n1
      col_var_production(i1) = -1 * h0(i1) * (u0(i1) ** 2)
    enddo
  else ! if n1 > n0
    do i1 = 1, n0
      col_var_production(i1) = -1 * h0(i1) * (u0(i1) ** 2)
    enddo
    do i1 = n0+1, n1
      col_var_production(i1) = 0.0
    enddo
  endif

  do i1 = 1, n1
    col_var_production(i1) = col_var_production(i1) + (h1(i1) * (u1(i1) ** 2))
      do j2 = isrc_start(i1), isrc_end(i1)
        if (j2 < itgt_start(i1) .or. j2 > itgt_end(i1)) then
          col_var_production(i1) = col_var_production(i1) + (h_sub(j2) * (u_sub(j2)**2))
        endif
      enddo
      do j2 = itgt_start(i1), itgt_end(i1)
        if (j2 < isrc_start(i1) .or. j2 > isrc_end(i1)) then
          col_var_production(i1) = col_var_production(i1) - (h_sub(j2) * (u_sub(j2)**2))
        endif
      enddo
  enddo

end subroutine remapping_variance_production

end module MOM_tracer_remap_variance_production
