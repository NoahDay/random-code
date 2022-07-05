program open_nc
  use netcdf
implicit none
  integer :: loop_w, ncid, status
  ! This is a comment line; it is ignored by the compiler
  print *, 'Opening netcdf'
  status = nf90_open(trim('/Users/noahday/GitHub/cice-dirs/input/CICE_data/grid/om2_1deg/icegrid_nonc.nc'), NF90_NOWRITE, ncid)
  if (status /= nf90_noerr) call handle_err(status)

end program open_nc
