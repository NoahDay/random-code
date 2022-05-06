program hello
implicit none
  integer :: loop_w
  ! This is a comment line; it is ignored by the compiler
  print *, 'Hello, World!'

  do loop_w=3,30,2
    print *, loop_w
  end do
end program hello