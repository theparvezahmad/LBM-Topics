program try
  implicit none

  double precision::A(5)
  integer::m(1)
  A = [10.5, 0.5, 2.3, 50.0, 1.2]
  m = maxloc(A)
  write (*, *) m
end program try
