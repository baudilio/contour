module common
  implicit none

  integer, parameter :: SP = Selected_Real_kind(p=6)
  integer, parameter :: DP = Selected_Real_Kind(p=15)

  real(KIND=SP), dimension(:), allocatable :: x, y
  real(KIND=SP), dimension(:), allocatable :: z ! array of contours
  real(KIND=SP), dimension(:,:), allocatable :: d

end module common
