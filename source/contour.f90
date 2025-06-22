program demo
  use common, only : d, x, y, z
  implicit none
  !

  ! Coordinates of the vertices of the grid
  real, parameter :: pxmin =  -1.5, pymin =  -1.5
  real, parameter :: pxmax = 1.5, pymax = 1.5

  integer, parameter :: nx = 300, ny = 250 ! number of grid points, nx_by_ny
  integer, parameter :: nc = 10  ! number of countours

  real ::  x1,y1,x2,y2
  real :: zmax,zmin
  integer :: i,j

  real :: a, b, c  ! The slopes of the linear interpolation for x and y

  ! ---

  allocate(d(0:nx, 0:ny))
  allocate(x(0:nx))
  allocate(y(0:ny))
  allocate(z(0:nc))


  !
  !     Set coordinates in Y array suitable for
  !     automatic plotting on the graphics screen
  !
  b = (pymax - pymin) / real(ny)
  do j=0, ny
     y(j) = j * b  + pymin
  end do

  !
  !     Set coordinates in X array suitable for
  !     automatic plotting on the graphics screen
  !
  a = (pxmax - pxmin) / real(nx)
  do i=0, nx
     x(i) = i * a  + pxmin
  end do


  !
  !     Create an artificial data surface and calculate the
  !     surface bounds for choosing the contour levels.
  !
  !     The function z = f(x,y), or the data read in from a file.
  !
  do i=0,nx
     do j=0,ny
        d(i,j) = x(i) + y(j)  ! A simple 3D ( in ℝ3 ) plane: x + y - z = 0
     end do
  end do

#ifdef DEBUG
  print '("zMin,zMax ", 2F8.3)', minval(d), maxval(d)
#endif

  zmin= minval(d)
  zmax= maxval(d)

  !
  !     Set a full contingent of contour levels
  !
  c = (zmax - zmin) / nc
  do i=0, nc
     z(i) = i * c  + zmin
  end do

#ifdef DEBUG
  print *, "The countour values"
  print '(F6.2)', z
#endif


  !
  !     Draw a border around the contour plot
  !
  x1 = pxmin
  y1 = pymin
  x2 = pxmax
  y2 = pymax
  call vecout(x1,y1,x1,y2,0.0)
  call vecout(x1,y2,x2,y2,0.0)
  call vecout(x2,y2,x2,y1,0.0)
  call vecout(x2,y1,x1,y1,0.0)

  !
  !     Call the contouring routine
  !
#ifndef DEBUG
  call conrec(d,0,nx,0,ny,x,y,nc,z)
#endif

  stop 0
end program demo
