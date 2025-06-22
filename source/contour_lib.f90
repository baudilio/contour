! subroutine conrec(d,ilb,iub,jlb,jub,x,y,nc,z)
subroutine conrec()
  use common, only : d, &  ! the data (values) in the grid, d := f(x,y)
       x, y, & ! the coordinates of the data points
       z ! the array of the contours
  implicit none

  ! -- Local varaibles
  integer :: ilb,iub, jlb,jub
  integer :: nc
  ! real, dimension(1:nc) :: z

  !
  ! -- Local varaibles
  !
  real, Dimension(0:4) :: h
  integer, dimension(0:4) :: sh
  real, dimension(0:4) :: xh, yh
  real :: x1,x2,y1,y2

  !    Data
  integer, parameter, dimension(1:4) :: im = [0, 1, 1, 0]
  integer, parameter, dimension(1:4) :: jm = [0, 0, 1, 1]
  integer, parameter, dimension(-1:1,-1:1,-1:1) :: castab = reshape([ &
       0,0,9,  0,1,5,  7,4,8, &
       0,3,6,  2,3,2,  6,3,0, &
       8,4,7,  5,1,0,  9,0,0 ], shape(castab))

  !
  !     Use statement functions for the line intersections
  !     BTA: OK - but obsolete prog. style. Funcs. to move into a module.
  integer :: p1,p2
  real :: xsect, ysect
  xsect(p1,p2) = (h(p2)*xh(p1) - h(p1)*xh(p2)) / (h(p2)-h(p1))
  ysect(p1,p2) = (h(p2)*yh(p1)-h(p1)*yh(p2))/(h(p2)-h(p1))

  !
  integer :: case
  integer :: i,j,k,m, m1,m2,m3
  real :: dmin, dmax


  ! ====

  ilb = lbound(d, dim=1)
  iub = ubound(d, dim=1)
  jlb = lbound(d, dim=2)
  jub = ubound(d, dim=2)
  nc = ubound(z, dim=1)

#ifdef DEBUG
  print *, "Arrays indices bounds:"
  print '(I0,":",I0,", ",I0,":",I0, "and z: ", I0)', ilb,iub,jlb,jub,nc
#endif

  !
  !     Scan the arrays, top down, left to right within rows
  !
  RY: do j=jub-1,jlb,-1
     RX: do i=ilb,iub-1

        dmin = min(d(i,j),d(i,j+1),d(i+1,j),d(i+1,j+1))
        dmax = max(d(i,j),d(i,j+1),d(i+1,j),d(i+1,j+1))

        ! BTA: if (dmax.ge.z(1) .and. dmin.le.z(nc)) then
        if (dmax < z(1) .OR. dmin > z(nc)) CYCLE

        CONTOUR: do k=1,nc
           if (z(k) < dmin .or. z(k) > dmax) cycle CONTOUR
           do m=4,0,-1
              if (m.gt.0) then
                 h(m)=d(i+im(m),j+jm(m))-z(k)
                 xh(m)=x(i+im(m))
                 yh(m)=y(j+jm(m))
              else
                 h(0)=0.25*(h(1)+h(2)+h(3)+h(4))
                 xh(0)=0.5*(x(i)+x(i+1))
                 yh(0)=0.5*(y(j)+y(j+1))
              endif
              if (h(m).gt.0.0) then
                 sh(m)=+1
              else if (h(m).lt.0.0) then
                 sh(m)=-1
              else
                 sh(m)=0
              endif
           end do

           !
           !     Note: at this stage the relative heights of the corners and the
           !     centre are in the h array, and the corresponding coordinates are
           !     in the xh and yh arrays. The centre of the box is indexed by 0
           !     and the 4 corners by 1 to 4 as shown below.
           !     Each triangle is then indexed by the parameter m, and the 3
           !     vertices of each triangle are indexed by parameters m1,m2,and m3.
           !     It is assumed that the centre of the box is always vertex 2 though
           !     this isimportant only when all 3 vertices lie exactly on the same
           !     contour level, in which case only the side of the box is drawn.
           !
           !
           !           vertex 4 +-------------------+ vertex 3
           !                    | \               / |
           !                    |   \    m-3    /   |
           !                    |     \       /     |
           !                    |       \   /       |
           !                    |  m=2    X   m=2   |       the centre is vertex 0
           !                    |       /   \       |
           !                    |     /       \     |
           !                    |   /    m=1    \   |
           !                    | /               \ |
           !           vertex 1 +-------------------+ vertex 2
           !
           !
           !
           !                    Scan each triangle in the box
           !
           DO60:  do m=1,4
              m1=m
              m2=0
              m3=1
              if (m .ne. 4) m3 = m + 1

              case = castab(sh(m1),sh(m2),sh(m3))

              select case (case)

              case(0)
                 cycle DO60

              case (1) ! Case 1 - Line between vertices 1 and 2
                 x1=xh(m1)
                 y1=yh(m1)
                 x2=xh(m2)
                 y2=yh(m2)

              case (2) ! Case 2 - Line between vertices 2 and 3
                 x1=xh(m2)
                 y1=yh(m2)
                 x2=xh(m3)
                 y2=yh(m3)

              case (3) ! Case 3 - Line between vertices 3 and 1
                 x1=xh(m3)
                 y1=yh(m3)
                 x2=xh(m1)
                 y2=yh(m1)

              case (4) ! Case 4 - Line between vertex 1 and side 2-3
                 x1=xh(m1)
                 y1=yh(m1)
                 x2=xsect(m2,m3)
                 y2=ysect(m2,m3)

              case (5) ! Case 5 - Line between vertex 2 and side 3-1
                 x1=xh(m2)
                 y1=yh(m2)
                 x2=xsect(m3,m1)
                 y2=ysect(m3,m1)

              case (6) ! Case 6 - Line between vertex 3 and side 1-2
                 x1=xh(m3)
                 y1=yh(m3)
                 x2=xsect(m1,m2)
                 y2=ysect(m1,m2)

              case (7) ! Case 7 - Line between sides 1-2 and 2-3
                 x1=xsect(m1,m2)
                 y1=ysect(m1,m2)
                 x2=xsect(m2,m3)
                 y2=ysect(m2,m3)

              case (8) ! Case 8 - Line between sides 2-3 and 3-1
                 x1=xsect(m2,m3)
                 y1=ysect(m2,m3)
                 x2=xsect(m3,m1)
                 y2=ysect(m3,m1)

              case (9) ! Case 9 - Line between sides 3-1 and 1-2
                 x1=xsect(m3,m1)
                 y1=ysect(m3,m1)
                 x2=xsect(m1,m2)
                 y2=ysect(m1,m2)

              case default
                 print *, " Something went wrong - IMPOSIBLE CASE: ", case
                 call exit(1)

              end select

              call vecout(x1,y1,x2,y2,z(k))

           end do DO60  ! 60 continue
        end do CONTOUR
     end do RX ! 90  continue
  end do RY ! 100 continue

  return

end subroutine conrec


!======================================================================
!
!     This is a sample vector output routine. For a local environment
!     either replace the VECOUT call in the main line, or better
!     replace the contents of this subroutine between the *'s shown.
!
!     There is often the requirement to distinguish each contour
!     line with a different colour or a different line style. This
!     can be done in many ways using the contour values z for a
!     particular line segment.
!
subroutine vecout(x1,y1,x2,y2,z)
  implicit none
  real, intent(in) :: x1,y1,x2,y2,z

  !BTA just write the args:
  write(*,100) x1,y1,x2,y2,z
100 format(4F12.4,F14.3)
100 format(4F12.4,F14.3)

end subroutine vecout
