      Subroutine getmap(old,new,nold,nnew,map)

      IMPLICIT None
      INTEGER, intent(in) :: nold,nnew
      REAL(8), DIMENSION(nold,3), intent(in) :: old
      REAL(8), DIMENSION(nnew,3), intent(in) :: new
      REAL(8), DIMENSION(nold,2), intent(out) :: map
      REAL(8), DIMENSION(3) :: rij
      REAL(8) :: r, thresh
      INTEGER :: i,j

      thresh=1.d-3
      do i=1,nold
          do j=1,nnew
              rij=(old(i,:)-new(j,:))
              r=dsqrt(rij(1)**2.0+rij(2)**2.0+rij(3)**2.0)

              if (r < thresh) then
                  map(i,1)=i
                  map(i,2)=j
              end if

          end do
      end do

      end Subroutine getmap
