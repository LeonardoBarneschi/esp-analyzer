      Subroutine getpot(atc,atq,gridc,nat,ngrid,pot)

      IMPLICIT None
      INTEGER, intent(in) :: nat,ngrid
      REAL(8), DIMENSION(nat), intent(in) :: atq
      REAL(8), DIMENSION(nat,3), intent(in) :: atc
      REAL(8), DIMENSION(ngrid,3), intent(in) :: gridc
      REAL(8), DIMENSION(ngrid), intent(out) :: pot
      REAL(8), DIMENSION(3) :: rij
      REAL(8) :: r, ang2au
      INTEGER :: i,j

      ang2au=1/0.52917721
      pot=0.d0
      do i=1,nat
          do j=1,ngrid
              rij=(atc(i,:)-gridc(j,:))
              r=dsqrt(rij(1)**2.0+rij(2)**2.0+rij(3)**2.0)
              r=r*ang2au
              pot(j)=pot(j)+atq(i)/r
          end do
      end do

      end Subroutine getpot
