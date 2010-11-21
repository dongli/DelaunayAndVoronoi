! *************************************************************************** !
! CommonUtils module                                                          !
!                                                                             !
! Description:                                                                !
!                                                                             !
! Author:                                                                     !
!                                                                             !
! *************************************************************************** !

module CommonUtils

    implicit none

    interface CommonUtils_Sort
        module procedure SortInteger
    end interface CommonUtils_Sort

contains

    subroutine SortInteger(array)
        integer, intent(inout) :: array(:)

        integer dimSize(1), temp, i, j

        dimSize = shape(array)

        do i = 1, dimSize(1)-1
            do j = i+1, dimSize(1)
                if (array(i) > array(j)) then
                    temp = array(i)
                    array(i) = array(j)
                    array(j) = temp
                end if
            end do
        end do

    end subroutine SortInteger

    subroutine CrossProduct(x1, y1, z1, x2, y2, z2, x, y, z)
        real(RealKind), intent(in) :: x1, y1, z1, x2, y2, z2
        real(RealKind), intent(out) :: x, y, z

        x = y1*z2-z1*y2
        y = z1*x2-x1*z2
        z = x1*y2-y1*x2

    end subroutine CrossProduct

end module CommonUtils
