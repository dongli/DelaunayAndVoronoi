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

end module CommonUtils

