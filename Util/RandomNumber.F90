! *************************************************************************** !
! RandomNumber module                                                         !
!                                                                             !
! Description:                                                                !
!                                                                             !
! Author:                                                                     !
!                                                                             !
! *************************************************************************** !

module RandomNumber

    implicit none

    private

    public RandomNumber_Start
    public RandomNumber_Get

    integer n
    integer, allocatable :: seed(:)

    interface RandomNumber_Get
        module procedure RandomNumber_Get_Double
    end interface RandomNumber_Get

contains

    subroutine RandomNumber_Start
        integer clock, i

        call random_seed(size=n)
        allocate(seed(n))

        call system_clock(count=clock)
        seed = clock+37*[(i-1, i = 1, n)]

        call random_seed(put=seed)

    end subroutine RandomNumber_Start

    subroutine RandomNumber_Get_Double(a, b, r)
        real(8), intent(in) :: a, b
        real(8), intent(out) :: r

        call random_number(r)
        r = (r*(b-a))+a
    
    end subroutine RandomNumber_Get_Double
    
end module RandomNumber

