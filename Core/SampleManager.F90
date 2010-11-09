! *************************************************************************** !
! SampleManager module                                                        !
!                                                                             !
! Description:                                                                !
!   This module manages sample data structure and provides interface for      !
!   other modules.
!                                                                             !
! Author:                                                                     !
!   DONG Li, dongli@lasg.iap.ac.cn                                            !
!                                                                             !
! *************************************************************************** !

module SampleManager

    use SphereService

    implicit none

    integer :: numSample = 0

    type Sample
        real(8) lon, lat
        real(8) x, y, z ! transformed Cartesian coordinates
        type(Sample), pointer :: prev => null()
        type(Sample), pointer :: next => null()
    end type

    type(Sample), pointer :: SMPHead

contains
    
    ! ************************************************************************ !
    ! SampleManager_Init
    ! Purpose:
    !   Create the doubled linked list of samples.
    ! ************************************************************************ !

    subroutine SampleManager_Init(n)
        integer, intent(in) :: n
    
        integer i
        type(Sample), pointer :: SMPPtr1, SMPPtr2

        numSample = n

        allocate(SMPHead)
        SMPPtr1 => SMPHead
        SMPPtr2 => SMPHead
        do i = 2, numSample
            allocate(SMPPtr1%next)
            SMPPtr1 => SMPPtr1%next
            SMPPtr1%prev => SMPPtr2
            SMPPtr2%next => SMPPtr1
            SMPPtr2 => SMPPtr1
        end do

    end subroutine SampleManager_Init

    subroutine SampleManager_Set(lon, lat)
        real(8), intent(in) :: lon(:), lat(:)

        integer i, dimSize(1)
        type(Sample), pointer :: SMPPtr1

        dimSize = shape(lon)
        if (dimSize(1) /= numSample) then
            ! Complain
        end if
        dimSize = shape(lat)
        if (dimSize(1) /= numSample) then
            ! Complain
        end if

        SMPPtr1 => SMPHead
        do i = 1, numSample
            SMPPtr1%lon = lon(i)
            SMPPtr1%lat = lat(i)
            call CartesianTransformOnUnitSphere(SMPPtr1%lon, SMPPtr1%lat, &
                SMPPtr1%x, SMPPtr1%y, SMPPtr1%z)
            SMPPtr1 => SMPPtr1%next
        end do
   

    end subroutine SampleManager_Set
    

end module SampleManager

