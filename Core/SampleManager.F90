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

    private

    public SampleManager_Init
    public SampleManager_Set

    public Sample

    public numSample
    public SMPHead

    integer :: numSample = 0

    type Sample
        real(RealKind) lon, lat
        real(RealKind) x, y, z ! transformed Cartesian coordinates
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
        type(Sample), pointer :: SMP1, SMP2

        numSample = n

        allocate(SMPHead)
        SMP1 => SMPHead
        SMP2 => SMPHead
        do i = 2, numSample
            allocate(SMP1%next)
            SMP1 => SMP1%next
            SMP1%prev => SMP2
            SMP2%next => SMP1
            SMP2 => SMP1
        end do

    end subroutine SampleManager_Init

    subroutine SampleManager_Set(lon, lat)
        real(RealKind), intent(in) :: lon(:), lat(:)

        integer i, dimSize(1)
        type(Sample), pointer :: SMP1

        dimSize = shape(lon)
        if (dimSize(1) /= numSample) then
            ! Complain
        end if
        dimSize = shape(lat)
        if (dimSize(1) /= numSample) then
            ! Complain
        end if

        SMP1 => SMPHead
        do i = 1, numSample
            SMP1%lon = lon(i)
            SMP1%lat = lat(i)
            call CartesianTransformOnUnitSphere(SMP1%lon, SMP1%lat, &
                SMP1%x, SMP1%y, SMP1%z)
            SMP1 => SMP1%next
        end do

    end subroutine SampleManager_Set

end module SampleManager

