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

    use RunManager
    use NFWrap
#if (defined TTS_Online)
    use MeshManager
#endif
    use SphereService

    implicit none

    private

    public SampleManager_Init
    public SampleManager_Set
    public SampleManager_NewDimVar
    public SampleManager_Output
    public SampleManager_SaveSpatialBound

    public Point, Sample

    public numSample
    public smpHead

    ! ======================================================================= !
    type Point
        real(RealKind) lon, lat
        real(RealKind) x, y, z
    end type Point

    type SpatialBound
        integer :: numVtx = 0
        type(Point), allocatable :: vtx(:)
        real(RealKind) area
    contains
        procedure :: dump => SpatialBound_dump
        procedure :: clone => SpatialBound_clone
        procedure :: getVertex => SpatialBound_getVertex
        procedure :: setVertex => SpatialBound_setVertex
    end type SpatialBound

    type, extends(Point) :: Sample
        integer :: id = -1
#if (defined TTS_Online)
        ! Location
        type(Location) loc ! The location where sample is at the mesh
#endif
        ! Spatial bound
        type(SpatialBound) bnd, oldBnd
        !
        type(Sample), pointer :: prev => null()
        type(Sample), pointer :: next => null()
    contains
        procedure :: getCoordinate => Sample_getCoordinate
        procedure :: setCoordinate => Sample_setCoordinate
        procedure :: CartesianTransform => Sample_CartesianTransform
#if (defined TTS_Online)
        procedure :: getLocation => Sample_getLocation
        procedure :: setLocation => Sample_setLocation
#endif
        procedure :: dump => Sample_dump
    end type

    integer :: numSample = 0
    type(Sample), pointer :: smpHead => null()

    integer, parameter :: maxNumVtx = 20

contains
    
    ! ************************************************************************ !
    ! SampleManager_Init
    ! Purpose:
    !   Create the doubled linked list of samples.
    ! ************************************************************************ !

    subroutine SampleManager_Init(n)
        integer, intent(in) :: n
    
        integer i
        type(Sample), pointer :: smp1, smp2

#if (!defined FullSpeed)
        call MsgManager_RecordSpeaker("SampleManager_Init")
#endif

        numSample = n

        do i = 1, numSample
            if (i == 1) then
                allocate(smpHead)
                smp1 => smpHead
            else
                smp2 => smp1
                allocate(smp1%next)
                smp1 => smp1%next
                smp1%prev => smp2
                smp2%next => smp1
            end if
            smp1%id = i
            smp1%bnd%numVtx = 1
            allocate(smp1%bnd%vtx(1))
            smp1%oldBnd%numVtx = 1
            allocate(smp1%oldBnd%vtx(1))
        end do

#if (!defined FullSpeed)
        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker
#endif

    end subroutine SampleManager_Init

    ! ************************************************************************ !
    ! SampleManager_Set                                                        !
    ! Purpose:                                                                 !
    !   Set the locations of each sample.                                      !
    ! ************************************************************************ !

    subroutine SampleManager_Set(lon, lat)
        real(RealKind), intent(in) :: lon(:), lat(:)

        integer i, dimSize(1)
        type(Sample), pointer :: smp1

#if (!defined FullSpeed)
        call MsgManager_RecordSpeaker("SampleManager_Set")
#endif

        dimSize = shape(lon)
        if (dimSize(1) /= numSample) then
            ! Complain
        end if
        dimSize = shape(lat)
        if (dimSize(1) /= numSample) then
            ! Complain
        end if

        smp1 => smpHead
        do i = 1, numSample
            smp1%lon = lon(i)
            smp1%lat = lat(i)
#if (defined TTS_Online)
            call MeshManager_LocationCheck([lon(i),lat(i)], smp1%loc)
#endif
            call smp1%CartesianTransform
            smp1 => smp1%next
        end do

#if (!defined FullSpeed)
        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker
#endif

    end subroutine SampleManager_Set

    subroutine SampleManager_SaveSpatialBound
        type(Sample), pointer :: smp
        integer i

        smp => smpHead
        do i = 1, numSample
            call smp%oldBnd%clone(smp%bnd)
            smp => smp%next
        end do
            
    end subroutine SampleManager_SaveSpatialBound

    subroutine SpatialBound_getVertex(bnd, i, x)
        class(SpatialBound), intent(in) :: bnd
        integer, intent(in) :: i
        real(RealKind), intent(out) :: x(2)

        x = [bnd%vtx(i)%lon,bnd%vtx(i)%lat]

    end subroutine SpatialBound_getVertex

    subroutine SpatialBound_setVertex(bnd, i, x)
        class(SpatialBound), intent(inout) :: bnd
        integer, intent(in) :: i
        real(RealKind), intent(in) :: x(2)

        bnd%vtx(i)%lon = x(1)
        bnd%vtx(i)%lat = x(2)
    
    end subroutine SpatialBound_setVertex
    
    subroutine Sample_getCoordinate(smp, x)
        class(Sample), intent(in) :: smp
        real(RealKind), intent(out) :: x(2)

        x = [smp%lon,smp%lat]
    
    end subroutine Sample_getCoordinate

    subroutine Sample_setCoordinate(smp, x)
        class(Sample), intent(inout) :: smp
        real(RealKind), intent(in) :: x(2)

        smp%lon = x(1)
        smp%lat = x(2)
        
    end subroutine Sample_setCoordinate

    subroutine Sample_CartesianTransform(smp)
        class(Sample), intent(inout) :: smp

        call CartesianTransformOnUnitSphere( &
            smp%lon, smp%lat, smp%x, smp%y, smp%z)
    
    end subroutine Sample_CartesianTransform
    
#if (defined TTS_Online)
    subroutine Sample_getLocation(smp, loc)
        class(Sample), intent(in) :: smp
        type(Location), intent(out) :: loc
    
        loc = smp%loc

    end subroutine Sample_getLocation

    subroutine Sample_setLocation(smp, loc)
        class(Sample), intent(inout) :: smp
        type(Location), intent(in) :: loc

        smp%loc = loc

    end subroutine Sample_setLocation
#endif

    subroutine SpatialBound_dump(bnd)
        class(SpatialBound), intent(in) :: bnd

        integer i

        write(*, "('Spatial bound')")
        write(*, "('  Area: ', F30.10)") bnd%area
        write(*, "('  Vertices:')")
        do i = 1, bnd%numVtx
            write(*, "('  ', 2F20.15)") &
                bnd%vtx(i)%lon*Rad2Deg, bnd%vtx(i)%lat*Rad2Deg
        end do
    
    end subroutine SpatialBound_dump
    
    subroutine SpatialBound_clone(a, b)
        class(SpatialBound), intent(inout) :: a
        type(SpatialBound), intent(in) :: b

        if (a%numVtx /= b%numVtx) then
            deallocate(a%vtx)
            a%numVtx = b%numVtx
            allocate(a%vtx(a%numVtx))
        end if
        a%vtx%lon = b%vtx%lon
        a%vtx%lat = b%vtx%lat
        a%area = b%area
    
    end subroutine SpatialBound_clone
    
    subroutine SampleManager_NewDimVar(fcard)
        type(FileCard), intent(inout) :: fcard
    
        call NFWrap_NewDim(fcard, "numSample", numSample)
        call NFWrap_NewDim(fcard, "maxNumVtx", maxNumVtx)
        call NFWrap_New1DVar(fcard, "lon", "double", "numSample", &
            "Sample longitude", "radian", .true.)
        call NFWrap_New1DVar(fcard, "lat", "double", "numSample", &
            "Sample latitude", "radian", .true.)
        call NFWrap_New1DVar(fcard, "numBndVtx", "integer", "numSample", &
            "Sample bound vertex number", "", .true.)
        call NFWrap_New2DVar(fcard, "lonBndVtx", "double", "maxNumVtx", "numSample", &
            "Sample bound vertex longitude", "radian", .true.)
        call NFWrap_New2DVar(fcard, "latBndVtx", "double", "maxNumVtx", "numSample", &
            "Sample bound vertex latitude", "radian", .true.)
        call NFWrap_New1DVar(fcard, "areaBnd", "double", "numSample", &
            "Sample bound area", "m2", .true.)
        call NFWrap_New1DVar(fcard, "numVoroVtx", "integer", "numSample", &
            "Sample Voronoi cell vertex number", "", .true.)
        call NFWrap_New2DVar(fcard, "lonVoroVtx", "double", "maxNumVtx", "numSample", &
            "Sample Voronoi cell vertex longitude", "radian", .true.)
        call NFWrap_New2DVar(fcard, "latVoroVtx", "double", "maxNumVtx", "numSample", &
            "Sample Voronoi cell vertex latitude", "radian", .true.)
        call NFWrap_New1DVar(fcard, "areaVoro", "double", "numSample", &
            "Sample Voronoi cell area", "m2", .true.)

    end subroutine SampleManager_NewDimVar
    
    subroutine SampleManager_Output(fcard)
        type(FileCard), intent(inout) :: fcard

        type(Sample), pointer :: smp
        real(RealKind) lon(numSample), lat(numSample)
        integer numBndVtx(numSample)
        real(RealKind) lonBndVtx(maxNumVtx,numSample)
        real(RealKind) latBndVtx(maxNumVtx,numSample)
        real(RealKind) areaBnd(numSample)
        integer numVoroVtx(numSample)
        real(RealKind) lonVoroVtx(maxNumVtx,numSample)
        real(RealKind) latVoroVtx(maxNumVtx,numSample)
        real(RealKind) areaVoro(numSample)
        integer i, j

        smp => smpHead
        do i = 1, numSample
            lon(i) = smp%lon
            lat(i) = smp%lat
            numBndVtx(i) = smp%oldBnd%numVtx
            areaBnd(i) = smp%oldBnd%area
            do j = 1, smp%oldBnd%numVtx
                lonBndVtx(j,i) = smp%oldBnd%vtx(j)%lon
                latBndVtx(j,i) = smp%oldBnd%vtx(j)%lat
            end do
            if (smp%bnd%numVtx > maxNumVtx) then
                call MsgManager_Speak(Error, "maxNumVtx has been exceeded.")
                print *, "The dirty vertices:"
                do j = 1, smp%bnd%numVtx
                    print *, smp%bnd%vtx(j)%lon*Rad2Deg, smp%bnd%vtx(j)%lat*Rad2Deg
                end do
                call RunManager_EndRun
            end if
            numVoroVtx(i) = smp%bnd%numVtx
            areaVoro(i) = smp%bnd%area
            do j = 1, smp%bnd%numVtx
                lonVoroVtx(j,i) = smp%bnd%vtx(j)%lon
                latVoroVtx(j,i) = smp%bnd%vtx(j)%lat
            end do
            smp => smp%next
        end do

        call NFWrap_Output1DVar(fcard, "lon", lon)
        call NFWrap_Output1DVar(fcard, "lat", lat)
        call NFWrap_Output1DVar(fcard, "numBndVtx", numBndVtx)
        call NFWrap_Output2DVar(fcard, "lonBndVtx", lonBndVtx)
        call NFWrap_Output2DVar(fcard, "latBndVtx", latBndVtx)
        call NFWrap_Output1DVar(fcard, "areaBnd", areaBnd)
        call NFWrap_Output1DVar(fcard, "numVoroVtx", numVoroVtx)
        call NFWrap_Output2DVar(fcard, "lonVoroVtx", lonVoroVtx)
        call NFWrap_Output2DVar(fcard, "latVoroVtx", latVoroVtx)
        call NFWrap_Output1DVar(fcard, "areaVoro", areaVoro)

    end subroutine SampleManager_Output
 
    subroutine Sample_dump(smp)
        class(Sample), intent(in) :: smp

        write(*, "('Sample ID - ')", advance="no")
        write(*, *) smp%id
        write(*, "('  Spherical coordinate: ', 2F20.15)") &
            smp%lon*Rad2Deg, smp%lat*Rad2Deg
        write(*, "('  Cartesian coordinate: ', 3F20.15)") &
            smp%x, smp%y, smp%z

    end subroutine Sample_dump

end module SampleManager

