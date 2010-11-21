! *************************************************************************** !
! DelaunayAndVoronoi module                                                   !
!                                                                             !
! Description:                                                                !
!   This module implements the construction of Delaunay triangulation (DT)    !
!   and Voronoi diagram (VD) on the sphere, with the requirement that all the !
!   vertices can move and DT/VD must be updated efficiently.                  !
!                                                                             !
!   Incremental algorithm is chosen for constructing "dynamic" Delaunay       !
!   triangulation. The efficiency is on the top of the priority list.         !
!                                                                             !
!   "Triangle-based" data structure is chosen because it is more efficient in !
!   storage and calculation than quad-edge.                                   !
!                                                                             !
!   Three "virtual vertices" are added for initially creating the triangles.  !
!   If the real vertices are only covered less than one hemisphere, the       !
!   virtual vertices will not be deleted due to the crossing of the edges,    !
!   where the  intersection points are exactly the virtual vertices.          !
!                                                                             !
!   An instantaneous history of triangle construction is been kept as         !
!   "obsolete" triangles and "temporal" triangles. The former one is a list   !
!   of existing triangles which is no longer valid (unsatisfying Delaunay     !
!   criterion), and the latter one is a list of temporally created triangles  !
!   that has not been recorded into main triangle list but are also invalid.  !
!   With "obsolete" triangle list, the update of point-in-triangle relation   !
!   can be restricted locally.                                                !
!                                                                             !
!   In addition to the pointer of "contained triangle" for each uninserted    !
!   vertex, a double linked lists of "included vertices" for each existing    !
!   triangle are also created and maintained to assist in restricting the     !
!   area of aforementioned update.                                            !
!                                                                             !
! Author:                                                                     !
!   DONG Li, dongli@lasg.iap.ac.cn                                            !
!                                                                             !
! Development log:                                                            !
!   - Nov 1 2010 ~ Nov 4 2010                                                 !
!     Implement the basic functionality of Delaunay triangulation.            !
!                                                                             !
! *************************************************************************** !

module DelaunayAndVoronoi

    use RunManager
    use MsgManager
    use FloatingPoint
    use CommonUtils
    use RandomNumber
    use SphereService
    use SampleManager
    use NFWrap

    implicit none

    private

    public DelaunayAndVoronoi_LinkSample
    public DelaunayAndVoronoi_CleanAll
    public ConstructDelaunayTriangulation
    public ExtractVoronoiDiagram
    public DelaunayAndVoronoi_Report
    public DelaunayAndVoronoi_Output
    public DelaunayAndVoronoi_OutputDelaunay

    ! ======================================================================= !
    ! Derived data types for constructing complicated data structure :)
    type DelaunayVertexPointer
        type(DelaunayVertex), pointer :: ptr => null()
    end type DelaunayVertexPointer

    type DelaunayTrianglePointer
        type(DelaunayTriangle), pointer :: ptr => null()
    end type DelaunayTrianglePointer

    type DelaunayVertexPointerList
        type(DelaunayVertex), pointer :: ptr => null()
        type(DelaunayVertexPointerList), pointer :: prev => null()
        type(DelaunayVertexPointerList), pointer :: next => null()
    end type DelaunayVertexPointerList

    type DelaunayTrianglePointerList
        type(DelaunayTriangle), pointer :: ptr => null()
        type(DelaunayTrianglePointerList), pointer :: prev => null()
        type(DelaunayTrianglePointerList), pointer :: next => null()
    end type DelaunayTrianglePointerList

    ! ======================================================================= !
    ! Delaunay triangulation data structure with double linked list
    type DelaunayVertex
        integer :: id = -1
        ! For geometry
        type(Sample), pointer :: smp => null()
        ! For topology
        integer :: numIcdDT = 0 ! Incident triangles
        type(DelaunayTrianglePointerList), pointer :: icdDTHead => null()
        integer :: numLinkDVT = 0 ! Link vertices
        type(DelaunayVertexPointerList), pointer :: linkDVTHead => null()
        ! For point-in-triangle relation
        type(DelaunayTriangle), pointer :: cntDT1
        type(DelaunayVertexPointerList), pointer :: stub1 ! For fast deletion
        type(DelaunayTriangle), pointer :: cntDT2
        type(DelaunayVertexPointerList), pointer :: stub2 ! For fast deletion
        integer :: edgeIdx = -1 ! The shared edge index if InTriangle return >0.
        ! For double link
        type(DelaunayVertex), pointer :: prev => null()
        type(DelaunayVertex), pointer :: next => null()
    end type DelaunayVertex

    type DelaunayTriangle
        integer :: id = -1
#if (defined NextWork)
        ! For geometry
        real(RealKind) edgeArc(3)
        real(RealKind) cornerAngle(3)
#endif
        ! For topology
        type(DelaunayVertexPointer) DVT(3)
        type(DelaunayTrianglePointer) adjDT(3) ! Three adjacent DTs
        type(Point) cnt ! Center of circumcircle
        real(RealKind) radius ! Radius of circumcircle in radian
        ! For validation
        integer :: numSubDT = 0 ! the number of subdivided triangles
        type(DelaunayTrianglePointer) subDT(3) ! Subdivided DTs (at most three)
        ! For local update of point-in-triangle relation
        integer :: numIncDVT = 0 ! the number of included vertices 
        type(DelaunayVertexPointerList), pointer :: incDVTHead => null()
        type(DelaunayVertexPointerList), pointer :: incDVTCurr => null()
        ! For double link
        type(DelaunayTriangle), pointer :: prev => null()
        type(DelaunayTriangle), pointer :: next => null()
    end type DelaunayTriangle

    ! ======================================================================= !
    ! Voronoi diagram data structure
    type VoronoiCell
        integer :: id = -1
        type(Sample), pointer :: smp => null()
        type(VoronoiCell), pointer :: prev => null()
        type(VoronoiCell), pointer :: next => null()
    end type VoronoiCell

    ! Data
    ! Delaunay vertex
    integer :: numDVT = 0
    type(DelaunayVertex), pointer :: DVTHead
    integer :: numVirtualDVT = 0
    type(DelaunayVertex), pointer :: VirtualDVTHead => null()
    ! Delaunay triangle
    integer :: numDT = 0, numTotalDT = 0
    type(DelaunayTriangle), pointer :: DTHead, DTCurr
    ! Obsolete DT record list
    ! Purpose: 1) Record the obsolete DTs for deletion;
    !          2) For "locally" updating of point-in-triangle relation
    integer :: numObsDT = 0
    type(DelaunayTrianglePointerList), pointer :: obsDTHead
    ! Temporal DT record list
    ! Purpose: 1) Record the temporal DTs for deletion;
    !          2) For updating point-in-triangle relation
    integer :: numTmpDT = 0 ! Record list of temporal DT
    type(DelaunayTrianglePointerList), pointer :: tmpDTHead
    ! Voronoi cell
    integer :: numVC = 0
    type(VoronoiCell), pointer :: VCHead

    ! ======================================================================= !
    ! Constant
    integer, parameter :: OrientLeft = 1
    integer, parameter :: OrientRight = 2
    integer, parameter :: OrientOn = 3
    integer, parameter :: InsideTriangle = 4
    integer, parameter :: OutsideTriangle = -4
    integer, parameter :: InsideCircle = 1
    integer, parameter :: OutsideCircle = 2
    integer, parameter :: OnTheCircle = 3

    ! ======================================================================== !
    ! Auxiliary data
    type(FileCard) fcard
    integer maxNumVVT
    integer, parameter :: im1(3) = [3,1,2], ip1(3) = [2,3,1]

    interface InCircle
        module procedure InCircle1, InCircle2
    end interface InCircle

#if (defined DelaunayAndVoronoi_Debug)
    public debug, doRandom, idx1, idx2, idx3
    logical :: debug
    logical :: printDet = .false. ! Switch on/off debug print of determinant.
    logical :: doRandom = .true.  ! Switch on/off randomized indices
    integer :: idx1 = 833         ! Fix the first inserted DVTs for debugging.
    integer :: idx2 = 1204
    integer :: idx3 = 1784
#endif

contains

    ! ######################################################################## !
    ! ########################### Public interface ########################### !
    ! ######################################################################## !

    ! ************************************************************************ !
    ! DelaunayAndVoronoi_LinkSample                                            !
    ! Purpose:                                                                 !
    !   The coordinates of the vertices are stored in samples of SampleManager !
    !   module, thus different modules can access the samples. This subroutine !
    !   links the vertices with the samples.                                   !
    ! ************************************************************************ !

    subroutine DelaunayAndVoronoi_LinkSample
        type(Sample), pointer :: smp
        type(DelaunayVertex), pointer :: DVT1, DVT2
        type(VoronoiCell), pointer :: VC1, VC2
        integer i

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_RecordSpeaker("DelaunayAndVoronoi_LinkSample")
#endif

        ! Initialize the Delaunay vertices,
        ! and link them with corresponding samples ...................... STEP 1
        numDVT = numSample
        smp => smpHead
        do i = 1, numDVT
            if (i == 1) then
                allocate(DVTHead)
                DVT1 => DVTHead
            else
                allocate(DVT1%next)
                DVT2 => DVT1
                DVT1 => DVT1%next
                DVT1%prev => DVT2
                DVT2%next => DVT1
            end if
            DVT1%id = i
            DVT1%smp => smp
            DVT1%numIcdDT = 1
            allocate(DVT1%icdDTHead)
            smp => smp%next
        end do

        ! Initialize the Voronoi cells .................................. STEP 2
        numVC = numSample
        smp => smpHead
        do i = 1, numVC
            if (i == 1) then
                allocate(VCHead)
                VC1 => VCHead
            else
                VC2 => VC1
                allocate(VC1%next)
                VC1 => VC1%next
                VC1%prev => VC2
                VC2%next => VC1
            end if
            VC1%id = i
            VC1%smp => smp
            smp => smp%next
        end do

        call RandomNumber_Start

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker
#endif

    end subroutine DelaunayAndVoronoi_LinkSample

    ! ************************************************************************ !
    ! ConstructDelaunayTriangulation                                           !
    ! Purpose:                                                                 !
    !   This subroutine is only called at the initial time step. When the      !
    !   samples are ready, construct Delaunay triangulation from them.         !
    ! ************************************************************************ !


    subroutine ConstructDelaunayTriangulation
        type(DelaunayVertex), pointer :: DVT1, DVT2, DVT3
        type(DelaunayTriangle), pointer :: DT1, DT2
        type(DelaunayVertexPointerList), pointer :: incDVT
        type(DelaunayTrianglePointerList), pointer :: obsDT, icdDT
        logical :: found, inserted(numDVT)
        integer i, j, k, l, ret, idx(3)

#if (!defined DelaunayAndVoronoi_FullSpeed || defined DelaunayAndVoronoi_Debug)
        call MsgManager_RecordSpeaker("ConstructDelaunayTriangulation")
#endif

        inserted = .false.

        ! Randomly choose the first three vertices to insert ............ STEP 0
        call GetTheFirstThreeDVTIndex(idx)
        DVT1 => DVTHead
        do i = 1, idx(1)-1
            DVT1 => DVT1%next
        end do
        DVT2 => DVT1%next
        do i = idx(1)+1, idx(2)-1
            DVT2 => DVT2%next
        end do
        DVT3 => DVT2%next
        do i = idx(2)+1, idx(3)-1
            DVT3 => DVT3%next
        end do

        ! Initialize the triangles comprised of the first three vertices and
        ! their "antipodal" counterparts ................................ STEP 1
        ! Note: On the sphere, the boundary condition is periodic, so there
        !       are seven other virtual triangles accompanying the real one.
        call InitDelaunayTriangle(DVT1, DVT2, DVT3)
        inserted(idx(1)) = .true.
        inserted(idx(2)) = .true.
        inserted(idx(3)) = .true.

        ! Initialize the point-in-triangle relation between the rest vertices
        ! and the created triangles ..................................... STEP 2
        DVT1 => DVTHead
        do i = 1, numDVT
            if (inserted(i)) then
                DVT1 => DVT1%next
                cycle
            end if
            DT1 => DTHead
            do j = 1, numDT
                ret = InTriangleRelaxed(DT1, DVT1)
                if (ret == InsideTriangle) then
                    ! DVT is included in some DT ........................ CASE 1
                    DVT1%cntDT1 => DT1
                    DVT1%cntDT2 => null()
                    ! Also record the DVT into the list of included vertices
                    ! of DT
                    call RecordIncludedVertex(DVT1, DT1)
                    exit
                else if (ret == OutsideTriangle) then
                    DT1 => DT1%next
                    cycle
                else if (ret > 0) then
                    ! DVT is on some edge of DT ......................... CASE 2
                    DT2 => DT1%adjDT(ret)%ptr
                    DVT1%cntDT1 => DT1
                    DVT1%cntDT2 => DT2
                    DVT1%edgeIdx = ret
                    call RecordIncludedVertex(DVT1, DT1, DT2)
                    exit
                else if (ret < 0) then
                    ! DVT coincides with some vertex of DT .............. CASE 3
                    DVT2 => DT1%DVT(-ret)%ptr
                    if (DVT2%id < 0) then
                        ! The vertex is virtual ....................... CASE 3-1
#if (!defined DelaunayAndVoronoi_FullSpeed || defined DelaunayAndVoronoi_Debug)
                        call MsgManager_Speak(Notice, "Virtual DVT "// &
                            trim(int2str(DVT2%id))//" coincides with DVT "// &
                            trim(int2str(DVT1%id))//", so make replacement.")
#endif
                        ! It is ok, just replace it with DVT1.
                        call ExtractIncidentDTAndLinkDVT(DVT2)
                        DVT1%icdDTHead%ptr => DVT2%icdDTHead%ptr
                        icdDT => DVT2%icdDTHead
                        do k = 1, DVT2%numIcdDT
                            do l = 1, 3
                                if (associated(icdDT%ptr%DVT(l)%ptr, DVT2)) exit
                            end do
                            icdDT%ptr%DVT(l)%ptr => DVT1
                            icdDT => icdDT%next
                        end do
                        inserted(i) = .true.
                        ! Delete it from virtual DVT list.
                        DVT3 => VirtualDVTHead
                        do k = 1, numVirtualDVT
                            if (associated(DVT2, DVT3)) then
                                if (associated(DVT3%prev)) then
                                    DVT3%prev%next => DVT3%next
                                else
                                    VirtualDVTHead => VirtualDVTHead%next
                                end if
                                if (associated(DVT3%next)) then
                                    DVT3%next%prev => DVT3%prev
                                end if
                                deallocate(DVT3)
                                numVirtualDVT = numVirtualDVT-1
                                exit
                            end if
                            DVT3 => DVT3%next
                        end do
                        exit
                    else
                        ! The vertex is real .......................... CASE 3-2
                        ! Complain this degenerate case.
                        call MsgManager_Speak(Error, &
                            "DVT "//trim(int2str(DVT1%id))//" coincides to "// &
                            "DVT "//trim(int2str(-ret))//" of "// &
                            "DT "//trim(int2str(DT1%id))//".")
                        call RunManager_EndRun
                    end if
                end if
            end do
            DVT1 => DVT1%next
        end do

#if (!defined DelaunayAndVoronoi_FullSpeed || defined DelaunayAndVoronoi_Debug)
        call MsgManager_Speak(Notice, "Initial point-in-triangle relation is updated.")
#endif

        ! Insert the rest vertices one at a time ........................ STEP 3
        DVT1 => DVTHead
        do i = 1, numDVT
            if (inserted(i)) then
                ! If DVT has replaced some virtual DVT, just skip it.
                DVT1 => DVT1%next
                cycle
            end if
            ! Update the Delauny triangulation ........................ STEP 3.1
            call InsertVertex(DVT1)
#if (!defined DelaunayAndVoronoi_FullSpeed)
#if (defined DelaunayAndVoronoi_Debug && defined DelaunayAndVoronoi_Verbose)
            call MsgManager_Speak(Notice, "DVT "// &
                trim(int2str(DVT1%id))//" has been inserted.")
#endif
#endif
            ! Update the point-in-triangle relation "locally" ......... STEP 3.2
            obsDT => obsDTHead
            do j = 1, numObsDT
                do while (associated(obsDT%ptr%incDVTHead))
                    call UpdatePointInTriangle(obsDT%ptr%incDVTHead%ptr, obsDT%ptr, found)
#if (defined DelaunayAndVoronoi_Debug)
                    if (.not. found) then
                        print *, "Point", obsDT%ptr%incDVTHead%ptr%id, "does not found its new DT."
                        call PrintTriangleTopology(obsDT%ptr)
                        do k = 1, 3
                            call obsDT%ptr%DVT(k)%ptr%smp%dump
                        end do
                        call obsDT%ptr%incDVTHead%ptr%smp%dump
                        stop
                    end if
#endif
                end do
                obsDT => obsDT%next
            end do
            ! Delete obsolete and temporal triangles .................. STEP 3.3
            call DeleteObsoleteTriangle
            call DeleteTemporalTriangle
            ! Shift the current DVT pointer ........................... STEP 3.4
            DVT1 => DVT1%next
        end do

        ! Extract the full list of incident DTs and link DVTs ........... STEP 4
        DVT1 => DVTHead
        do i = 1, numDVT
            call ExtractIncidentDTAndLinkDVT(DVT1)
            DVT1 => DVT1%next
        end do
        DVT1 => VirtualDVTHead
        do i = 1, numVirtualDVT
            call ExtractIncidentDTAndLinkDVT(DVT1)
            DVT1 => DVT1%next
        end do

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_Speak(Notice, &
            "Full lists of incident DTs and link DVTs for each DVT are extraced.")
#endif

#if (!defined KeepVirtualDVT)
        DVT1 => VirtualDVTHead
        do while (associated(DVT1))
            call DeleteVertex(DVT1)
            numVirtualDVT = numVirtualDVT-1
            DVT2 => DVT1
            DVT1 => DVT1%next
            deallocate(DVT2)
        end do
#if (!defined DelaunayAndVoronoi_FullSpeed || defined DelaunayAndVoronoi_Debug)
        if (numVirtualDVT /= 0) then
            call MsgManager_Speak(Notice, trim(int2str(numVirtualDVT))// &
                " virtual DVTs are deleted.")
        end if
#endif
#endif

        ! Calculate the circumcenter and circumradius
        call CalcCircumcircle

#if (!defined DelaunayAndVoronoi_FullSpeed || defined DelaunayAndVoronoi_Debug)
        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker
#endif

    end subroutine ConstructDelaunayTriangulation

    ! ************************************************************************ !
    ! ExtractVoronoiDiagram                                                    !
    ! Purpose:                                                                 !
    !   Get the dual (Voronoi diagram) of Delaunay triangulation.              !
    ! ************************************************************************ !

    subroutine ExtractVoronoiDiagram
        type(DelaunayVertex), pointer :: DVT
        type(VoronoiCell), pointer :: VC
        type(DelaunayTrianglePointerList), pointer :: icdDT
#if (!defined NoFixNegativeArea)
        type(DelaunayTrianglePointerList), pointer :: icdDTNext
        type(DelaunayVertex), pointer :: linkDVT
        real(RealKind) dlon, dlat
        logical flag
#endif
        integer i, j
#if (defined DelaunayAndVoronoi_Debug)
        real(RealKind) :: totalArea = 0.0e0
#endif

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_RecordSpeaker("ExtractVoronoiDiagram")
#endif

        maxNumVVT = 0
#if (defined DelaunayAndVoronoi_Debug)
        totalArea = 0.0
#endif
        DVT => DVTHead
        VC => VCHead
        do i = 1, numDVT ! = numVC
#if (!defined NoFixNegativeArea)
            ! Handle the too-short-distance between circumcenters
            ! Strategy: Merge them!
            flag = .false.
            icdDT => DVT%icdDTHead
            do
                icdDTNext => icdDT%next
                dlon = icdDT%ptr%cnt%lon-icdDTNext%ptr%cnt%lon
                dlat = icdDT%ptr%cnt%lat-icdDTNext%ptr%cnt%lat
                if (abs(dlat) < 1.0e-6) then
                    if (PI05-abs(icdDT%ptr%cnt%lat) < 1.0e-6 .or. &
                        (abs(dlon) < 1.0e-6 .or. PI2-abs(dlon) < 1.0e-6)) then
                        do j = 1, 3
                            if (associated(icdDTNext%ptr%DVT(j)%ptr, DVT)) exit
                        end do
                        linkDVT => icdDTNext%ptr%DVT(ip1(j))%ptr
                        call MergeIncidentTriangle(DVT, icdDT%ptr, icdDTNext%ptr, icdDT%ptr)
                        call DeleteLinkVertex(DVT, linkDVT)
                    else
                        icdDT => icdDT%next
                        flag = .true.
                    end if
                else
                    icdDT => icdDT%next
                    flag = .true.
                end if
                if (flag .and. associated(icdDT, DVT%icdDTHead)) exit
            end do
#endif
            ! Record the maximum number of Voronoi vertices
            if (DVT%numIcdDT > maxNumVVT) maxNumVVT = DVT%numIcdDT
            ! Set the spatial bound
            if (VC%smp%bnd%numVtx /= DVT%numIcdDT) then
                deallocate(VC%smp%bnd%vtx)
                VC%smp%bnd%numVtx = DVT%numIcdDT
                allocate(VC%smp%bnd%vtx(VC%smp%bnd%numVtx))
            end if
            icdDT => DVT%icdDTHead
            do j = 1, DVT%numIcdDT
                VC%smp%bnd%vtx(j)%lon = icdDT%ptr%cnt%lon
                VC%smp%bnd%vtx(j)%lat = icdDT%ptr%cnt%lat
                icdDT => icdDT%next
            end do
            call calcArea(VC%smp%bnd%numVtx, &
                VC%smp%bnd%vtx(:)%lon, VC%smp%bnd%vtx(:)%lat, VC%smp%bnd%area)
            if (VC%smp%bnd%area < 0) then
                write(*, "('DVT', I5, ' has negative area!')") DVT%id
                call VC%smp%bnd%dump
            end if
#if (defined DelaunayAndVoronoi_Debug)
            totalArea = totalArea+VC%smp%bnd%area
#endif
            DVT => DVT%next
            VC => VC%next
        end do

#if (defined DelaunayAndVoronoi_Debug)
        print *, "PI is", totalArea/Re2/4.0, ", is that almostly right?"
#endif

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker
#endif

    end subroutine ExtractVoronoiDiagram

    ! ************************************************************************ !
    ! DelaunayAndVoronoi_CleanAll                                              !
    ! Purpose:                                                                 !
    !   Clean all the Delaunay triangles and reset the values for a new start. !
    ! ************************************************************************ !

    subroutine DelaunayAndVoronoi_CleanAll
        type(DelaunayTriangle), pointer :: DT
        type(DelaunayVertex), pointer :: DVT
        type(DelaunayTrianglePointerList), pointer :: icdDT1, icdDT2
        type(DelaunayVertexPointerList), pointer :: linkDVT
        integer i, j

        ! Delete triangles
        do while (associated(DTHead))
            DT => DTHead%next
            deallocate(DTHead)
            numDT = numDT-1
            DTHead => DT
        end do
        numTotalDT = 0

        ! Delete virtual vertices
        do while (associated(VirtualDVTHead))
            DVT => VirtualDVTHead%next
            deallocate(VirtualDVTHead)
            numVirtualDVT = numVirtualDVT-1
            VirtualDVTHead => DVT
        end do

        ! Delete incident triangle and link vertex lists
        DVT => DVTHead
        do i = 1, numDVT
            icdDT1 => DVT%icdDTHead%next
            do while (DVT%numIcdDT /= 1)
                icdDT2 => icdDT1%next
                deallocate(icdDT1)
                DVT%numIcdDT = DVT%numIcdDT-1
                icdDT1 => icdDT2
            end do
            nullify(DVT%icdDTHead%prev)
            do j = 1, DVT%numLinkDVT
                linkDVT => DVT%linkDVTHead%next
                deallocate(DVT%linkDVTHead)
                DVT%numLinkDVT = DVT%numLinkDVT-1
                DVT%linkDVTHead => linkDVT
            end do
            nullify(DVT%linkDVTHead)
            nullify(DVT%cntDT1)
            nullify(DVT%cntDT2)
            nullify(DVT%stub1)
            nullify(DVT%stub2)
            DVT%edgeIdx = -1
            DVT => DVT%next
        end do

    end subroutine DelaunayAndVoronoi_CleanAll
    
    ! ************************************************************************ !
    ! Output subroutines                                                       !
    ! ************************************************************************ !

    subroutine DelaunayAndVoronoi_Output(filePath)
        character(*), intent(in) :: filePath

        real(4), allocatable :: lon(:), lat(:)
        integer, allocatable :: triangle(:,:)
        type(DelaunayVertex), pointer :: DVT
        type(DelaunayTriangle), pointer :: DT
        real(4), allocatable :: lonVVT(:,:), latVVT(:,:)
        real(4), allocatable :: areaVC(:)
        integer, allocatable :: numVVT(:)
        type(VoronoiCell), pointer :: VC
        integer i, j, numPoint

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_RecordSpeaker("DelaunayAndVoronoi_Output")
#endif

        numPoint = numDVT+numVirtualDVT

        allocate(lon(numPoint))
        allocate(lat(numPoint))

        DVT => DVTHead
        do i = 1, numDVT
            lon(i) = real(DVT%smp%lon)
            lat(i) = real(DVT%smp%lat)
            DVT => DVT%next
        end do
        DVT => VirtualDVTHead
        do i = 1, numVirtualDVT
            lon(numDVT+i) = real(DVT%smp%lon)
            lat(numDVT+i) = real(DVT%smp%lat)
            DVT => DVT%next
        end do

        allocate(triangle(3,numDT))

        DT => DTHead
        do i = 1, numDT
            do j = 1, 3
                triangle(j,i) = DT%DVT(j)%ptr%id
            end do
            DT => DT%next
        end do

        allocate(lonVVT(maxNumVVT,numVC))
        allocate(latVVT(maxNumVVT,numVC))
        allocate(areaVC(numVC))
        allocate(numVVT(numVC))

        VC => VCHead
        do j = 1, numVC
            do i = 1, VC%smp%bnd%numVtx
                lonVVT(i,j) = real(VC%smp%bnd%vtx(i)%lon)
                latVVT(i,j) = real(VC%smp%bnd%vtx(i)%lat)
            end do
            areaVC(j) = VC%smp%bnd%area
            numVVT(j) = VC%smp%bnd%numVtx
            VC => VC%next
        end do

        call NFWrap_CreateIrregular(filePath, timeVariant=.false., card=fcard)
        call NFWrap_NewDim(fcard, "numPnt", numPoint)
        call NFWrap_NewDim(fcard, "numTri", numDT)
        call NFWrap_NewDim(fcard, "numTriVtx", 3)
        call NFWrap_NewDim(fcard, "maxNumVtx", maxNumVVT)
        call NFWrap_New1DVar(fcard, "lonPnt", "float", "numPnt", &
            "Point longitude", "radian", .false.)
        call NFWrap_New1DVar(fcard, "latPnt", "float", "numPnt", &
            "Point latitude", "radian", .false.)
        call NFWrap_New2DVar(fcard, "triIdx", "integer", &
            "numTriVtx", "numTri", "Triangle vertex index", "", .false.)
        call NFWrap_New1DVar(fcard, "numVtx", "integer", "numPnt", &
            "Real Voronoi cell vertex number", "", .false.)
        call NFWrap_New2DVar(fcard, "lonVtx", "float", &
            "maxNumVtx", "numPnt", "Voronoi cell vertex longitude", "radian", .false.)
        call NFWrap_New2DVar(fcard, "latVtx", "float", &
            "maxNumVtx", "numPnt", "Voronoi cell vertex latitude", "radian", .false.)
        call NFWrap_New1DVar(fcard, "area", "float", "numPnt", &
            "Voronoi cell area", "m2", .false.)
        call NFWrap_Output1DVar(fcard, "lonPnt", lon)
        call NFWrap_Output1DVar(fcard, "latPnt", lat)
        call NFWrap_Output2DVar(fcard, "triIdx", triangle)
        call NFWrap_Output1DVar(fcard, "numVtx", numVVT)
        call NFWrap_Output2DVar(fcard, "lonVtx", lonVVT)
        call NFWrap_Output2DVar(fcard, "latVtx", latVVT)
        call NFWrap_Output1DVar(fcard, "area", areaVC)
        call NFWrap_Close(fcard)

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_Speak(Notice, "File "//trim(filePath)//" is generated.")
        call MsgManager_DeleteSpeaker
#endif

    end subroutine DelaunayAndVoronoi_Output

    subroutine DelaunayAndVoronoi_OutputDelaunay(filePath)
        character(*), intent(in) :: filePath

        real(4), allocatable :: lon(:), lat(:)
        integer, allocatable :: triangle(:,:)
        type(DelaunayVertex), pointer :: DVT
        type(DelaunayTriangle), pointer :: DT
        integer i, j, numPoint

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_RecordSpeaker("DelaunayAndVoronoi_Output")
#endif

        numPoint = numDVT+numVirtualDVT

        allocate(lon(numPoint))
        allocate(lat(numPoint))

        DVT => DVTHead
        do i = 1, numDVT
            lon(i) = real(DVT%smp%lon)
            lat(i) = real(DVT%smp%lat)
            DVT => DVT%next
        end do
        DVT => VirtualDVTHead
        do i = 1, numVirtualDVT
            lon(numDVT+i) = real(DVT%smp%lon)
            lat(numDVT+i) = real(DVT%smp%lat)
            DVT => DVT%next
        end do

        allocate(triangle(3,numDT))

        DT => DTHead
        do i = 1, numDT
            do j = 1, 3
                triangle(j,i) = DT%DVT(j)%ptr%id
                ! This is problematic!!!
                if (triangle(j,i) < 0) then
                    triangle(j,i) = numDVT-triangle(j,i)
                end if
            end do
            DT => DT%next
        end do

        call NFWrap_CreateIrregular(filePath, .false., fcard)
        call NFWrap_NewDim(fcard, "numPnt", numPoint)
        call NFWrap_NewDim(fcard, "numTri", numDT)
        call NFWrap_NewDim(fcard, "numTriVtx", 3)
        call NFWrap_New1DVar(fcard, "lonPnt", "float", "numPnt", &
            "Point longitude", "radian", .false.)
        call NFWrap_New1DVar(fcard, "latPnt", "float", "numPnt", &
            "Point latitude", "radian", .false.)
        call NFWrap_New2DVar(fcard, "triIdx", "integer", "numTriVtx",&
            "numTri", "Triangle vertex index", "", .false.)
        call NFWrap_Output1DVar(fcard, "lonPnt", lon)
        call NFWrap_Output1DVar(fcard, "latPnt", lat)
        call NFWrap_Output2DVar(fcard, "triIdx", triangle)
        call NFWrap_Close(fcard)

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_Speak(Notice, "File "//trim(filePath)//" is generated.")
        call MsgManager_DeleteSpeaker
#endif

    end subroutine DelaunayAndVoronoi_OutputDelaunay

    ! ######################################################################## !
    ! ######################### Private procedures ########################### !
    ! ######################################################################## !

    subroutine GetTheFirstThreeDVTIndex(idx)
        integer, intent(out) :: idx(3)

        logical :: success = .false.

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_RecordSpeaker("GetTheFirstThreeDVTIndex")
#endif

        do while (.not. success)
            call RandomNumber_Get(1, numDVT, idx)
            if (idx(1) /= idx(2) .and. &
                idx(1) /= idx(3) .and. &
                idx(2) /= idx(3)) then
                success = .true.
            end if
        end do
        call CommonUtils_Sort(idx)
#if (defined DelaunayAndVoronoi_Debug)
        if (.not. doRandom) then
            idx = [idx1,idx2,idx3]
            print *, "The first vertices are fixed to "//trim(int2str(idx1))// &
                ", "//trim(int2str(idx2))//", "//trim(int2str(idx3))//" for debug!"
        end if
#endif
#if (!defined DelaunayAndVoronoi_FullSpeed || defined DelaunayAndVoronoi_Debug)
        call MsgManager_Speak(Notice, &
            "The randomly chosen vertices are: "// &
            trim(int2str(idx(1)))//", "//&
            trim(int2str(idx(2)))//", "//&
            trim(int2str(idx(3)))//".")
#endif

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_DeleteSpeaker
#endif

    end subroutine GetTheFirstThreeDVTIndex
    
    ! ************************************************************************ !
    ! InitDelaunayTriangle                                                     !
    ! Purpose:                                                                 !
    !   Initialize Delaunay triangulation by creating the first real triangle  !
    !   and the associated seven virtual triangles.                            !
    ! ************************************************************************ !

    subroutine InitDelaunayTriangle(DVT1_in, DVT2_in, DVT3_in)
        type(DelaunayVertex), intent(in), target :: DVT1_in, DVT2_in, DVT3_in

        type(DelaunayVertexPointer) DVT(6)  ! 3 real + "3" virtual vertices
        type(DelaunayTrianglePointer) DT(8) ! 1 real + "7" virtual triangles
        type(DelaunayVertex), pointer :: DVT1, DVT2
        integer map(24)
        integer i, j, ret

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_RecordSpeaker("InitDelaunayTriangle")
#endif

        ! Get the three real vertices and set up the three virtual
        ! vertices ...................................................... STEP 1
        ! Note: The virtual vertices are "antipodal" to the corresponding
        !       vertices of the first three inserted ones.
        ! ............................................................. STEP 1-1
        DVT(1)%ptr => DVT1_in
        DVT(2)%ptr => DVT2_in
        DVT(3)%ptr => DVT3_in
        ! ............................................................. STEP 1-2
        do i = 1, 3
            if (.not. associated(VirtualDVTHead)) then
                allocate(VirtualDVTHead)
                numVirtualDVT = 1
                DVT1 => VirtualDVTHead
            else
                allocate(DVT1%next)
                numVirtualDVT = numVirtualDVT+1
                DVT2 => DVT1
                DVT1 => DVT1%next
                DVT1%prev => DVT2
                DVT2%next => DVT1
            end if
            DVT1%id = -i
            allocate(DVT1%smp)
            DVT1%smp%id = -i
            DVT1%numIcdDT = 1
            allocate(DVT1%icdDTHead)
            DVT(i+3)%ptr => DVT1
            call InverseRotationTransform(DVT(i)%ptr%smp%lon, &
                DVT(i)%ptr%smp%lat, DVT1%smp%lon, DVT1%smp%lat, Zero, -PI05)
            if (DVT1%smp%lon < Zero) then
                DVT1%smp%lon = DVT1%smp%lon+PI2
            end if
            call CartesianTransformOnUnitSphere(DVT1%smp%lon, DVT1%smp%lat, &
                DVT1%smp%x, DVT1%smp%y, DVT1%smp%z)
        end do

        ! Create the first eight triangles .............................. STEP 2
        ! Note: One of them is real Delaunay triangle, and the rest are
        !       virtual to cover all the sphere.
        do i = 1, 8
            call NewDelaunayTriangle(DT(i)%ptr)
        end do

        ! Link the triangle with its adjacent triangles ................. STEP 3
        map = [5,2,4, 6,3,1, 7,4,2, 8,1,3, 1,8,6, 2,5,7, 3,6,8, 4,7,5]
        j = 1
        do i = 1, 8
            DT(i)%ptr%adjDT(1)%ptr => DT(map(j))%ptr; j = j+1
            DT(i)%ptr%adjDT(2)%ptr => DT(map(j))%ptr; j = j+1
            DT(i)%ptr%adjDT(3)%ptr => DT(map(j))%ptr; j = j+1
        end do

        ! Link the triangle with its vertices ........................... STEP 4
        ret = Orient(DVT(1)%ptr%smp%x, DVT(1)%ptr%smp%y, DVT(1)%ptr%smp%z, &
                     DVT(2)%ptr%smp%x, DVT(2)%ptr%smp%y, DVT(2)%ptr%smp%z, &
                     DVT(3)%ptr%smp%x, DVT(3)%ptr%smp%y, DVT(3)%ptr%smp%z)
        if (ret == OrientLeft) then ! Counter-clockwise
            map = [1,2,3, 1,3,5, 1,5,6, 1,6,2, 4,3,2, 4,5,3, 4,6,5, 4,2,6]
        else if (ret == OrientRight) then ! Clockwise
            map = [1,3,2, 1,2,6, 1,6,5, 1,5,3, 4,2,3, 4,6,2, 4,5,6, 4,3,5]
        else if (ret == OrientOn) then
            call MsgManager_Speak(Error, &
                "The first three points are collinear.", &
                "Change the order of points to ensure the non-collinearity "// &
                "of the first three ones. If you are sure they are not, "// &
                "this may be caused by the failure of floating-point"// &
                "calculation.")
            call RunManager_EndRun
        end if
        j = 1
        do i = 1, 8
            DT(i)%ptr%DVT(1)%ptr => DVT(map(j))%ptr; j = j+1
            DT(i)%ptr%DVT(2)%ptr => DVT(map(j))%ptr; j = j+1
            DT(i)%ptr%DVT(3)%ptr => DVT(map(j))%ptr; j = j+1
        end do

        ! Link the vertex with one of its incident DT ................... STEP 5
        ! Note: Any one is ok, the full list will be extracted at the end of
        !       the Delaunay triangulation.
        DVT(1)%ptr%icdDTHead%ptr => DT(1)%ptr
        DVT(2)%ptr%icdDTHead%ptr => DT(1)%ptr
        DVT(3)%ptr%icdDTHead%ptr => DT(1)%ptr
        DVT(4)%ptr%icdDTHead%ptr => DT(7)%ptr
        DVT(5)%ptr%icdDTHead%ptr => DT(7)%ptr
        DVT(6)%ptr%icdDTHead%ptr => DT(7)%ptr

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker
#endif

    end subroutine InitDelaunayTriangle
    
    ! ************************************************************************ !
    ! InsertVertex                                                             !
    ! Purpose:                                                                 !
    !   Insert one vertex into the Delaunay triangulation. It is the core of   !
    !   the incremental algorithm.                                             !
    ! ************************************************************************ !

    subroutine InsertVertex(DVT)
        type(DelaunayVertex), intent(inout), target :: DVT

        type(DelaunayVertex), pointer :: DVT1, DVT2, DVT3, DVT4
        type(DelaunayTriangle), pointer :: adjDT1, adjDT2, adjDT3, adjDT4
        type(DelaunayTriangle), pointer :: oldDT1, oldDT2, adjDT
        type(DelaunayTrianglePointer) newDT(4)
        integer i, j, idx

        if (associated(DVT%cntDT2)) then
            ! DVT is on the edge shared by two adjacent triangles.
            ! Note: In this case, there is axial dissymetry with respect to
            !       DVT topologically. Therefore, the indices of the vertices
            !       of oldDT1 does matter.
            ! Point for short-hand ...................................... STEP 0
            oldDT1 => DVT%cntDT1
            oldDT2 => DVT%cntDT2
            do idx = 1, 3
                if (associated(oldDT2%adjDT(idx)%ptr, oldDT1)) exit
            end do
            DVT1 => oldDT1%DVT(DVT%edgeIdx)%ptr
            DVT2 => oldDT1%DVT(ip1(DVT%edgeIdx))%ptr
            DVT3 => oldDT1%DVT(im1(DVT%edgeIdx))%ptr
            DVT4 => oldDT2%DVT(idx)%ptr
            adjDT1 => oldDT1%adjDT(ip1(DVT%edgeIdx))%ptr
            adjDT2 => oldDT1%adjDT(im1(DVT%edgeIdx))%ptr
            adjDT3 => oldDT2%adjDT(ip1(idx))%ptr
            adjDT4 => oldDT2%adjDT(im1(idx))%ptr
            ! Subdivide the two triangles into four new ones ............ STEP 1
            do i = 1, 4
                call NewDelaunayTriangle(newDT(i)%ptr)
            end do
            ! Set up the topology of the four triangles ................. STEP 2
            ! For oldDT1's newDTs
            ! ................................................... newDT(1)
            newDT(1)%ptr%DVT(1)%ptr => DVT3
            newDT(1)%ptr%DVT(2)%ptr => DVT1
            newDT(1)%ptr%DVT(3)%ptr => DVT
            newDT(1)%ptr%adjDT(1)%ptr => newDT(2)%ptr
            newDT(1)%ptr%adjDT(2)%ptr => newDT(4)%ptr
            newDT(1)%ptr%adjDT(3)%ptr => adjDT1
            ! ................................................... newDT(2)
            newDT(2)%ptr%DVT(1)%ptr => DVT1
            newDT(2)%ptr%DVT(2)%ptr => DVT2
            newDT(2)%ptr%DVT(3)%ptr => DVT
            newDT(2)%ptr%adjDT(1)%ptr => newDT(3)%ptr
            newDT(2)%ptr%adjDT(2)%ptr => newDT(1)%ptr
            newDT(2)%ptr%adjDT(3)%ptr => adjDT2
            ! For oldDT2's newDTs
            ! ................................................... newDT(3)
            newDT(3)%ptr%DVT(1)%ptr => DVT2
            newDT(3)%ptr%DVT(2)%ptr => DVT4
            newDT(3)%ptr%DVT(3)%ptr => DVT
            newDT(3)%ptr%adjDT(1)%ptr => newDT(4)%ptr
            newDT(3)%ptr%adjDT(2)%ptr => newDT(2)%ptr
            newDT(3)%ptr%adjDT(3)%ptr => adjDT3
            ! ................................................... newDT(4)
            newDT(4)%ptr%DVT(1)%ptr => DVT4
            newDT(4)%ptr%DVT(2)%ptr => DVT3
            newDT(4)%ptr%DVT(3)%ptr => DVT
            newDT(4)%ptr%adjDT(1)%ptr => newDT(1)%ptr
            newDT(4)%ptr%adjDT(2)%ptr => newDT(3)%ptr
            newDT(4)%ptr%adjDT(3)%ptr => adjDT4
            ! Link newly inserted DVT to its incident DTs ............... STEP 3
            DVT%icdDTHead%ptr => newDT(1)%ptr
            ! Make change to the old triangle's adjDT and DVT ........... STEP 4
            ! ...................................................... adjDT
            do i = 1, 3
                if (associated(adjDT1%adjDT(i)%ptr, oldDT1)) then
                    adjDT1%adjDT(i)%ptr => newDT(1)%ptr
                    exit
                end if
            end do
            do i = 1, 3
                if (associated(adjDT2%adjDT(i)%ptr, oldDT1)) then
                    adjDT2%adjDT(i)%ptr => newDT(2)%ptr
                    exit
                end if
            end do
            do i = 1, 3
                if (associated(adjDT3%adjDT(i)%ptr, oldDT2)) then
                    adjDT3%adjDT(i)%ptr => newDT(3)%ptr
                    exit
                end if
            end do
            do i = 1, 3
                if (associated(adjDT4%adjDT(i)%ptr, oldDT2)) then
                    adjDT4%adjDT(i)%ptr => newDT(4)%ptr
                    exit
                end if
            end do
            ! ......................................................... DVT
            if (associated(DVT1%icdDTHead%ptr, oldDT1)) &
                DVT1%icdDTHead%ptr => newDT(1)%ptr
            if (associated(DVT2%icdDTHead%ptr, oldDT1) .or. &
                associated(DVT2%icdDTHead%ptr, oldDT2)) &
                DVT2%icdDTHead%ptr => newDT(2)%ptr
            if (associated(DVT3%icdDTHead%ptr, oldDT1) .or. &
                associated(DVT3%icdDTHead%ptr, oldDT2)) &
                DVT3%icdDTHead%ptr => newDT(4)%ptr
            if (associated(DVT4%icdDTHead%ptr, oldDT2)) &
                DVT2%icdDTHead%ptr => newDT(3)%ptr
            ! Validate the three triangles .............................. STEP 5
            do i = 1, 4
                call ValidateNewTriangle(newDT(i)%ptr)
            end do
            ! Record the obsolete triangles ............................. STEP 6
            call RecordObsoleteTriangle(oldDT1)
            oldDT1%numSubDT = 2
            oldDT1%subDT(1)%ptr => newDT(1)%ptr
            oldDT1%subDT(2)%ptr => newDT(2)%ptr
            call RecordObsoleteTriangle(oldDT2)
            oldDT2%numSubDT = 2
            oldDT2%subDT(1)%ptr => newDT(3)%ptr
            oldDT2%subDT(2)%ptr => newDT(4)%ptr
            ! Remove DVT from point-in-triangle relation ................ STEP 7
            DVT%cntDT1 => null()
            DVT%cntDT2 => null()
            call DeleteIncludedVertex(DVT, oldDT1, oldDT2)
        else
            ! DVT is only contained by one triangle.
            ! Note: In this case, there is no axial dissymetry with respect to
            !       DVT topologically. Therefore, the indices of the vertices
            !       of oldDT1 does not matter.
            ! For short-hand ............................................ STEP 0
            oldDT1 => DVT%cntDT1
            DVT1 => oldDT1%DVT(1)%ptr
            DVT2 => oldDT1%DVT(2)%ptr
            DVT3 => oldDT1%DVT(3)%ptr
            adjDT1 => oldDT1%adjDT(3)%ptr
            adjDT2 => oldDT1%adjDT(1)%ptr
            adjDT3 => oldDT1%adjDT(2)%ptr
            ! Subdivide the old triangle into three new ones ............ STEP 1
            do i = 1, 3
                call NewDelaunayTriangle(newDT(i)%ptr)
            end do
            ! Set up the topology of the three triangles ................ STEP 2
            ! ................................................... newDT(1)
            newDT(1)%ptr%DVT(1)%ptr => DVT1
            newDT(1)%ptr%DVT(2)%ptr => DVT2
            newDT(1)%ptr%DVT(3)%ptr => DVT
            newDT(1)%ptr%adjDT(1)%ptr => newDT(2)%ptr
            newDT(1)%ptr%adjDT(2)%ptr => newDT(3)%ptr
            newDT(1)%ptr%adjDT(3)%ptr => adjDT1
            ! ................................................... newDT(1)
            newDT(2)%ptr%DVT(1)%ptr => DVT2
            newDT(2)%ptr%DVT(2)%ptr => DVT3
            newDT(2)%ptr%DVT(3)%ptr => DVT
            newDT(2)%ptr%adjDT(1)%ptr => newDT(3)%ptr
            newDT(2)%ptr%adjDT(2)%ptr => newDT(1)%ptr
            newDT(2)%ptr%adjDT(3)%ptr => adjDT2
            ! ................................................... newDT(1)
            newDT(3)%ptr%DVT(1)%ptr => DVT3
            newDT(3)%ptr%DVT(2)%ptr => DVT1
            newDT(3)%ptr%DVT(3)%ptr => DVT
            newDT(3)%ptr%adjDT(1)%ptr => newDT(1)%ptr
            newDT(3)%ptr%adjDT(2)%ptr => newDT(2)%ptr
            newDT(3)%ptr%adjDT(3)%ptr => adjDT3
            ! Link the newly inserted DVT to one of its incident DT ..... STEP 3
            DVT%icdDTHead%ptr => newDT(1)%ptr
            ! Make change to the old triangle's adjDT ....... ........... STEP 4
            ! Note: There is an awkful situation as that we know
            !       the adjacent triangle of oldDT1, but we do not
            !       know what is the index of oldDT1 in that triangle.
            !       So we do not know which adjDT of that triangle
            !       to point to the newDT directly, but have to check it.
            do i = 1, 3
                if (associated(adjDT1%adjDT(i)%ptr, oldDT1)) then
                    adjDT1%adjDT(i)%ptr => newDT(1)%ptr
                    exit
                end if
            end do
            do i = 1, 3
                if (associated(adjDT2%adjDT(i)%ptr, oldDT1)) then
                    adjDT2%adjDT(i)%ptr => newDT(2)%ptr
                    exit
                end if
            end do
            do i = 1, 3
                if (associated(adjDT3%adjDT(i)%ptr, oldDT1)) then
                    adjDT3%adjDT(i)%ptr => newDT(3)%ptr
                    exit
                end if
            end do
            if (associated(DVT1%icdDTHead%ptr, oldDT1)) &
                DVT1%icdDTHead%ptr => newDT(1)%ptr
            if (associated(DVT2%icdDTHead%ptr, oldDT1)) &
                DVT2%icdDTHead%ptr => newDT(2)%ptr
            if (associated(DVT3%icdDTHead%ptr, oldDT1)) &
                DVT3%icdDTHead%ptr => newDT(3)%ptr
            ! Validate the three triangles .............................. STEP 5
            do i = 1, 3
                call ValidateNewTriangle(newDT(i)%ptr)
            end do
            ! Record the obsolete triangle .............................. STEP 6
            call RecordObsoleteTriangle(oldDT1)
            oldDT1%numSubDT = 3
            oldDT1%subDT(1)%ptr => newDT(1)%ptr
            oldDT1%subDT(2)%ptr => newDT(2)%ptr
            oldDT1%subDT(3)%ptr => newDT(3)%ptr
            ! Remove DVT from point-in-triangle relation ................ STEP 7
            DVT%cntDT1 => null()
            call DeleteIncludedVertex(DVT, oldDT1)
        end if

    end subroutine InsertVertex

    ! ************************************************************************ !
    ! DeleteVertex                                                             !
    ! Purpose:                                                                 !
    !   The inverse of InsertVertex. It is used for deleting virtual vertices. !
    ! ************************************************************************ !

    subroutine DeleteVertex(DVT)
        type(DelaunayVertex), intent(inout), target :: DVT

        type(DelaunayTrianglePointerList), pointer :: icdDT
        type(DelaunayVertexPointerList), pointer :: linkDVT, restDVT
        type(DelaunayTriangle), pointer :: DT1, DT2, DT3, newDT1, newDT2
        type(DelaunayVertex), pointer :: DVT1, DVT2, DVT3
        logical empty
        integer i, j, k, ret

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_RecordSpeaker("DeleteVertex")
#endif

#if (defined DelaunayAndVoronoi_Debug)
        print *, "Deleting vertex", DVT%id, "..."
        !call DVT%smp%dump
        !linkDVT => DVT%linkDVTHead
        !do i = 1, DVT%numLinkDVT
        !    call linkDVT%ptr%smp%dump
        !    linkDVT => linkDVT%next
        !end do
        !read(*, *)
#endif

        icdDT => DVT%icdDTHead
        linkDVT => DVT%linkDVTHead
        do while (DVT%numIcdDT > 3)
            DT1 => icdDT%ptr
            DT2 => icdDT%next%ptr
            DVT1 => linkDVT%ptr
            DVT2 => linkDVT%next%ptr
            DVT3 => linkDVT%next%next%ptr
            ! Check 1: Is potential DT convex?
            ret = Orient(DVT1%smp%x, DVT1%smp%y, DVT1%smp%z, &
                         DVT2%smp%x, DVT2%smp%y, DVT2%smp%z, &
                         DVT3%smp%x, DVT3%smp%y, DVT3%smp%z)
            if (ret == OrientRight) then
                icdDT => icdDT%next
                linkDVT => linkDVT%next
                cycle ! NOT PASS
            end if
            ! Check 2: Does potential DT contain DVT?
            ret = Orient(DVT3%smp%x, DVT3%smp%y, DVT3%smp%z, &
                         DVT1%smp%x, DVT1%smp%y, DVT1%smp%z, &
                         DVT%smp%x,  DVT%smp%y,  DVT%smp%z)
            if (ret == OrientLeft) then
                icdDT => icdDT%next
                linkDVT => linkDVT%next
                cycle ! NOT PASS
            end if
            ! Check 3: Does potential DT satisfy empty-circumcircle rule?
            empty = .true.
            restDVT => linkDVT%next%next%next
            do i = 1, DVT%numIcdDT-3
                ret = InCircle(DVT1%smp%x, DVT1%smp%y, DVT1%smp%z, &
                               DVT2%smp%x, DVT2%smp%y, DVT2%smp%z, &
                               DVT3%smp%x, DVT3%smp%y, DVT3%smp%z, &
                               restDVT%ptr%smp%x, restDVT%ptr%smp%y, &
                               restDVT%ptr%smp%z)
                if (ret == InsideCircle) then
                    empty = .false.
                    exit
                else if (ret == OnTheCircle) then
#if (!defined DelaunayAndVoronoi_FullSpeed)
                    call MsgManager_Speak(Warning, "Encounter cocircular vertices.")
#endif
                end if
                restDVT => restDVT%next
            end do
            if (.not. empty) then
                icdDT => icdDT%next
                linkDVT => linkDVT%next
                cycle ! NOT PASS
            end if
            ! So far, the potential DT is a real DT.
            do i = 1, 3
                if (associated(DT1%DVT(i)%ptr, DVT)) exit
            end do
            do j = 1, 3
                if (associated(DT2%DVT(j)%ptr, DVT)) exit
            end do
            call Flip22(DT1, DT2, newDT1, newDT2, [im1(i),i,ip1(i),im1(j)])
            call DeleteDelaunayTriangle(DT1)
            call DeleteDelaunayTriangle(DT2)
        end do

        DT1 => DVT%icdDTHead%ptr
        DT2 => DVT%icdDTHead%next%ptr
        DT3 => DVT%icdDTHead%next%next%ptr
        DVT1 => DVT%linkDVTHead%ptr
        DVT2 => DVT%linkDVTHead%next%ptr
        DVT3 => DVT%linkDVTHead%next%next%ptr
        do i = 1, 3
            if (associated(DT1%DVT(i)%ptr, DVT)) exit
        end do
        do j = 1, 3
            if (associated(DT2%DVT(j)%ptr, DVT)) exit
        end do
        do k = 1, 3
            if (associated(DT3%DVT(k)%ptr, DVT)) exit
        end do
        call Flip31(DT1, DT2, DT3, newDT1, [ip1(i),im1(i),i,j,k])
        call DeleteDelaunayTriangle(DT1)
        call DeleteDelaunayTriangle(DT2)
        call DeleteDelaunayTriangle(DT3)

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_Speak(Notice, &
            "DVT "//trim(int2str(DVT%id))//" has been deleted from the list.")
        call MsgManager_DeleteSpeaker
#endif

    end subroutine DeleteVertex

    ! ************************************************************************ !
    ! Flip??                                                                   !
    ! Purpose:                                                                 !
    !   These flips are the basic operations in the action of insertion and    !
    !   deletion of vertices.                                                  !
    ! Note:                                                                    !
    !   Flip24 and Flip 13 have not been extracted from InsertVertex yet.      !
    ! ************************************************************************ !

    subroutine Flip22(oldDT1, oldDT2, newDT1, newDT2, idxMap)
        type(DelaunayTriangle), intent(in), target :: oldDT1, oldDT2
        type(DelaunayTriangle), intent(out), pointer :: newDT1, newDT2
        integer, intent(in) :: idxMap(4)

        type(DelaunayVertex), pointer :: DVT1, DVT2, DVT3, DVT4
        type(DelaunayTriangle), pointer :: adjDT1, adjDT2, adjDT3, adjDT4
        integer i

        ! Point for short-hand .......................................... STEP 0
        DVT1 => oldDT1%DVT(idxMap(1))%ptr
        DVT2 => oldDT1%DVT(idxMap(2))%ptr
        DVT3 => oldDT1%DVT(idxMap(3))%ptr
        DVT4 => oldDT2%DVT(idxMap(4))%ptr
        adjDT1 => oldDT1%adjDT(idxMap(1))%ptr
        adjDT2 => oldDT1%adjDT(idxMap(2))%ptr
        adjDT3 => oldDT2%adjDT(ip1(idxMap(4)))%ptr
        adjDT4 => oldDT2%adjDT(im1(idxMap(4)))%ptr
        ! Create two new DTs and set up their topology .................. STEP 1
        call NewDelaunayTriangle(newDT1)
        call NewDelaunayTriangle(newDT2)
        newDT1%DVT(1)%ptr => DVT1
        newDT1%DVT(2)%ptr => DVT4
        newDT1%DVT(3)%ptr => DVT3
        newDT1%adjDT(1)%ptr => newDT2
        newDT1%adjDT(2)%ptr => adjDT2
        newDT1%adjDT(3)%ptr => adjDT3 ! This is correct
        newDT2%DVT(1)%ptr => DVT4
        newDT2%DVT(2)%ptr => DVT2
        newDT2%DVT(3)%ptr => DVT3
        newDT2%adjDT(1)%ptr => adjDT1
        newDT2%adjDT(2)%ptr => newDT1
        newDT2%adjDT(3)%ptr => adjDT4
        ! Make change to the old triangle's adjDTs and DVTs ............. STEP 2
        do i = 1, 3
            if (associated(adjDT1%adjDT(i)%ptr, oldDT1)) then
                adjDT1%adjDT(i)%ptr => newDT2
                exit
            end if
        end do
        do i = 1, 3
            if (associated(adjDT2%adjDT(i)%ptr, oldDT1)) then
                adjDT2%adjDT(i)%ptr => newDT1
                exit
            end if
        end do
        do i = 1, 3
            if (associated(adjDT3%adjDT(i)%ptr, oldDT2)) then
                adjDT3%adjDT(i)%ptr => newDT1
                exit
            end if
        end do
        do i = 1, 3
            if (associated(adjDT4%adjDT(i)%ptr, oldDT2)) then
                adjDT4%adjDT(i)%ptr => newDT2
                exit
            end if
        end do
        if (associated(DVT1%icdDTHead%prev)) then
            ! The full list of incident DTs has been extracted.
            ! This may be called from DeleteVertex.
            call MergeIncidentTriangle(DVT1, oldDT1, oldDT2, newDT1)
            call MergeIncidentTriangle(DVT2, oldDT1, oldDT2, newDT2)
            call SplitIncidentTriangle(DVT3, oldDT1, newDT1, newDT2)
            call SplitIncidentTriangle(DVT4, oldDT2, newDT2, newDT1)
            call DeleteLinkVertex(DVT1, DVT2)
            call DeleteLinkVertex(DVT2, DVT1)
            call AddLinkVertex(DVT3, DVT1, DVT2, DVT4)
            call AddLinkVertex(DVT4, DVT2, DVT1, DVT3)
        else
            ! The full list of incident DTs has not been extracted.
            ! This may be called from ValidateNewTriangle.
            DVT1%icdDTHead%ptr => newDT1
            DVT2%icdDTHead%ptr => newDT2
            DVT3%icdDTHead%ptr => newDT2
            DVT4%icdDTHead%ptr => newDT1
        end if

    end subroutine Flip22

    subroutine Flip31(oldDT1, oldDT2, oldDT3, newDT, idxMap)
        type(DelaunayTriangle), intent(in), target :: oldDT1, oldDT2, oldDT3
        type(DelaunayTriangle), intent(out), pointer :: newDT
        integer, intent(in) :: idxMap(5)

        type(DelaunayVertex), pointer :: DVT1, DVT2, DVT3, DVT4
        type(DelaunayTriangle), pointer :: adjDT1, adjDT2, adjDT3
        integer i

        ! Point for short-hand .......................................... STEP 0
        DVT1 => oldDT1%DVT(idxMap(1))%ptr
        DVT2 => oldDT1%DVT(idxMap(2))%ptr
        DVT3 => oldDT1%DVT(idxMap(3))%ptr
        DVT4 => oldDT2%DVT(im1(idxMap(4)))%ptr
        adjDT1 => oldDT1%adjDT(idxMap(3))%ptr
        adjDT2 => oldDT2%adjDT(idxMap(4))%ptr
        adjDT3 => oldDT3%adjDT(idxMap(5))%ptr
        ! Create two new DTs and set up their topology .................. STEP 1
        call NewDelaunayTriangle(newDT)
        newDT%DVT(1)%ptr => DVT1
        newDT%DVT(2)%ptr => DVT2
        newDT%DVT(3)%ptr => DVT4
        newDT%adjDT(1)%ptr => adjDT2
        newDT%adjDT(2)%ptr => adjDT3
        newDT%adjDT(3)%ptr => adjDT1
        ! Make change to the old triangle's adjDTs and DVTs ............. STEP 2
        do i = 1, 3
            if (associated(adjDT1%adjDT(i)%ptr, oldDT1)) then
                adjDT1%adjDT(i)%ptr => newDT
            end if
        end do
        do i = 1, 3
            if (associated(adjDT2%adjDT(i)%ptr, oldDT2)) then
                adjDT2%adjDT(i)%ptr => newDT
            end if
        end do
        do i = 1, 3
            if (associated(adjDT3%adjDT(i)%ptr, oldDT1)) then
                adjDT3%adjDT(i)%ptr => newDT
            end if
        end do
        ! This should be called after the full list of incident DTs
        ! has been extracted. (Called in DeleteVertex)
        call MergeIncidentTriangle(DVT1, oldDT1, oldDT3, newDT)
        call MergeIncidentTriangle(DVT2, oldDT2, oldDT1, newDT)
        call MergeIncidentTriangle(DVT4, oldDT3, oldDT2, newDT)
        call DeleteLinkVertex(DVT1, DVT3)
        call DeleteLinkVertex(DVT2, DVT3)
        call DeleteLinkVertex(DVT4, DVT3)

    end subroutine Flip31

    ! ************************************************************************ !
    ! ValidateNewTriangle                                                      !
    ! Purpose:                                                                 !
    !   Validate the new triangle with its third vertex as the newly inserted  !
    !   vertex.                                                                !
    ! ************************************************************************ !

    recursive subroutine ValidateNewTriangle(newDT)
        type(DelaunayTriangle), intent(inout), target :: newDT
    
        type(DelaunayTriangle), pointer :: opsDT, newerDT1, newerDT2
        type(DelaunayVertex), pointer :: opsDVT
        integer i, idx, ret

        opsDT => newDT%adjDT(3)%ptr
        ! Get the opposite vertex to the newly inserted vertex with respect to
        ! the shared edge by newDT and opsDT.
        do idx = 1, 3
            if (associated(opsDT%adjDT(idx)%ptr, newDT)) then
                opsDVT => opsDT%DVT(idx)%ptr
                exit
            end if
        end do
        ! Check if the opposite vertex is in the circumcircle of DT.
        ret = InCircle(newDT, opsDVT)
        if (ret == OutsideCircle) then
            return ! THE ONLY EXIT
        else if (ret == InsideCircle) then
            call Flip22(newDT, opsDT, newerDT1, newerDT2, [1,2,3,idx])
            ! * Recursively validate the two new DT.
            call ValidateNewTriangle(newerDT1)
            call ValidateNewTriangle(newerDT2)
            ! * Record the temporal triangle
            call RecordTemporalTriangle(newDT)
            newDT%numSubDT = 2
            newDT%subDT(1)%ptr => newerDT1
            newDT%subDT(2)%ptr => newerDT2
            ! * Record the obsolete triangle
            call RecordObsoleteTriangle(opsDT)
            opsDT%numSubDT = 2
            opsDT%subDT(1)%ptr => newerDT1
            opsDT%subDT(2)%ptr => newerDT2
        else if (ret == OnTheCircle) then
            ! How to deal with this cocircle degenerate case?
#if (!defined DelaunayAndVoronoi_FullSpeed)
            call MsgManager_Speak(Warning, "Encounter cocircular vertices.")
#endif
        end if

    end subroutine ValidateNewTriangle

    ! ************************************************************************ !
    ! UpdatePointInTriangle                                                    !
    ! Purpose:                                                                 !
    !   After the insertion of one vertex, the Delaunay triangulation will     !
    !   change, so is the point-in-triangle relation. Update the relation by   !
    !   testing in which subdivided triangle the included vertices of the      !
    !   obsolete triangles are. These tests is limited locally by first        !
    !   recording the included vertices in each Delaunay triangle, and then    !
    !   recording the obsolete triangle pointers and its subdivided triangles. !
    ! ************************************************************************ !

    recursive subroutine UpdatePointInTriangle(DVT, DT, found)
        type(DelaunayVertex), intent(inout) :: DVT
        type(DelaunayTriangle), intent(inout), target :: DT
        logical, intent(inout) :: found
    
        integer i, ret, onPlane(2)

        if (DT%numSubDT /= 0) then
            ! DT has been subdivided, go to next level ................. ENTRY 1
            do i = 1, DT%numSubDT
                call UpdatePointInTriangle(DVT, DT%subDT(i)%ptr, found)
                if (found) return
            end do
        else
            ! DT is at the last level without being subdivided ......... ENTRY 2
            onPlane = 0
            ! Check edge 3->1
            ret = Orient(DT%DVT(3)%ptr%smp%x, DT%DVT(3)%ptr%smp%y, &
                         DT%DVT(3)%ptr%smp%z, DT%DVT(1)%ptr%smp%x, &
                         DT%DVT(1)%ptr%smp%y, DT%DVT(1)%ptr%smp%z, &
                         DVT%smp%x, DVT%smp%y, DVT%smp%z)
            if (ret == OrientRight) then
                found = .false.
                return
            else if (ret == OrientOn) then
                onPlane(1) = 2
            end if
            ! Check edge 2->3
            ret = Orient(DT%DVT(2)%ptr%smp%x, DT%DVT(2)%ptr%smp%y, &
                         DT%DVT(2)%ptr%smp%z, DT%DVT(3)%ptr%smp%x, &
                         DT%DVT(3)%ptr%smp%y, DT%DVT(3)%ptr%smp%z, &
                         DVT%smp%x, DVT%smp%y, DVT%smp%z)
            if (ret == OrientRight) then
                found = .false.
                return
            else if (ret == OrientOn) then
                onPlane(2) = 1
            end if
            ! Summary
            found = .true.
            ! Note: Delete DVT from its contained DTs to avoid duplicate update.
            if (associated(DVT%cntDT2)) then
                call DeleteIncludedVertex(DVT, DVT%cntDT1, DVT%cntDT2)
            else
                call DeleteIncludedVertex(DVT, DVT%cntDT1)
            end if
            if (onPlane(1) == 0 .and. onPlane(2) == 0) then
                DVT%cntDT1 => DT
                DVT%cntDT2 => null()
                call RecordIncludedVertex(DVT, DT)
            else if (onPlane(1) == 2 .and. onPlane(2) == 0) then
                DVT%cntDT1 => DT
                DVT%cntDT2 => DT%adjDT(2)%ptr
                DVT%edgeIdx = onPlane(1)
                call RecordIncludedVertex(DVT, DT, DT%adjDT(2)%ptr)
            else if (onPlane(1) == 0 .and. onPlane(2) == 1) then
                DVT%cntDT1 => DT
                DVT%cntDT2 => DT%adjDT(1)%ptr
                DVT%edgeIdx = onPlane(2)
                call RecordIncludedVertex(DVT, DT, DT%adjDT(1)%ptr)
            else
                call MsgManager_Speak(Error, "DVT "//trim(int2str(DVT%id))// &
                    " coincides with vertex 3 of DT "//trim(int2str(DT%id)))
                call PrintTriangleTopology(DT)
                print *, DVT%smp%x, DVT%smp%y, DVT%smp%z
                print *, DT%DVT(3)%ptr%smp%x, DT%DVT(3)%ptr%smp%y, DT%DVT(3)%ptr%smp%z
                call RunManager_EndRun
            end if
        end if

    end subroutine UpdatePointInTriangle

    ! ************************************************************************ !
    ! ExtractIncidentDTAndLinkDVT                                              !
    ! Purpose:                                                                 !
    !   During the Delaunay triangulation, only one incident triangle of each  !
    !   vertex is recorded. This subroutine extract all the incident triangles !
    !   and link vertices to form a ring in counter-clockwise order.           !
    ! ************************************************************************ !

    subroutine ExtractIncidentDTAndLinkDVT(DVT)
        type(DelaunayVertex), intent(inout), target :: DVT

        type(DelaunayTrianglePointerList), pointer :: icdDT1, icdDT2
        type(DelaunayVertexPointerList), pointer :: linkDVT1, linkDVT2

        integer i

        icdDT1 => DVT%icdDTHead
        IncidentDTLoop: do
            ! Get the link DVT
            do i = 1, 3
                if (associated(icdDT1%ptr%DVT(i)%ptr, DVT)) then
                    if (.not. associated(DVT%linkDVTHead)) then
                        allocate(DVT%linkDVTHead)
                        linkDVT1 => DVT%linkDVTHead
                    else
                        linkDVT2 => linkDVT1
                        allocate(linkDVT1%next)
                        linkDVT1 => linkDVT1%next
                        linkDVT1%prev => linkDVT2
                        linkDVT2%next => linkDVT1
                    end if
                    DVT%numLinkDVT = DVT%numLinkDVT+1
                    linkDVT1%ptr => icdDT1%ptr%DVT(ip1(i))%ptr
                    if (associated(linkDVT1%ptr, DVT%linkDVTHead%ptr) .and. &
                        DVT%numLinkDVT > 1) then
                        ! The ring has formed.
                        linkDVT2%next => DVT%linkDVTHead
                        DVT%linkDVTHead%prev => linkDVT2
                        icdDT2%next => DVT%icdDTHead
                        DVT%icdDTHead%prev => icdDT2
                        DVT%numIcdDT = DVT%numIcdDT-1
                        DVT%numLinkDVT = DVT%numLinkDVT-1
                        exit IncidentDTLoop
                    end if
                    exit
                end if
            end do
            ! Shift to next incident DT
            icdDT2 => icdDT1
            allocate(icdDT1%next)
            icdDT1 => icdDT1%next
            icdDT1%prev => icdDT2
            icdDT2%next => icdDT1
            DVT%numIcdDT = DVT%numIcdDT+1
            icdDT1%ptr => icdDT2%ptr%adjDT(ip1(i))%ptr
        end do IncidentDTLoop

    end subroutine ExtractIncidentDTAndLinkDVT

    ! ************************************************************************ !
    ! Orient                                                                   !
    ! Purpose:                                                                 !
    !   This is the one of two basic geometric tests to give the relation      !
    !   between a point (x0,y0,z0) and a plane containing (x1,y1,z1),          !
    !   (x2,y2,z2), and (0,0,0). There are three types of relation, i.e. left, !
    !   right, and on, where left is defined relative to an observer at        !
    !   (x1,y1,z1) facing (x2,y2,z2).                                          !
    ! Return value:                                                            !
    !   OrientLeft                                                             !
    !   OrientRight                                                            !
    !   OrientOn                                                               !
    ! ************************************************************************ !

    integer function Orient(x1, y1, z1, x2, y2, z2, x0, y0, z0)
        real(RealKind), intent(in) :: x1, y1, z1, x2, y2, z2, x0, y0, z0

        real(RealKind) det

        det = x0*(y1*z2-y2*z1)-y0*(x1*z2-x2*z1)+z0*(x1*y2-x2*y1)

#if (defined DelaunayAndVoronoi_Debug)
        if (printDet) then
            print *, "x1 =", x1, "y1 =", y1, "z1 =", z1
            print *, "x2 =", x2, "y2 =", y2, "z2 =", z2
            print *, "x0 =", x0, "y0 =", y0, "z0 =", z0
            print *, "In OrientRelaxed, det =", det
            read(*, *)
        end if
#endif

        if (det > eps) then
            Orient = OrientLeft
        else if (-det > eps) then
            Orient = OrientRight
        else
            Orient = OrientOn
        end if

    end function Orient

    ! ************************************************************************ !
    ! InTriangle(Relaxed)                                                      !
    ! Purpose:                                                                 !
    !   This is a derived geometric test from "Orient". It gives the relation  !
    !   between a point (DVT) and a triangle (DT). There are four types of     !
    !   relation, i.e. inside, outside, on one edge and on one vertex.         !
    !   For relaxed InTriangle, the overlap with virtual DVTs situation will   !
    !   be checked.                                                            !
    ! Return value:                                                            !
    !   InsideTriangle                                                         !
    !   OutsideTriangle                                                        !
    !    1 ~  3 for on one edge, where the value indicates the edge index      !
    !   -1 ~ -3 for on one vertex, where the abs(value) indicates the vertex   !
    !           index                                                          !
    ! ************************************************************************ !

    integer function InTriangle(DT, DVT)
        type(DelaunayTriangle), intent(in) :: DT
        type(DelaunayVertex), intent(in) :: DVT

        real(RealKind) x(3), y(3), z(3), x0, y0, z0
        integer i, k, ret, onPlane(2)

        ! Copy for short-hand
        do i = 1, 3
            x(i) = DT%DVT(i)%ptr%smp%x
            y(i) = DT%DVT(i)%ptr%smp%y
            z(i) = DT%DVT(i)%ptr%smp%z
        end do
        x0 = DVT%smp%x
        y0 = DVT%smp%y
        z0 = DVT%smp%z

        k = 0
        do i = 1, 3
            ret = Orient(x(i),      y(i),      z(i),    &
                         x(ip1(i)), y(ip1(i)), z(ip1(i)), &
                         x0,        y0,        z0)
            if (ret == OrientRight) then
                InTriangle = OutsideTriangle
                return
            else if (ret == OrientOn) then
                k = k+1; onPlane(k) = im1(i)
            end if
        end do

        if (k == 0) then
            InTriangle = InsideTriangle
        else if (k == 1) then ! on the edge
            InTriangle = onPlane(k)
        else if (k == 2) then ! on the vertex
            if (onPlane(1) == 1 .and. onPlane(2) == 2) then
                InTriangle = -3
            else if (onPlane(1) == 2 .and. onPlane(2) == 1) then
                InTriangle = -3
            else if (onPlane(1) == 2 .and. onPlane(2) == 3) then
                InTriangle = -1
            else if (onPlane(1) == 3 .and. onPlane(2) == 2) then
                InTriangle = -1
            else if (onPlane(1) == 1 .and. onPlane(2) == 3) then
                InTriangle = -2
            else if (onPlane(1) == 3 .and. onPlane(2) == 1) then
                InTriangle = -2
            end if
        end if

    end function InTriangle

    integer function InTriangleRelaxed(DT, DVT)
        type(DelaunayTriangle), intent(in) :: DT
        type(DelaunayVertex), intent(in) :: DVT

        real(RealKind), parameter :: eps = 1.0e-4
        real(RealKind) x(3), y(3), z(3), x0, y0, z0
        integer i, j, k, ret, onPlane(2)

        ! Copy for short-hand
        do i = 1, 3
            if (DT%DVT(i)%ptr%id < 0) then
                ! For virtual DVT, first judge whether DVT overlaps it.
                if (abs(DT%DVT(i)%ptr%smp%lon-DVT%smp%lon) < eps .and. &
                    abs(DT%DVT(i)%ptr%smp%lat-DVT%smp%lat) < eps) then
                    InTriangleRelaxed = -i
                    return
                end if
            end if
            x(i) = DT%DVT(i)%ptr%smp%x
            y(i) = DT%DVT(i)%ptr%smp%y
            z(i) = DT%DVT(i)%ptr%smp%z
        end do
        x0 = DVT%smp%x
        y0 = DVT%smp%y
        z0 = DVT%smp%z

        k = 0
        do i = 1, 3
            j = ip1(i)
            ret = Orient(x(i), y(i), z(i), x(j), y(j), z(j), x0, y0, z0)
            if (ret == OrientRight) then
                InTriangleRelaxed = OutsideTriangle
                return
            else if (ret == OrientOn) then
                k = k+1; onPlane(k) = im1(i)
            end if
        end do

        if (k == 0) then
            InTriangleRelaxed = InsideTriangle
        else if (k == 1) then ! on the edge
            InTriangleRelaxed = onPlane(k)
        else if (k == 2) then ! on the vertex
            if (onPlane(1) == 1 .and. onPlane(2) == 2) then
                InTriangleRelaxed = -3
            else if (onPlane(1) == 2 .and. onPlane(2) == 1) then
                InTriangleRelaxed = -3
            else if (onPlane(1) == 2 .and. onPlane(2) == 3) then
                InTriangleRelaxed = -1
            else if (onPlane(1) == 3 .and. onPlane(2) == 2) then
                InTriangleRelaxed = -1
            else if (onPlane(1) == 1 .and. onPlane(2) == 3) then
                InTriangleRelaxed = -2
            else if (onPlane(1) == 3 .and. onPlane(2) == 1) then
                InTriangleRelaxed = -2
            end if
        end if

    end function InTriangleRelaxed

    ! ************************************************************************ !
    ! InCircle                                                                 !
    ! Purpose:                                                                 !
    !   This is the one of two basic geometric tests to check whether or not a !
    !   point (DVT) is inside the circumcircle (here is spherical circle) of a !
    !   Delaunay triangle (DT).                                                !
    !   In Cartesian coordinate system, this test turns into checking whether  !
    !   or not DVT is above the plane defined by the three vertices of DT.     !
    ! ************************************************************************ !

    integer function InCircle1(DT, DVT)
        type(DelaunayTriangle), intent(in) :: DT
        type(DelaunayVertex), intent(in) :: DVT

        real(RealKind) dx(3), dy(3), dz(3), x0, y0, z0, det
        integer i

        x0 = DVT%smp%x
        y0 = DVT%smp%y
        z0 = DVT%smp%z
        do i = 1, 3
            dx(i) = DT%DVT(i)%ptr%smp%x-x0
            dy(i) = DT%DVT(i)%ptr%smp%y-y0
            dz(i) = DT%DVT(i)%ptr%smp%z-z0
        end do

        det = dx(3)*(dy(2)*dz(1)-dy(1)*dz(2)) &
             -dy(3)*(dx(2)*dz(1)-dx(1)*dz(2)) &
             +dz(3)*(dx(2)*dy(1)-dx(1)*dy(2))

        if (det > eps) then
            InCircle1 = InsideCircle
        else if (-det > eps) then
            InCircle1 = OutsideCircle
        else
            InCircle1 = OnTheCircle
        end if

    end function InCircle1

    integer function InCircle2(x1, y1, z1, x2, y2, z2, x3, y3, z3, x0, y0, z0)
        real(RealKind), intent(in) :: x1, y1, z1, x2, y2, z2, x3, y3, z3
        real(RealKind), intent(in) :: x0, y0, z0

        real(RealKind) dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, det

        det = dx3*(dy2*dz1-dy1*dz2) &
             -dy3*(dx2*dz1-dx1*dz2) &
             +dz3*(dx2*dy1-dx1*dy2)

        if (det > eps) then
            InCircle2 = InsideCircle
        else if (-det > eps) then
            InCircle2 = OutsideCircle
        else
            InCircle2 = OnTheCircle
        end if
 

    end function InCircle2

    ! ************************************************************************ !
    ! CalcCircumcircle                                                         !
    ! Purpose:                                                                 !
    !   Calculate the circumcenter and circumradius of DT.                     !
    ! ************************************************************************ !

    subroutine CalcCircumcircle
        type(DelaunayTriangle), pointer :: DT
        real(RealKind) E2(3), E3(3), N(3), L, tmp
        integer i

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_RecordSpeaker("CalcCircumcircle")
#endif

        DT => DTHead
        do i = 1, numDT
            E3 = [DT%DVT(2)%ptr%smp%x-DT%DVT(1)%ptr%smp%x, &
                  DT%DVT(2)%ptr%smp%y-DT%DVT(1)%ptr%smp%y, &
                  DT%DVT(2)%ptr%smp%z-DT%DVT(1)%ptr%smp%z]
            E2 = [DT%DVT(3)%ptr%smp%x-DT%DVT(1)%ptr%smp%x, &
                  DT%DVT(3)%ptr%smp%y-DT%DVT(1)%ptr%smp%y, &
                  DT%DVT(3)%ptr%smp%z-DT%DVT(1)%ptr%smp%z]
            N = [E3(2)*E2(3)-E3(3)*E2(2), &
                 E3(3)*E2(1)-E3(1)*E2(3), &
                 E3(1)*E2(2)-E3(2)*E2(1)]
            L = N(1)**2+N(2)**2+N(3)**2
            if (L < eps) then
                call MsgManager_Speak(Error, "Vertices "// &
                    trim(int2str(DT%DVT(1)%ptr%id))//", "// &
                    trim(int2str(DT%DVT(2)%ptr%id))//", "// &
                    trim(int2str(DT%DVT(3)%ptr%id))//" are collinear.")
                call RunManager_EndRun
            end if
            L = sqrt(L)
            DT%cnt%x = N(1)/L
            DT%cnt%y = N(2)/L
            DT%cnt%z = N(3)/L
            call InverseCartesianTransformOnUnitSphere( &
                DT%cnt%lon, DT%cnt%lat, DT%cnt%x, DT%cnt%y, DT%cnt%z)
            tmp = (DT%DVT(1)%ptr%smp%x*DT%cnt%x+ &
                   DT%DVT(1)%ptr%smp%y*DT%cnt%y+ &
                   DT%DVT(1)%ptr%smp%z*DT%cnt%z)
            DT%radius = acos(tmp)
            DT => DT%next
        end do

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker
#endif

    end subroutine CalcCircumcircle

    ! ************************************************************************ !
    ! NewDelaunayTriangle                                                      !
    ! Purpose:                                                                 !
    !   Create the storage of a new Delaunay triangle in the double linked     !
    !   list.                                                                  !
    ! ************************************************************************ !

    subroutine NewDelaunayTriangle(DT)
        type(DelaunayTriangle), intent(out), pointer :: DT

        type(DelaunayTriangle), pointer, save :: DTCurr => null()

        numDT = numDT+1
        numTotalDT = numTotalDT+1
        if (numDT == 1) then
            allocate(DTHead)
            DTCurr => DTHead
        else
            DT => DTCurr
            allocate(DTCurr%next)
            DTCurr => DTCurr%next
            DTCurr%prev => DT
            DT%next => DTCurr
        end if
        DTCurr%id = numTotalDT
        DT => DTCurr

#if (!defined DelaunayAndVoronoi_FullSpeed)
#if (defined DelaunayAndVoronoi_Debug && defined DelaunayAndVoronoi_Verbose)
        call MsgManager_Speak(Notice, "DT "//trim(int2str(DT%id))//" has been newed.")
#endif
#endif

    end subroutine NewDelaunayTriangle
 
    ! ************************************************************************ !
    ! DeleteDelaunayTriangle                                                   !
    ! Purpose:                                                                 !
    !   Delete the storage of a Delaunay triangle from the double linked list. !
    ! ************************************************************************ !

    subroutine DeleteDelaunayTriangle(DT)
        type(DelaunayTriangle), intent(inout), pointer :: DT

#if (!defined DelaunayAndVoronoi_FullSpeed)
#if (defined DelaunayAndVoronoi_Debug && defined DelaunayAndVoronoi_Verbose)
        call MsgManager_Speak(Notice, "DT "//trim(int2str(DT%id))//" has been deleted.")
#endif
#endif

        numDT = numDT-1
        if (associated(DT%prev)) then
            DT%prev%next => DT%next
        else
            ! Delete the first triangle on the double linked list
            DTHead => DTHead%next
        end if
        if (associated(DT%next)) then
            DT%next%prev => DT%prev
        end if
        deallocate(DT)
    
    end subroutine DeleteDelaunayTriangle

    subroutine SplitIncidentTriangle(DVT, oldDT, newDT1, newDT2)
        type(DelaunayVertex), intent(inout), target :: DVT
        type(DelaunayTriangle), intent(in), target :: oldDT, newDT1, newDT2

        type(DelaunayTrianglePointerList), pointer :: icdDT1, icdDT2
        integer i

#if (!defined DelaunayAndVoronoi_FullSpeed)
#if (defined DelaunayAndVoronoi_Debug && DVT_ID)
        if (DVT%id == DVT_ID) then
            print *, "DVT", DVT%id, "is in SplitIncidentTriangle."
            print *, "  oldDT", oldDT%id, "--> newDT1", newDT1%id, "newDT2", newDT2%id
            print *, "Before split:"
            call PrintVertexTopology(DVT)
        end if
#endif
#endif

        icdDT1 => DVT%icdDTHead
        do i = 1, DVT%numIcdDT
            if (associated(icdDT1%ptr, oldDT)) then
                icdDT1%ptr => newDT1
                allocate(icdDT2)
                icdDT2%ptr => newDT2
                icdDT2%prev => icdDT1
                icdDT2%next => icdDT1%next
                if (associated(icdDT1%next)) then
                    icdDT1%next%prev => icdDT2
                end if
                icdDT1%next => icdDT2
                DVT%numIcdDT = DVT%numIcdDT+1
                exit
            end if
            icdDT1 => icdDT1%next
        end do

#if (!defined DelaunayAndVoronoi_FullSpeed)
#if (defined DelaunayAndVoronoi_Debug && DVT_ID)
        if (DVT%id == DVT_ID) then
            print *, "After split:"
            call PrintVertexTopology(DVT)
        end if
#endif
#endif

    end subroutine SplitIncidentTriangle

    subroutine MergeIncidentTriangle(DVT, oldDT1, oldDT2, newDT)
        type(DelaunayVertex), intent(inout), target :: DVT
        type(DelaunayTriangle), intent(in), target :: oldDT1, oldDT2, newDT

        type(DelaunayTrianglePointerList), pointer :: icdDT, tmp
        integer i

#if (!defined DelaunayAndVoronoi_FullSpeed)
#if (defined DelaunayAndVoronoi_Debug && DVT_ID)
        if (DVT%id == DVT_ID) then
            print *, "DVT", DVT%id, "is in MergeIncidentTriangle."
            print *, "  oldDT1", oldDT1%id, "oldDT2", oldDT2%id, "--> newDT", newDT%id
            print *, "Before merge:"
            call PrintVertexTopology(DVT)
        end if
#endif
#endif

        icdDT => DVT%icdDTHead
        do i = 1, DVT%numIcdDT
            if (associated(icdDT%ptr, oldDT1)) then
                icdDT%ptr => newDT
                exit
            end if
            icdDT => icdDT%next
        end do
        icdDT => DVT%icdDTHead
        do i = 1, DVT%numIcdDT
            if (associated(icdDT%ptr, oldDT2)) then
                icdDT%prev%next => icdDT%next
                icdDT%next%prev => icdDT%prev
                if (associated(icdDT, DVT%icdDTHead)) then
                    DVT%icdDTHead => icdDT%next
                end if
                deallocate(icdDT)
                exit
            end if
            icdDT => icdDT%next
        end do
        DVT%numIcdDT = DVT%numIcdDT-1

#if (!defined DelaunayAndVoronoi_FullSpeed)
#if (defined DelaunayAndVoronoi_Debug && DVT_ID)
        if (DVT%id == DVT_ID) then
            print *, "After merge:"
            call PrintVertexTopology(DVT)
        end if
#endif
#endif

    end subroutine MergeIncidentTriangle

    subroutine DeleteLinkVertex(DVT, delDVT)
        type(DelaunayVertex), intent(inout), target :: DVT
        type(DelaunayVertex), intent(in), target :: delDVT

        type(DelaunayVertexPointerList), pointer :: linkDVT
#if (defined DelaunayAndVoronoi_Debug)
        integer count
        count = 0
#endif

#if (!defined DelaunayAndVoronoi_FullSpeed)
        call MsgManager_RecordSpeaker("DeleteLinkVertex")

#if (defined DelaunayAndVoronoi_Debug && DVT_ID)
        if (DVT%id == DVT_ID) then
            print *, "DVT", DVT%id, "is in DeleteLinkVertex to delete DVT", delDVT%id
            print *, "Before delete:"
            call PrintVertexTopology(DVT)
        end if
#endif
#endif

        linkDVT => DVT%linkDVTHead
        do
            if (associated(linkDVT%ptr, delDVT)) then
                linkDVT%prev%next => linkDVT%next
                linkDVT%next%prev => linkDVT%prev
                if (associated(linkDVT, DVT%linkDVTHead)) then
                    DVT%linkDVTHead => linkDVT%next
                end if
                deallocate(linkDVT)
                DVT%numLinkDVT = DVT%numLinkDVT-1
#if (!defined DelaunayAndVoronoi_FullSpeed)
#if (defined DelaunayAndVoronoi_Debug && DVT_ID)
                if (DVT%id == DVT_ID) then
                    print *, "After delete:"
                    call PrintVertexTopology(DVT)
                end if
#endif
                call MsgManager_DeleteSpeaker
#endif
                return
            end if
#if (!defined DelaunayAndVoronoi_FullSpeed)
#if (defined DelaunayAndVoronoi_Debug)
            count = count+1
            if (count > DVT%numLinkDVT) then
                call MsgManager_Speak(Error, "DVT "// &
                    trim(int2str(delDVT%id))//" is not a link DVT to DVT "// &
                    trim(int2str(DVT%id)))
                call RunManager_EndRun
            end if
#endif
#endif
            linkDVT => linkDVT%next
        end do

        call MsgManager_Speak(Error, "DeleteLinkVertex error!")
        call RunManager_EndRun

    end subroutine DeleteLinkVertex

    subroutine AddLinkVertex(DVT, DVT1, DVT2, addDVT)
        type(DelaunayVertex), intent(inout), target :: DVT
        type(DelaunayVertex), intent(in), target :: DVT1, DVT2, addDVT

        type(DelaunayVertexPointerList), pointer :: linkDVT1, linkDVT2, linkDVT3

#if (!defined DelaunayAndVoronoi_FullSpeed)
#if (defined DelaunayAndVoronoi_Debug && DVT_ID)
        if (DVT%id == DVT_ID) then
            print *, "DVT", DVT%id, "is in AddLinkVertex to add DVT", addDVT%id
            print *, "Before add:"
            call PrintVertexTopology(DVT)
        end if
#endif
#endif

        linkDVT1 => DVT%linkDVTHead
        do
            linkDVT2 => linkDVT1%next
            if (associated(linkDVT1%ptr, DVT1) .and. &
                associated(linkDVT2%ptr, DVT2)) then
                allocate(linkDVT3)
                DVT%numLinkDVT = DVT%numLinkDVT+1
                linkDVT3%ptr => addDVT
                linkDVT1%next => linkDVT3
                linkDVT3%prev => linkDVT1
                linkDVT2%prev => linkDVT3
                linkDVT3%next => linkDVT2
#if (!defined DelaunayAndVoronoi_FullSpeed)
#if (defined DelaunayAndVoronoi_Debug && DVT_ID)
                if (DVT%id == DVT_ID) then
                    print *, "After add:"
                    call PrintVertexTopology(DVT)
                end if
#endif
#endif
                return
            end if
            linkDVT1 => linkDVT2
        end do

        call MsgManager_Speak(Error, "AddLinkVertex error!")
        call RunManager_EndRun

    end subroutine AddLinkVertex

    ! ************************************************************************ !
    ! RecordIncludedVertex                                                     !
    ! Purpose:                                                                 !
    !   Record the included vertex into the Delaunay triangle for speeding up  !
    !   the update of point-in-triangle relation by limiting the checking      !
    !   area.                                                                  !
    ! ************************************************************************ !

    subroutine RecordIncludedVertex(DVT, DT1, DT2)
        type(DelaunayVertex), intent(inout), target :: DVT
        type(DelaunayTriangle), intent(inout) :: DT1
        type(DelaunayTriangle), intent(inout), optional :: DT2

        type(DelaunayVertexPointerList), pointer :: tmp

        if (DT1%numIncDVT == 0) then
            allocate(DT1%incDVTHead)
            DT1%incDVTCurr => DT1%incDVTHead
        else
            allocate(DT1%incDVTCurr%next)
            tmp => DT1%incDVTCurr
            DT1%incDVTCurr => DT1%incDVTCurr%next
            DT1%incDVTCurr%prev => tmp
            tmp%next => DT1%incDVTCurr
        end if
        DT1%numIncDVT = DT1%numIncDVT+1
        DT1%incDVTCurr%ptr => DVT
        DVT%stub1 => DT1%incDVTCurr ! Keep the pointer for fast deletion.

        if (present(DT2)) then
            if (DT2%numIncDVT == 0) then
                allocate(DT2%incDVTHead)
                DT2%incDVTCurr => DT2%incDVTHead
            else
                allocate(DT2%incDVTCurr%next)
                tmp => DT2%incDVTCurr
                DT2%incDVTCurr => DT2%incDVTCurr%next
                DT2%incDVTCurr%prev => tmp
                tmp%next => DT2%incDVTCurr
            end if
            DT2%numIncDVT = DT2%numIncDVT+1
            DT2%incDVTCurr%ptr => DVT
            DVT%stub2 => DT2%incDVTCurr ! Keep the pointer for fast deletion.
        end if

    end subroutine RecordIncludedVertex

    ! ************************************************************************ !
    ! DeleteIncludedVertex                                                     !
    ! Purpose:                                                                 !
    !   After the insertion of the vertex, delete it from the included vertex  !
    !   list of the triangle.                                                  !
    ! ************************************************************************ !

    subroutine DeleteIncludedVertex(DVT, DT1, DT2)
        type(DelaunayVertex), intent(inout) :: DVT
        type(DelaunayTriangle), intent(inout) :: DT1
        type(DelaunayTriangle), intent(inout), optional :: DT2

        if (associated(DVT%stub1%prev)) then
            DVT%stub1%prev%next => DVT%stub1%next
        else
            DT1%incDVTHead => DT1%incDVTHead%next
        end if
        if (associated(DVT%stub1%next)) then
            DVT%stub1%next%prev => DVT%stub1%prev
        end if
        deallocate(DVT%stub1)
        DT1%numIncDVT = DT1%numIncDVT-1

        if (present(DT2)) then
            if (associated(DVT%stub2%prev)) then
                DVT%stub2%prev%next => DVT%stub2%next
            else
                DT2%incDVTHead => DT2%incDVTHead%next
            end if
            if (associated(DVT%stub2%next)) then
                DVT%stub2%next%prev => DVT%stub2%prev
            end if
            deallocate(DVT%stub2)
            DT2%numIncDVT = DT2%numIncDVT-1
        end if

    end subroutine DeleteIncludedVertex
 
    ! ************************************************************************ !
    ! RecordObsoleteTriangle                                                   !
    ! Purpose:                                                                 !
    !   Record the obsolete triangles for later deletion.                      !
    ! ************************************************************************ !

    subroutine RecordObsoleteTriangle(obsDT)
        type(DelaunayTriangle), intent(in), target :: obsDT

        type(DelaunayTrianglePointerList), pointer, save :: obsDTCurr => null()

        numObsDT = numObsDT+1
        if (numObsDT == 1) then
            allocate(obsDTHead)
            obsDTCurr => obsDTHead
        else
            allocate(obsDTCurr%next)
            obsDTCurr => obsDTCurr%next
        end if
        obsDTCurr%ptr => obsDT

    end subroutine RecordObsoleteTriangle

    subroutine DeleteObsoleteTriangle
        type(DelaunayTrianglePointerList), pointer :: obsDT
        integer i

        do i = 1, numObsDT
            obsDT => obsDTHead
            call DeleteDelaunayTriangle(obsDT%ptr)
            obsDTHead => obsDT%next
            deallocate(obsDT)
        end do
        numObsDT = 0
    
    end subroutine DeleteObsoleteTriangle
 
    ! ************************************************************************ !
    ! RecordTemporalTriangle                                                   !
    ! Purpose:                                                                 !
    !   Record the temporal triangles for later deletion.                      !
    ! ************************************************************************ !

    subroutine RecordTemporalTriangle(tmpDT)
        type(DelaunayTriangle), intent(in), target :: tmpDT

        type(DelaunayTrianglePointerList), pointer, save :: tmpDTCurr => null()

        numTmpDT = numTmpDT+1
        if (numTmpDT == 1) then
            allocate(tmpDTHead)
            tmpDTCurr => tmpDTHead
        else
            allocate(tmpDTCurr%next)
            tmpDTCurr => tmpDTCurr%next
        end if
        tmpDTCurr%ptr => tmpDT

    end subroutine RecordTemporalTriangle

    subroutine DeleteTemporalTriangle
        type(DelaunayTrianglePointerList), pointer :: tmpDT
        integer i

        do i = 1, numTmpDT
            tmpDT => tmpDTHead
            call DeleteDelaunayTriangle(tmpDT%ptr)
            tmpDTHead => tmpDT%next
            deallocate(tmpDT)
        end do
        numTmpDT = 0

    end subroutine DeleteTemporalTriangle

    subroutine DelaunayAndVoronoi_Report
        type(DelaunayTriangle), pointer :: DT
        type(DelaunayVertex), pointer :: DVT
        type(VoronoiCell), pointer :: VC
        integer i

        DT => DTHead
        do i = 1, numDT
            call PrintTriangleTopology(DT)
            DT => DT%next
        end do
        print *, "-------------------------"
        DVT => DVTHead
        do i = 1, numDVT
            call PrintVertexTopology(DVT)
            DVT => DVT%next
        end do
        print *, "-------------------------"
        VC => VCHead
        do i = 1, numVC
            call PrintVoronoiCell(VC)
            VC => VC%next
        end do
    
    end subroutine DelaunayAndVoronoi_Report
    
    subroutine PrintTriangleTopology(DT)
        type(DelaunayTriangle), intent(in) :: DT

        write(*, "('Triangle ID - ')", advance="no")
        write(*, *) trim(int2str(DT%id))
        write(*, "('  Vertex ID: ', 3I8)") &
            DT%DVT(1)%ptr%id, DT%DVT(2)%ptr%id, DT%DVT(3)%ptr%id
        write(*, "('  Adjacent triangle ID: ', 3I8)") &
            DT%adjDT(1)%ptr%id, DT%adjDT(2)%ptr%id, DT%adjDT(3)%ptr%id

    end subroutine PrintTriangleTopology

    subroutine PrintVertexTopology(DVT)
        type(DelaunayVertex), intent(in), target :: DVT

        type(DelaunayTrianglePointerList), pointer :: icdDT
        type(DelaunayVertexPointerList), pointer :: linkDVT
        integer i

        write(*, "('Vertex ID - ')", advance="no")
        write(*, *) trim(int2str(DVT%id))
        write(*, "('  Incident triangle ID: ')", advance="no")
        icdDT => DVT%icdDTHead
        do i = 1, DVT%numIcdDT
            write(*, "(I8)", advance="no") icdDT%ptr%id
            icdDT => icdDT%next
        end do
        write(*, *)
        if (associated(DVT%linkDVTHead)) then
            write(*, "('  Link vertex ID: ')", advance="no")
            linkDVT => DVT%linkDVTHead
            do i = 1, DVT%numLinkDVT
                write(*, "(I8)", advance="no") linkDVT%ptr%id
                linkDVT => linkDVT%next
            end do
            write(*, *)
        end if
    
    end subroutine PrintVertexTopology

    subroutine PrintVoronoiCell(VC)
        type(VoronoiCell), intent(in) :: VC
        integer i

        write(*, "('Voronoi cell ID - ')", advance="no")
        write(*, *) trim(int2str(VC%id))
        write(*, "('  Vertex coordinate (in degree): ')")
        do i = 1, VC%smp%bnd%numVtx
            write(*, "(2F10.3)") &
                VC%smp%bnd%vtx(i)%lon*Rad2Deg, VC%smp%bnd%vtx(i)%lat*Rad2Deg
        end do
    
    end subroutine PrintVoronoiCell
    
    subroutine PrintTriangleContainedVertex(DT)
        type(DelaunayTriangle), intent(in) :: DT

        type(DelaunayVertexPointerList), pointer :: incDVT
        integer i

        write(*, "('Triangle ID - ')", advance="no")
        write(*, *) trim(int2str(DT%id))
        write(*, "('  Contained vertex IDs ')", advance="no")
        write(*, *) trim(int2str(DT%numIncDVT))
        incDVT => DT%incDVTHead
        do i = 1, DT%numIncDVT
            write(*, "(I6)", advance="no") incDVT%ptr%id
            incDVT => incDVT%next
        end do
        write(*, *)

    end subroutine PrintTriangleContainedVertex

    subroutine PrintObsoleteTriangle
        type(DelaunayTrianglePointerList), pointer :: obsDT
        integer i

        write(*, "('Obsolete triangles ID:')")
        obsDT => obsDTHead
        do i = 1, numObsDT
            write(*, "(I8)", advance="no") obsDT%ptr%id
        end do
        write(*, *)

    end subroutine PrintObsoleteTriangle
 
    subroutine PrintTemporalTriangle
        type(DelaunayTrianglePointerList), pointer :: tmpDT
        integer i

        write(*, "('Temporal triangles ID:')")
        tmpDT => tmpDTHead
        do i = 1, numTmpDT
            write(*, "(I8)", advance="no") tmpDT%ptr%id
        end do
        write(*, *)

    end subroutine PrintTemporalTriangle

end module DelaunayAndVoronoi
 
