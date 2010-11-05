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
    use SphereService
    use SampleManager
    use NFWrap

    implicit none

    private

    public DelaunayAndVoronoi_LinkSample
    public ConstructDelaunayTriangulation
    public DelaunayAndVoronoi_Output
    public ExtractVoronoiDiagram
    public DelaunayAndVoronoi_Report

    ! ======================================================================= !
    ! Derived data types
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
        ! For debug
        integer :: id = -1
        ! For geometry
        type(Sample), pointer :: SMP => null()
        ! For topology
        integer :: numIcdDT = 0 ! Incident triangles
        type(DelaunayTrianglePointerList), pointer :: icdDTHead => null()
        type(DelaunayTrianglePointerList), pointer :: icdDTCurr => null()
        ! For point-in-triangle relation
        type(DelaunayTrianglePointer) cntDT(2)
        integer :: idx = -1 ! The shared edge index if InTriangle return >0.
        type(DelaunayVertexPointerList), pointer :: stub1, stub2
        ! For double link
        type(DelaunayVertex), pointer :: prev => null()
        type(DelaunayVertex), pointer :: next => null()
    end type DelaunayVertex

    type DelaunayTriangle
        ! For debug
        integer :: id = -1
        ! For topology
        type(DelaunayVertexPointer) DVT(3)
        type(DelaunayTrianglePointer) adjDT(3) ! Three adjacent DTs
        real(8) lonCirc, latCirc ! Center of circumcircle
        real(8) radius ! Radius of circumcircle in radian
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
        ! For debug
        integer :: id = -1
        ! For geometry
        integer :: numVVT
        real(8), allocatable :: lonVVT(:), latVVT(:)
        type(VoronoiCell), pointer :: prev => null()
        type(VoronoiCell), pointer :: next => null()
    end type VoronoiCell

    ! Data
    ! Delaunay vertex
    integer :: numDVT = 0
    type(DelaunayVertex), pointer :: DVTHead, DVTCurr
    type(DelaunayVertexPointer) VirtualDVT(3)
    ! Delaunay triangle
    integer :: numDT = 0, numTotalDT = 0
    type(DelaunayTriangle), pointer :: DTHead, DTCurr
    ! Obsolete DT record list
    ! Purpose: 1) Record the obsolete DTs for deletion;
    !          2) For "locally" updating of point-in-triangle relation
    integer :: numObsDT = 0
    type(DelaunayTrianglePointerList), pointer :: obsDTHead, obsDTCurr
    ! Temporal DT record list
    ! Purpose: 1) Record the temporal DTs for deletion;
    !          2) For updating point-in-triangle relation
    integer :: numTmpDT = 0 ! Record list of temporal DT
    type(DelaunayTrianglePointerList), pointer :: tmpDTHead, tmpDTCurr
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

contains

    subroutine DelaunayAndVoronoi_LinkSample
        type(Sample), pointer :: SMP
        type(DelaunayVertex), pointer :: DVT1, DVT2
        type(VoronoiCell), pointer :: VC1, VC2
        integer i

        call MsgManager_RecordSpeaker("DelaunayAndVoronoi_LinkSample")

        numDVT = numSample
        SMP => SMPHead
        allocate(DVTHead)
        DVTHead%id = 1
        DVTHead%SMP => SMP
        SMP => SMP%next
        DVT1 => DVTHead
        DVT2 => DVTHead
        do i = 2, numDVT
            allocate(DVT1%next)
            DVT1 => DVT1%next
            DVT1%id = i
            DVT1%SMP => SMP
            SMP => SMP%next
            DVT1%prev => DVT2
            DVT2%next => DVT1
            DVT2 => DVT1
        end do
        DVTCurr => DVTHead

        numVC = numSample
        allocate(VCHead)
        VCHead%id = 1
        VCHead%numVVT = 1
        allocate(VCHead%lonVVT(1))
        allocate(VCHead%latVVT(1))
        VC1 => VCHead
        VC2 => VCHead
        do i = 2, numVC
            allocate(VC1%next)
            VC1 => VC1%next
            VC1%id = i
            VC1%numVVT = 1
            allocate(VC1%lonVVT(1))
            allocate(VC1%latVVT(1))
            VC1%prev => VC2
            VC2%next => VC1
            VC2 => VC1
        end do

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker

    end subroutine DelaunayAndVoronoi_LinkSample
    
    ! ************************************************************************ !
    ! InitDelaunayTriangle                                                     !
    ! Purpose:                                                                 !
    !   Initialize Delaunay triangulation by creating the first real triangle  !
    !   and the associated seven virtual triangles.                            !
    ! ************************************************************************ !

    subroutine InitDelaunayTriangle()
        type(DelaunayVertexPointer) DVT(6)  ! 3 real + "3" virtual vertices
        type(DelaunayTrianglePointer) DT(8) ! 1 real + "7" virtual triangles
        integer map(24)
        integer i, j, ret

        call MsgManager_RecordSpeaker("InitDelaunayTriangle")

        ! Get the first three real vertices and set up the three virtual 
        ! vertices ...................................................... STEP 1
        ! Note: The virtual vertices are "antipodal" to the corresponding
        !       vertices of the first three inserted ones.
        do i = 1, 3
            DVT(i)%ptr => DVTCurr
            DVTCurr => DVTCurr%next
            allocate(VirtualDVT(i)%ptr)
            VirtualDVT(i)%ptr%id = -i
            allocate(VirtualDVT(i)%ptr%SMP)
            DVT(i+3)%ptr => VirtualDVT(i)%ptr ! Point to the virtual vertex
            call InverseRotationTransform(DVT(i)%ptr%SMP%lon, &
                DVT(i)%ptr%SMP%lat, VirtualDVT(i)%ptr%SMP%lon, &
                VirtualDVT(i)%ptr%SMP%lat, 0.0d0, -PI05)
            call CartesianTransform(VirtualDVT(i)%ptr%SMP%lon, &
                VirtualDVT(i)%ptr%SMP%lat, VirtualDVT(i)%ptr%SMP%x, &
                VirtualDVT(i)%ptr%SMP%y, VirtualDVT(i)%ptr%SMP%z)
        end do

        ! Create the first eight triangles .............................. STEP 2
        ! Note: One of them is real Delaunay triangle, and the rest are
        !       virtual to cover all the sphere.
        numDT = 1
        numTotalDT = 1
        allocate(DTHead)
        DTHead%id = 1
        DTCurr => DTHead
        DT(1)%ptr => DTHead
        do i = 2, 8
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
        ret = Orient(DVT(1)%ptr%SMP%x, DVT(1)%ptr%SMP%y, DVT(1)%ptr%SMP%z, &
                     DVT(2)%ptr%SMP%x, DVT(2)%ptr%SMP%y, DVT(2)%ptr%SMP%z, &
                     DVT(3)%ptr%SMP%x, DVT(3)%ptr%SMP%y, DVT(3)%ptr%SMP%z)
        if (ret == OrientLeft) then ! Counter-clockwise
            map = [1,2,3, 1,3,5, 1,5,6, 1,6,2, 4,3,2, 4,5,3, 4,6,5, 4,2,6]
        else if (ret == OrientRight) then ! Clockwise
            map = [1,3,2, 1,2,6, 1,6,5, 1,5,3, 4,2,3, 4,6,2, 4,5,6, 4,3,5]
        else if (ret == OrientOn) then
            ! The first three vertices are colinear.
        end if
        j = 1
        do i = 1, 8
            DT(i)%ptr%DVT(1)%ptr => DVT(map(j))%ptr; j = j+1
            DT(i)%ptr%DVT(2)%ptr => DVT(map(j))%ptr; j = j+1
            DT(i)%ptr%DVT(3)%ptr => DVT(map(j))%ptr; j = j+1
        end do

        ! Link the vertex with its incident triangles ................... STEP 5
        ! Note: The order is no matter.
        if (ret == OrientLeft) then ! Counter-clockwise
            map = [1,2,3,4, 1,4,8,5, 1,5,6,2, 5,8,7,6, 2,6,7,3, 3,7,8,4]
        else if (ret == OrientRight) then ! Clockwise
            map = [1,2,3,4, 1,5,6,2, 1,4,8,5, 5,8,7,6, 3,7,8,4, 2,6,7,3]
        else if (ret == OrientOn) then
            call MsgManager_Speak(Error, "The first vertices are colinear")
            call RunManager_EndRun
        end if
        do i = 1, 6
            do j = 1, 4
                call RecordIncidentTriangle(DVT(i)%ptr, DT(map((i-1)*4+j))%ptr)
            end do
        end do

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker

    end subroutine InitDelaunayTriangle
    
    ! ************************************************************************ !
    ! InsertVertex                                                             !
    ! Purpose:                                                                 !
    !   Insert one vertex into the Delaunay triangulation, and do some flip    !
    !   operations if there are "topological events".                          !
    ! ************************************************************************ !

    subroutine InsertVertex(DVT)
        type(DelaunayVertex), intent(inout), target :: DVT

        type(DelaunayVertex), pointer :: DVT1, DVT2, DVT3, DVT4
        type(DelaunayTriangle), pointer :: adjDT1, adjDT2, adjDT3, adjDT4
        type(DelaunayTriangle), pointer :: oldDT1, oldDT2, adjDT
        type(DelaunayTrianglePointer) newDT(4)
        integer, save :: ip1(3) = [2,3,1],   im1(3) = [3,1,2]
        integer i, j, idx

        if (associated(DVT%cntDT(2)%ptr)) then
            ! DVT is on the edge shared by two adjacent triangles.
            ! Note: In this case, there is axial dissymetry with respect to
            !       DVT topologically. Therefore, the indices of the vertices
            !       of oldDT1 does matter.
            ! For short-hand ............................................ STEP 0
            oldDT1 => DVT%cntDT(1)%ptr
            oldDT2 => DVT%cntDT(2)%ptr
            do idx = 1, 3
                if (associated(oldDT2%adjDT(idx)%ptr, oldDT1)) exit
            end do
            DVT1 => oldDT1%DVT(DVT%idx)%ptr          
            DVT2 => oldDT1%DVT(ip1(DVT%idx))%ptr     
            DVT3 => oldDT1%DVT(im1(DVT%idx))%ptr     
            DVT4 => oldDT2%DVT(idx)%ptr              
            adjDT1 => oldDT1%adjDT(ip1(DVT%idx))%ptr 
            adjDT2 => oldDT1%adjDT(im1(DVT%idx))%ptr 
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
            ! Change existing vertices' incident triangle
            ! ................................................... newDT(4)
            newDT(4)%ptr%DVT(1)%ptr => DVT4
            newDT(4)%ptr%DVT(2)%ptr => DVT3
            newDT(4)%ptr%DVT(3)%ptr => DVT
            newDT(4)%ptr%adjDT(1)%ptr => newDT(1)%ptr
            newDT(4)%ptr%adjDT(2)%ptr => newDT(3)%ptr
            newDT(4)%ptr%adjDT(3)%ptr => adjDT4
            ! Link newly inserted DVT to its incident DTs ............... STEP 3
            do i = 1, 4
                call RecordIncidentTriangle(DVT, newDT(i)%ptr)
            end do
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
            call ReplaceIncidentTriangle(DVT2, oldDT1, newDT(2)%ptr)
            call ReplaceIncidentTriangle(DVT2, oldDT2, newDT(3)%ptr)
            call ReplaceIncidentTriangle(DVT3, oldDT1, newDT(1)%ptr)
            call ReplaceIncidentTriangle(DVT3, oldDT2, newDT(4)%ptr)
            call SplitIncidentTriangle(DVT1, oldDT1, newDT(2)%ptr, newDT(1)%ptr)
            call SplitIncidentTriangle(DVT4, oldDT2, newDT(4)%ptr, newDT(3)%ptr)
#if (defined DEBUG && defined OUTPUTSTAGE)
            ! Testing output ...
            call DelaunayAndVoronoi_Output("stage1-"//trim(int2str(DVT%id))//".nc")
#endif
            ! Validate the three triangles .............................. STEP 5
            do i = 1, 4
                call ValidateNewTriangle(newDT(i)%ptr)
            end do
#if (defined DEBUG && defined OUTPUTSTAGE)
            ! Testing output ...
            call DelaunayAndVoronoi_Output("stage2-"//trim(int2str(DVT%id))//".nc")
#endif
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
            DVT%cntDT(1)%ptr => null()
            DVT%cntDT(2)%ptr => null()
            call DeleteIncludedVertex(DVT, oldDT1, oldDT2)
        else
            ! DVT is only contained by one triangle.
            ! Note: In this case, there is no axial dissymetry with respect to
            !       DVT topologically. Therefore, the indices of the vertices
            !       of oldDT1 does not matter.
            ! For short-hand ............................................ STEP 0
            oldDT1 => DVT%cntDT(1)%ptr
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
            ! Link newly inserted DVT to its incident DTs ............... STEP 3
            do i = 1, 3
                call RecordIncidentTriangle(DVT, newDT(i)%ptr)
            end do
            ! Make change to the old triangle's adjDT and DVT ........... STEP 4
            ! Note: There is an awkful situation as that we know
            !       the adjacent triangle of oldDT1, but we do not
            !       know what is the index of oldDT1 in that triangle.
            !       So we do not know to change which adjDT of that
            !       triangle to the newDT directly, but have to check it.
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
            call SplitIncidentTriangle(DVT1, oldDT1, newDT(1)%ptr, newDT(3)%ptr)
            call SplitIncidentTriangle(DVT2, oldDT1, newDT(2)%ptr, newDT(1)%ptr)
            call SplitIncidentTriangle(DVT3, oldDT1, newDT(3)%ptr, newDT(2)%ptr)
#if (defined DEBUG && defined OUTPUTSTAGE)
            ! Testing output ...
            call DelaunayAndVoronoi_Output("stage1-"//trim(int2str(DVT%id))//".nc")
#endif
            ! Validate the three triangles .............................. STEP 5
            do i = 1, 3
                call ValidateNewTriangle(newDT(i)%ptr)
            end do
#if (defined DEBUG && defined OUTPUTSTAGE)
            ! Testing output ...
            call DelaunayAndVoronoi_Output("stage2-"trim(int2str(DVT%id))//".nc")
#endif
            ! Record the obsolete triangle .............................. STEP 6
            call RecordObsoleteTriangle(oldDT1)
            oldDT1%numSubDT = 3
            oldDT1%subDT(1)%ptr => newDT(1)%ptr
            oldDT1%subDT(2)%ptr => newDT(2)%ptr
            oldDT1%subDT(3)%ptr => newDT(3)%ptr
            ! Remove DVT from point-in-triangle relation ................ STEP 7
            DVT%cntDT(1)%ptr => null()
            call DeleteIncludedVertex(DVT, oldDT1)
        end if

    end subroutine InsertVertex
    
    ! ************************************************************************ !
    ! ValidateNewTriangle                                                      !
    ! Purpose:                                                                 !
    !   Validate the new triangle with its DVT(3) as the newly inserted vertex !
    ! ************************************************************************ !

    recursive subroutine ValidateNewTriangle(newDT)
        type(DelaunayTriangle), intent(inout), target :: newDT
    
        type(DelaunayTriangle), pointer :: opsDT, newerDT1, newerDT2
        type(DelaunayVertex), pointer :: opsDVT
        type(DelaunayTriangle), pointer :: adjDT1, adjDT2, adjDT3, adjDT4
        type(DelaunayVertex), pointer :: DVT1, DVT2, DVT3, DVT4
        integer :: im1(3) = [3,1,2], ip1(3) = [2,3,1]
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
            return
        else if (ret == InsideCircle) then
            ! Flip the shared edge
            ! * Create two newer DT to replace newDT and opsDT.
            call NewDelaunayTriangle(newerDT1)
            call NewDelaunayTriangle(newerDT2)
            newerDT1%DVT(1)%ptr => newDT%DVT(1)%ptr
            newerDT1%DVT(2)%ptr => opsDVT
            newerDT1%DVT(3)%ptr => newDT%DVT(3)%ptr
            newerDT1%adjDT(1)%ptr => newerDT2
            newerDT1%adjDT(2)%ptr => newDT%adjDT(2)%ptr
            newerDT1%adjDT(3)%ptr => opsDT%adjDT(ip1(idx))%ptr
            newerDT2%DVT(1)%ptr => opsDVT
            newerDT2%DVT(2)%ptr => newDT%DVT(2)%ptr
            newerDT2%DVT(3)%ptr => newDT%DVT(3)%ptr
            newerDT2%adjDT(1)%ptr => newDT%adjDT(1)%ptr
            newerDT2%adjDT(2)%ptr => newerDT1
            newerDT2%adjDT(3)%ptr => opsDT%adjDT(im1(idx))%ptr
            ! * Make change to the old triangle's adjDT and DVT
            ! For short-hand
            DVT1 => newDT%DVT(1)%ptr
            DVT2 => newDT%DVT(2)%ptr
            DVT3 => newDT%DVT(3)%ptr
            DVT4 => opsDT%DVT(idx)%ptr
            adjDT1 => newDT%adjDT(1)%ptr
            adjDT2 => newDT%adjDT(2)%ptr
            adjDT3 => opsDT%adjDT(ip1(idx))%ptr
            adjDT4 => opsDT%adjDT(im1(idx))%ptr
            do i = 1, 3
                if (associated(adjDT1%adjDT(i)%ptr, newDT)) then
                    adjDT1%adjDT(i)%ptr => newerDT2
                    exit
                end if
            end do
            do i = 1, 3
                if (associated(adjDT2%adjDT(i)%ptr, newDT)) then
                    adjDT2%adjDT(i)%ptr => newerDT1
                    exit
                end if
            end do
            do i = 1, 3
                if (associated(adjDT3%adjDT(i)%ptr, opsDT)) then
                    adjDT3%adjDT(i)%ptr => newerDT1
                    exit
                end if
            end do
            do i = 1, 3
                if (associated(adjDT4%adjDT(i)%ptr, opsDT)) then
                    adjDT4%adjDT(i)%ptr => newerDT2
                    exit
                end if
            end do
            call MergeIncidentTriangle(DVT1, newDT, opsDT, newerDT1)
            call MergeIncidentTriangle(DVT2, newDT, opsDT, newerDT2)
            call SplitIncidentTriangle(DVT3, newDT, newerDT1, newerDT2)
            call SplitIncidentTriangle(DVT4, opsDT, newerDT2, newerDT1)
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
            call MsgManager_Speak(Warning, "Encounter cocircle vertices.")
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
    !   recording the obsolete triangle pointers.                              !
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
            ret = Orient(DT%DVT(3)%ptr%SMP%x, DT%DVT(3)%ptr%SMP%y, &
                         DT%DVT(3)%ptr%SMP%z, DT%DVT(1)%ptr%SMP%x, &
                         DT%DVT(1)%ptr%SMP%y, DT%DVT(1)%ptr%SMP%z, &
                         DVT%SMP%x, DVT%SMP%y, DVT%SMP%z)
            if (ret == OrientRight) then
                found = .false.
                return
            else if (ret == OrientOn) then
                onPlane(1) = 2
            end if
            ! Check edge 2->3
            ret = Orient(DT%DVT(2)%ptr%SMP%x, DT%DVT(2)%ptr%SMP%y, &
                         DT%DVT(2)%ptr%SMP%z, DT%DVT(3)%ptr%SMP%x, &
                         DT%DVT(3)%ptr%SMP%y, DT%DVT(3)%ptr%SMP%z, &
                         DVT%SMP%x, DVT%SMP%y, DVT%SMP%z)
            if (ret == OrientRight) then
                found = .false.
                return
            else if (ret == OrientOn) then
                onPlane(2) = 1
            end if
            ! Summary
            found = .true.
            if (onPlane(1) == 0 .and. onPlane(2) == 0) then
                DVT%cntDT(1)%ptr => DT
                DVT%cntDT(2)%ptr => null()
                call RecordIncludedVertex(DVT, DT)
            else if (onPlane(1) == 2 .and. onPlane(2) == 0) then
                DVT%cntDT(1)%ptr => DT
                DVT%cntDT(2)%ptr => DT%adjDT(2)%ptr
                DVT%idx = onPlane(1)
                call RecordIncludedVertex(DVT, DT, DT%adjDT(2)%ptr)
            else if (onPlane(1) == 0 .and. onPlane(2) == 1) then
                DVT%cntDT(1)%ptr => DT
                DVT%cntDT(2)%ptr => DT%adjDT(1)%ptr
                DVT%idx = onPlane(2)
                call RecordIncludedVertex(DVT, DT, DT%adjDT(1)%ptr)
            else
                call MsgManager_Speak(Error, "DVT "//trim(int2str(DVT%id))// &
                    " is on vertex 3 of DT "//trim(int2str(DT%id)))
                call RunManager_EndRun
            end if
        end if

    end subroutine UpdatePointInTriangle

    ! ************************************************************************ !
    ! ConstructDelaunayTriangulation                                           !
    ! Purpose:                                                                 !
    !   This subroutine is only called at the initial time step. When the      !
    !   samples are ready, construct Delaunay triangulation from them.         !
    ! ************************************************************************ !

    subroutine ConstructDelaunayTriangulation
        type(DelaunayVertex), pointer :: DVT
        type(DelaunayTriangle), pointer :: DT1, DT2
        logical found
        integer i, j, k, ret
    
        call MsgManager_RecordSpeaker("ConstructDelaunayTriangulation")

        ! Maybe add random permutation of the vertices.

        ! Initialize the triangles comprised of the first three vertices and
        ! their antipodal counterparts, and shift DVTCurr ............... STEP 1
        ! Note: On the sphere, the boundary condition is periodic, so there
        !       are seven other virtual triangles accompanying the one real 
        !       triangle.
        call InitDelaunayTriangle
#if (defined DEBUG && defined OUTPUTSTAGE)
        ! Testing output ...
        call DelaunayAndVoronoi_Output("init.nc")
#endif

        ! Initialize the point-in-triangle relation between the rest vertices
        ! and the created triangles ..................................... STEP 2
        DVT => DVTCurr
        do i = 4, numDVT
            DT1 => DTHead
            do j = 1, numDT
                ret = InTriangle(DT1, DVT)
                if (ret == InsideTriangle) then
                    ! Point to the triangle that contains DVT
                    DVT%cntDT(1)%ptr => DT1
                    DVT%cntDT(2)%ptr => null()
                    ! Also record the DVT into the list of included vertices
                    ! of DT
                    call RecordIncludedVertex(DVT, DT1)
                    exit
                else if (ret == OutsideTriangle) then
                    DT1 => DT1%next
                    cycle
                else if (ret > 0) then
                    ! DVT is on some edge of DT.
                    DT2 => DT1%adjDT(ret)%ptr
                    DVT%cntDT(1)%ptr => DT1
                    DVT%cntDT(2)%ptr => DT2
                    DVT%idx = ret
                    call RecordIncludedVertex(DVT, DT1, DT2)
                    exit
                else if (ret < 0) then
                    ! DVT coincides with some vertex of DT.
                    ! Complain this degenerate case.
                    call MsgManager_Speak(Error, "DVT "//trim(int2str(DVT%id))// &
                        " coincides to vertex "//trim(int2str(-ret))// &
                        " of DT "//trim(int2str(DT1%id))//".")
                    call RunManager_EndRun
                end if
            end do
            DVT => DVT%next
        end do

        !call MsgManager_Speak(Notice, "Initialization of "// &
        !    "point-in-triangle is finished.")

        ! Insert the rest vertices one at a time ........................ STEP 3
        do i = 4, numDVT
            ! Update the Delauny triangulation ........................ STEP 3.1
            call InsertVertex(DVTCurr)
            !call MsgManager_Speak(Notice, "DVT "// &
            !    trim(int2str(DVTCurr%id))//" has been inserted.")
            ! Update the point-in-triangle relation "locally" ......... STEP 3.2
            obsDTCurr => obsDTHead
            do j = 1, numObsDT
                obsDTCurr%ptr%incDVTCurr => obsDTCurr%ptr%incDVTHead
                do k = 1, obsDTCurr%ptr%numIncDVT
                    call UpdatePointInTriangle( &
                        obsDTCurr%ptr%incDVTCurr%ptr, obsDTCurr%ptr, found)
                    obsDTCurr%ptr%incDVTCurr => obsDTCurr%ptr%incDVTCurr%next
                end do
                obsDTCurr => obsDTCurr%next
            end do
            ! Delete obsolete and temporal triangles .................. STEP 3.3
            call DeleteObsoleteTriangle
            call DeleteTemporalTriangle
#if (defined DEBUG && defined OUTPUTSTAGE)
            ! Testing output ...
            call DelaunayAndVoronoi_Output("final"//trim(int2str(i))//".nc")
#endif
            ! Shift the current DVT pointer ........................... STEP 3.4
            DVTCurr => DVTCurr%next
        end do

        ! Calculate the circumcenter and circumradius
        call CalcCircumcircle

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
        integer i, j

        maxNumVVT = 0
        DVT => DVTHead
        VC => VCHead
        do i = 1, numDVT ! = numVC
            if (DVT%numIcdDT > maxNumVVT) maxNumVVT = DVT%numIcdDT
            if (VC%numVVT /= DVT%numIcdDT) then
                deallocate(VC%lonVVT)
                deallocate(VC%latVVT)
                VC%numVVT = DVT%numIcdDT
                allocate(VC%lonVVT(VC%numVVT))
                allocate(VC%latVVT(VC%numVVT))
            end if
            icdDT => DVT%icdDTHead
            do j = 1, DVT%numIcdDT
                VC%lonVVT(j) = icdDT%ptr%lonCirc
                VC%latVVT(j) = icdDT%ptr%latCirc
                icdDT => icdDT%next
            end do
            DVT => DVT%next
            VC => VC%next
        end do

    end subroutine ExtractVoronoiDiagram
    
    ! ************************************************************************ !
    ! Orient                                                                   !
    ! Purpose:                                                                 !
    !   This is the one of two basic geometric tests to give the relation      !
    !   between a point (x0,y0,z0) and a plane containing (x1,y1,z1),          !
    !   (x2,y2,z2), and (0,0,0). There are three types of relation, i.e. left, !
    !   right, and on, where left is defined relative to an observer at        !
    !   (x1,y1,z1) facing (x2,y2,z2).                                          !
    ! ************************************************************************ !

    integer function Orient(x1, y1, z1, x2, y2, z2, x0, y0, z0)
        real(8), intent(in) :: x1, y1, z1, x2, y2, z2, x0, y0, z0
    
        real(8) det

        det = x0*(y1*z2-y2*z1)-y0*(x1*z2-x2*z1)+z0*(x1*y2-x2*y1)

        if (det > eps) then
            Orient = OrientLeft
        else if (-det > eps) then
            Orient = OrientRight
        else
            Orient = OrientOn
        end if

    end function Orient

    ! ************************************************************************ !
    ! InTriangle                                                               !
    ! Purpose:                                                                 !
    !   This is a derived geometric test from "Orient". It gives the relation  !
    !   between a point (DVT) and a triangle (DT). There are four types of     !
    !   relation, i.e. inside, outside, on the edge and on the vertex          !
    ! Return value:                                                            !
    !   InsideTriangle                                                         !
    !   OutsideTriangle                                                        !
    !    1 ~  3 for on the edge, where the value indicates the edge index      !
    !   -1 ~ -3 for on the vertex, where the abs(value) indicates the vertex   !
    !           index                                                          !
    ! ************************************************************************ !

    integer function InTriangle(DT, DVT)
        type(DelaunayTriangle), intent(in) :: DT
        type(DelaunayVertex), intent(in) :: DVT

        real(8) x(3), y(3), z(3), x0, y0, z0
        integer, save :: j(3) = [2,3,1], l(3) = [3,1,2]
        integer i, k, ret, onPlane(2)

        do i = 1, 3
            x(i) = DT%DVT(i)%ptr%SMP%x
            y(i) = DT%DVT(i)%ptr%SMP%y
            z(i) = DT%DVT(i)%ptr%SMP%z
        end do
        x0 = DVT%SMP%x
        y0 = DVT%SMP%y
        z0 = DVT%SMP%z

        k = 0
        do i = 1, 3
            ret = Orient(x(i),    y(i),    z(i),    &
                         x(j(i)), y(j(i)), z(j(i)), &
                         x0,      y0,      z0)
            if (ret == OrientRight) then
                InTriangle = OutsideTriangle
                return
            else if (ret == OrientOn) then
                k = k+1; onPlane(k) = l(i)
            end if
        end do

        if (k == 0) then
            InTriangle = InsideTriangle
        else if (k == 1) then ! on the edge
            InTriangle = onPlane(k)
        else if (k == 2) then ! on the vertex
            if (onPlane(1) == 1 .and. onPlane(2) == 2) then
                InTriangle = -3
            else if (onPlane(1) == 2 .and. onPlane(2) == 3) then
                InTriangle = -1
            else if (onPlane(1) == 1 .and. onPlane(2) == 3) then
                InTriangle = -2
            end if
        end if

    end function InTriangle

    ! ************************************************************************ !
    ! InCircle                                                                 !
    ! Purpose:                                                                 !
    !   This is the one of two basic geometric tests to check whether or not a !
    !   point (DVT) is inside the circumcircle (here is spherical circle) of a !
    !   Delaunay triangle (DT).                                                !
    !   In Cartesian coordinate system, this test turns into checking whether  !
    !   or not DVT is above the plane defined by the three vertices of DT.     !
    ! ************************************************************************ !

    integer function InCircle(DT, DVT)
        type(DelaunayTriangle), intent(in) :: DT
        type(DelaunayVertex), intent(in) :: DVT

        real(8) dx(3), dy(3), dz(3), x0, y0, z0, det
        integer i

        x0 = DVT%SMP%x
        y0 = DVT%SMP%y
        z0 = DVT%SMP%z
        do i = 1, 3
            dx(i) = DT%DVT(i)%ptr%SMP%x-x0
            dy(i) = DT%DVT(i)%ptr%SMP%y-y0
            dz(i) = DT%DVT(i)%ptr%SMP%z-z0
        end do

        det = dx(3)*(dy(2)*dz(1)-dy(1)*dz(2)) &
             -dy(3)*(dx(2)*dz(1)-dx(1)*dz(2)) &
             +dz(3)*(dx(2)*dy(1)-dx(1)*dy(2))

        if (det > eps) then
            InCircle = InsideCircle
        else if (-det > eps) then
            InCircle = OutsideCircle
        else
            InCircle = OnTheCircle
        end if

    end function InCircle

    subroutine CalcCircumcircle
        type(DelaunayTriangle), pointer :: DT
        real(8) E2(3), E3(3), N(3), L, C(3), tmp
        integer i

        DT => DTHead
        do i = 1, numDT
            E3 = [DT%DVT(2)%ptr%SMP%x-DT%DVT(1)%ptr%SMP%x, &
                  DT%DVT(2)%ptr%SMP%y-DT%DVT(1)%ptr%SMP%y, &
                  DT%DVT(2)%ptr%SMP%z-DT%DVT(1)%ptr%SMP%z]
            E2 = [DT%DVT(3)%ptr%SMP%x-DT%DVT(1)%ptr%SMP%x, &
                  DT%DVT(3)%ptr%SMP%y-DT%DVT(1)%ptr%SMP%y, &
                  DT%DVT(3)%ptr%SMP%z-DT%DVT(1)%ptr%SMP%z]
            N = [E3(2)*E2(3)-E3(3)*E2(2), &
                 E3(3)*E2(1)-E3(1)*E2(3), &
                 E3(1)*E2(2)-E3(2)*E2(1)]
            L = N(1)**2+N(2)**2+N(3)**2
            if (L == 0.0d0) then
                call MsgManager_Speak(Error, "Vertices are colinear")
                call RunManager_EndRun
            end if
            L = sqrt(L)
            C = N/L*Re
            call InverseCartesianTransform(DT%lonCirc, DT%latCirc, C(1), C(2), C(3))
            tmp = (DT%DVT(1)%ptr%SMP%x*C(1)+ &
                   DT%DVT(1)%ptr%SMP%y*C(2)+ &
                   DT%DVT(1)%ptr%SMP%z*C(3))/Re**2
            DT%radius = acos(tmp)
            DT => DT%next
        end do

    end subroutine CalcCircumcircle
    
    ! ************************************************************************ !
    ! NewDelaunayTriangle                                                      !
    ! Purpose:                                                                 !
    !   Create the storage of a new Delaunay triangle in the double linked     !
    !   list.                                                                  !
    ! ************************************************************************ !

    subroutine NewDelaunayTriangle(DT)
        type(DelaunayTriangle), intent(out), pointer :: DT

        type(DelaunayTriangle), pointer :: DT1

        numDT = numDT+1
        numTotalDT = numTotalDT+1
        DT1 => DTCurr
        allocate(DT1%next)
        DT1 => DT1%next
        DT1%id = numTotalDT
        DT1%prev => DTCurr
        DTCurr%next => DT1
        DTCurr => DT1
        DT => DT1
 
        !call MsgManager_Speak(Notice, "DT "//trim(int2str(DT%id))// &
        !    " has been newed.")

    end subroutine NewDelaunayTriangle
 
    ! ************************************************************************ !
    ! DeleteDelaunayTriangle                                                   !
    ! Purpose:                                                                 !
    !   Delete the storage of a Delaunay triangle from the double linked list. !
    ! ************************************************************************ !

    subroutine DeleteDelaunayTriangle(DT)
        type(DelaunayTriangle), intent(inout), pointer :: DT

        !call MsgManager_Speak(Notice, "DT "//trim(int2str(DT%id))// &
        !    " has been deleted.")

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
    
    subroutine RecordIncidentTriangle(DVT, DT)
        type(DelaunayVertex), intent(inout) :: DVT
        type(DelaunayTriangle), intent(in), target :: DT

        type(DelaunayTrianglePointerList), pointer :: tmp

        if (DVT%numIcdDT == 0) then
            ! Initialize the list head.
            allocate(DVT%icdDTHead)
            DVT%icdDTCurr => DVT%icdDTHead
            DVT%icdDTCurr%ptr => DT
            DVT%numIcdDT = 1
        else
            ! Append an element to the end.
            allocate(DVT%icdDTCurr%next)
            tmp => DVT%icdDTCurr
            DVT%icdDTCurr => DVT%icdDTCurr%next
            DVT%icdDTCurr%prev => tmp
            tmp%next => DVT%icdDTCurr
            DVT%icdDTCurr%ptr => DT
            DVT%numIcdDT = DVT%numIcdDT+1
        end if
    
    end subroutine RecordIncidentTriangle

    subroutine ReplaceIncidentTriangle(DVT, oldDT, newDT)
        type(DelaunayVertex), intent(inout) :: DVT
        type(DelaunayTriangle), intent(in), target :: oldDT, newDT

        type(DelaunayTrianglePointerList), pointer :: icdDT
        integer i

        icdDT => DVT%icdDTHead
        do i = 1, DVT%numIcdDT
            if (associated(icdDT%ptr, oldDT)) then
                icdDT%ptr => newDT
                exit
            end if
            icdDT => icdDT%next
        end do

    end subroutine ReplaceIncidentTriangle

    subroutine SplitIncidentTriangle(DVT, oldDT, newDT1, newDT2)
        type(DelaunayVertex), intent(inout) :: DVT
        type(DelaunayTriangle), intent(in), target :: oldDT, newDT1, newDT2

        type(DelaunayTrianglePointerList), pointer :: icdDT1, icdDT2
        integer i

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

    end subroutine SplitIncidentTriangle

    subroutine MergeIncidentTriangle(DVT, oldDT1, oldDT2, newDT)
        type(DelaunayVertex), intent(inout) :: DVT
        type(DelaunayTriangle), intent(in), target :: oldDT1, oldDT2, newDT

        type(DelaunayTrianglePointerList), pointer :: icdDT, tmp
        logical flag
        integer i

        !if (DVT%id == 3) then
        !    print *, "Merge oldDT1", oldDT1%id, "and oldDT2", oldDT2%id, "to newDT", newDT%id
        !    print *, "Before merge:"
        !    icdDT => DVT%icdDTHead
        !    do i = 1, DVT%numIcdDT
        !        print *, icdDT%ptr%id, &
        !            merge(icdDT%prev%ptr%id, -1, associated(icdDT%prev)), &
        !            merge(icdDT%next%ptr%id, -1, associated(icdDT%next))
        !        icdDT => icdDT%next
        !    end do
        !end if

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
                if (associated(icdDT%prev)) then
                    icdDT%prev%next => icdDT%next
                else
                    DVT%icdDTHead => DVT%icdDTHead%next
                end if
                if (associated(icdDT%next)) then
                    icdDT%next%prev => icdDT%prev
                end if
                deallocate(icdDT)
                exit
            end if
            icdDT => icdDT%next
        end do
        DVT%numIcdDT = DVT%numIcdDT-1

        !flag = .false.
        !icdDT => DVT%icdDTHead
        !do i = 1, DVT%numIcdDT
        !    if (associated(icdDT%ptr, oldDT1)) then
        !        icdDT%ptr => newDT
        !        if (flag .eqv. .false.) then
        !            flag = .true.
        !            icdDT => icdDT%next
        !            cycle
        !        else
        !            exit
        !        end if
        !    else if (associated(icdDT%ptr, oldDT2)) then
        !        ! Delete icdDT
        !        if (associated(icdDT%prev)) then
        !            icdDT%prev%next => icdDT%next
        !            if (associated(icdDT%next)) then
        !                icdDT%next%prev => icdDT%prev
        !            end if
        !        else
        !            DVT%icdDTHead => icdDT%next
        !        end if
        !        if (flag .eqv. .false.) then
        !            flag = .true.
        !            tmp => icdDT
        !            icdDT => icdDT%next
        !            deallocate(tmp) !!!
        !            cycle
        !        else
        !            exit
        !        end if
        !    end if
        !    icdDT => icdDT%next
        !end do
        !DVT%numIcdDT = DVT%numIcdDT-1

        !if (DVT%id == 3) then
        !    print *, "After merge:"
        !    icdDT => DVT%icdDTHead
        !    do i = 1, DVT%numIcdDT
        !        print *, icdDT%ptr%id, &
        !            merge(icdDT%prev%ptr%id, -1, associated(icdDT%prev)), &
        !            merge(icdDT%next%ptr%id, -1, associated(icdDT%next))
        !        icdDT => icdDT%next
        !    end do
        !end if

    end subroutine MergeIncidentTriangle
    
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

        ! Testing output ...
        !if (DVT%id == 7) then
        !    print *, "Before record DVT", DVT%id
        !    call PrintTriangleContainedVertex(DT1)
        !    if (present(DT2)) call PrintTriangleContainedVertex(DT2)
        !end if

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

        ! Testing output ...
        !if (DVT%id == 7) then
        !    print *, "After record DVT", DVT%id
        !    call PrintTriangleContainedVertex(DT1)
        !    if (present(DT2)) call PrintTriangleContainedVertex(DT2)
        !end if

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

        ! Testing output ...
        !if (DVT%id == 7) then
        !    print *, "Before delete DVT", DVT%id
        !    call PrintTriangleContainedVertex(DT1)
        !    if (present(DT2)) call PrintTriangleContainedVertex(DT2)
        !end if

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

        ! Testing output ...
        !if (DVT%id == 7) then
        !    print *, "After delete DVT", DVT%id
        !    call PrintTriangleContainedVertex(DT1)
        !    if (present(DT2)) call PrintTriangleContainedVertex(DT2)
        !end if

    end subroutine DeleteIncludedVertex
 
    ! ************************************************************************ !
    ! RecordObsoleteTriangle                                                   !
    ! Purpose:                                                                 !
    !   Record the obsolete triangles for later deletion.                      !
    ! ************************************************************************ !

    subroutine RecordObsoleteTriangle(obsDT)
        type(DelaunayTriangle), intent(in), target :: obsDT

        if (numObsDT == 0) then
            allocate(obsDTHead)
            obsDTCurr => obsDTHead
            obsDTCurr%ptr => obsDT
            numObsDT = 1
        else
            allocate(obsDTCurr%next)
            obsDTCurr => obsDTCurr%next
            obsDTCurr%ptr => obsDT
            numObsDT = numObsDT+1
        end if

    end subroutine RecordObsoleteTriangle

    subroutine DeleteObsoleteTriangle
        integer i

        do i = 1, numObsDT
            obsDTCurr => obsDTHead
            call DeleteDelaunayTriangle(obsDTCurr%ptr)
            obsDTHead => obsDTCurr%next
            deallocate(obsDTCurr)
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

        if (numTmpDT == 0) then
            allocate(tmpDTHead)
            tmpDTCurr => tmpDTHead
            tmpDTCurr%ptr => tmpDT
            numTmpDT = 1
        else
            allocate(tmpDTCurr%next)
            tmpDTCurr => tmpDTCurr%next
            tmpDTCurr%ptr => tmpDT
            numTmpDT = numTmpDT+1
        end if

    end subroutine RecordTemporalTriangle

    subroutine DeleteTemporalTriangle
        integer i

        do i = 1, numTmpDT
            tmpDTCurr => tmpDTHead
            call DeleteDelaunayTriangle(tmpDTCurr%ptr)
            tmpDTHead => tmpDTCurr%next
            deallocate(tmpDTCurr)
        end do
        numTmpDT = 0

    end subroutine DeleteTemporalTriangle

    subroutine DelaunayAndVoronoi_Report()
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
        !print *, "-------------------------"
        !VC => VCHead
        !do i = 1, numVC
        !    call PrintVoronoiCell(VC)
        !    VC => VC%next
        !end do
    
    end subroutine DelaunayAndVoronoi_Report
    
    subroutine PrintTriangleTopology(DT)
        type(DelaunayTriangle), intent(in) :: DT

        write(*, "('Triangle (ID -', I3, ')')") DT%id
        write(*, "('  Vertex ID: ', 3I3)") &
            DT%DVT(1)%ptr%id, DT%DVT(2)%ptr%id, DT%DVT(3)%ptr%id
        write(*, "('  Adjacent triangle ID: ', 3I3)") &
            DT%adjDT(1)%ptr%id, DT%adjDT(2)%ptr%id, DT%adjDT(3)%ptr%id

    end subroutine PrintTriangleTopology

    subroutine PrintVertexTopology(DVT)
        type(DelaunayVertex), intent(in), target :: DVT

        type(DelaunayTrianglePointerList), pointer :: icdDT
        integer i

        write(*, "('Vertex (ID -', I3, ')')") DVT%id
        write(*, "('  Incident triangle ID: ')", advance="no")
        icdDT => DVT%icdDTHead
        do i = 1, DVT%numIcdDT
            write(*, "(I3)", advance="no") icdDT%ptr%id
            icdDT => icdDT%next
        end do
        write(*, *)
    
    end subroutine PrintVertexTopology

    subroutine PrintVoronoiCell(VC)
        type(VoronoiCell), intent(in) :: VC
        integer i

        write(*, "('Voronoi cell (ID -', I3, ')')") VC%id
        write(*, "('  Vertex coordinate (in degree): ')")
        do i = 1, VC%numVVT
            write(*, "(2F10.3)") VC%lonVVT(i)*Rad2Deg, VC%latVVT(i)*Rad2Deg
        end do
    
    end subroutine PrintVoronoiCell
    
    subroutine PrintTriangleContainedVertex(DT)
        type(DelaunayTriangle), intent(in) :: DT

        type(DelaunayVertexPointerList), pointer :: incDVT
        integer i

        write(*, "('Triangle (ID -', I3, ')')") DT%id
        write(*, "('  Contained vertex ID: ')")
        incDVT => DT%incDVTHead
        do i = 1, DT%numIncDVT
            write(*, "(I3)", advance="no") incDVT%ptr%id
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
            write(*, "(I3)", advance="no") obsDT%ptr%id
        end do
        write(*, *)

    end subroutine PrintObsoleteTriangle
 
    subroutine PrintTemporalTriangle
        type(DelaunayTrianglePointerList), pointer :: tmpDT
        integer i

        write(*, "('Temporal triangles ID:')")
        tmpDT => tmpDTHead
        do i = 1, numTmpDT
            write(*, "(I3)", advance="no") tmpDT%ptr%id
        end do
        write(*, *)

    end subroutine PrintTemporalTriangle

    subroutine DelaunayAndVoronoi_Output(filePath)
        character(*), intent(in) :: filePath

        real(4), allocatable :: lon(:), lat(:)
        real(4), allocatable :: lonCirc(:), latCirc(:), radius(:)
        integer, allocatable :: triangle(:,:)
        real(4), allocatable :: lonVVT(:,:), latVVT(:,:)
        integer, allocatable :: numVVT(:)
        type(DelaunayVertex), pointer :: DVT
        type(DelaunayTriangle), pointer :: DT
        type(VoronoiCell), pointer :: VC
        integer, save :: tag = 0
        integer i, j

        call MsgManager_RecordSpeaker("DelaunayAndVoronoi_Output")

        allocate(lon(numDVT))
        allocate(lat(numDVT))

        DVT => DVTHead
        do i = 1, numDVT
            lon(i) = real(DVT%SMP%lon*Rad2Deg)
            lat(i) = real(DVT%SMP%lat*Rad2Deg)
            DVT => DVT%next
        end do

        allocate(lonCirc(numDT))
        allocate(latCirc(numDT))
        allocate(radius(numDT))

        allocate(triangle(3,numDT))

        DT => DTHead
        do i = 1, numDT
            do j = 1, 3
                triangle(j,i) = DT%DVT(j)%ptr%id
            end do
            lonCirc(i) = real(DT%lonCirc*Rad2Deg)
            latCirc(i) = real(DT%latCirc*Rad2Deg)
            radius(i) = real(DT%radius*Rad2Deg)
            DT => DT%next
        end do

        allocate(lonVVT(maxNumVVT,numVC))
        allocate(latVVT(maxNumVVT,numVC))
        allocate(numVVT(numVC))

        VC => VCHead
        do j = 1, numVC
            do i = 1, VC%numVVT
                lonVVT(i,j) = real(VC%lonVVT(i)*Rad2Deg)
                latVVT(i,j) = real(VC%latVVT(i)*Rad2Deg)
            end do
            numVVT(j) = VC%numVVT
            VC => VC%next
        end do

        tag = tag+1
        call NFWrap_CreateIrregular(trim(int2str(tag))//"_"//filePath, &
            timeVariant=.false., card=fcard)
        call NFWrap_NewDim(fcard, "numPoint", numDVT)
        call NFWrap_NewDim(fcard, "numTriangle", numDT)
        call NFWrap_NewDim(fcard, "numTriangleVertex", 3)
        call NFWrap_NewDim(fcard, "maxNumVoronoiCellVertex", maxNumVVT)
        call NFWrap_New1DVar(fcard, &
            varName="lonPoint", dataType="float", &
            dimName="numPoint", &
            longName="Point longitude", &
            unitName="degree_east", timeVariant=.false.)
        call NFWrap_New1DVar(fcard, &
            varName="latPoint", dataType="float", &
            dimName="numPoint", &
            longName="Point latitude", &
            unitName="degree_north", timeVariant=.false.)
        call NFWrap_New1DVar(fcard, &
            varName="lonCirc", dataType="float", &
            dimName="numTriangle", &
            longName="Circumcenter longitude", &
            unitName="degree_east", timeVariant=.false.)
        call NFWrap_New1DVar(fcard, &
            varName="latCirc", dataType="float", &
            dimName="numTriangle", &
            longName="Circumcenter latitude", &
            unitName="degree_north", timeVariant=.false.)
        call NFWrap_New1DVar(fcard, &
            varName="radius", dataType="float", &
            dimName="numTriangle", &
            longName="Circumradius", &
            unitName="degree", timeVariant=.false.)
        call NFWrap_New2DVar(fcard, &
            varName="triangle", dataType="integer", &
            dimName1="numTriangleVertex", dimName2="numTriangle", &
            longName="Triangle vertex index", &
            unitName="", &
            timeVariant=.false.)
        call NFWrap_New1DVar(fcard, &
            varName="numVVT", dataType="integer", &
            dimName="numPoint", &
            longName="Real Voronoi cell vertex number", &
            unitName="", timeVariant=.false.)
        call NFWrap_New2DVar(fcard, &
            varName="lonVVT", dataType="float", &
            dimName1="maxNumVoronoiCellVertex", dimName2="numPoint", &
            longName="Voronoi cell vertex longitude", &
            unitName="degree_east", timeVariant=.false.)
        call NFWrap_New2DVar(fcard, &
            varName="latVVT", dataType="float", &
            dimName1="maxNumVoronoiCellVertex", dimName2="numPoint", &
            longName="Voronoi cell vertex latitude", &
            unitName="degree_north", timeVariant=.false.)
        call NFWrap_Output1DVar(fcard, "lonPoint", lon)
        call NFWrap_Output1DVar(fcard, "latPoint", lat)
        call NFWrap_Output1DVar(fcard, "lonCirc", lonCirc)
        call NFWrap_Output1DVar(fcard, "latCirc", latCirc)
        call NFWrap_Output1DVar(fcard, "radius", radius)
        call NFWrap_Output2DVar(fcard, "triangle", triangle)
        call NFWrap_Output1DVar(fcard, "numVVT", numVVT)
        call NFWrap_Output2DVar(fcard, "lonVVT", lonVVT)
        call NFWrap_Output2DVar(fcard, "latVVT", latVVT)
        call NFWrap_Close(fcard)

        call MsgManager_Speak(Notice, "File "//trim(filePath)//" is generated.")
        call MsgManager_DeleteSpeaker

    end subroutine DelaunayAndVoronoi_Output
 
end module DelaunayAndVoronoi
 
