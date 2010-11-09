module NFWrap

    use MsgManager
    use RunManager
    use netcdf

    implicit none

    integer, parameter :: maxNumDim = 100
    integer, parameter :: maxNumVar = 100

    type Dimension
        integer id
        character(50) name
        integer size
    end type Dimension

    type Variable
        integer id
        character(50) name
        integer numDim
        integer, allocatable :: dimSize(:)
        integer, allocatable :: dimId(:)
        logical :: timeVariant = .false.
    end type Variable

    type FileCard
        integer id
        ! dimensions
        integer timeDimId ! short hand
        integer :: numDim = 0
        type(Dimension) dim(maxNumDim)
        ! variables
        integer timeVarId ! short hand
        integer :: numVar = 0
        type(Variable) var(maxNumVar)
        ! others
        integer :: timeStep = 0
    contains
        procedure :: addDim => FileCard_addDim
        procedure :: getDim => FileCard_getDim
        procedure :: addVar => FileCard_addVar
        procedure :: getVar => FileCard_getVar
    end type FileCard

    interface NFWrap_Output1DVar
        module procedure NFWrap_Output1DIntVar
        module procedure NFWrap_Output1DFloatVar
        module procedure NFWrap_Output1DDoubleVar
    end interface NFWrap_Output1DVar

    interface NFWrap_Output2DVar
        module procedure NFWrap_Output2DIntVar
        module procedure NFWrap_Output2DFloatVar
        module procedure NFWrap_Output2DDoubleVar
    end interface NFWrap_Output2DVar

contains

    ! ************************************************************************ !
    !                         File access interface                            !
    ! ************************************************************************ !

    subroutine NFWrap_CreateSpherical2D(filePath, numLon, numLat, card)
        character(*), intent(in) :: filePath
        integer, intent(in) :: numLon, numLat
        type(FileCard), intent(out) :: card

        integer ierr
        type(Dimension), pointer :: dim
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_CreateSpherical2D")

        ierr = nf90_create(filePath, nf90_clobber, card%id)
        call NFWrap_HandleError(ierr)
        
        ! =============================================================================== !
        ! time dimension
        call card%addDim("time", 0)
        call card%getDim("time", dim)
        ierr = nf90_def_dim(card%id, "time", nf90_unlimited, dim%id)
        call NFWrap_HandleError(ierr)
        call card%addVar("time", 0, timeVariant=.true.)
        call card%getVar("time", var)
        ierr = nf90_def_var(card%id, "time", nf90_double, [dim%id], var%id)
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "long_name", "seconds since 0")
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "units", "seconds")
        call NFWrap_HandleError(ierr)
        card%timeDimId = dim%id ! short hand
        card%timeVarId = var%id ! short hand

        ! =============================================================================== !
        ! longitude dimension
        call card%addDim("lon", numLon)
        call card%getDim("lon", dim)
        ierr = nf90_def_dim(card%id, "lon", numLon, dim%id)
        call NFWrap_HandleError(ierr)
        call card%addVar("lon", 1, [numLon], timeVariant=.false.)
        call card%getVar("lon", var)
        ierr = nf90_def_var(card%id, "lon", nf90_double, [dim%id], var%id)
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "long_name", "longitude")
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "units", "degree_E")
        call NFWrap_HandleError(ierr)

        ! =============================================================================== !
        ! latitude dimension
        call card%addDim("lat", numLat)
        call card%getDim("lat", dim)
        ierr = nf90_def_dim(card%id, "lat", numLat, dim%id)
        call NFWrap_HandleError(ierr)
        call card%addVar("lat", 1, [numLat], timeVariant=.false.)
        call card%getVar("lat", var)
        ierr = nf90_def_var(card%id, "lat", nf90_double, [dim%id], var%id)
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "long_name", "latitude")
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "units", "degree_N")
        call NFWrap_HandleError(ierr)
        
        ierr = nf90_enddef(card%id)
        call NFWrap_HandleError(ierr)

        call MsgManager_Speak(Notice, "File """//trim(filePath)//""" is created.")
        call MsgManager_DeleteSpeaker

    end subroutine NFWrap_CreateSpherical2D

    subroutine NFWrap_CreateIrregular(filePath, timeVariant, card)
        character(*), intent(in) :: filePath
        logical, intent(in) :: timeVariant
        type(FileCard), intent(out) :: card
    
        integer ierr
        type(Dimension), pointer :: dim
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_CreateIrregular")

        ierr = nf90_create(filePath, nf90_clobber, card%id)
        call NFWrap_HandleError(ierr)

        ! =============================================================================== !
        ! time dimension
        if (timeVariant) then
            call card%addDim("time", 0)
            call card%getDim("time", dim)
            ierr = nf90_def_dim(card%id, "time", nf90_unlimited, dim%id)
            call NFWrap_HandleError(ierr)
            call card%addVar("time", 0, timeVariant=.true.)
            call card%getVar("time", var)
            ierr = nf90_def_var(card%id, "time", nf90_double, [dim%id], var%id)
            call NFWrap_HandleError(ierr)
            ierr = nf90_put_att(card%id, var%id, "long_name", "seconds since 0")
            call NFWrap_HandleError(ierr)
            ierr = nf90_put_att(card%id, var%id, "units", "seconds")
            call NFWrap_HandleError(ierr)
            card%timeDimId = dim%id ! short hand
            card%timeVarId = var%id ! short hand
        end if

        ierr = nf90_enddef(card%id)
        call NFWrap_HandleError(ierr)

        call MsgManager_Speak(Notice, "File """//trim(filePath)//""" is created.")
        call MsgManager_DeleteSpeaker

    end subroutine NFWrap_CreateIrregular

    subroutine NFWrap_OpenForRead(filePath, card)
        character(*), intent(in) :: filePath
        type(FileCard), intent(out) :: card
    
        integer i, j, k, ierr

        call MsgManager_RecordSpeaker("NFWrap_OpenForRead")

        ierr = nf90_open(filePath, nf90_nowrite, card%id)
        call NFWrap_HandleError(ierr)

        ierr = nf90_inquire(card%id, card%numDim, card%numVar)
        call NFWrap_HandleError(ierr)

        ierr = nf90_inquire(card%id, unlimitedDimID=card%timeDimId)
        call NFWrap_HandleError(ierr)

        ! inquire dimensions
        do i = 1, card%numDim
            card%dim(i)%id = i
            ierr = nf90_inquire_dimension(card%id, card%dim(i)%id, &
                name=card%dim(i)%name, len=card%dim(i)%size)
            call NFWrap_HandleError(ierr)
        end do

        ! inquire variables
        do i = 1, card%numVar
            card%var(i)%id = i
            ierr = nf90_inquire_variable(card%id, card%var(i)%id, &
                name=card%var(i)%name, ndims=card%var(i)%numDim)
            call NFWrap_HandleError(ierr)
            allocate(card%var(i)%dimSize(card%var(i)%numDim))
            allocate(card%var(i)%dimId(card%var(i)%numDim))
            ierr = nf90_inquire_variable(card%id, card%var(i)%id, &
                dimids=card%var(i)%dimId)
            call NFWrap_HandleError(ierr)
            do j = 1, card%var(i)%numDim
                k = card%var(i)%dimId(j)
                card%var(i)%dimSize(j) = card%dim(k)%size
                if (card%var(i)%dimId(j) == card%timeDimId) then
                    card%var(i)%timeVariant = .true.
                end if
            end do
        end do

        write(*, "('**********************************************')")
        write(*, "('NetCDF file ', A)") trim(filePath)
        write(*, "('  Dimensions:')")
        write(*, "('     ', A10, ' ', A10)") "Name", "Size"
        do i = 1, card%numDim
            write(*, "('  ', I1, '. ', A10, ' ', I10)") i, &
                trim(card%dim(i)%name), card%dim(i)%size
        end do
        write(*, "('  Variables:')")
        write(*, "('     ', A10, ' ', A10)") "Name", "Dimension"
        do i = 1, card%numVar
            write(*, "('  ', I1, '. ', A10)", advance="no") i, &
                trim(card%var(i)%name)
            do j = 1, card%var(i)%numDim
                k = card%var(i)%dimId(j)
                write(*, "(' ', A6)", advance="no") &
                    trim(card%dim(k)%name)
            end do
            write(*, *)
        end do
        write(*, "('**********************************************')")

        call MsgManager_DeleteSpeaker

    end subroutine NFWrap_OpenForRead

    subroutine NFWrap_Close(card)
        type(FileCard), intent(in) :: card

        integer ierr

        ierr = nf90_close(card%id)

    end subroutine NFWrap_Close
    
    ! ************************************************************************ !
    !                           Dimension interface                            !
    ! ************************************************************************ !

    subroutine NFWrap_NewDim(card, dimName, dimSize)
        type(FileCard), intent(inout) :: card
        character(*), intent(in) :: dimName
        integer, intent(in) :: dimSize

        type(Dimension), pointer :: dim
        integer ierr

        call MsgManager_RecordSpeaker("NFWrap_NewDim")

        call card%addDim(dimName, dimSize)
        call card%getDim(dimName, dim)
        
        ierr = nf90_redef(card%id)
        call NFWrap_HandleError(ierr)
        ierr = nf90_def_dim(card%id, dimName, dimSize, dim%id)
        call NFWrap_HandleError(ierr)

        ierr = nf90_enddef(card%id)
        call NFWrap_HandleError(ierr)

        call MsgManager_DeleteSpeaker
    
    end subroutine NFWrap_NewDim

    subroutine NFWrap_GetDimSize(card, dimName, dimSize)
        type(FileCard), intent(inout) :: card
        character(*), intent(in) :: dimName
        integer, intent(out) :: dimSize

        type(Dimension), pointer :: dim
        integer ierr
 
        call card%getDim(dimName, dim)
        dimSize = dim%size

    end subroutine NFWrap_GetDimSize
 
    ! ************************************************************************ !
    !                           Variable interface                             !
    ! ************************************************************************ !

    subroutine NFWrap_NewScalar(card, varName, longName, unitName)
        type(FileCard), intent(inout), target :: card
        character(*), intent(in) :: varName, longName, unitName
    
        integer ierr
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_NewScalar")
    
        call card%addVar(varName, 0, timeVariant=.true.)
        call card%getVar(varName, var)

        ierr = nf90_redef(card%id)
        call NFWrap_HandleError(ierr)
        ierr = nf90_def_var(card%id, varName, nf90_double, &
            [card%timeDimId], var%id)
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "long_name", longName)
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "units", unitName)
        call NFWrap_HandleError(ierr)

        ierr = nf90_enddef(card%id)
        call NFWrap_HandleError(ierr)

        call MsgManager_DeleteSpeaker
    
    end subroutine NFWrap_NewScalar
    
    subroutine NFWrap_New1DVar(card, varName, dataType, dimName, &
        longName, unitName, timeVariant)
        type(FileCard), intent(inout), target :: card
        character(*), intent(in) :: varName, dimName, longName, unitName
        character(*), intent(in), optional :: dataType
        logical, intent(in) :: timeVariant

        integer ierr, xtype
        type(Dimension), pointer :: dim
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_New1DVar")

        call card%getDim(dimName, dim)

        call card%addVar(varName, 1, [dim%size], timeVariant)
        call card%getVar(varName, var)

        if (present(dataType)) then
            select case (dataType)
            case ("float")
                xtype = nf90_float
            case ("double")
                xtype = nf90_double
            case ("integer")
                xtype = nf90_int
            end select
        else
            xtype = nf90_double
        end if

        ierr = nf90_redef(card%id)
        call NFWrap_HandleError(ierr)
        if (timeVariant) then
            ierr = nf90_def_var(card%id, varName, xtype, &
                [dim%id,card%timeDimId], var%id)
        else
            ierr = nf90_def_var(card%id, varName, xtype, &
                [dim%id], var%id)
        end if
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "long_name", longName)
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "units", unitName)
        call NFWrap_HandleError(ierr)

        ierr = nf90_enddef(card%id)
        call NFWrap_HandleError(ierr)

        call MsgManager_DeleteSpeaker
    
    end subroutine NFWrap_New1DVar

    subroutine NFWrap_New2DVar(card, varName, dataType, dimName1, dimName2, &
        longName, unitName, timeVariant)
        type(FileCard), intent(inout), target :: card
        character(*), intent(in) :: varName, dimName1, dimName2
        character(*), intent(in), optional :: dataType
        character(*), intent(in) :: longName, unitName
        logical, intent(in) :: timeVariant

        integer ierr, xtype
        type(Dimension), pointer :: dim1, dim2
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_New2DVar")

        call card%getDim(dimName1, dim1)
        call card%getDim(dimName2, dim2)

        call card%addVar(varName, 2, [dim1%size,dim2%size], timeVariant)
        call card%getVar(varName, var)

        if (present(dataType)) then
            select case (dataType)
            case ("float")
                xtype = nf90_float
            case ("double")
                xtype = nf90_double
            case ("integer")
                xtype = nf90_int
            end select
        else
            xtype = nf90_double
        end if

        ierr = nf90_redef(card%id)
        call NFWrap_HandleError(ierr)
        if (timeVariant) then
            ierr = nf90_def_var(card%id, varName, xtype, &
                [dim1%id,dim2%id,card%timeDimId], var%id)
        else
            ierr = nf90_def_var(card%id, varName, xtype, &
                [dim1%id,dim2%id], var%id)
        end if
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "long_name", longName)
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "units", unitName)
        call NFWrap_HandleError(ierr)

        ierr = nf90_enddef(card%id)
        call NFWrap_HandleError(ierr)

        call MsgManager_DeleteSpeaker

    end subroutine NFWrap_New2DVar

    subroutine NFWrap_Advance(card, timeStep, time)
        type(FileCard), intent(inout) :: card
        integer, intent(in) :: timeStep
        real(8), intent(in) :: time

        integer ierr

        call MsgManager_RecordSpeaker("NFWrap_Advance")
    
        card%timeStep = timeStep+1

        ierr = nf90_put_var(card%id, card%timeVarId, time, [card%timeStep])
    
        call MsgManager_DeleteSpeaker
    
    end subroutine NFWrap_Advance

    subroutine NFWrap_OutputScalar(card, varName, varValue)
        type(FileCard), intent(inout), target :: card
        character(*), intent(in) :: varName
        real(8), intent(in) :: varValue

        integer dimSize(1), ierr
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_OutputScalar")

        call card%getVar(varName, var)

        ierr = nf90_put_var(card%id, var%id, varValue, [card%timeStep])
        call NFWrap_HandleError(ierr)

        call MsgManager_DeleteSpeaker
    
    end subroutine NFWrap_OutputScalar

    subroutine NFWrap_Output1DIntVar(card, varName, varValue)
        type(FileCard), intent(inout), target :: card
        character(*), intent(in) :: varName
        integer, intent(in) :: varValue(:)

        integer dimSize(1), ierr
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_Output1DVar")

        dimSize = shape(varValue)

        call card%getVar(varName, var)

        if (dimSize(1) /= var%dimSize(1)) then
            call MsgManager_Speak(Error, &
                "Output variable """//trim(varName)// &
                """: Dimension not match! ("//trim(int2str(dimSize(1)))// &
                " /= "//trim(int2str(var%dimSize(1)))//")")
            call RunManager_EndRun
        end if

        if (var%timeVariant) then
            ierr = nf90_put_var(card%id, var%id, varValue, [1,card%timeStep])
        else
            ierr = nf90_put_var(card%id, var%id, varValue)
        end if
        call NFWrap_HandleError(ierr)

#if (defined DEBUG)
        ierr = nf90_sync(card%id)
        call NFWrap_HandleError(ierr)
#endif

        call MsgManager_DeleteSpeaker
    
    end subroutine NFWrap_Output1DIntVar

    subroutine NFWrap_Output1DFloatVar(card, varName, varValue)
        type(FileCard), intent(inout), target :: card
        character(*), intent(in) :: varName
        real(4), intent(in) :: varValue(:)

        integer dimSize(1), ierr
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_Output1DVar")

        dimSize = shape(varValue)

        call card%getVar(varName, var)

        if (dimSize(1) /= var%dimSize(1)) then
            call MsgManager_Speak(Error, &
                "Output variable """//trim(varName)// &
                """: Dimension not match! ("//trim(int2str(dimSize(1)))// &
                " /= "//trim(int2str(var%dimSize(1)))//")")
            call RunManager_EndRun
        end if

        if (var%timeVariant) then
            ierr = nf90_put_var(card%id, var%id, varValue, [1,card%timeStep])
        else
            ierr = nf90_put_var(card%id, var%id, varValue)
        end if
        call NFWrap_HandleError(ierr)

#if (defined DEBUG)
        ierr = nf90_sync(card%id)
        call NFWrap_HandleError(ierr)
#endif

        call MsgManager_DeleteSpeaker
    
    end subroutine NFWrap_Output1DFloatVar

    subroutine NFWrap_Output1DDoubleVar(card, varName, varValue)
        type(FileCard), intent(inout), target :: card
        character(*), intent(in) :: varName
        real(8), intent(in) :: varValue(:)

        integer dimSize(1), ierr
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_Output1DVar")

        dimSize = shape(varValue)

        call card%getVar(varName, var)

        if (dimSize(1) /= var%dimSize(1)) then
            call MsgManager_Speak(Error, &
                "Output variable """//trim(varName)// &
                """: Dimension not match! ("//trim(int2str(dimSize(1)))// &
                " /= "//trim(int2str(var%dimSize(1)))//")")
            call RunManager_EndRun
        end if

        if (var%timeVariant) then
            ierr = nf90_put_var(card%id, var%id, varValue, [1,card%timeStep])
        else
            ierr = nf90_put_var(card%id, var%id, varValue)
        end if
        call NFWrap_HandleError(ierr)

#if (defined DEBUG)
        ierr = nf90_sync(card%id)
        call NFWrap_HandleError(ierr)
#endif

        call MsgManager_DeleteSpeaker
    
    end subroutine NFWrap_Output1DDoubleVar
    
    subroutine NFWrap_Output2DIntVar(card, varName, varValue)
        type(FileCard), intent(inout), target :: card
        character(*), intent(in) :: varName
        integer, intent(in) :: varValue(:,:)

        integer dimSize(2), ierr
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_Output2DVar")
    
        dimSize = shape(varValue)

        call card%getVar(varName, var)

        if (dimSize(1) /= var%dimSize(1)) then
            call MsgManager_Speak(Error, &
                "Output variable """//trim(varName)//""": First dimension not match ("//&
                trim(int2str(dimSize(1)))//" /= "//trim(int2str(var%dimSize(1)))//")!")
            call RunManager_EndRun
        else if (dimSize(2) /= var%dimSize(2)) then
            call MsgManager_Speak(Error, &
                "Output variable """//trim(varName)//""": Second dimension not match ("//&
                trim(int2str(dimSize(2)))//" /= "//trim(int2str(var%dimSize(2)))//")!")
            call RunManager_EndRun
        end if

        if (var%timeVariant) then
            ierr = nf90_put_var(card%id, var%id, varValue, [1,1,card%timeStep])
        else
            ierr = nf90_put_var(card%id, var%id, varValue)
        end if
        call NFWrap_HandleError(ierr)
 
#if (defined DEBUG)
        ierr = nf90_sync(card%id)
        call NFWrap_HandleError(ierr)
#endif

        call MsgManager_DeleteSpeaker
    
    end subroutine NFWrap_Output2DIntVar

    subroutine NFWrap_Output2DFloatVar(card, varName, varValue)
        type(FileCard), intent(inout), target :: card
        character(*), intent(in) :: varName
        real(4), intent(in) :: varValue(:,:)

        integer dimSize(2), ierr
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_Output2DVar")
    
        dimSize = shape(varValue)

        call card%getVar(varName, var)

        if (dimSize(1) /= var%dimSize(1)) then
            call MsgManager_Speak(Error, &
                "Output variable """//trim(varName)//""": First dimension not match ("//&
                trim(int2str(dimSize(1)))//" /= "//trim(int2str(var%dimSize(1)))//")!")
            call RunManager_EndRun
        else if (dimSize(2) /= var%dimSize(2)) then
            call MsgManager_Speak(Error, &
                "Output variable """//trim(varName)//""": Second dimension not match ("//&
                trim(int2str(dimSize(2)))//" /= "//trim(int2str(var%dimSize(2)))//")!")
            call RunManager_EndRun
        end if

        if (var%timeVariant) then
            ierr = nf90_put_var(card%id, var%id, varValue, [1,1,card%timeStep])
        else
            ierr = nf90_put_var(card%id, var%id, varValue)
        end if
        call NFWrap_HandleError(ierr)
 
#if (defined DEBUG)
        ierr = nf90_sync(card%id)
        call NFWrap_HandleError(ierr)
#endif

        call MsgManager_DeleteSpeaker
    
    end subroutine NFWrap_Output2DFloatVar

    subroutine NFWrap_Output2DDoubleVar(card, varName, varValue)
        type(FileCard), intent(inout), target :: card
        character(*), intent(in) :: varName
        real(8), intent(in) :: varValue(:,:)

        integer dimSize(2), ierr
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_Output2DVar")
    
        dimSize = shape(varValue)

        call card%getVar(varName, var)

        if (dimSize(1) /= var%dimSize(1)) then
            call MsgManager_Speak(Error, &
                "Output variable """//trim(varName)//""": First dimension not match ("//&
                trim(int2str(dimSize(1)))//" /= "//trim(int2str(var%dimSize(1)))//")!")
            call RunManager_EndRun
        else if (dimSize(2) /= var%dimSize(2)) then
            call MsgManager_Speak(Error, &
                "Output variable """//trim(varName)//""": Second dimension not match ("//&
                trim(int2str(dimSize(2)))//" /= "//trim(int2str(var%dimSize(2)))//")!")
            call RunManager_EndRun
        end if

        if (var%timeVariant) then
            ierr = nf90_put_var(card%id, var%id, varValue, [1,1,card%timeStep])
        else
            ierr = nf90_put_var(card%id, var%id, varValue)
        end if
        call NFWrap_HandleError(ierr)
 
#if (defined DEBUG)
        ierr = nf90_sync(card%id)
        call NFWrap_HandleError(ierr)
#endif

        call MsgManager_DeleteSpeaker
    
    end subroutine NFWrap_Output2DDoubleVar

    subroutine NFWrap_Input1DVar(card, varName, varValue, timeStep)
        type(FileCard), intent(in), target :: card
        character(*), intent(in) :: varName
        real(8), intent(out) :: varValue(:)
        integer, intent(in), optional :: timeStep

        integer dimSize(1), ierr
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_Input1DVar")

        dimSize = shape(varValue)

        call card%getVar(varName, var)

        if (dimSize(1) /= var%dimSize(1)) then
            call MsgManager_Speak(Error, &
                "Output variable """//trim(varName)//""": Dimension not match!")
            call RunManager_EndRun
        end if

        if (var%timeVariant .and. present(timeStep)) then
            ierr = nf90_get_var(card%id, var%id, varValue, [1,timeStep], [dimSize(1),1])
        else
            ierr = nf90_get_var(card%id, var%id, varValue)
        end if

        call MsgManager_DeleteSpeaker

    end subroutine NFWrap_Input1DVar

    subroutine NFWrap_Input2DVar(card, varName, varValue, timeStep)
        type(FileCard), intent(in), target :: card
        character(*), intent(in) :: varName
        real(8), intent(out) :: varValue(:,:)
        integer, intent(in), optional :: timeStep

        integer dimSize(2), ierr
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_Input1DVar")
 
        dimSize = shape(varValue)

        call card%getVar(varName, var)

        if (dimSize(1) /= var%dimSize(1) .or. dimSize(2) /= var%dimSize(2)) then
            call MsgManager_Speak(Error, &
                "Output variable """//trim(varName)//""": Dimension not match!")
            call RunManager_EndRun
        end if

        if (var%timeVariant .and. present(timeStep)) then
            ierr = nf90_get_var(card%id, var%id, varValue, [1,1,timeStep], [dimSize(1),dimSize(2),1])
        else
            ierr = nf90_get_var(card%id, var%id, varValue)
        end if

        call MsgManager_DeleteSpeaker
   
    end subroutine NFWrap_Input2DVar

    ! ************************************************************************ !
    !                           Utilities                                      !
    ! ************************************************************************ !

    subroutine NFWrap_HandleError(ierr)
        integer, intent(in) :: ierr

        if (ierr /= nf90_noerr) then
            call MsgManager_Speak(Error, &
                "NetCDF error """//trim(nf90_strerror(ierr))//""".")
            stop
        end if

    end subroutine NFWrap_HandleError

    ! ************************************************************************ !
    !                     Type-bounded procedures                              !
    ! ************************************************************************ !

    subroutine FileCard_addDim(card, name, size)
        class(FileCard), intent(inout) :: card
        character(*), intent(in) :: name
        integer, intent(in) :: size

        integer i

        call MsgManager_RecordSpeaker("FileCard_addDim")

        card%numDim = card%numDim+1
        if (card%numDim > maxNumDim) then
            call MsgManager_Speak(Error, &
                "Exceed maximum number of dimensions!")
            call RunManager_EndRun
        end if
        i = card%numDim
        card%dim(i)%name = name
        card%dim(i)%size = size

        call MsgManager_DeleteSpeaker

    end subroutine FileCard_addDim

    subroutine FileCard_getDim(card, dimName, dimPtr)
        class(FileCard), intent(inout), target :: card
        character(*), intent(in) :: dimName
        type(Dimension), intent(out), pointer :: dimPtr

        integer i

        call MsgManager_RecordSpeaker("FileCard_getDim")

        do i = 1, card%numDim
            if (dimName == card%dim(i)%name) then
                dimPtr => card%dim(i)
                call MsgManager_DeleteSpeaker
                return
            end if
        end do

        call MsgManager_Speak(Error, &
            "Unknown dimension """//trim(dimName)//""".")
        call RunManager_EndRun
    
    end subroutine FileCard_getDim

    subroutine FileCard_addVar(card, name, numDim, dimSize, timeVariant)
        class(FileCard), intent(inout) :: card 
        character(*), intent(in) :: name
        integer, intent(in) :: numDim
        integer, intent(in), optional :: dimSize(numDim)
        logical, intent(in) :: timeVariant

        integer i

        call MsgManager_RecordSpeaker("FileCard_addVar")

        card%numVar = card%numVar+1
        if (card%numVar > maxNumVar) then
            call MsgManager_Speak(Error, &
                "Exceed maximum number of variables!")
            call RunManager_EndRun
        end if
        i = card%numVar
        card%var(i)%name = name
        card%var(i)%numDim = numDim
        if (numDim /= 0 .and. present(dimSize)) then
            allocate(card%var(i)%dimSize(numDim))
            card%var(i)%dimSize = dimSize
        end if
        card%var(i)%timeVariant = timeVariant

        call MsgManager_DeleteSpeaker

    end subroutine FileCard_addVar
    
    subroutine FileCard_getVar(card, varName, varPtr)
        class(FileCard), intent(in), target :: card
        character(*), intent(in) :: varName
        type(Variable), intent(out), pointer :: varPtr

        integer i

        call MsgManager_RecordSpeaker("FileCard_getVar")

        do i = 1, card%numVar
            if (varName == card%var(i)%name) then
                varPtr => card%var(i)
                call MsgManager_DeleteSpeaker
                return
            end if
        end do

        call MsgManager_Speak(Error, &
            "Unknown variable """//trim(varName)//""".")
        call RunManager_EndRun

    end subroutine FileCard_getVar
    
end module NFWrap
