module NFWrap

    use MsgManager
    use RunManager
    use netcdf

    implicit none

    integer, parameter :: maxNumDim = 100
    integer, parameter :: maxNumVar = 100
    integer, parameter :: maxNumAtt = 100

    type Dimension
        integer id
        character(50) name
        integer size
    end type Dimension

    type Attribute
        character(50) name
        character(50) type
        character(50) valueStr
        integer valueInt
        real(4) valueReal
        real(8) valueDouble
    end type Attribute

    type Variable
        integer id
        character(50) name
        integer numDim
        integer, allocatable :: dimSize(:)
        integer, allocatable :: dimId(:)
        integer numAtt
        type(Attribute), allocatable :: att(:)
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

    interface NFWrap_Output3DVar
        module procedure NFWrap_Output3DDoubleVar
    end interface NFWrap_Output3DVar

    interface NFWrap_Input1DVar
        module procedure NFWrap_Input1DFloatVar
        module procedure NFWrap_Input1DDoubleVar
    end interface NFWrap_Input1DVar

    interface NFWrap_Input3DVar
        module procedure NFWrap_Input3DDoubleVar
    end interface NFWrap_Input3DVar

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
        ierr = nf90_put_att(card%id, var%id, "units", "degree_east")
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
        ierr = nf90_put_att(card%id, var%id, "units", "degree_north")
        call NFWrap_HandleError(ierr)
        
        ierr = nf90_enddef(card%id)
        call NFWrap_HandleError(ierr)

        call MsgManager_Speak(Notice, "File """//trim(filePath)//""" is created.")
        call MsgManager_DeleteSpeaker

    end subroutine NFWrap_CreateSpherical2D

    subroutine NFWrap_CreateSpherical3D(filePath, numLon, numLat, numLev, card)
        character(*), intent(in) :: filePath
        integer, intent(in) :: numLon, numLat, numLev
        type(FileCard), intent(out) :: card

        integer ierr
        type(Dimension), pointer :: dim
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_CreateSpherical3D")

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
        ierr = nf90_put_att(card%id, var%id, "units", "degree_east")
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
        ierr = nf90_put_att(card%id, var%id, "units", "degree_north")
        call NFWrap_HandleError(ierr)

        ! =============================================================================== !
        ! level dimension
        call card%addDim("lev", numLev)
        call card%getDim("lev", dim)
        ierr = nf90_def_dim(card%id, "lev", numLev, dim%id)
        call NFWrap_HandleError(ierr)
        call card%addVar("lev", 1, [numLev], timeVariant=.false.)
        call card%getVar("lev", var)
        ierr = nf90_def_var(card%id, "lev", nf90_double, [dim%id], var%id)
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "long_name", "level")
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "units", "hPa")
        call NFWrap_HandleError(ierr)

        ierr = nf90_enddef(card%id)
        call NFWrap_HandleError(ierr)

        call MsgManager_Speak(Notice, "File """//trim(filePath)//""" is created.")
        call MsgManager_DeleteSpeaker

    end subroutine NFWrap_CreateSpherical3D

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
            ierr = nf90_def_var(card%id, "time", nf90_float, [dim%id], var%id)
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
            ! dimension
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
            ! attribute
            ierr = nf90_inquire_variable(card%id, card%var(i)%id, &
                nAtts=card%var(i)%numAtt)
            call NFWrap_HandleError(ierr)
            allocate(card%var(i)%att(card%var(i)%numAtt))
            do j = 1, card%var(i)%numAtt
                ierr = nf90_inq_attname(card%id, card%var(i)%id, j, &
                    card%var(i)%att(j)%name)
                call NFWrap_HandleError(ierr)
                ierr = nf90_inquire_attribute(card%id, card%var(i)%id, &
                    card%var(i)%att(j)%name, k)
                call NFWrap_HandleError(ierr)
                select case (k)
                case (nf90_char)
                    card%var(i)%att(j)%type = "string"
                    ierr = nf90_get_att(card%id, card%var(i)%id, &
                        card%var(i)%att(j)%name, card%var(i)%att(j)%valueStr)
                case (nf90_int)
                    card%var(i)%att(j)%type = "integer"
                    ierr = nf90_get_att(card%id, card%var(i)%id, &
                        card%var(i)%att(j)%name, card%var(i)%att(j)%valueInt)
                case (nf90_float)
                    card%var(i)%att(j)%type = "float"
                    ierr = nf90_get_att(card%id, card%var(i)%id, &
                        card%var(i)%att(j)%name, card%var(i)%att(j)%valueReal)
                case (nf90_double)
                    card%var(i)%att(j)%type = "double"
                    ierr = nf90_get_att(card%id, card%var(i)%id, &
                        card%var(i)%att(j)%name, card%var(i)%att(j)%valueDouble)
                end select
            end do
        end do

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

    subroutine NFWrap_New3DVar(card, varName, dataType, dimName1, dimName2, &
        dimName3, longName, unitName, timeVariant)
        type(FileCard), intent(inout), target :: card
        character(*), intent(in) :: varName, dimName1, dimName2, dimName3
        character(*), intent(in), optional :: dataType
        character(*), intent(in) :: longName, unitName
        logical, intent(in) :: timeVariant

        integer ierr, xtype
        type(Dimension), pointer :: dim1, dim2, dim3
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_New2DVar")

        call card%getDim(dimName1, dim1)
        call card%getDim(dimName2, dim2)
        call card%getDim(dimName3, dim3)

        call card%addVar(varName, 3, &
            [dim1%size,dim2%size,dim3%size], timeVariant)
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
                [dim1%id,dim2%id,dim3%id,card%timeDimId], var%id)
        else
            ierr = nf90_def_var(card%id, varName, xtype, &
                [dim1%id,dim2%id,dim3%id], var%id)
        end if
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "long_name", longName)
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(card%id, var%id, "units", unitName)
        call NFWrap_HandleError(ierr)

        ierr = nf90_enddef(card%id)
        call NFWrap_HandleError(ierr)

        call MsgManager_DeleteSpeaker

    end subroutine NFWrap_New3DVar

    subroutine NFWrap_Advance(card, timeStep, time)
        type(FileCard), intent(inout) :: card
        integer, intent(in) :: timeStep
        real(4), intent(in) :: time

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

    subroutine NFWrap_Output3DDoubleVar(card, varName, varValue)
        type(FileCard), intent(inout), target :: card
        character(*), intent(in) :: varName
        real(8), intent(in) :: varValue(:,:,:)

        integer dimSize(3), ierr
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_Output3DVar")
    
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
        else if (dimSize(3) /= var%dimSize(3)) then
            call MsgManager_Speak(Error, &
                "Output variable """//trim(varName)//""": Second dimension not match ("//&
                trim(int2str(dimSize(3)))//" /= "//trim(int2str(var%dimSize(3)))//")!")
            call RunManager_EndRun
        end if

        if (var%timeVariant) then
            ierr = nf90_put_var(card%id, var%id, varValue, [1,1,1,card%timeStep])
        else
            ierr = nf90_put_var(card%id, var%id, varValue)
        end if
        call NFWrap_HandleError(ierr)
 
#if (defined DEBUG)
        ierr = nf90_sync(card%id)
        call NFWrap_HandleError(ierr)
#endif

        call MsgManager_DeleteSpeaker
    
    end subroutine NFWrap_Output3DDoubleVar

    subroutine NFWrap_Input1DFloatVar(card, varName, varValue, timeStep)
        type(FileCard), intent(in), target :: card
        character(*), intent(in) :: varName
        real(4), intent(out) :: varValue(:)
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

    end subroutine NFWrap_Input1DFloatVar

    subroutine NFWrap_Input1DDoubleVar(card, varName, varValue, timeStep)
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

    end subroutine NFWrap_Input1DDoubleVar

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

    subroutine NFWrap_Input3DDoubleVar(card, varName, varValue, timeStep)
        type(FileCard), intent(in), target :: card
        character(*), intent(in) :: varName
        real(8), intent(out) :: varValue(:,:,:)
        integer, intent(in), optional :: timeStep

        integer dimSize(3), ierr
        type(Variable), pointer :: var

        call MsgManager_RecordSpeaker("NFWrap_Input3DVar")

        dimSize = shape(varValue)

        call card%getVar(varName, var)

        if (dimSize(1) /= var%dimSize(1) .or. &
            dimSize(2) /= var%dimSize(2) .or. &
            dimSize(3) /= var%dimSize(3)) then
            call MsgManager_Speak(Error, &
                "Output variable """//trim(varName)//""": Dimension not match!")
            call RunManager_EndRun
        end if

        if (var%timeVariant .and. present(timeStep)) then
            ierr = nf90_get_var(card%id, var%id, varValue, &
                [1,1,1,timeStep], [dimSize(1),dimSize(2),dimSize(3),1])
        else
            ierr = nf90_get_var(card%id, var%id, varValue)
        end if

        call MsgManager_DeleteSpeaker

    end subroutine NFWrap_Input3DDoubleVar
    
    ! ************************************************************************ !
    !                        Attribute interface                               !
    ! ************************************************************************ !

    ! To be continued ...

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
