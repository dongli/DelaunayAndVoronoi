module CommonTypes

    implicit none

    ! ======================================================================== !
    !                                     1
    type SubHandle
        character(20) :: modName = "N/A"
        character(20) :: subName = "N/A"
        procedure(), nopass, pointer :: handle
        type(SubHandle), pointer :: next
    end type SubHandle

    type OperationList
        character(20) :: name = "N/A"
        integer :: num = 0
        type(SubHandle), pointer :: head
    contains
        procedure :: register => OperationList_register
        procedure :: dump => OperationList_dump
    end type OperationList

    ! ======================================================================== !
    !                                     2
    type RealPtr
        real(RealKind), pointer :: ptr
    end type RealPtr

    type TwoTimeLevel
        type(RealPtr) value(2)
    contains
        procedure :: getNew => TwoTimeLevel_getNew
        procedure :: setNew => TwoTimeLevel_setNew
        procedure :: getOld => TwoTimeLevel_getOld
        procedure :: init => TwoTimeLevel_init
        procedure :: link => TwoTimeLevel_link
        procedure :: mirror => TwoTimeLevel_mirror
        procedure :: save => TwoTimeLevel_save
    end type TwoTimeLevel

contains

    subroutine OperationList_register(a, modName, subName, handle)
        class(OperationList), intent(inout) :: a
        character(*), intent(in) :: modName, subName
        procedure(), intent(in), pointer :: handle

        type(SubHandle), pointer :: sh
        integer i

        a%num = a%num+1
        if (a%num == 1) then
            allocate(a%head)
            nullify(a%head%next)
            sh => a%head
        else
            sh => a%head
            do i = 2, a%num-1
                sh => sh%next
            end do
            allocate(sh%next)
            sh => sh%next
            nullify(sh%next)
        end if

        sh%modName = modName
        sh%subName = subName
        sh%handle => handle

    end subroutine OperationList_register

    subroutine OperationList_dump(a)
        class(OperationList), intent(in) :: a
    
        type(SubHandle), pointer :: sh
        integer i

        write(*, "('Operation list of ', A)") trim(a%name)
        sh => a%head
        do i = 1, a%num
            write(*, "(' handle ', I2, ': ', A, '::', A)") i, trim(sh%modName), trim(sh%subName)
            sh => sh%next
        end do

    end subroutine OperationList_dump
    
    real(RealKind) function TwoTimeLevel_getNew(a) result(res)
        class(TwoTimeLevel), intent(in) :: a

        res = a%value(1)%ptr

    end function TwoTimeLevel_getNew

    subroutine TwoTimeLevel_setNew(a, value)
        class(TwoTimeLevel), intent(out) :: a
        real(RealKind), intent(in) :: value
    
        a%value(1)%ptr = value

    end subroutine TwoTimeLevel_setNew
    
    real(RealKind) function TwoTimeLevel_getOld(a) result(res)
        class(TwoTimeLevel), intent(in) :: a

        res = a%value(2)%ptr

    end function TwoTimeLevel_getOld

    subroutine TwoTimeLevel_init(a)
        class(TwoTimeLevel), intent(out) :: a

        integer i

        allocate(a%value(1)%ptr)
        allocate(a%value(2)%ptr)
        
    end subroutine TwoTimeLevel_init

    subroutine TwoTimeLevel_link(a, value)
        class(TwoTimeLevel), intent(out) :: a
        real(RealKind), intent(in), target :: value

        integer i        

        a%value(1)%ptr => value

        allocate(a%value(2)%ptr)

    end subroutine TwoTimeLevel_link
    
    subroutine TwoTimeLevel_mirror(a, b)
        class(TwoTimeLevel), intent(out) :: a
        class(TwoTimeLevel), intent(in) :: b

        integer i

        a%value(1)%ptr => b%value(1)%ptr
        a%value(2)%ptr => b%value(2)%ptr
        
    end subroutine TwoTimeLevel_mirror
    
    subroutine TwoTimeLevel_save(a)
        class(TwoTimeLevel), intent(inout) :: a

        integer i

        a%value(2)%ptr = a%value(1)%ptr
        
    end subroutine TwoTimeLevel_save

end module CommonTypes
