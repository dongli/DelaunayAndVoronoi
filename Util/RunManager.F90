module RunManager

    use MsgManager
    use CommonTypes

    implicit none

    private

    public RunManager_RegisterOperation
    public RunManager_EndRun

    type(OperationList), target :: EndRun

contains

    subroutine RunManager_RegisterOperation(hostName, modName, subName, handle)
        character(*), intent(in) :: hostName, modName, subName
        procedure(), intent(in), pointer :: handle

        call MsgManager_RecordSpeaker("RunManager_RegisterOperation")

        select case (hostName)
        case ("EndRun")
            call EndRun%register(modName, subName, handle)
        case default
            call MsgManager_Speak(Error, &
                "No host subroutine "// &
                trim(hostName)//" for registration.")
            stop
        end select
#if (defined VERBOSE)
        call MsgManager_Speak(Notice, &
            "Add operation """// &
            trim(msg%content)//trim(modName)//"::" &
            trim(msg%content)//trim(subName)// &
            """ to """//trim(hostName)//""".")
#endif
        call MsgManager_DeleteSpeaker

    end subroutine RunManager_RegisterOperation

    subroutine RunManager_EndRun
        type(SubHandle), pointer :: sh
        integer i

        sh => EndRun%head
        do i = 1, EndRun%num
            call sh%handle
            sh => sh%next
        end do

        stop

    end subroutine RunManager_EndRun

end module RunManager
