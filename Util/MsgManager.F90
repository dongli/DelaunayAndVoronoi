module MsgManager

    implicit none

    integer, parameter :: Notice  = 1
    integer, parameter :: Error   = 2
    integer, parameter :: Warning = 3

    integer, parameter :: maxStrLen = 256

    ! ======================================================================== !

    type Speaker
        character(maxStrLen) :: name = ""
    end type Speaker

    integer, parameter :: maxNumSpeaker = 20

    type SpeakerStack
        integer :: n = 0
        type(Speaker) s(maxNumSpeaker)
    contains
        procedure :: printCallFlow => SpeakerStack_printCallFlow
    end type SpeakerStack

    type(SpeakerStack), target :: speakers
    type(Speaker), pointer :: currSpeaker

    ! ======================================================================== !

    type StringPair
        character(maxStrLen) key
        character(maxStrLen) value
    end type StringPair

    type Config
        character(30) modName
        type(StringPair) pair
    end type Config

    integer, parameter :: maxNumConfig = 50

    type ConfigStack
        integer :: n = 0
        type(Config) s(maxNumConfig)
    end type ConfigStack

    type(ConfigStack), target :: configs

    ! ======================================================================== !

    interface real2str
        module procedure db2str
    end interface real2str

contains

    subroutine MsgManager_RecordSpeaker(name)
        character(*), intent(in) :: name

        speakers%n = speakers%n+1
        if (speakers%n > maxNumSpeaker) then
            call MsgManager_Speak(Error, &
                "MsgManager_RecordSpeaker encounters an internal error!")
            stop
        end if
        currSpeaker => speakers%s(speakers%n)
        currSpeaker%name = name

    end subroutine MsgManager_RecordSpeaker

    subroutine MsgManager_DeleteSpeaker

        if (speakers%n == 0) then
            call MsgManager_Speak(Error, &
                "MsgManager_DeleteSpeaker encounters an internal error!")
            stop
        end if
        speakers%n = speakers%n-1
        if (speakers%n /= 0) then
            currSpeaker => speakers%s(speakers%n)
        else
            nullify(currSpeaker)
        end if
    
    end subroutine MsgManager_DeleteSpeaker
    
    subroutine MsgManager_Speak(type, content, comment)
        integer, intent(in) :: type
        character(*), intent(in) :: content
        character(*), intent(in), optional :: comment

        if (speakers%n == 0) then
            write(*, "('Error: MsgManager_Speak: Internal error!')")
            stop
        end if

        select case (type)
        case (Notice)
            write(*, "('Notice: ')", advance="no")
        case (Error)
            write(*, "('******')")
            write(*, "('Error: ')", advance="no")
        case (Warning)
            write(*, "('Warning: ')", advance="no")
        end select
        call speakers%printCallFlow
        if (content /= "") then
            write(*, "(A)") trim(content)
        else
            write(*, "('Some one want to say, but give no content.')")
        end if
        if (comment /= "") then
            write(*, "('Comment: ', A)") trim(comment)
        end if

    end subroutine MsgManager_Speak

    subroutine MsgManager_AddConfig(modName, key, value)
        character(*), intent(in) :: modName, key, value
    
        configs%n = configs%n+1
        if (configs%n > maxNumConfig) then
            call MsgManager_Speak(Error, &
                "MsgManager_AddConfig encounters an internal error!")
            stop
        end if
        configs%s(configs%n)%modName = modName
        configs%s(configs%n)%pair%key = key
        configs%s(configs%n)%pair%value = value

    end subroutine MsgManager_AddConfig

    subroutine MsgManager_ShowConfig
        integer i
        character(30) :: currModName = "Nothing"

        write(*, "('****************CONFIGURATION*****************')")
        do i = 1, configs%n
            if (currModName /= configs%s(i)%modName) then
                if (currModName /= "Nothing") write(*, *)
                currModName = configs%s(i)%modName
                write(*, "('Module ', A, ': ')") trim(currModName)
            end if
            write(*, "(2x, A, ': ', A)") trim(configs%s(i)%pair%key), trim(configs%s(i)%pair%value)
        end do
        write(*, "('**********************************************')")

    end subroutine MsgManager_ShowConfig
    
    function int2str(i) result(s)
        integer, intent(in) :: i
        character(maxStrLen) s

        character(20) formatStr
        integer n

        n = 1
        do while(i/(10**n) /= 0)
            n = n+1
        end do
        if (i < 0) n = n+1
        write(formatStr, "('(I', I10, ')')") n
        write(s, formatStr) i

    end function int2str

    function db2str(d, m, n) result(s)
        real(8), intent(in) :: d
        integer, intent(in) :: m, n
        character(maxStrLen) :: s

        character(maxStrLen) formatStr

        write(formatStr, "('(F', I10, '.', I5, ')')") m, n
        write(s, formatStr) d

    end function db2str

    subroutine SpeakerStack_printCallFlow(a)
        class(SpeakerStack), intent(in) :: a

        integer i
    
        do i = 1, speakers%n
            if (speakers%s(i)%name /= "") then
                write(*, "(A, ': ')", advance="no") trim(speakers%s(i)%name)
            end if
        end do
    
    end subroutine SpeakerStack_printCallFlow
    

end module MsgManager
