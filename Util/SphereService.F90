module SphereService

    use MsgManager
    use RunManager
    use FloatingPoint

    implicit none

    real(RealKind), parameter :: Zero  = 0.0
    real(RealKind), parameter :: PI    = 4.0*atan(1.0)
    real(RealKind), parameter :: PI2   = 2.0*PI
    real(RealKind), parameter :: PI05  = 0.5*PI
    real(RealKind), parameter :: PI025 = 0.25*PI
    real(RealKind), parameter :: PI15  = 1.5*PI
    real(RealKind), parameter :: Re    = 6371.229
    real(RealKind), parameter :: Re2   = Re**2.0
    real(RealKind), parameter :: Rad2Deg = 180.0/PI

    interface calcArea
        module procedure calcArea1, calcArea2
    end interface calcArea

contains

    subroutine RotationTransform(lonP, latP, lonO, latO, lonR, latR)
        real(RealKind), intent(in) :: lonP, latP ! rotated pole coordinate
        real(RealKind), intent(in) :: lonO, latO ! original coordinate
        real(RealKind), intent(out), optional :: lonR, latR ! rotated coordinate

        real(RealKind) temp1, temp2, temp3, dlon

        call MsgManager_RecordSpeaker("RotationTransform")

        dlon = lonO-lonP
        if (present(lonR)) then
            temp1 = cos(latO)*sin(dlon)
            temp2 = cos(latO)*sin(latP)*cos(dlon)-cos(latP)*sin(latO)
            lonR = atan2(temp1, temp2)
            if (lonR < 0.0d0) lonR = PI2+lonR
        end if
        if (present(latR)) then
            temp1 = sin(latO)*sin(latP)
            temp2 = cos(latO)*cos(latP)*cos(dlon)
            temp3 = temp1+temp2
            temp3 = min(max(temp3, -1.0d0), 1.0d0)
            latR = asin(temp3)
        end if

        call MsgManager_DeleteSpeaker

    end subroutine RotationTransform

    subroutine InverseRotationTransform(lonP, latP, lonO, latO, lonR, latR)
        real(RealKind), intent(in) :: lonP, latP  ! rotated pole coordinate
        real(RealKind), intent(out) :: lonO, latO ! original coordinate
        real(RealKind), intent(in) :: lonR, latR  ! rotated coordinate

        real(RealKind) temp1, temp2, temp3

        call MsgManager_RecordSpeaker("InverseRotationTransform")

        temp1 = cos(latR)*sin(lonR)
        temp2 = sin(latR)*cos(latP)+cos(latR)*cos(lonR)*sin(latP)
        if (abs(temp2) < eps) temp2 = 0.0d0 ! DONG Li: This trick is due to the inaccuracy of trigonometry calculation.
        lonO = atan2(temp1, temp2)
        lonO = lonP+lonO
#if (defined DEBUG)
        if (.not. FloatingPoint_Check(lonO, "lonO")) then
            call RunManager_EndRun
        end if
#endif
        if (lonO > PI2) lonO = lonO-PI2
        temp1 = sin(latR)*sin(latP)
        temp2 = cos(latR)*cos(latP)*cos(lonR)
        temp3 = temp1-temp2
        temp3 = min(max(temp3, -1.0d0), 1.0d0)
        latO = asin(temp3)
#if (defined DEBUG)
        if (.not. FloatingPoint_Check(latO, "latO")) then
            call RunManager_EndRun
        end if
#endif

        call MsgManager_DeleteSpeaker

    end subroutine InverseRotationTransform

    subroutine CartesianTransform(lon, lat, x, y, z)
        real(RealKind), intent(in) :: lon, lat
        real(RealKind), intent(out) :: x, y, z

        real(RealKind) ReCosLat

        ReCosLat = Re*cos(lat)
        x = ReCosLat*cos(lon)
        y = ReCosLat*sin(lon)
        z = Re*sin(lat)
    
    end subroutine CartesianTransform

    subroutine InverseCartesianTransform(lon, lat, x, y, z)
        real(RealKind), intent(out) :: lon, lat
        real(RealKind), intent(in) :: x, y, z

        lon = atan2(y, x)
        lat = asin(z/Re)
    
    end subroutine InverseCartesianTransform

    subroutine CartesianTransformOnUnitSphere(lon, lat, x, y, z)
        real(RealKind), intent(in) :: lon, lat
        real(RealKind), intent(out) :: x, y, z

        real(RealKind) CosLat

        CosLat = cos(lat)
        x = CosLat*cos(lon)
        y = CosLat*sin(lon)
        z = sin(lat)

    end subroutine CartesianTransformOnUnitSphere

    subroutine InverseCartesianTransformOnUnitSphere(lon, lat, x, y, z)
        real(RealKind), intent(out) :: lon, lat
        real(RealKind), intent(in) :: x, y, z

        lon = atan2(y, x)
        lat = asin(z)

    end subroutine InverseCartesianTransformOnUnitSphere

    ! ************************************************************************ !
    ! Area                                                                     !
    ! Purpose:                                                                 !
    !   Calculate the area of a spherical polygon.                             !
    ! ************************************************************************ !

    subroutine calcArea1(n, x, y, z, area)
        integer, intent(in) :: n
        real(RealKind), intent(in) :: x(n), y(n), z(n)
        real(RealKind), intent(out) :: area

        real(RealKind) norm(3,n), len, dot, angle(n)
        real(RealKind) excess ! The spherical excess

        integer i, j

        do i = 1, n
            j = merge(i-1, n, i /= 1)
            norm(1,i) = y(i)*z(j)-z(i)*y(j)
            norm(2,i) = x(i)*z(j)-z(i)*x(j)
            norm(3,i) = x(i)*y(j)-y(i)*x(j)
            len = sqrt(norm(1,i)**2+norm(2,i)**2+norm(3,i)**2)
            norm(:,i) = norm(:,i)/len
        end do

        do i = 1, n
            j = merge(i+1, 1, i /= n)
            dot = norm(1,i)*norm(1,j)+norm(2,i)*norm(2,j)+norm(3,i)*norm(3,j)
            angle(i) = acos(-dot)
        end do

        excess = sum(angle)-(n-2)*PI

        area = excess*Re2

    end subroutine calcArea1

    subroutine calcArea2(n, lon, lat, area)
        integer, intent(in) :: n
        real(RealKind), intent(in) :: lon(n), lat(n)
        real(RealKind), intent(out) :: area

        real(RealKind) x(n), y(n), z(n)

        integer i

        do i = 1, n
            call CartesianTransformOnUnitSphere(lon(i), lat(i), &
                x(i), y(i), z(i))
        end do
        call calcArea1(n, x, y, z, area)

    end subroutine calcArea2

end module SphereService
