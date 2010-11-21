module SphereService

    use MsgManager
    use RunManager
    use FloatingPoint
    use CommonUtils

    implicit none

    real(RealKind), parameter :: Zero  = 0.0
    real(RealKind), parameter :: Equator = 0.0
    real(RealKind), parameter :: PI    = 4.0*atan(1.0)
    real(RealKind), parameter :: PI2   = 2.0*PI
    real(RealKind), parameter :: PI05  = 0.5*PI
    real(RealKind), parameter :: PI025 = 0.25*PI
    real(RealKind), parameter :: PI15  = 1.5*PI
    real(RealKind), parameter :: Re    = 6371.229
    real(RealKind), parameter :: Re2   = Re**2.0
    real(RealKind), parameter :: Rad2Deg = 180.0/PI

    interface GreatCircleArcIntersection
        module procedure GreatCircleArcIntersection1
        module procedure GreatCircleArcIntersection2
    end interface GreatCircleArcIntersection

    interface calcArea
        module procedure calcArea1, calcArea2
    end interface calcArea

contains

    ! ************************************************************************ !
    ! Rotation transform                                                       !
    ! Purpose:                                                                 !
    !   Calculate the rotating transformation and its inverse of the original  !
    !   coordinate system (lonO,latO) to the rotated one (lonR, latR) with     !
    !   the north pole (lonP,latP) defined at the original coordinate system.  !
    ! ************************************************************************ !

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

        real(RealKind) sinLonR, cosLonR, sinLatR, cosLatR, sinLatP, cosLatP
        real(RealKind) temp1, temp2, temp3

        call MsgManager_RecordSpeaker("InverseRotationTransform")

        sinLonR = sin(lonR)
        cosLonR = cos(lonR)
        sinLatR = sin(latR)
        cosLatR = cos(latR)
        sinLatP = sin(latP)
        cosLatP = cos(LatP)

        temp1 = cosLatR*sinLonR
        temp2 = sinLatR*cosLatP+cosLatR*cosLonR*sinLatP
        ! This trick is due to the inaccuracy of trigonometry calculation.
        if (abs(temp2) < eps) temp2 = 0.0d0
        lonO = atan2(temp1, temp2)
        lonO = lonP+lonO
#if (defined DEBUG)
        if (.not. FloatingPoint_Check(lonO, "lonO")) then
            call RunManager_EndRun
        end if
#endif
        if (lonO > PI2) lonO = lonO-PI2
        temp1 = sinLatR*sinLatP
        temp2 = cosLatR*cosLatP*cosLonR
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

    ! ************************************************************************ !
    ! Cartesian transform                                                      !
    ! Purpose:                                                                 !
    !   Calculate the transformation and its inverse of spherical coordinate   !
    !   (lon,lat) to the Cartesian coordinate (x,y,z).                         !
    ! ************************************************************************ !

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

        if (lon < Zero) then
            lon = lon+PI2
        end if

    end subroutine InverseCartesianTransformOnUnitSphere

    ! ************************************************************************ !
    ! GreatCircleArcIntersection                                               !
    ! Purpose:                                                                 !
    !   Calculate the intersection between two great circle arcs.              !
    ! ************************************************************************ !

    subroutine GreatCircleArcIntersection1(x1, y1, z1, x2, y2, z2, x3, y3, z3, &
                                          x4, y4, z4, xi, yi, zi)
        real(RealKind), intent(in) :: x1, y1, z1, x2, y2, z2 ! The first arc
        real(RealKind), intent(in) :: x3, y3, z3, x4, y4, z4 ! The second arc
        real(RealKind), intent(out) :: xi, yi, zi ! The intersection point

        real(RealKind) n1(3), n2(3), len

        call CrossProduct(x1, y1, z1, x2, y2, z2, n1(1), n1(2), n1(3))
        call CrossProduct(x3, y3, z3, x4, y4, z4, n2(1), n2(2), n2(3))
        call CrossProduct(n1(1), n1(2), n1(3), n2(1), n2(2), n2(3), xi, yi, zi)

        len = sqrt(xi**2+yi**2+zi**2)
        if (abs(len) < eps) then
            call MsgManager_Speak(Error, &
                "Encounter problematic intersection situation, "// &
                "the vector's length is too small.")
            call RunManager_EndRun
        end if
        xi = xi/len
        yi = yi/len
        zi = zi/len

    end subroutine GreatCircleArcIntersection1

    subroutine GreatCircleArcIntersection2(lon1, lat1, lon2, lat2, lon3, lat3, &
                                           lon4, lat4, loni, lati)
        real(RealKind), intent(in) :: lon1, lat1, lon2, lat2 ! The first arc
        real(RealKind), intent(in) :: lon3, lat3, lon4, lat4 ! The second arc
        real(RealKind), intent(out) :: loni, lati ! The intersection point

        real(RealKind) x1, y1, z1, x2, y2, z2
        real(RealKind) x3, y3, z3, x4, y4, z4
        real(RealKind) xi, yi, zi

        call CartesianTransformOnUnitSphere(lon1, lat1, x1, y1, z1)
        call CartesianTransformOnUnitSphere(lon2, lat2, x2, y2, z2)
        call CartesianTransformOnUnitSphere(lon3, lat3, x3, y3, z3)
        call CartesianTransformOnUnitSphere(lon4, lat4, x4, y4, z4)

        call GreatCircleArcIntersection1(x1, y1, z1, x2, y2, z2, x3, y3, z3, &
                                         x4, y4, z4, xi, yi, zi)

        call InverseCartesianTransformOnUnitSphere(loni, lati, xi, yi, zi)

        if (loni < Zero) then
            loni = loni+PI2
        end if

    end subroutine GreatCircleArcIntersection2

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
            call CrossProduct(x(i), y(i), z(i), x(j), y(j), z(j), &
                norm(1,i), norm(2,i), norm(3,i))
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

        real(RealKind), parameter :: eps = 1.0e-6
        real(RealKind) x(n), y(n), z(n)

        integer i

        do i = 1, n
            call CartesianTransformOnUnitSphere(lon(i),lat(i), x(i),y(i),z(i))
        end do
        call calcArea1(n, x, y, z, area)

    end subroutine calcArea2

end module SphereService
