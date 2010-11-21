program TestIntersection

    use SphereService

    implicit none

    real(RealKind) lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4
    real(RealKind) lon, lat

    write(*, "('Input the first great circle arc:')")
    write(*, "('lon1,lat1 > ')", advance="no")
    read(*, *) lon1, lat1
    write(*, "('lon2,lat2 > ')", advance="no")
    read(*, *) lon2, lat2

    write(*, "('Input the second great circle arc:')")
    write(*, "('lon3,lat3 > ')", advance="no")
    read(*, *) lon3, lat3
    write(*, "('lon4,lat4 > ')", advance="no")
    read(*, *) lon4, lat4

    lon1 = lon1/Rad2Deg; lat1 = lat1/Rad2Deg
    lon2 = lon2/Rad2Deg; lat2 = lat2/Rad2Deg
    lon3 = lon3/Rad2Deg; lat3 = lat3/Rad2Deg
    lon4 = lon4/Rad2Deg; lat4 = lat4/Rad2Deg

    call GreatCircleArcIntersection(lon1, lat1, lon2, lat2, lon3, lat3, &
                                    lon4, lat4, lon, lat)

    lon = lon*Rad2Deg; lat = lat*Rad2Deg

    write(*, "('The intersection is:')")
    write(*, *) lon, lat

end program TestIntersection
