#!/bin/bash

for file in $(ls *.nc)
do
    figure=$(basename ${file} .nc)
    echo -n "Plotting ${file} ... "
    ncl 1> /dev/null <<-EOF
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

    begin

        f = addfile("${file}", "r")

        numPoint = dimsizes(f->lonPoint)

        numTriangle = dimsizes(f->triangle(:,0))

        wks = gsn_open_wks("pdf", "${figure}")

        ; Plot map
        mapRes = True
        mapRes@gsnFrame = False
        mapRes@mpGreatCircleLinesOn = True
        mapRes@mpGridAndLimbOn = True
        mapRes@mpGridLineColor = "Background"
        mapRes@mpProjection = "Satellite"
        mapRes@mpCenterLonF = 60.0
        mapRes@mpCenterLatF = 45.0
        ;mapRes@mpProjection = "Stereographic"

        map = gsn_csm_map(wks, mapRes)

        if (True) then
            ; Plot Delaunay triangle edges
            edgeRes = True
            edgeRes@gsLineThicknessF = 2.
            edgeRes@gsLineColor = "blue"

            lon = new(4, "float")
            lat = new(4, "float")
            do i = 0, numTriangle-1
                numVertex = 0
                do j = 0, 2
                    k = f->triangle(i,j)-1
                    if (k .ge. 0) then
                        lon(numVertex) = f->lonPoint(k)
                        lat(numVertex) = f->latPoint(k)
                        numVertex = numVertex+1
                    end if
                end do
                if (numVertex .gt. 1) then
                    lon(numVertex) = lon(0)
                    lat(numVertex) = lat(0)
                    gsn_polyline(wks, map, lon(0:numVertex), lat(0:numVertex), edgeRes)
                end if
            end do

            delete(lon)
            delete(lat)
        end if

        if (False) then
            ; Plot circumcirlces
            circRes = True
            circRes@gsLineThicknessF = 2.
            circRes@gsLineColor = "green"

            arclon = new(50, float)
            arclat = new(50, float)
            do i = 0, numTriangle-1
                nggcog(f->latCirc(i), f->lonCirc(i), f->radius(i), arclat, arclon)
                gsn_polyline(wks, map, arclon, arclat, circRes)
            end do
        end if

        if (True) then
            ; Plot Voronoi cell edges
            edgeRes = True
            edgeRes@gsLineThicknessF = 2.
            edgeRes@gsLineColor = "yellow"

            do i = 0, numPoint-1
                lon = new(f->numVVT(i)+1, "float")
                lat = new(f->numVVT(i)+1, "float")
                do j = 0, f->numVVT(i)-1
                    lon(j) = f->lonVVT(i,j)
                    lat(j) = f->latVVT(i,j)
                end do
                lon(f->numVVT(i)) = lon(0)
                lat(f->numVVT(i)) = lat(0)
                gsn_polyline(wks, map, lon, lat, edgeRes)
                delete(lon)
                delete(lat)
            end do
        end if

        ; Plot vertices
        vertexRes = True
        vertexRes@gsMarkerIndex = 1
        vertexRes@gsMarkerSizeF = 0.02
        vertexRes@gsMarkerColor = "red"

        gsn_polymarker(wks, map, f->lonPoint, f->latPoint, vertexRes)

        frame(wks)

    end
EOF
    if [ -f "${figure}.pdf" ]; then
        echo "${figure}.pdf generated."
    else
        echo "NCL error!"
    fi
done
