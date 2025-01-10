;----------------------------------------------------------------------
; circ.ncl 
;
; Concepts illustrated:
;   - Drawing a black-and-white streamline plot over a map
;   - Changing the length of streamlines
;   - Setting the spacing for streamline arrowheads
;----------------------------------------------------------------------
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   

begin
  yr = (/1979, 2017/)
  yr0 = 1948
  rg = (/-30,130,-10,70/)
  rect_lon = (/75, 75, 103, 103, 75/)
  rect_lat = (/26.5, 39.5, 39.5, 26.5, 26.5/)
  
  fu = addfile("~/RD/uwnd.mon.mean.nc","r")
  fv = addfile("~/RD/vwnd.mon.mean.nc","r")
  fq = addfile("~/RD/shum.mon.mean.nc","r")
  lev = fu->level
  IS = minind(abs(lev-300))
  u    = fu->uwnd(:,0:IS,:,:)    ; read in u data
  v    = fv->vwnd(:,0:IS,:,:)    ; read in v data
  q    = fq->shum(:,0:IS,:,:)/1000    ; read in q data
  
  lv  = lev(0:IS)
  lat = fu->lat
  lon = fu->lon

;---Create plots
  wks  = gsn_open_wks("pdf","qtp_iwvf_ncep")               ; send graphics to EPS file
  
  mpid = new(1,graphic)
  
  mth = (/6, 7, 8/)
  xmth = dimsizes(mth)
  nyr = (yr(1)-yr(0)+1)
  nt = (nyr-1)*xmth
  nlev = 8
  nlat = 73
  nlon = 144
  
  ky = 0
  t = new((/(yr(1)-yr(0))*xmth/),"integer")
  do iy = 0,nyr-2
    do im = 0,(xmth-1)
      y = yr(0)-yr0+iy
      m = mth(im)
      t(ky) = y*12+m-1
      ky = ky+1
    end do
  end do
  
  qu = new((/nt,nlev,nlat,nlon/),float)
  qv = new((/nt,nlev,nlat,nlon/),float)
  qu = q(t,:,:,:)*u(t,:,:,:)/9.81
  qv = q(t,:,:,:)*v(t,:,:,:)/9.81
  qu!0         = "time"
  qu!1         = "level"
  qu!2         = "lat"
  qu!3         = "lon"
  qv!0         = "time"
  qv!1         = "level"
  qv!2         = "lat"
  qv!3         = "lon"
  
  psf = new((/nt,nlat,nlon/),float)
  psf!0         = "time"
  psf!1         = "lat"
  psf!2         = "lon"
  psf(:,:,:)    = 1000;
  iwfu = new((/nt,nlat,nlon/),float)
  iwfv = new((/nt,nlat,nlon/),float)
  iwfu = vibeta(lv,qu(time|:,lat|:,lon|:,level|:),1,psf,1000,300)*100
  iwfv = vibeta(lv,qv(time|:,lat|:,lon|:,level|:),1,psf,1000,300)*100

  u1 = dim_avg_n_Wrap(iwfu,0)
  v1 = dim_avg_n_Wrap(iwfv,0)
  delete(t)
;  print(u1(1,:))
  
  lon = fspan(-177.5,180,144)
  u2 = new((/73,144/),float)
  v2 = new((/73,144/),float)
  do i = 0,72
    do j = 0,143
      if(j .lt. 72) then
        u2(i,j+72) = u1(i,j)
        v2(i,j+72) = v1(i,j)
      else
        u2(i,j-72) = u1(i,j)
        v2(i,j-72) = v1(i,j)
      end if
    end do
  end do
  u2!0         = "lat"
  u2!1         = "lon"
  u2&lon       = lon
  u2&lat       = lat
  
  v2!0         = "lat"
  v2!1         = "lon"
  v2&lon       = lon
  v2&lat       = lat
  
  res           = True         ; plot mods desired
  res@gsnDraw   = False        ; don't draw instantly
  res@gsnFrame  = False
  vcres = res
  
  res@mpMinLonF = rg(0)
  res@mpMaxLonF = rg(1)
  res@mpMinLatF = rg(2)
  res@mpMaxLatF = rg(3)
  res@mpCenterRotF = 0.0
  res@mpFillOn  = True
  res@mpOutlineOn            = True
  res@mpLimbLineThicknessF   = 0.5
  res@mpLandFillColor        = (/1, 1, 1/)              ; turn off map fill
  res@mpInlandWaterFillColor = (/0.557, 0.898, 0.933/)  ; fill inlandwaters
  res@mpOceanFillColor       = (/0.557, 0.898, 0.933/)  ; fill oceans
  
  res@gsnLeftString = ""
  res@gsnLeftStringFontHeightF = 0.015
  res@gsnRightString = ""
  res@gsnRightStringFontHeightF = 0.015
  
  res@tmXBMinorOn = False
  res@tmXBTickSpacingF = 30
  res@tmXBLabelFontHeightF = 0.015
  res@tmXBLabels = (/"-30~S~o~N~E","0~S~o~N~E","30~S~o~N~E","60~S~o~N~E","90~S~o~N~E","120~S~o~N~E"/)
  res@tmYLMinorOn = False
  res@tmYLTickSpacingF = 10
  res@tmYLLabelFontHeightF = 0.015
  res@tmYLLabels = (/"10~S~o~N~S","EQ","10~S~o~N~N","20~S~o~N~N","30~S~o~N~N","40~S~o~N~N","50~S~o~N~N","60~S~o~N~N","70~S~o~N~N"/)
  
  vcres@vcGlyphStyle             = "CurlyVector"         ; select wind barbs
  vcres@vcLineArrowColor         = "black"
  vcres@vcLineArrowThicknessF    = 1.2
  vcres@vcLineArrowHeadMaxSizeF  = 0.008
  vcres@vcRefMagnitudeF          = 300.0
  vcres@vcRefLengthF             = 0.08
  vcres@vcMinDistanceF           = 0.015
  vcres@vcRefAnnoOrthogonalPosF  = 0.15
  vcres@vcRefAnnoParallelPosF    = 0.999
  vcres@vcRefAnnoString1On       = "False"
  vcres@vcRefAnnoString2On       = "True"
  vcres@vcRefAnnoString2         = "300 kg/m/s"
  vcres@gsnLeftString = ""
  vcres@gsnRightString = ""

;---适合中国的国界地图 
  shpfn1="~/RD/nb_bou1_4m/bou1_4l.shp"
  lnres        = True
  lnres@minlon = rg(0)
  lnres@maxlon = rg(1)
  lnres@minlat = rg(2)
  lnres@maxlat = rg(3)
  lnres@gsLineColor = "black"
  lnres@gsLineThicknessF = 1
  
  vcid = gsn_csm_vector(wks,u2,v2,vcres)
  mpid = gsn_csm_map(wks,res)
  shpp = gsn_add_shapefile_polylines(wks,mpid,shpfn1,lnres)
  overlay(mpid,vcid)

  rect = new(4,graphic)
  plres                  = True
  plres@gsLineColor      = "black"
  plres@gsLineThicknessF = 2
  do i = 0,3
    rect(i) = gsn_add_polyline(wks,mpid,rect_lon(i:i+1),rect_lat(i:i+1),plres)
  end do
  
  draw(mpid)
  frame(wks)

end
