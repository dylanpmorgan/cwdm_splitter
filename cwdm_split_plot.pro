FUNCTION cwdm_split_plot,light=light

COMMON share1

cleanplot,/silent

pyps,filename=out_eps_path+STRTRIM(outfile_id,2)+'.eps',/color,/encapsulated
   DEVICE, SET_FONT='Palatino-Roman',/ISOLATIN1

   erase
   loadct,39
   angstrom = string(197B)

;; Do a few calculations to try and find the right yaxis-plotting.
;; ind = where(COM.obswave ge 3600. AND COM.obswave le 5500 or COM.obswave ge 5700)
;; minmaxflux = minmax([COM.obsflux[ind],COM.obsflux[ind]-COM.wdflux[ind],COM.dmflux[ind]])

;; ind2 = where(COM.obswave ge 6000 and COM.obswave le 8000)
;; testxyout = minmax([COM.obsflux[ind2],COM.obsflux[ind2]-COM.wdflux[ind2],COM.dmflux[ind2]])
;; testxyout2 = where([COM.obsflux[ind2],COM.obsflux[ind2]-COM.wdflux[ind2],COM.dmflux[ind2]] gt minmaxflux[1]-((minmaxflux[1]-minmaxflux[0])*0.20),cnt)

;; if cnt gt 1 then minmaxflux[1] += minmaxflux[1]*0.15

ind = where(com.obswave ge 3800)
fluxes = [com.obsflux[ind],com.wdflux[ind],com.dmflux[ind], $
          com.obsflux[ind]-com.dmflux[ind],com.obsflux[ind]-com.wdflux[ind]]
meanclip,fluxes,msigma,3.,subs=subs
perfranges = minmax(fluxes[subs])
perfranges[1] += perfranges[1]*0.15

IF KEYWORD_SET(light) THEN BEGIN
   !x.margin = [5,3]
   !y.margin = [3,2]
   multiplot,[1,2]
   ;; Plot WD
   plot,COM.obswave,COM.obsflux, $
        ;yrange=[minmaxflux[0],minmaxflux[1]],/ystyle, $
        yrange=perfranges,/ystyle, $
        xrange=minmax(COM.obswave),/xstyle, $
        charsize=0.9
   oplot,COM.obswave,COM.obsflux-COM.dmflux,color=80
   oplot,COM.obswave,COM.wdflux,color=50

   ;; xyouts dM info
   xyouts,0.55,0.48,'Resulting dM Spectrum (WD+dM - WD)', $
          charsize=0.8,charthick=1.1,color=210,/normal
   xyouts,0.55,0.45,'dM Best Fit: '+string(COM.dm_spt[COM.lap],f='(F3.1)')+$
                    " ["+COM.dm_act[COM.lap]+"]", $
          charsize=0.8,charthick=1.1,color=250,/normal

   multiplot
   ;; Plot dM
   plot,COM.obswave,COM.obsflux,/noerase, $
        xtit=textoidl(' \lambda')+' ('+angstrom+')', $
        ;yrange=[minmaxflux[0],minmaxflux[1]],/ystyle, $
        yrange=perfranges,/ystyle, $
        xrange=minmax(COM.obswave),/xstyle, $
        charsize=0.9
   oplot,COM.obswave,COM.obsflux-COM.wdflux,color=210
   oplot,COM.obswave,COM.dmflux,color=250

   ;; xyouts WD info
   xyouts,0.55,.89,'Resulting WD Spectrum (WD+dM - dM)',charsize=0.8,charthick=1.2,color=80,/normal
   xyouts,0.55,.86,'WD Best Fit: '+COM.wdtype[COM.lap]+', '+textoidl('T_{eff}: ')+string(COM.wdteff[COM.lap],f='(F7.1)')+textoidl(', log(g): ')+string(COM.wdg[COM.lap],f='(F4.2)'),charsize=0.8,charthick=1.2,color=50,/normal

   xyouts,0.03,0.35,textoidl('Flux (10^{-17} erg s^{-1} cm^{-2}'+angstrom+'^{-1})'),charsize=1.0,orientation=90,/normal

   tit = STRSPLIT(raw_spec_fname,'.',/EXTRACT)
   xyouts,0.33,0.94,tit[0],charsize=1.5,/normal
ENDIF ELSE BEGIN

   ;; Plot WD
   plot,COM.obswave,COM.obsflux,position=[.08,.50,.50,.90], $
        xtit=textoidl(' \lambda')+' ('+angstrom+')', $
        yrange=perfranges,/ystyle, $
        xrange=[min(COM.obswave)-20.,max(COM.obswave)+20.],/xstyle, $
        charsize=0.9
   oplot,COM.obswave,COM.obsflux-COM.dmflux,color=80
   oplot,COM.obswave,COM.wdflux,color=50

   ;; xyouts dM info
   xyouts,0.65,.86,'Resulting dM Spectrum (WD+dM - WD)',charsize=0.6,charthick=1.1,color=210,/normal
   xyouts,0.65,.84,'dM Best Fit: '+string(COM.dm_spt[COM.lap],f='(F3.1)')+" ["+COM.dm_act[COM.lap]+"]",charsize=0.6,charthick=1.1,color=250,/normal

   ;; Plot dM
   plot,COM.obswave,COM.obsflux,position=[.50,.50,.92,.90],/noerase, $
        xtit=textoidl(' \lambda')+' ('+angstrom+')', $
        yrange=perfranges,ystyle=12, $
        xrange=[min(COM.obswave)-20.,max(COM.obswave)+20.],/xstyle, $
        charsize=0.9
   axis,yaxis=1,charsize=0.8
   oplot,COM.obswave,COM.obsflux-COM.wdflux,color=210
   oplot,COM.obswave,COM.dmflux,color=250

   ;; xyouts WD info
   xyouts,0.22,.86,'Resulting WD Spectrum (WD+dM - dM)',charsize=0.6,charthick=1.2,color=80,/normal
   xyouts,0.22,.84,'WD Best Fit: '+COM.wdtype[COM.lap]+', '+textoidl('T_{eff}: ')+string(COM.wdteff[COM.lap],f='(F7.1)')+textoidl(', log(g): ')+string(COM.wdg[COM.lap],f='(F4.2)'),charsize=0.6,charthick=1.2,color=50,/normal

;; ind = where(obswave le 5500 or obswave ge 5700)
;; minmaxflux = minmax([obsflux[ind],obsflux[ind]-COM.dmflux[ind],COM.wdflux[ind]])
;; if cnt gt 1 then minmaxflux[1] += minmaxflux[1]*0.20

   ;; Plot WD RV fitting
   plot,COM.obswave,COM.obsflux,position=[.08,.08,.36,.40],/noerase, $
        xrange=[4800,4920],/xstyle, $
        xtit=textoidl(' \lambda')+' ('+angstrom+')', $
        charsize=0.6
   oplot,COM.obswave,COM.wd_restflux,color=100

   lines=[4861.36,4340.47,4101.74]
   airtovac,lines
   wrv = COM.wdrv[COM.lap]
   wrv_err = COM.wdrv_err[COM.lap]

   oplot,[1,1]*lines[0],[-1e6,1e6]
   oplot,[1,1]*lines[0]+(lines[0]*(wrv/299792.458)),[-1e6,1e6],color=150
   oplot,[1,1]*lines[0]+(lines[0]*((wrv-wrv_err)/299792.458)),[-1e6,1e6],color=150,line=3
   oplot,[1,1]*lines[0]+(lines[0]*((wrv+wrv_err)/299792.458)),[-1e6,1e6],color=150,line=3

   xyouts,.10,.10,strtrim(string(wrv,f='(F7.1)'),2)+textoidl('\pm')+strtrim(string(wrv_err,f='(F5.1)'),2),charsize=0.7,/normal

   ;; Plot dM RV fitting
   plot,COM.obswave,COM.obsflux,position=[.64,.08,.92,.40],/noerase, $
        xrange=[8150,8230],/xstyle, $
        ytickformat='(A1)', $
        xtit=textoidl(' \lambda')+' ('+angstrom+')', $
        charsize=0.6
   axis,yaxis=1,charsize=0.7
   oplot,COM.obswave,COM.dm_restflux,color=50

   lines = [8183.255,8194.790]
   airtovac,lines
   drv = COM.dmrv[COM.lap]
   drv_err = COM.dmrv_err[COM.lap]

   oplot,[1,1]*lines[0],[-1e6,1e6]
   oplot,[1,1]*lines[0]+(lines[0]*(drv/299792.458)),[-1e6,1e6],color=250
   oplot,[1,1]*lines[0]+(lines[0]*((drv-drv_err)/299792.458)),[-1e6,1e6],color=250,line=3
   oplot,[1,1]*lines[0]+(lines[0]*((drv+drv_err)/299792.458)),[-1e6,1e6],color=250,line=3

   oplot,[1,1]*lines[1],[-1e6,1e6]
   oplot,[1,1]*lines[1]+(lines[1]*(drv/299792.458)),[-1e6,1e6],color=50
   oplot,[1,1]*lines[1]+(lines[1]*((drv-drv_err)/299792.458)),[-1e6,1e6],color=50,line=3
   oplot,[1,1]*lines[1]+(lines[1]*((drv+drv_err)/299792.458)),[-1e6,1e6],color=50,line=3

   xyouts,.82,.10,string(drv,f='(F7.1)')+textoidl('\pm')+strtrim(string(drv_err,f='(F5.1)'),2),charsize=0.7,/normal

   ;; Plot SDSS gri snapshot
   ;; First if gri isn't already downloaded, grab gri snapshot.
   IF file_test(out_jpeg_path+outfile_id+'.jpeg') EQ 0 THEN BEGIN
      link = "'"+'http://skyservice.pha.jhu.edu/dr8/ImgCutout/getjpeg.aspx?'+ $
             'ra='+strtrim(ra,2)+$
             '&dec='+strtrim(dec,2)+$
             '&scale=0.40&width=100&height=100&opt='+"'"

      spawn_str = 'wget -O '+out_jpeg_path+outfile_id+'.jpeg '+link
      spawn,spawn_str

      read_jpeg,out_jpeg_path+outfile_id+'.jpeg',gri_snap
   ENDIF ELSE read_jpeg,out_jpeg_path+outfile_id+'.jpeg',gri_snap

   tvimage,gri_snap,position=[.36,.08,.64,.40],/keep_aspect_ratio

  ;; xyouts top top title
   IF dec LT 0 THEN radecstr = strtrim(string(ra),2)+strtrim(string(dec),2) ELSE $
      radecstr =strtrim(string(ra),2)+'+'+strtrim(string(dec),2) 
   xyouts,0.33,0.96,radecstr,charsize=1.5,/normal

   ;; xyouts bottom top title
   xyouts,.30,.92,'Plate:'+STRTRIM(STRING(spec_plate),2)+' MJD:'+STRTRIM(STRING(spec_mjd),2)+' Fiber:'+STRTRIM(STRING(spec_fiberid),2),charsize=1.3,/normal

   ;; xyouts top ylabel
   xyouts,0.03,0.28,textoidl('Flux (10^{-17} erg s^{-1} cm^{-2}'+angstrom+'^{-1})'),charsize=1.0,orientation=90,/normal
ENDELSE

pyps,/close

END

