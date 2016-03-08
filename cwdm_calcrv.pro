FUNCTION cwdm_calcrv,flux,norm_temp_restflux,wd=wd,dm=dm

COMMON share1

IF keyword_set(dm) THEN BEGIN
   rvct=0L
   allrvs = fltarr(10000)*0.-9999.

   feat=[6250.,7150.,7700.,8200.]
   offset = 150.
   nudge = set_range(0,205,5)-100.

   FOR mm=0,3 DO BEGIN
      FOR nn=0,N_ELEMENTS(nudge)-1 DO BEGIN    
         rvind = WHERE(COM.obswave GE feat[mm]-offset+nudge[nn] AND COM.obswave LE feat[mm]+offset+nudge[nn],rvind_ct)
         IF rvind_ct GT 40. THEN BEGIN
            xcorl,flux[rvind],norm_temp_restflux[rvind],15,dmft,ff
            IF ff EQ 1 THEN allrvs[rvct] = -9999. ELSE allrvs[rvct] = dmft*velspacing
            rvct++
         ENDIF
      ENDFOR
   ENDFOR

   allrvs_ind = where(allrvs ne -9999.,allrvs_cnt)
   IF allrvs_cnt gt 0 then allrvs = allrvs[where(allrvs ne -9999.)]

   median = median(allrvs)
   mad = mad(allrvs)

   sigcut = where(abs(allrvs-median) le 6.*mad)

   median = median(allrvs[sigcut])
   mad = mad(allrvs[sigcut])

   COM.dmrv[COM.lap] = median
   COM.dmrv_err[COM.lap] = sqrt(mad^2 + 7.^2)

   ;; pyps,filename='/Users/dpmorg/Desktop/research/cwdm_flares/plots/split_rvs/'+strtrim(outfname,2)+'_dmrvs.eps',/color,/encapsulated
   ;; DEVICE, SET_FONT='Palatino-Roman',/ISOLATIN1
   ;; !x.margin = [7,4]
   ;; !y.margin = [5,3]
   ;; angstrom = string(197B)
   
   ;;    plot,COM.obswave,COM.obsflux,xrange=[8100,8300]
   ;;    oplot,COM.obswave,norm_temp_restflux,color=50
   ;;    lines = [8183.255,8194.790]
   ;;    airtovac,lines
   ;;    oplot,[1,1]*lines[0],[-1e6,1e6]
   ;;    oplot,[1,1]*lines[0]+(lines[0]*(COM.dmrv/299792.458)),[-1e6,1e6],color=250
   ;;    oplot,[1,1]*lines[0]+(lines[0]*((COM.dmrv-COM.dmrv_err)/299792.458)),[-1e6,1e6],color=250,line=3
   ;;    oplot,[1,1]*lines[0]+(lines[0]*((COM.dmrv+COM.dmrv_err)/299792.458)),[-1e6,1e6],color=250,line=3

   ;;    oplot,[1,1]*lines[1],[-1e6,1e6]
   ;;    oplot,[1,1]*lines[1]+(lines[1]*(COM.dmrv/299792.458)),[-1e6,1e6],color=50
   ;;    oplot,[1,1]*lines[1]+(lines[1]*((COM.dmrv-COM.dmrv_err)/299792.458)),[-1e6,1e6],color=50,line=3
   ;;    oplot,[1,1]*lines[1]+(lines[1]*((COM.dmrv+COM.dmrv_err)/299792.458)),[-1e6,1e6],color=50,line=3

   ;; pyps,/close

ENDIF

IF keyword_set(wd) THEN BEGIN
   rvct=0L
   allrvs = fltarr(10000)*0.-9999.

   lines=[4861.36,4340.47,4101.74,3970.075]
   airtovac,lines
   offset = 75.
   nudge = set_range(0,202.,2)-100.

   FOR mm=0,3 DO BEGIN
      FOR nn=0,N_ELEMENTS(nudge)-1 DO BEGIN    
         rvind = WHERE(COM.obswave GE lines[mm]-offset+nudge[nn] AND COM.obswave LE lines[mm]+offset+nudge[nn],rvind_ct)
         IF rvind_ct GT 40. THEN BEGIN
            xcorl,flux[rvind],norm_temp_restflux[rvind],10,wdft,ff
            IF ff EQ 1 THEN allrvs[rvct] = -9999. ELSE allrvs[rvct] = wdft*velspacing
            rvct++
         ENDIF
      ENDFOR
   ENDFOR

   allrvs_ind = where(allrvs ne -9999.,allrvs_cnt)
   IF allrvs_cnt gt 0 then allrvs = allrvs[where(allrvs ne -9999.)]

   median = median(allrvs)
   mad = mad(allrvs)

   sigcut = where(abs(allrvs-median) le 6.*mad)

   median = median(allrvs[sigcut])
   mad = mad(allrvs[sigcut])

   COM.wdrv[COM.lap] = median
   COM.wdrv_err[COM.lap] = sqrt(mad^2 + 7.^2)

   ;; pyps,filename='/Users/dpmorg/Desktop/research/cwdm_flares/plots/split_rvs/'+strtrim(outfname,2)+'_wdrvs.eps',/color,/encapsulated
   ;; DEVICE, SET_FONT='Palatino-Roman',/ISOLATIN1
   ;; !x.margin = [7,4]
   ;; !y.margin = [5,3]
   ;; angstrom = string(197B)

   ;;    plot,COM.obswave,COM.obsflux,xrange=[4250,4400]
   ;;    oplot,COM.obswave,norm_temp_restflux,color=100
   ;;    lines=[4861.36,4340.47,4101.74]
   ;;    airtovac,lines
   ;;    oplot,[1,1]*lines[1],[-1e6,1e6]
   ;;    oplot,[1,1]*lines[1]+(lines[1]*(COM.wdrv/299792.458)),[-1e6,1e6],color=150
   ;;    oplot,[1,1]*lines[1]+(lines[1]*((COM.wdrv-COM.wdrv_err)/299792.458)),[-1e6,1e6],color=150,line=3
   ;;    oplot,[1,1]*lines[1]+(lines[1]*((COM.wdrv+COM.wdrv_err)/299792.458)),[-1e6,1e6],color=150,line=3

   ;; pyps,/close
ENDIF

END
