PRO cwdm_calcmultirv,cube,allspecs
loadct,39
COMMON share1,obswave,obsflux,dmflux,wdflux,splined_norm_dmtemp,splined_norm_wdmod,rvstr,i,j,str_ct,samp_ra,samp_dec,specfnames

; Procedure for calculating RV of dMs and WDs using xcorrelation
; technique. 
; When using the SDSS spectrum, we use the each individual composite
; spectra in order to determine any variability in the RV calculations
; and thus try to observe orbital motion.


; Read in the S82 fits file containing multiple spectra for each
; object.

;; Fits table of the sample without duplicates. Output structure will
;; be based off of this
rvstr_file = '/Users/dpmorg/Desktop/research/cwdm_flares/data/s82_cwdm_rvs.fits'

samp = cube

;; allspecs - Fits table with all the duplicates.

;; samp = allspecs
;; n = n_elements(samp.ra)-1
;; gcirc,2,samp[0:n].ra,samp[0:n].dec,samp[1:*].ra,samp[1:*].dec,dist
;; dup_ind = where(dist le 3.0)
;; remove,dup_ind,samp

; Read in template information. Need rest templates
readcol,'/Users/dpmorg/Google Drive/idlpros/wdmsplitter/dependencies/template.list',f='A',wdtemp_id,/silent 
wd_m=mrdfits('/Users/dpmorg/Google Drive/idlpros/wdmsplitter/dependencies/koester_data_rebin.fits',1,/SILENT)

readcol,'/Users/dpmorg/Google Drive/idlpros/wdmsplitter/dependencies/bochamp_temp.list',f='A',dmtemp_id,/silent 
dm_t=mrdfits('/Users/dpmorg/Google Drive/idlpros/wdmsplitter/dependencies/bochanski_data_rebin.fits',1,/SILENT)

; End Input Variables
; ~~~~~~~~~~~~~~

; ~~~~~~~~~~~~~~
; Create structure to contain RV information

rvstr = replicate( {ra:-999.,dec:-999.,filenames:replicate('-',250),dm_snr:replicate(-999.,250),wd_snr:replicate(-999.,250),dmrv: replicate(-999.,250), dmrv_err: replicate(-999.,250), wdrv: replicate(-999.,250), wdrv_err: replicate(-999.,250), mjd: replicate(-999L,250), tai_beg: replicate(-999., 250), tai_end: replicate(-999., 250)}, n_elements(samp.ra))

FOR i=0,n_elements(samp.ra)-1 DO BEGIN

   ;; Find spectra of the same object
   gcirc,2,samp[i].ra,samp[i].dec,allspecs.ra,allspecs.dec,dist
   dup_ind = where(dist le 3.0,dupct)
   if dupct eq 0 then continue
   samp_ra = strtrim(string(samp[i].ra,f='(f10.5)'),2)
   samp_dec = strtrim(string(samp[i].dec,f='(f12.6)'),2)
   rvstr[i].ra = samp[i].ra
   rvstr[i].dec = samp[i].dec

   ; For some reason there the same spectrum is reported many
   ; times. Match plate,mjd, and fiber
   specfnames = allspecs[dup_ind].specfile
   splitfnames = allspecs[dup_ind].outfile_id

   ;; plates = allspecs[dup_ind].plate
   ;; mjds = allspecs[dup_ind].mjd
   ;; fibers = allspecs[dup_ind].fiberid
   ;; detectors= strlowcase(strtrim(allspecs[dup_ind].instrument,2))
   ;; runs = strtrim(allspecs[dup_ind].run2d,2)

   ;; sort_mjds = sort(mjds)
   ;; plates = plates[sort_mjds]
   ;; mjds = mjds[sort_mjds]
   ;; fibers = fibers[sort_mjds]
   ;; detectors = detectors[sort_mjds]
   ;; runs = runs[sort_mjds]

   ;; plates_uniq = string(plates[uniq(mjds)],f='(i04)')
   ;; mjds_uniq = strtrim(mjds[uniq(mjds)],2)
   ;; fibers_uniq = string(fibers[uniq(mjds)],f='(i04)')
   ;; detectors_uniq = detectors[uniq(mjds)]
   ;; runs_uniq = runs[uniq(mjds)]
   ;; specfnames = 'spec-'+plates_uniq+'-'+mjds_uniq+'-'+fibers_uniq+'.fits'
   ;; splitfnames = strtrim(.fits'

   str_ct = 0L
   FOR j=0,N_ELEMENTS(specfnames)-1 DO BEGIN
      rvstr[i].filenames[str_ct] = specfnames[j]
                                
      ;; Download spectrum if necessary
      IF file_test('/Users/dpmorg/Desktop/research/cwdm_flares/data/raw_spectra/'+specfnames[j]) THEN BEGIN
      ENDIF ELSE BEGIN
         IF strtrim(detectors_uniq[j],2) EQ 'segue1' THEN detector = 'sdss' ELSE detector = detectors_uniq[j]
         spawnstr = 'wget --content-disposition http://data.sdss3.org/sas/dr12/'+detector+'/spectro/redux/'+runs_uniq[j]+'/spectra/'+plates_uniq[j]+'/'+specfnames[j]
         spawn,spawnstr
      ENDELSE 

      ;; Read in the split fits file. This is where we'll grab
      ;; our best fit template.
      FOR m=0,N_ELEMENTS(splitfnames)-1 DO BEGIN
         IF file_test('/Users/dpmorg/Desktop/research/cwdm_flares/data/split_spectra/'+strtrim(splitfnames[m],2)+'.fits') THEN split=mrdfits('/Users/dpmorg/Desktop/research/cwdm_flares/data/split_spectra/'+strtrim(splitfnames[m],2)+'.fits',1)
      ENDFOR

      ct = 1
      repeat begin
         spec= mrdfits('/Users/dpmorg/Desktop/research/cwdm_flares/data/raw_spectra/'+strtrim(specfnames[j],2), ct, hdr, status=status,/silent)
         ;; print,status
         IF status LT 0 THEN BEGIN
            GOTO, JUMP1
         ENDIF

         ;; grab MJD, TAI-BEG, TAI-END from header
         if ct gt 1 then begin
            rvstr[i].filenames[str_ct] = specfnames[j]
            rvstr[i].mjd[str_ct] = sxpar(hdr,'MJD')
            rvstr[i].tai_beg[str_ct] = sxpar(hdr,'TAI-BEG')
            rvstr[i].tai_end[str_ct] = sxpar(hdr,'TAI-END')
         endif

         obswave=10^spec.loglam
         obsflux = spec.flux
         obsnoise = sqrt(1/spec.ivar)

         ;; Trim obswave and split.wave to the same start and end
         ;; wavelength values.
         st = max([min(obswave),min(split.wave)])
         ed = min([max(obswave),max(split.wave)]) 

         obsind = where(obswave GE st AND obswave LE ed)
         obswave = obswave[obsind]
         obsflux = obsflux[obsind]
         obsnoise = obsnoise[obsind]

         split_trim = split[where(split.wave GE st AND split.wave LE ed)]

        ; ~~~~~~~~~~~~~~
        ; Calculate RVs for the dM
        ; ~~~~~~~~~~~~~~
         ;; IF max(obswave) LT 8300 THEN BEGIN
         ;;    stop
         ;;    GOTO,JUMP_DM
         ;; ENDIF

         dind = where( (obswave ge 6100 and obswave le 6200) or (obswave ge 8140 and obswave le 8180) or (obswave ge 8220 and obswave le 8260))
         rvstr[i].dm_snr[str_ct] = median(obsflux[dind]/obsnoise[dind])
    
         ; Normalize rest dM template
         ;; dmtempflux = split_trim.dm_temp_flux
         ;; dmtempwave = split_trim.wavelength
         ;; norm_ind = where(obswave ge 7450 and obswave le 7550)
         ;; norm_indt = where(dmtempwave ge 7450 and dmtempwave le 7550)
         ;; splined_norm_dmtemp = (dmtempflux/median(dmtempflux[norm_indt]))*median(obsflux[norm_ind])
         splined_norm_dmtemp = interpol(split.dm_temp_flux,split.wave,obswave,/spline)
         dmflux = obsflux; - splined_norm_wdmod
         

         upwave = max(obswave)
         lowave = 6100.

         ; Calculate RVs
         rvct=0L
         allrvs = fltarr(10000)*0.-9999.

         feat=[6250.,7150.,7700.,8200.]
         offset = 150.
         nudge = set_range(0,205,5)-100.

         !p.multi=0
         FOR mm=0,3 DO BEGIN
            FOR nn=0,N_ELEMENTS(nudge)-1 DO BEGIN    
               rvind = WHERE(obswave GE feat[mm]-offset+nudge[nn] AND obswave LE feat[mm]+offset+nudge[nn],rvind_ct)
               IF rvind_ct GT 30. THEN BEGIN
                  ;; print,rvind_ct
                  xcorl,dmflux[rvind],splined_norm_dmtemp[rvind],15,dmft,ff
                  IF ff EQ 1 THEN allrvs[rvct] = -9999. ELSE allrvs[rvct] = dmft*69.09
                  rvct++
                  ;; plot,obswave,obsflux,xrange=[feat[mm]-400.,feat[mm]+400.]
                  ;; oplot,obswave[rvind],obsflux[rvind],color=190.
                  ;; oplot,obswave,splined_norm_dmtemp,color=50
                  ;; oplot,obswave[rvind],splined_norm_dmtemp[rvind],color=90
                  ;; print,allrvs[rvct-1]
                  ;; wait,0.1
               ENDIF
            ENDFOR
         ENDFOR

         allrvs_ind = where(allrvs ne -9999.,allrvs_cnt)
         IF allrvs_cnt gt 0 then allrvs = allrvs[where(allrvs ne -9999.)]

         ;; median = median(allrvs)
         ;; mom = moment(allrvs,mdev=mad)

         median = median(allrvs)
         mad = mad(allrvs)

         sigcut = where(abs(allrvs-median) le 6.*mad)

         median = median(allrvs[sigcut])
         mad = mad(allrvs[sigcut])

         ;; median = robust_mean(allrvs,3,mad) 

         rvstr[i].dmrv[str_ct] = median
         rvstr[i].dmrv_err[str_ct] = sqrt(mad^2 + 7.^2)
        ; ~~~~~~~~~~~~~~
        ; END dM
        ; ~~~~~~~~~~~~~~
         JUMP_DM:print,'Wave OOR for dM, Skipped'

        ; ~~~~~~~~~~~~~~
        ; Calculate RVs for the WD
        ; ~~~~~~~~~~~~~~
         IF min(obswave) gt 4000. THEN BEGIN
            ;; stop
            GOTO,JUMP_WD
         ENDIF

         wind = where(obswave ge 4190 and obswave le 4210 or obswave ge 4500 and obswave le 4600)
         rvstr[i].wd_snr[str_ct] = median(obsflux[wind]/obsnoise[wind])

         ;Normalize rest WD template
         ;; wdmodflux = split_trim.wd_temp_flux
         ;; wdmodwave = split_trim.wavelength         
         ;; norm_ind = where(obswave ge 4550 and obswave le 4650)
         ;; norm_indt = where(wdmodwave ge 4550 and wdmodwave le 4650)
         ;; splined_norm_wdmod= (wdmodflux/median(wdmodflux[norm_indt]))*median(obsflux[norm_ind])
         splined_norm_wdmod = interpol(split.wd_temp_flux,split.wave,obswave,/spline)
         wdflux = obsflux; - splined_norm_dmtemp

         upwave = 6000.
         lowave = min(obswave)

         ; Calculate RVs
         rvct=0L
         allrvs = fltarr(10000)*0.-9999.

         lines=[4861.36,4340.47,4101.74,3970.075]
         airtovac,lines
         offset = 75.
         nudge = set_range(0,202.,2)-100.

         FOR mm=0,3 DO BEGIN
            FOR nn=0,N_ELEMENTS(nudge)-1 DO BEGIN    
               rvind = WHERE(obswave GE lines[mm]-offset+nudge[nn] AND obswave LE lines[mm]+offset+nudge[nn],rvind_ct)
               IF rvind_ct GT 20. THEN BEGIN
                  ;; print,rvind_ct
                  xcorl,wdflux[rvind],splined_norm_wdmod[rvind],10,wdft,ff
                  IF ff EQ 1 THEN allrvs[rvct] = -9999. ELSE allrvs[rvct] = wdft*69.09
                  rvct++
                  ;; plot,obswave,obsflux,xrange=[lines[mm]-300.,lines[mm]+300.]
                  ;; oplot,obswave[rvind],obsflux[rvind],color=190.
                  ;; oplot,obswave,splined_norm_wdmod,color=50
                  ;; oplot,obswave[rvind],splined_norm_wdmod[rvind],color=90
                  ;; print,allrvs[rvct-1]
                  ;; wait,0.1
               ENDIF
            ENDFOR
         ENDFOR

         allrvs_ind = where(allrvs ne -9999.,allrvs_cnt)
         IF allrvs_cnt gt 0 then allrvs = allrvs[where(allrvs ne -9999.)]

         median = median(allrvs)
         mom = moment(allrvs,mdev=mad)

         ;; median = robust_mean(allrvs,3,mad) 

         rvstr[i].wdrv[str_ct] = median
         rvstr[i].wdrv_err[str_ct] = sqrt(mad^2 + 7.^2)
        ; ~~~~~~~~~~~~~~
        ; End WD
        ; ~~~~~~~~~~~~~~

         JUMP_WD:print,'Wave OOR for WD, Skipped'

         plotmultirv

         JUMP1: print,'EOF; Next Spectrum'

         print,'i: ',i
         print,'j: ',j
         print,'ct: ',ct
         print,'str_ct: ',str_ct
         print,'Object: ',samp_ra,samp_dec,f='(a,x,f10.6,x,f11.6)'
         print,'File: ',specfnames[j]

         if ct eq 1 then ct = 3
         ct++
         str_ct++

         wait,0.5
;         stop
      endrep until status lt 0

   ENDFOR
ENDFOR

mwrfits,rvstr,rvstr_file,/CREATE

END

PRO plotmultirv

COMMON share1

!p.multi=[0,2,2]
angstrom = string(197B)

ind = where(obswave le 5500 or obswave ge 5700)
minmaxflux = minmax([obsflux[ind],obsflux[ind]-wdflux[ind],dmflux[ind]])
ind2 = where(obswave ge 6000 and obswave le 8000)
testxyout = minmax([obsflux[ind2],obsflux[ind2]-wdflux[ind2],dmflux[ind2]])
testxyout2 = where([obsflux[ind2],obsflux[ind2]-wdflux[ind2],dmflux[ind2]] gt minmaxflux[1]-((minmaxflux[1]-minmaxflux[0])*0.20),cnt)

if cnt gt 1 then minmaxflux[1] += minmaxflux[1]*0.15
   plot,obswave,obsflux, $
      tit='OBJ: '+samp_ra+' '+samp_dec+' File: '+strtrim(specfnames[j],2), $
      yrange=[minmaxflux[0],minmaxflux[1]],/ystyle, $
      xrange=[3700,9300],/xstyle;, $
   oplot,obswave,dmflux,color=210
   oplot,obswave,splined_norm_dmtemp,color=250

ind = where(obswave le 5500 or obswave ge 5700)
minmaxflux = minmax([obsflux[ind],obsflux[ind]-dmflux[ind],wdflux[ind]])

if cnt gt 1 then minmaxflux[1] += minmaxflux[1]*0.20
   plot,obswave,obsflux, $
      xtit=textoidl(' \lambda')+' ('+angstrom+')', $
;      ytit='Flux', $
      yrange=[minmaxflux[0],minmaxflux[1]],/ystyle, $
      xrange=[3700,9300],/xstyle
   oplot,obswave,wdflux,color=80
   oplot,obswave,splined_norm_wdmod,color=50

   plot,obswave,obsflux,xrange=[8150,8250],xtit='Med RV: '+strtrim(rvstr[i].dmrv[str_ct] ,2)+', Mad RV: '+strtrim(rvstr[i].dmrv_err[str_ct] ,2),xcharsize=1.3
   oplot,obswave,splined_norm_dmtemp,color=100
   lines = [8183.255,8194.790]
   airtovac,lines
   oplot,[1,1]*lines[0],[-1e6,1e6]
   oplot,[1,1]*lines[0]+(lines[0]*(rvstr[i].dmrv[str_ct])/299792.458),[-1e6,1e6],color=250
   oplot,[1,1]*lines[0]+(lines[0]*((rvstr[i].dmrv[str_ct]-rvstr[i].dmrv_err[str_ct])/299792.458)),[-1e6,1e6],color=250,line=3
   oplot,[1,1]*lines[0]+(lines[0]*((rvstr[i].dmrv[str_ct]+rvstr[i].dmrv_err[str_ct])/299792.458)),[-1e6,1e6],color=250,line=3

   oplot,[1,1]*lines[1],[-1e6,1e6]
   oplot,[1,1]*lines[1]+(lines[1]*(rvstr[i].dmrv[str_ct]/299792.458)),[-1e6,1e6],color=50
   oplot,[1,1]*lines[1]+(lines[1]*((rvstr[i].dmrv[str_ct]-rvstr[i].dmrv_err[str_ct])/299792.458)),[-1e6,1e6],color=50,line=3
   oplot,[1,1]*lines[1]+(lines[1]*((rvstr[i].dmrv[str_ct]-rvstr[i].dmrv_err[str_ct])/299792.458)),[-1e6,1e6],color=50,line=3

   plot,obswave,obsflux,xrange=[4250,4400],xtit='Med RV: '+strtrim(rvstr[i].wdrv[str_ct],2)+', Mad RV: '+strtrim(rvstr[i].wdrv_err[str_ct],2)+' S2N: ',xcharsize=1.3
   oplot,obswave,splined_norm_wdmod,color=100
   lines=[4861.36,4340.47,4101.74]
   airtovac,lines
   ;; oplot,[1,1]*lines[0],[-1e6,1e6]
   ;; oplot,[1,1]*lines[0]+(lines[0]*(rvstr[i].wdrv[str_ct])/299792.458),[-1e6,1e6],color=250
   ;; oplot,[1,1]*lines[0]+(lines[0]*((rvstr[i].wdrv[str_ct]-rvstr[i].wdrv_err[str_ct])/299792.458)),[-1e6,1e6],color=250,line=3
   ;; oplot,[1,1]*lines[0]+(lines[0]*((rvstr[i].wdrv[str_ct]-rvstr[i].wdrv_err[str_ct])/299792.458)),[-1e6,1e6],color=250,line=3

   oplot,[1,1]*lines[1],[-1e6,1e6]
   oplot,[1,1]*lines[1]+(lines[1]*(rvstr[i].wdrv[str_ct]/299792.458)),[-1e6,1e6],color=150
   oplot,[1,1]*lines[1]+(lines[1]*((rvstr[i].wdrv[str_ct]-rvstr[i].wdrv_err[str_ct])/299792.458)),[-1e6,1e6],color=150,line=3
   oplot,[1,1]*lines[1]+(lines[1]*((rvstr[i].wdrv[str_ct]+rvstr[i].wdrv_err[str_ct])/299792.458)),[-1e6,1e6],color=150,line=3

   ;; oplot,[1,1]*lines[2],[-1e6,1e6]
   ;; oplot,[1,1]*lines[2]+(lines[2]*(rvstr[i].wdrv[str_ct]/299792.458)),[-1e6,1e6],color=50
   ;; oplot,[1,1]*lines[2]+(lines[2]*((rvstr[i].wdrv[str_ct]-rvstr[i].wdrv_err[str_ct])/299792.458)),[-1e6,1e6],color=50,line=3
   ;; oplot,[1,1]*lines[2]+(lines[2]*((rvstr[i].wdrv[str_ct]-rvstr[i].wdrv_err[str_ct])/299792.458)),[-1e6,1e6],color=50,line=3

END

