PRO cwdm_splitter,cube,outfile,beg=beg,fin=fin,fast=fast,new=new

; Create a COMMON block, this allows us to access variables in
; different functions and procedures used by wdmsplitter.
COMMON share1,out_file,out_fits_path,out_eps_path,out_jpeg_path,spec_path,wd_m,splined_wdtemps,wdtemp_id,splined_dmtemps,dmtemp_id,COM,raw_spec_fname,outfile_id,ra,dec,spec_sn,spec_plate,spec_mjd,spec_fiberid,dm_t_ind,wd_m_ind,int_label,velspacing

;; Setting paths -- May need to change
out_fits_path = '/Users/dpmorg/gdrive/research/cwdm_flares/new/data/split_fits/allspec/'
out_eps_path = '/Users/dpmorg/gdrive/research/cwdm_flares/new/plots/split_eps/allspec/'
out_jpeg_path = '/Users/dpmorg/gdrive/research/cwdm_flares/new/plots/gri_snapshots/'
spec_path = '/Users/dpmorg/gdrive/research/cwdm_flares/new/data/spectra/'
temp_path = '/Users/dpmorg/gdrive/research/cwdm_splitter/data/'

;; Input fits structure parameters.
;; REQ: Spec file names -- set location below with spec_path
;;      Outfile ID
;;      RA
;;      DEC
;; OPT: snMedian -- Just for plot labeling purposes
;;      PLATE -- Just for plot labeling purposes
;;      MJD -- Just for plot labeling purposes
;;      FIBERID -- Just for plot labeling purposes
IF TAG_EXIST(cube, 'specfile') EQ 1 THEN BEGIN
   raw_spec_fnames = cube.specfile
ENDIF ELSE BEGIN
   raw_spec_fnames = replicate('', n_elements(cube.ra))
   FOR i=0,n_elements(cube.ra)-1 DO BEGIN
      raw_spec_fnames[i] = 'spec-'+string(cube[i].plate,f='(i04)')+'-'+strtrim(cube[i].mjd,2)+'-'+string(cube[i].fiberid,f='(i04)')+'.fits'
   ENDFOR
ENDELSE

IF TAG_EXIST(cube, 'outfile_id') EQ 1 THEN BEGIN
   out_ids = cube.outfile_id
ENDIF ELSE BEGIN
   out_ids = replicate('', n_elements(cube.ra))
   FOR i=0,n_elements(cube.ra)-1 DO BEGIN
      ;; Create our naming scheme for the files. Format is
      ;; s82cwdm_RA_DEC_PLATE_MJD_FIBER Add .eps or .fits to this name
      rastr = strtrim(string(cube[i].ra,f='(f009.5)'),2)
      decstr =  strtrim(string(abs(cube[i].dec),f='(f8.5)'),2)

      IF cube[i].dec LT 0 THEN radecstr = rastr+'-'+decstr ELSE radecstr = rastr+'+'+decstr 
      out_ids[i] = 's82cwdm_'+radecstr
   ENDFOR
ENDELSE 

ras = cube.ra
decs = cube.dec

IF TAG_EXIST(cube,'snMedian') EQ 1 THEN spec_sns = cube.snMedian $
   ELSE spec_sns = replicate('',n_elements(cube.ra))
IF TAG_EXIST(cube,'PLATE') EQ 1 THEN spec_plates = cube.plate $
   ELSE spec_plates = replicate('',n_elements(cube.ra))
IF TAG_EXIST(cube,'MJD') EQ 1 THEN spec_mjds = cube.mjd $
   ELSE spec_mjds = replicate('',n_elements(cube.ra))
IF TAG_EXIST(cube,'FIBERID') EQ 1 THEN spec_fiberids = cube.fiberid $
   ELSE spec_fiberids = replicate('',n_elements(cube.ra))

;; If the input structure doesn't have the output tags, add the tags.
IF TAG_EXIST(cube,'DM_TEMP') EQ 0 AND $
   TAG_EXIST(cube,'DM') EQ 0 AND $
   TAG_EXIST(cube,'WD_TEMP') EQ 0 AND $
   TAG_EXIST(cube,'WD_TYPE') EQ 0 AND $
   TAG_EXIST(cube,'WD_TEFF') EQ 0 AND $
   TAG_EXIST(cube,'WD_LOGG') EQ 0 AND $
   TAG_EXIST(cube,'DM_RV_RAW') EQ 0 AND $
   TAG_EXIST(cube,'DM_RVERR_RAW') EQ 0 AND $
   TAG_EXIST(cube,'WD_RV_RAW') EQ 0 AND $
   TAG_EXIST(cube,'WD_RVERR_RAW') EQ 0 THEN BEGIN
    
   FLOAT = FLOAT(-9999)
   DBL = DOUBLE(-9999)
   STRING = '-'
   split = REPLICATE({SPECFILE: STRING, $
                      OUTFILE_ID: STRING, $
                      DM_TEMP: STRING, $
                      DM: FLOAT, $
                      WD_TEMP: STRING, $
                      WD_TYPE: STRING, $
                      WD_TEFF: FLOAT, $
                      WD_LOGG: FLOAT, $
                      DM_RV_RAW: FLOAT, $
                      DM_RVERR_RAW: FLOAT, $
                      WD_RV_RAW: FLOAT, $
                      WD_RVERR_RAW: FLOAT},N_ELEMENTS(cube.ra))

   ;split.outfile_id = out_ids
   ;split.specfile = raw_spec_fnames

   cube = struct_addtags(cube,split)
ENDIF
                    
;; Load WD models (Koester) and dM templates (Bochanski)
;; To save run time on splining template wavelengths to sloan wavelengths
;; (which are always the same), the _rebin fits file has the templates
;; and models pre-splined. 
readcol,temp_path+'template.list',f='A',wdtemp_id,/silent 
wd_m=mrdfits(temp_path+'koester_data_rebin.fits',1,/SILENT)
;wd_m=mrdfits(temp_path+'koester_data.fits',1,/SILENT)

readcol,temp_path+'bochamp_temp.list',f='A',dmtemp_id,/silent 
dm_t=mrdfits(temp_path+'bochanski_data_rebin.fits',1,/SILENT)
;dm_t=mrdfits(temp_path+'bochanski_data.fits',1,/SILENT)

;; Begin to loop through the WD+dM and begin the splitting process.
for i=beg,fin do begin
   ;; Start timer
   st = systime(/seconds)

   ;; Tracking
   print,'Begin bisecting spectrum '+strtrim(string(i),2)

   ;; Needed for file naming and plot labeling - passed through COMMON block
   int_label = STRING(i,f='(I003)')
   raw_spec_fname=STRTRIM(raw_spec_fnames[i],2)
   outfile_id = int_label+'_'+STRTRIM(out_ids[i],2)
   ra = ras[i]
   dec = decs[i]
   spec_sn = spec_sns[i]
   spec_plate = spec_plates[i]
   spec_mjd = spec_mjds[i]
   spec_fiberid = spec_fiberids[i]

   ;; Check to make sure the spectra file is available.
   IF file_test(spec_path+raw_spec_fname) THEN BEGIN
      ;; spec = mrdfits(spec_path+raw_spec_fname,1,/silent)
      spec = mrdfits(spec_path+raw_spec_fname,1,hdr,/silent)
   ENDIF ELSE BEGIN
      print,'Spectra File Missing - Downloading'
      download_sdss_spectra,cube[i]
      IF file_test(spec_path+raw_spec_fname) THEN BEGIN
         spec = mrdfits(spec_path+raw_spec_fname,1,hdr,/silent)
      ENDIF ELSE BEGIN
         print, 'Failed to download - STOPPING'
         continue
      ENDELSE
   ENDELSE 

   ;; Set necessary parameters for 
   obswave=10^spec.loglam
   obsflux = spec.flux
   obsnoise = sqrt(1/spec.ivar)
   noise0 = where(obsnoise eq 0.0,noisect)
   if noisect ne 0 then obsnoise[where(obsnoise eq 0.0)] = 1.0
   testfile = STRSPLIT(raw_spec_fname,'/',/EXTRACT)
   ;; IF testfile[0] EQ 'coadded' THEN BEGIN
   ;;    obswave = TRANSPOSE(spec[0,*])
   ;;    obsflux = TRANSPOSE(spec[1,*])
   ;;    obsnoise = TRANSPOSE(spec[2,*])
   ;;    noise0 = where(obsnoise eq 0.0,noisect)
   ;;    if noisect ne 0 then obsnoise[where(obsnoise eq 0.0)] = median(obsnoise)
   ;; ENDIF ELSE BEGIN
   ;;    crval1 = SXPAR(hdr,'CRVAL1')
   ;;    naxis1 = SXPAR(hdr,'NAXIS1')

   ;;    obswave = findgen(naxis1)+crval1
   ;;    obsflux = spec[*,0,0]*1d17 ;10^-17 erg/s/cm^2/A
   ;;    obsnoise = spec[*,0,3]*1d17 ;10^-17 erg/s/cm^2/A
   ;;    noise0 = where(obsnoise eq 0.0,noisect)
   ;;    if noisect ne 0 then obsnoise[where(obsnoise eq 0.0)] = median(obsnoise)
   ;; ENDELSE

   ;; Put on a logarithmic wavelength spacing, even in velocity
   ;; space. We want it to be as close to the velocity spacing of the
   ;; Halpha line as possible. This is since we're not doing the
   ;; wavelength rebinning perfectly and will not be conserving flux.
   ;; nwave = N_ELEMENTS(obswave)-1
   ;; delwave = obswave[0:nwave-1]-obswave[1:nwave-2]
   ;; halpha_loc = where(obswave gt 6555 and obswave lt 6565)
   ;; halpha_vel = mean((delwave[halpha_loc]/obswave[halpha_loc]))
   ;; velspacing = abs(halpha_vel*2.99792458d5)
   
   ;; lo = alog(min(obswave))
   ;; hi = alog(max(obswave))
   ;; new_nwave = round((hi-lo)/abs(halpha_vel))
   ;; xst = (hi-lo)/(new_nwave)
   ;; newwave = exp(lo+xst*findgen(new_nwave))

   ;; ;; Interpolate flux over new wavelength spacings.
   ;; obsflux = interpol(obsflux,obswave,newwave,/spline)
   ;; obsnoise = interpol(obsnoise,obswave,newwave,/spline)
   ;; obswave = newwave

   ;; ;; Convert vactoair
   ;; IF testfile[1] EQ 'hilt' THEN BEGIN
   ;;    airtovac,obswave
   ;; ENDIF

   ;; Check to make sure the raw sdss spectra actually has flux measurements.
   flux_chk = where(abs(obsflux) gt 0.,flux_chk0)
   if flux_chk0 eq 0 then begin
      print,'Flux information missing - Skipping'
      continue
   endif

   ;; Structure for storing and passing information between the
   ;; splitting procedures.
   COM = {obswave:obswave,obsflux:obsflux,obsnoise:obsnoise, $
          dmflux:fltarr(n_elements(obsflux)),wdflux:fltarr(n_elements(obsflux)), $
          dm_restflux:fltarr(n_elements(obsflux)),wd_restflux:fltarr(n_elements(obsflux)), $
          dm_spt:strarr(10),dm_act:strarr(10), $
          wdtype:strarr(10),wdg:fltarr(10),wdteff:fltarr(10), $
          dm_normind:fltarr(2,10)*0.-9999.,wd_normind:fltarr(2,10)*0.-9999., $
          dmrv:fltarr(10)*0.,dmrv_err:fltarr(10)*0., $
          wdrv:fltarr(10)*0.,wdrv_err:fltarr(10)*0.,$ 
          wdsn:0.,dmsn:0.,lap:0L}

   ;; Spline the wd models and dm templates
   ;; Create placeholder array for splined WD templates.
   ;; splined_wdtemps = dblarr(n_elements(wd_m.wave[*,0]),n_elements(obswave))
   ;; FOR w=0,n_elements(wd_m.wave[*,0])-1 DO BEGIN
   ;;    ;; Remove any erroneous data points.
   ;;    dn = where(finite(wd_m.wave[w,*]))
   ;;    ;; Spline the WD template fluxes to the SDSS wavelength.
   ;;    spf = interpol(wd_m.flux[w,dn],wd_m.wave[w,dn],obswave,/spline)
   ;;    ;; Insert into placeholder array.
   ;;    splined_wdtemps[w,0:n_elements(spf)-1] = spf
   ;; ENDFOR

   ;; splined_dmtemps = dblarr(n_elements(dm_t.wave[*,0]),n_elements(obswave))
   ;; for m=0,n_elements(dm_t.wave[*,0])-1 do begin
   ;;    dn = where(finite(dm_t.wave[m,*]))
   ;;    spf = interpol(dm_t.flux[m,dn],dm_t.wave[m,dn],obswave,/spline)
   ;;    missing_data_chk = where(obswave le min(dm_t.wave[m,dn]) or obswave ge max(dm_t.wave[m,dn]),missing_data_cnt)
   ;;    IF missing_data_cnt GT 0 THEN spf[missing_data_chk] = 0.0
   ;;    splined_dmtemps[m,0:n_elements(spf)-1] = spf
   ;; endfor

   ;; Trim WD models and dM templates to the observed wavelength ranges.
   wd_m_ind = where(wd_m.wave[0,*] ge min(obswave) and wd_m.wave[0,*] le max(obswave))
   splined_wdtemps = wd_m.flux[*,wd_m_ind]
   dm_t_ind = where(dm_t.wave[0,*] ge min(obswave) and dm_t.wave[0,*] le max(obswave))
   splined_dmtemps = dm_t.flux[*,dm_t_ind]

   ;; Start the fitting process. Start with WD, then dM, then
   ;; iterate on WD and dM until all the parameters agree or
   ;; we reach 10 iterations.
   cwdm_split_wd
   cwdm_split_dm
   repeat begin
      COM.lap++
      cwdm_split_wd
      cwdm_split_dm
   endrep until (COM.dm_spt[COM.lap-1] eq COM.dm_spt[COM.lap]) and (COM.wdg[COM.lap-1] eq COM.wdg[COM.lap]) and (COM.wdteff[COM.lap-1] eq COM.wdteff[COM.lap]) or COM.lap eq 9

   ;; Print the following information, useful as a counter.
   print,'# of laps: ',strtrim(COM.lap,2)
   print,'dM types: ',string(COM.dm_spt[0:COM.lap],f='(F3.1)')
   print,'WD Log(g): ',string(COM.wdg[0:COM.lap],f='(F4.2)')
   print,'WD Teff: ',string(COM.wdteff[0:COM.lap],f='(F7.1)')

   ;; Save each flux component in sloan file format.
   cwdm_make_file

   ;; Make plots
   makeplots=cwdm_split_plot()

   ;; Store the pertinent information in a structure that will be
   ;; saved as fits file.
   j = i
   cube[j].specfile = strtrim(raw_spec_fname,2)
   cube[j].outfile_id = strtrim(outfile_id,2)
   ;; cube[j].ra = ra
   ;; cube[j].dec = dec
   cube[j].dm_temp = strtrim(dmtemp_id[dm_t_ind],2)
   cube[j].wd_temp = strtrim(wdtemp_id[wd_m_ind],2)
   cube[j].dm = float(COM.dm_spt[COM.lap])
   cube[j].wd_type = strtrim(COM.wdtype[COM.lap],2)
   cube[j].wd_teff = float(strtrim(COM.wdteff[COM.lap],2))
   cube[j].wd_logg = float(strtrim(COM.wdg[COM.lap],2))
   cube[j].dm_rv_raw = COM.dmrv[COM.lap]
   cube[j].dm_rverr_raw = COM.dmrv_err[COM.lap]
   cube[j].wd_rv_raw = COM.wdrv[COM.lap]
   cube[j].wd_rverr_raw = COM.wdrv_err[COM.lap]
   IF KEYWORD_SET(fast) THEN BEGIN
      IF j MOD 50 EQ 0 THEN mwrfits,cube,outfile,/create
   ENDIF ELSE mwrfits,cube,outfile,/create

;; End timer, print elapsed time.
et = systime(/seconds) & print,'Time taken to complete: ',et-st,' seconds'
endfor

;; For the lols
print,'Master Computer: Procedure successfully executed.  Awaiting further commands.'

END
