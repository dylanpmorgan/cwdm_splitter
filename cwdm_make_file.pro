PRO cwdm_make_file

COMMON share1

dn = where(finite(COM.obsflux))
splitspec = REPLICATE({wave:0.,flux:0.,noise:0.,dm_flux:0.,wd_flux:0.,wd_wave_rvcorr:0.,dm_wave_rvcorr:0.,wd_temp_flux:0.,dm_temp_flux:0.},N_ELEMENTS(dn))

splitspec.wave = COM.obswave[dn]
splitspec.flux = COM.obsflux[dn]
splitspec.noise = COM.obsnoise[dn]
splitspec.wd_flux = COM.obsflux[dn]-COM.dmflux[dn]
splitspec.dm_flux = COM.obsflux[dn]-COM.wdflux[dn]

;; Apply the wavelength correction calculated before
wd_wc = COM.obswave*((COM.wdrv[COM.lap])/299792.458)
splitspec.wd_wave_rvcorr = COM.obswave+wd_wc

dm_wc = COM.obswave*((COM.dmrv[COM.lap])/299792.458)
splitspec.dm_wave_rvcorr = COM.obswave+dm_wc

splitspec.dm_temp_flux = COM.dmflux[dn]
splitspec.wd_temp_flux = COM.wdflux[dn]

mwrfits,splitspec,out_fits_path+strtrim(outfile_id,2)+'.fits',/create

;;   spec_path= path+'spectra_sample1/'
;;   naxis1 = n_elements(wave[dn])
;;   naxis2 = 2
;;   basewave = alog10(wave[dn[0]])
;;   step = (alog10(wave[dn]) - basewave) / findgen(naxis1)
;;   hdr = headfits(path+'launchpad/'+fname)

;;   sxaddpar,hdr,'COEFF0',basewave
;;   sxaddpar,hdr,'COEFF1',step[10]
;;   sxaddpar,hdr,'NAXIS1',naxis1
;;   sxaddpar,hdr,'NAXIS2',naxis2
;;   sxaddpar,hdr,'CRVAL1',basewave
;;   sxaddpar,hdr,'CD1_1',step[10]
;; ;  sxaddpar,hdr,'TEMPLATE','WD: '+temp_id

 ;;  mwrfits,[[flux[dn]],[noise[dn]]],path+'WD+dM_extras/'+strtrim(num,2)+'_wdspec.fits',hdr,/create

;; ;dM
;;   wave = sdsswave
;;   flux = sdssflux - wdflux
;;   noise = sdssnoise
;;   num = snum
;; ; spec_id = spec_ids[snum]
;;   temp_id = dmtemp_id[dmnum]
;;   dn = where(finite(flux))

;;   naxis1 = n_elements(wave[dn])
;;   naxis2 = 2
;;   basewave = alog10(wave[dn[0]])
;;   step = (alog10(wave[dn]) - basewave) / findgen(naxis1)
;;   hdr = headfits(path+'launchpad/'+fname)
;;   sxaddpar,hdr,'COEFF0',basewave
;;   sxaddpar,hdr,'COEFF1',step[10]
;;   sxaddpar,hdr,'NAXIS1',naxis1
;;   sxaddpar,hdr,'NAXIS2',naxis2
;;   sxaddpar,hdr,'CRVAL1',basewave
;;   sxaddpar,hdr,'CD1_1',step[10]
;; ;  sxaddpar,hdr,'TEMPLATE','dM: '+temp_id

;;   mwrfits,[[flux[dn]],[noise[dn]]],path+'WD+dM_extras/'+strtrim(num,2)+'_dmspec.fits',hdr,/create


END
