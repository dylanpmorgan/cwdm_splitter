PRO cwdm_split_wd

COMMON share1

;; Flux of the WD (SDSS - Best fit dM)
flux = COM.obsflux - COM.dmflux

;; Wavelength correction using the previous best wdrv
IF COM.lap GT 0 THEN wc = COM.obswave*((COM.wdrv[COM.lap-1])/299792.458) ELSE wc = 0.

;Applying wavelength correction to the template wavelength array

wave = COM.obswave+wc

;Hydrogen balmer lines
lines=[6562.85,4861.36,4340.47,4101.74,3970.075,3889.064]
;Convert lines from air wavelengths to vacuum wavelengths
airtovac,lines
;It may be necessary to mask out the Hydrogen balmer emission from the
;active dM in the WD part of the spectrum.  Here we calculate the
;approximate emission positions using the shift due to the dM rv.
IF COM.lap GT 0 THEN lines2 = lines+lines*(COM.dmrv[COM.lap-1]/299792.458) ELSE lines2 = lines

;Create a placeholder for [Best lower range, Best higher range, Best
;chisq, best template number]
bests = [0.,0.,9.9d20,0]
;It's not clear to me which normalization range works the best,
;so I opted to try multiple ranges and pick the one with the best
;chisq.
wd_nind = [4400,4450,4600,4650,6350,6400]
for l=0,116 do begin
;;Define normalization region.
  wind = where(COM.obswave ge 3800+(l*25) and COM.obswave le 3800+((l+1)*25),windct)
;; for l=0,4,2 do begin
;;    wind = where(COM.obswave ge wd_nind[l] and COM.obswave le wd_nind[l+1],windct)
                                ;Create a normalizing array with mean
                                ;values of the above defined regions
                                ;for every template.
  if windct lt 5 then continue
   mean_array = (total(splined_wdtemps[*,wind],2)/n_elements(wind)) # (fltarr(n_elements(COM.obswave))*0.0+1.0)
                                ;Now normalize all the templates simultaneously
   nf = (splined_wdtemps/mean_array)*mean(flux[wind])

;----Define regions to weight----
   cw = 3.
   sigma = COM.obsnoise
   ;~~~~Wavelengths less than 6600 angstroms~~~~
   ;; chir = where(abs(flux) gt 0. and $
   ;;              finite(nf) and $
   ;;              COM.obswave le 7500)

   cc = 3. ;7
   beta = 120.
   delta = 100.
   gamma = 65.
   eps = 45.
   zeta = 30.

   wi = 50. ;30 ;og 50
   cc = 3. ;5 ;og 3
   ex = 80. ;50 ;ex 80

   chir = where($;COM.obswave ge 3500. and COM.obswave le lines2[5]-wi or $
                COM.obswave ge lines2[5]-wi and COM.obswave le lines2[5]-cc or $
                COM.obswave ge lines2[5]+cc and COM.obswave le lines2[5]+wi or $
                COM.obswave ge lines2[4]-wi and COM.obswave le lines2[4]-cc or $
                COM.obswave ge lines2[4]+cc and COM.obswave le lines2[4]+wi or $
                COM.obswave ge lines2[3]-wi and COM.obswave le lines2[3]-cc or $
                COM.obswave ge lines2[3]+cc and COM.obswave le lines2[3]+wi or $
                COM.obswave ge lines2[2]-ex and COM.obswave le lines2[2]-cc or $
                COM.obswave ge lines2[2]+cc and COM.obswave le lines2[2]+ex); or $
                ;COM.obswave ge lines2[1]-ex and COM.obswave le lines2[1]-cc or $
                ;COM.obswave ge lines2[1]+cc and COM.obswave le lines2[1]+ex)

   k = 0
   IF COM.lap gt 0 THEN BEGIN
      IF COM.dm_act[COM.lap-1] eq 'active' then ha_em = [4,6] else ha_em = [2,3]
   ENDIF ELSE ha_em = [2,3]

   ;; chir = where(COM.obswave ge lines2[5]-zeta[k] and COM.obswave le lines2[5]-ha_em[0] or $
   ;;              COM.obswave ge lines2[5]+ha_em[0] and COM.obswave le lines2[5]+zeta[k] or $
   ;;              COM.obswave ge lines2[4]-eps[k] and COM.obswave le lines2[4]-ha_em[0] or $
   ;;              COM.obswave ge lines2[4]+ha_em[0] and COM.obswave le lines2[4]+eps[k] or $
   ;;              COM.obswave ge lines2[3]-gamma[k] and COM.obswave le lines2[3]-ha_em[0] or $
   ;;              COM.obswave ge lines2[3]+ha_em[0] and COM.obswave le lines2[3]+gamma[k] or $
   ;;              COM.obswave ge lines2[2]-delta[k] and COM.obswave le lines2[2]-ha_em[1] or $
   ;;              COM.obswave ge lines2[2]+ha_em[1] and COM.obswave le lines2[2]+delta[k] or $
   ;;              COM.obswave ge 4550. and COM.obswave le 4650. or $
   ;;              COM.obswave ge lines2[1]-beta[k] and COM.obswave le lines2[1]-ha_em[1] or $
   ;;              COM.obswave ge lines2[1]+ha_em[1] and COM.obswave le lines2[1]+beta[k] or $
   ;;              COM.obswave ge 5500. and COM.obswave le 5600. and finite(COM.obsnoise))
   

   ;; chir = where(COM.obswave ge min(COM.obswave)+10. and COM.obswave le lines2[4]-cc or $
   ;;                  COM.obswave ge lines2[4]+cc and COM.obswave le lines2[4]+wi or $
   ;;                  COM.obswave ge lines2[3]-wi and COM.obswave le lines2[3]-cc or $
   ;;                  COM.obswave ge lines2[3]+cc and COM.obswave le lines2[3]+wi or $
   ;;                  COM.obswave ge lines2[2]-ex and COM.obswave le lines2[2]-cc or $
   ;;                  COM.obswave ge lines2[2]+cc and COM.obswave le lines2[1]-cc or $
   ;;                  COM.obswave ge lines2[1]+cc and COM.obswave le lines2[1]+ex) 
   
;   print,'WD S/N',median(COM.obsflux[chir]/COM.obsnoise[chir])

   sigma_big = sigma[chir]##(fltarr(n_elements(nf[*,0]),1)*0.0+1.0)
   flux_big = flux[chir]##(fltarr(n_elements(nf[*,0]),1)*0.0+1.0)
   ranarr = randomu(seed,n_elements(chir))*2.-1.
   ranarr_big = ranarr##(fltarr(n_elements(nf[*,0]),1)*0.0+1.0)
   nf_noise = nf[*,chir];-(sigma_big*ranarr_big)
   wdchisq_big = total(((flux_big - nf_noise)^2)/(sigma_big^2),2)/(n_elements(flux[chir])-n_elements(splined_wdtemps[*,0])-1.)
   ;; sigma_big = sigma##(fltarr(n_elements(nf[*,0]),1)*0.0+1.0)
   ;; flux_big = smooth(flux,5)##(fltarr(n_elements(nf[*,0]),1)*0.0+1.0)
   ;; wdchisq_big = total(((flux_big - nf)^2)/(sigma_big^2),2)/(n_elements(flux)-n_elements(splined_wdtemps[*,0])-1.)
   wdnum0 = where(wdchisq_big eq min(wdchisq_big),wdnum0_ct)
   IF wdnum0_ct gt 1 then wdnum0 = wdnum0[0]

   if wdchisq_big[wdnum0] lt bests[2] then begin
      ;; bests = [wd_nind[l],wd_nind[l+1],wdchisq_big[wdnum0],wdnum0]
      bests = [3800+(l*25),3800+((l+1)*25),wdchisq_big[wdnum0],wdnum0]
      COM.wdflux = nf[wdnum0,*]
      COM.wdchisq[COM.lap] = wdchisq_big[wdnum0]
      ;COM.wdchisq_all[COM.lap,*] = wdchisq_big
      wdsplt = strsplit(wdtemp_id[wdnum0[0]],'/,_,.',/extract)
      COM.wdtype[COM.lap] = wdsplt[0]
      COM.wdg[COM.lap] = FLOAT(strmid(wdsplt[2],0,1)+'.'+strmid(wdsplt[2],1,3))
      COM.wdteff[COM.lap] = FLOAT(strmid(wdsplt[1],2))
      COM.wd_normind[*,COM.lap] = [bests[0],bests[1]]
      wd_m_ind = wdnum0
   endif
endfor

; Use rest best fit template/model and the original flux and
; wavelength to calculate the RV
;; dn = where(finite(wd_m.wave[bests[3],*]))
;; splined_rest_tempflux = interpol(wd_m.flux[bests[3],dn],wd_m.wave[bests[3],dn],COM.obswave,/spline)
;; stop
splined_rest_tempflux = splined_wdtemps[bests[3],*]
clean_flux = where(COM.obswave le 3900 or COM.obswave ge 9100,clean_fluxct)
IF clean_fluxct gt 0 then splined_rest_tempflux[clean_flux] = !VALUES.F_NAN
norm_ind = where(COM.obswave ge bests[0] and COM.obswave le bests[1])
norm_rest_tempflux = (splined_rest_tempflux/mean(splined_rest_tempflux[norm_ind]))*mean(flux[norm_ind])

xind = where(COM.obswave gt 3900 and COM.obswave le 6000)
xcorl,flux[xind],norm_rest_tempflux[xind],10,wdft,ff,/fine
if ff eq 1 then begin 
   wdrv = 0.
endif else begin
   wdrv = wdft*69.09
endelse

COM.wd_restflux = norm_rest_tempflux
wdrvs_calcrv = cwdm_calcrv(flux,norm_rest_tempflux,/wd)

;; COM.wdsn = median(COM.obsflux[chir]/COM.obsnoise[chir])

END
