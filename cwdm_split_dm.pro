PRO cwdm_split_dm

COMMON share1

flux = COM.obsflux - COM.wdflux

IF COM.lap GT 0 THEN wc = COM.obswave*((COM.dmrv[COM.lap-1])/299792.458) ELSE wc = 0
wave = COM.obswave + wc
sigma = COM.obsnoise

bests = [0.,0.,9.999d10,0]
; Normalization indices
dm_nind = [7350,7400,8050,8100]
for l=0,148 do begin
  dind = where(COM.obswave ge 5500.+(l*25) and COM.obswave le 5500.+((l+1)*25),dindct)
;; for l=0,1 do begin
;;    dind = where(COM.obswave ge dm_nind[l] and COM.obswave le dm_nind[l+1],dindct)
   if dindct lt 5 then continue
   mean_array = (total(smooth(splined_dmtemps[*,dind],5),2)/n_elements(dind)) # (fltarr(n_elements(COM.obswave))*0.0+1.0)
   nf = (smooth(splined_dmtemps,5)/mean_array)*mean(flux[dind])

   ;; chir = where($;COM.obswave ge 5100. and COM.obswave le 5250. or $
   ;;              $;COM.obswave ge 5400. and COM.obswave le 5500. or $
   ;;              $;COM.obswave ge 5800. and COM.obswave le 6100. or $
   ;;              COM.obswave ge 6150. and COM.obswave le 6350. or $
   ;;              COM.obswave ge 6500. and COM.obswave le 6600. or $
   ;;              COM.obswave ge 6610. and COM.obswave le 6990. or $
   ;;              ;COM.obswave ge 7000. and COM.obswave le 7300. or $
   ;;              COM.obswave ge 7570. and COM.obswave le 7800. or $
   ;;              COM.obswave ge 8150. and COM.obswave le 8250. or $
   ;;              COM.obswave ge 8400. and COM.obswave le 8750.)
   chir = where(COM.obswave ge 5000. and COM.obswave le 9000. and finite(COM.obsnoise))

   sigma_big = (sigma[chir])##(fltarr(57,1)*0.0+1.0)
   flux_big = flux[chir]##(fltarr(57,1)*0.0+1.0)
   ;; ranarr = randomu(seed,n_elements(chir))*2.-1.
   ;; ranarr_big = ranarr##(fltarr(57,1)*0.0+1.0)
   nf_noise = nf[*,chir];-(sigma_big*ranarr_big)
   dmchisq_big = total(((flux_big - nf_noise)^2)/(sigma_big^2),2)/(n_elements(flux[chir])-n_elements(splined_dmtemps[*,0])-1.)

   dmnum0 = where(dmchisq_big eq min(dmchisq_big),dmnum0_ct)
   IF dmnum0_ct gt 1 then dmnum0 = dmnum0[0]

   if dmchisq_big[dmnum0] lt bests[2] then begin

      bests = [5500.+(l*25),5500+((l+1)*25),dmchisq_big[dmnum0],dmnum0]
      COM.dmflux = nf[dmnum0,*]
      dmsplt = strsplit(dmtemp_id[dmnum0[0]],'.',/extract) 
      if dmsplt[1] eq '5' then begin 
         type = dmsplt[0]+"."+dmsplt[1]
         act = dmsplt[2]
      endif else begin 
         type = dmsplt[0]
         act = dmsplt[1]
      endelse
      COM.dm_spt[COM.lap] = FLOAT(strmid(type,1,3))
      COM.dm_act[COM.lap] = act
      COM.dmchisq[COM.lap] = dmchisq_big[dmnum0]
      ;; COM.dmchisq_all[COM.lap,*] = dmchisq_big
      COM.dm_normind[*,COM.lap] = [bests[0],bests[1]]
      dm_t_ind = dmnum0
   endif
endfor

;dn = where(finite(dm_t.wave[bests[3],*]))
splined_rest_tempflux = splined_dmtemps[bests[3],*]
;splined_rest_tempflux = interpol(dm_t.flux[bests[3],dn],dm_t.wave[bests[3],dn],COM.obswave,/spline)
clean_flux = where(COM.obswave le 3900 or COM.obswave ge 9100,clean_fluxct)
IF clean_fluxct GT 0 then splined_rest_tempflux[clean_flux] = !VALUES.F_NAN
norm_ind = where(COM.obswave ge bests[0] and COM.obswave le bests[1])
norm_rest_tempflux = (splined_rest_tempflux/mean(splined_rest_tempflux[norm_ind]))*mean(flux[norm_ind])
xind=where(COM.obswave ge 7000. and COM.obswave le 8500.)
xcorl,flux[xind],norm_rest_tempflux[xind],10,dmft,ff,/fine

if ff eq 1 then begin 
   dmrv_flag='Fail' 
   dmrv = 0.
endif else begin
   dmrv_flag='Success'
   dmrv = dmft*69.09
endelse

COM.dm_restflux = norm_rest_tempflux
wdrvs_calcrv = cwdm_calcrv(flux,norm_rest_tempflux,/dm)

;; COM.dmsn = median(COM.obsflux[chir]/COM.obsnoise[chir])
END
