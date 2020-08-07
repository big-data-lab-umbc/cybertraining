; IDL procedure  -  read_ckd.pro
;
;****************************************************************************
;
; NAME:
;
;      read_ckd
;
; PURPOSE:
;
;      Read an SHDOM-related correlated-k distribution file.
;
; CATEGORY:
;
;      SHDOM utility
;
; CALLING SEQUENCE:
;
;      read_ckd, filename
;
; INPUTS:
;
;      filename:  a string with the directory path and file name
;
; OPTIONAL KEYWORD OUTPUTS:
;
;      NUMB:    number of spectral bands
;      NUMG:    number of g's in each spectral band
;      WAVENO:  starting and ending wavenumbers for each spectral band
;      SOLF:    solar flux [W/m2] for each spectral band
;      DELG:    fractional weight for each g in each band
;      NUMZ:    number of z levels
;      CO2:     concentration of CO2 [ppmv]
;      CH4:     concentration of CH4 [ppmv]
;      N2O:     concentration of N2O [ppmv]
;      Z:       profile heights [km]
;      PRES:    profile pressures [mb]
;      TEMP:    profile temperatures [K]
;      QH2O:    profile water vapor
;      QO3:     profile ozone
;      ABSORB:  absorption coefficients [km^-1] for each band, g, and level
;
; MODIFICATION HISTORY:
;
;      Written, Tim Benner, August 1999
;
; NOTES:
;
;      For use with CKD files created for the Spherical Harmonics Discrete
;        Ordinate Method (SHDOM) (Evans, 1998) radiative transfer model.
;
;****************************************************************************

PRO read_ckd, filename, NUMB=numb, NUMG=numg, WAVENO=waveno, SOLF=solf, $
              DELG=delg, NUMZ=numz, CO2=co2, CH4=ch4, N2O=n2o, $
              Z=z, PRES=pres, TEMP=temp, QH2O=qh2o, QO3=qo3, ABSORB=absorb

;****************************************************************************
; Open the file.
;****************************************************************************

ss = ' '
openr, unit, filename, /get_lun
readf, unit, ss

;****************************************************************************
; Read the number of bands.
;****************************************************************************

readf, unit, ss
i = strpos(ss,'!')
ss = strtrim(strmid(ss,0,i),2)
numb = long(ss)

;****************************************************************************
; Read the number of g's.
;****************************************************************************

readf, unit, ss
maxg = 0
numg = fltarr(numb)
for b=0,numb-1 do begin
    readf, unit, n, wn0, wn1, sf, ng
    numg(b) = ng
  endfor
maxg = fix(max(numg))

;****************************************************************************
; Close the file.
;****************************************************************************

close, unit
free_lun, unit

;****************************************************************************
; Open the file again.
;****************************************************************************

ss = ' '
openr, unit, filename, /get_lun
for i=0,2 do readf, unit, ss

;****************************************************************************
; Read the band and g information.
;****************************************************************************

waveno = fltarr(numb,2)
solf = fltarr(numb)
delg = fltarr(numb,maxg)
for b=0,numb-1 do begin
    dg = fltarr(numg(b))
    readf, unit, n, wn0, wn1, sf, ng, dg
    waveno(b,0) = wn0
    waveno(b,1) = wn1
    solf(b) = sf
    numg(b) = ng
    delg(b,0:numg(b)-1) = dg
  endfor

;****************************************************************************
; Read the number of levels.
;****************************************************************************

readf, unit, ss
i = strpos(ss,'!')
ss = strtrim(strmid(ss,0,i),2)
numz = long(ss)

;****************************************************************************
; Read the gas concentration information.
;****************************************************************************

readf, unit, ss
i = strpos(ss,'!')
ss = strtrim(strmid(ss,0,i),2)
i = strpos(ss,' ')
co2 = float(strtrim(strmid(ss,0,i),2))
ss = strtrim(strmid(ss,i,strlen(ss)-i),2)
i = strpos(ss,' ')
ch4 = float(strtrim(strmid(ss,0,i),2))
n2o = float(strtrim(strmid(ss,i+1,strlen(ss)-i-1),2))

;****************************************************************************
; Read the profile information.
;****************************************************************************

z = fltarr(numz)
pres = fltarr(numz)
temp = fltarr(numz)
qh2o = fltarr(numz)
qo3 = fltarr(numz)
readf, unit, ss
for k=0,numz-1 do begin
    readf, unit, h, p, t, qh, qo
    z(k) = h
    pres(k) = p
    temp(k) = t
    qh2o(k) = qh
    qo3(k) = qo
  endfor

;****************************************************************************
; Read the absorption information.
;****************************************************************************

absorb = fltarr(numb,maxg,numz)
readf, unit, ss
for b=0,numb-1 do begin
    a = fltarr(numg(b))
    for k=0,numz-1 do begin
        readf, unit, n1, n2, a
        absorb(b,0:numg(b)-1,k) = a
      endfor
  endfor

;****************************************************************************
; Close the file.
;****************************************************************************

close, unit
free_lun, unit

END
