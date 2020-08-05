; IDL procedure  -  vertical_integral.pro
;
;****************************************************************************
;
; NAME:
;      vertical_integral
;
; PURPOSE:
;      Computes the vertical integral of a field
;
; CATEGORY:
;      SHDOM utility
;
; CALLING SEQUENCE:
;      intfield = vertical_integral(zgrid, field, dim)
;
; INPUTS:
;      zgrid:     An array of nx heights at which the field is defined.
;      field:     An array (2D or 3D)
;      dim:       Integer giving the dimension that is vertical (e.g. 1,2,3)
;
; OUTPUTS:
;      intfield:  The output vertically integrated field (1D or 2D)
;
; MODIFICATION HISTORY:
;      Written by:   Frank Evans, May 1999
;
; NOTES:
;      Uses trapezoidal integration.
;
;****************************************************************************

FUNCTION vertical_integral, zgrid, field, dim

;****************************************************************************

fsize = size(field)
ndim = fsize(0)
if (ndim ne 2 and ndim ne 3) then $
  message, 'Input field must be 2 or 3 dimensions.'
n1 = fsize(1)
n2 = fsize(2)
n3 = fsize(3)

if (ndim eq 2) then begin
  if (dim eq 1) then begin
    intfield = reform(0.0*field(0,*),n2)
    for i = 0, n_elements(zgrid)-2 do $
      intfield=intfield + (zgrid(i+1)-zgrid(i))*0.5*(field(i,*)+field(i+1,*))
  endif else begin
    intfield = reform(0.0*field(*,0),n1)
    for i = 0, n_elements(zgrid)-2 do $
      intfield=intfield + (zgrid(i+1)-zgrid(i))*0.5*(field(*,i)+field(*,i+1))
  endelse
endif else begin
  if (dim eq 1) then begin
    intfield = reform(0.0*field(0,*,*),n2,n3)
    for i = 0, n_elements(zgrid)-2 do $
      intfield=intfield+ (zgrid(i+1)-zgrid(i))*0.5*(field(i,*,*)+field(i+1,*,*))
  endif else begin
    intfield = reform(0.0*field(*,*,0),n1,n2)
    for i = 0, n_elements(zgrid)-2 do $
      intfield=intfield+ (zgrid(i+1)-zgrid(i))*0.5*(field(*,*,i)+field(*,*,i+1))
  endelse
endelse

RETURN, intfield

END

