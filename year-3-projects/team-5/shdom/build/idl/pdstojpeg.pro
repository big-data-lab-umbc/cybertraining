;  IDL program for converting SHDOM PDS format images to JPEG.

   ; The base names of the input PDS images and the output JPEG images
infilebase = 'stcu067av'
outfilebase = 'stcu067av'

nl=240  & ns=320    ; number of lines and samples in PDS format images
nbytes=1            ; number of bytes per pixel (1 or 2)
swap=0              ; 1 for swapping bytes of integers, 0 otherwise
labrecs=7           ; number of records PDS header occupies (given in header)
pixscale=2          ; image magnify scale (integer)
bytescaling=1.0     ; scaling factor for integer to byte pixels
jpegquality=90      ; JPEG output quality
Nimage=100          ; number of images to process
Spds=1              ; starting PDS image number
Sjpeg=1             ; starting JPEG image number


if (nbytes eq 2) then begin
  buf = INTARR(ns,nl+labrecs)
endif else begin
  buf = BYTARR(ns,nl+labrecs)
endelse

nlo = pixscale*nl
nso = pixscale*ns


for i = Spds, Nimage-1+Spds do begin
  infile = infilebase+string(i,format='(I3.3)')+'.pds'
  j = i-Spds+Sjpeg
  outfile = outfilebase+string(j,format='(I3.3)')+'.jpeg'

  OPENR,1, infile
  READU,1, buf
  CLOSE,1

  img = buf[*,labrecs:nl+labrecs-1]
  if (swap eq 1) then img = SWAP_ENDIAN(img)
  img = REBIN(img,nso,nlo)
  bytimg = BYTE((bytescaling*img) < 255)

  WRITE_JPEG, outfile, bytimg, quality=jpegquality, order=1
endfor

end
