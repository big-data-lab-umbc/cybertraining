;  IDL program for animating JPEG images.

   ; The base names of the input JPEG images
infilebase = 'stcu067av'

nlo=480  & nso=640  ; number of lines and samples in images
Nimage=100          ; number of images
Simage=1            ; starting image


XINTERANIMATE, SET=[nso,nlo,Nimage] , /SHOWLOAD

j=0

for i = Simage, Simage+Nimage-1 do begin
  infile = infilebase+string(i,format='(i3.3)')+'.jpeg'
  READ_JPEG, infile, bytimg

  XINTERANIMATE, FRAME=j, IMAGE=bytimg
  j = j + 1
endfor

XINTERANIMATE, /KEEP_PIXMAPS

end
