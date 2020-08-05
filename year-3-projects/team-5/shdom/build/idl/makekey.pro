PRO MAKEKEY,x0,y0,xsize,ysize,xlaboff,ylaboff, $
	COLORS=colors, LABELS=labels, ORIENTATION=orientation, $
	NOBORDER=noborder, BCOLOR=bcolor, THICK=thick, $
	PATTERN=pattern, LINESPACING=linespacing, LINEANGLE=lineangle, $
	ALIGNMENT=alignment, CHARSIZE=charsize, CHARTHICK=charthick
;+
; NAME:
; MAKEKEY
;
; PURPOSE:
; Create a color key for a plot.
;
; CATEGORY:
;
; CALLING SEQUENCE:
;
; INPUTS:
; X0:  		x position of lower left corner of key.
; Y0:  		y position of lower left corner of key.
; XSIZE:	x size of key.
; YSIZE:	y size of key.
; XLABOFF:	x label offset (relative to left side of box).
; YLABOFF:	y label offset (relative to bottom side of box).
;	Positions/sizes are in normalized coordinates
;
; OUTPUTS:
; None.
;
; OPTIONAL INPUT/OUTPUT PARAMETERS:
;
; COLORS: array of colors for each box 
; LABELS: array of strings for labels (can have 1 to Nbox+1 labels)
; ORIENTATION: orientation of key: 0=left to right(default), 1=top to bottom
; NOBORDER: don't put a border around each box
; BCOLOR: index (or scalar) of colors for border (default=!P.color)
; THICK: thickness of border lines (default=!P.thick)
; PATTERN: a MxMxN array of patterns, N=index, MxM is pattern
; LINESPACING: array of line spacings (in cm) for box fill lines 
; LINEANGLE: array of orientation angles of box fill lines
; ALIGNMENT: label justify(0=left justify,0.5=center(default),1=right)
; CHARSIZE: size of labels (default=!P.charsize)
; CHARTHICK: thickness of vector drawn characters (default=!P.charthick)
;	If LINESPACING or LINEANGLE are specified then line fills are done.
;	If PATTERN is specified then pattern fills are done.
;	Otherwise solid fills using COLORS are done.
;	Number of boxes is determined from COLORS, LINE*, or PATTERN.
;	The number of labels is usually one more than the number of
;	  boxes, but can be less to label boxes less frequently.
;
; EXAMPLE:
; 	MAKEKEY, 0.15,0.05, 0.75,0.04, 0.0,-0.02, $
;		LABELS=STRING(4*INDGEN(17)), COLORS=colora
;
;  
; COMMON BLOCKS:
; None.
;
; SIDE EFFECTS:
; None.
;
; MODIFICATION HISTORY:
; Written, C. Torrence, Nov. 9, 1993.
; Modified, Frank Evans, CU, 2/24/94

;-
;	ON_ERROR,2  ;return to caller if error
	x0=FLOAT(x0)
	y0=FLOAT(y0)
	xsize=FLOAT(xsize)
	ysize=FLOAT(ysize)


	nbox = N_ELEMENTS(colors)
	patflag = (N_ELEMENTS(pattern) GT 0)
	IF patflag THEN BEGIN
	    tmp = SIZE(pattern)
	    IF tmp(0) EQ 3 THEN nbox = tmp(3)
	ENDIF
        lineflag = (N_ELEMENTS(linespacing) GT 0) $
		OR (N_ELEMENTS(lineangle) GT 0)
	if lineflag THEN $
	     nbox = N_ELEMENTS(lineangle) < N_ELEMENTS(linespacing)

	IF nbox EQ 0 THEN BEGIN
	    PRINT, 'MAKEKEY: Zero boxes specified.'
	    RETURN
	ENDIF
	IF N_ELEMENTS(colors) EQ 0 THEN	colors=REPLICATE(!P.color,nbox)
	IF N_ELEMENTS(colors) EQ 1 THEN	colors=REPLICATE(colors(0),nbox)

	nlab = N_ELEMENTS(labels)
	IF nlab LE 0 THEN labels=REPLICATE(nbox,'') $
		ELSE labels=STRCOMPRESS(labels,/rem)
        steplab = FIX(nbox/(nlab-1)) > 1

	IF N_ELEMENTS(noborder) LE 0 THEN	noborder=0
	IF N_ELEMENTS(bcolor) LE 0 THEN		bcolor=!P.color
	IF N_ELEMENTS(bcolor) EQ 1 THEN		bcolor=BYTARR(nbox)+bcolor
	IF N_ELEMENTS(orientation) LE 0 THEN	orientation=0
	IF N_ELEMENTS(charsize) LE 0 THEN	charsize=!P.charsize
	IF N_ELEMENTS(alignment) LE 0 THEN	alignment=0.5
	IF N_ELEMENTS(thick) LE 0 THEN		thick=!P.thick
	IF N_ELEMENTS(charthick) LE 0 THEN	charthick=!P.charthick
	CASE orientation OF
	0:BEGIN
            xbox = xsize/nbox
            ybox = ysize
	    x = x0 + FINDGEN(nbox+1)*xbox
	    y = y0 + FLTARR(nbox+1)
	  END
	1:BEGIN
            xbox = xsize
            ybox = ysize/nbox
	    x = x0 + FLTARR(nbox+1)
	    y = y0 + FINDGEN(nbox+1)*ybox
	  END
	ENDCASE
	xl = x + xlaboff
	yl = y + ylaboff


;  Make the boxes and draw lines around them (if desired)
	FOR i = 0, nbox-1 DO BEGIN
	    IF patflag THEN $
		POLYFILL, [x(i),x(i)+xbox,x(i)+xbox,x(i)], $
		    [y(i),y(i),y(i)+ybox,y(i)+ybox], /NORMAL, $
		    PATTERN=pattern(*,*,i) $
	    ELSE IF lineflag THEN $
		POLYFILL, [x(i),x(i)+xbox,x(i)+xbox,x(i)], $
		    [y(i),y(i),y(i)+ybox,y(i)+ybox], /NORMAL, $
		    /LINE_FILL, SPACING=linespacing(i), $
		    ORIENTATION=lineangle(i), COLOR=colors(i) $
	    ELSE $
		POLYFILL, [x(i),x(i)+xbox,x(i)+xbox,x(i)], $
		    [y(i),y(i),y(i)+ybox,y(i)+ybox], /NORMAL, $
		    COLOR=colors(i)
	    IF noborder EQ 0 THEN $
		PLOTS,[x(i),x(i)+xbox,x(i)+xbox,x(i),x(i)],  $
		      [y(i),y(i),y(i)+ybox,y(i)+ybox,y(i)], /NORMAL, $
		    COLOR=bcolor(i), /NOCLIP,THICK=thick
	ENDFOR

;  Label the boxes
	FOR i = 0, nbox DO BEGIN
	    IF (i MOD steplab) EQ 0 THEN BEGIN 
		j = FIX(i/steplab)
		IF j LT nlab THEN $
		    XYOUTS,xl(i),yl(i),labels(j), /NORMAL, $
			ALIGNMENT=alignment, $
			CHARSIZE=charsize, CHARTHICK=charthick
	    ENDIF
	ENDFOR

END

