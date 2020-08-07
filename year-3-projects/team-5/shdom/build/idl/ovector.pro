PRO OVECTOR,U,V,XX,YY, $
	Missing=missing,Length=length,Standard=standard,Hsize=hsize, $
	Dots=dots,Zero=zero,color=color,Notail=notail,NoAspect=noAspect
;
;+
; NAME:
;	OVECTOR
;
; PURPOSE:
;	Produce a two-dimensional velocity field overplot.
;
;	A directed arrow is drawn at each point showing the direction and 
;	magnitude of the field.
;       
; CATEGORY:
;	Plotting, two-dimensional.
;
; CALLING SEQUENCE:
;	OVECTOR, U, V [, X, Y]
;
; INPUTS:
;	U:	The X component of the two-dimensional field.  
;		U must be a two-dimensional array.
;
;	V:	The Y component of the two dimensional field.  Y must have
;		the same dimensions as X.  The vector at point (i,j) has a 
;		magnitude of:
;
;			(U(i,j)^2 + V(i,j)^2)^0.5
;
;		and a direction of:
;
;			ATAN2(V(i,j),U(i,j)).
;
; OPTIONAL INPUT PARAMETERS:
; 	X:	Optional abcissae values. X must be a vector with a length 
;		equal to the first dimension of U and V.
;
;	Y:	Optional ordinate values. Y must be a vector with a length
;		equal to the first dimension of U and V.
;
;   MISSING: Missing data value. Vectors with a LENGTH greater
;		than MISSING are ignored.
;
;   LENGTH: The vector length as fraction of Max(magnitude)
;
;   STANDARD: The standard length of a vector, instead of Max(magnitude).
;             This is useful when plotting several graphs, each of
;             which should have vectors of the same standard length.
;
;	HSIZE:	Size of arrowhead, default is 0.2 of length.
;
;	COLOR:	The color used for the vectors.
;
; KEYWORD INPUT PARAMETERS:
;
;	DOTS:	Set this keyword to place a dot at each missing point. 
;           Has effect only IF MISSING is specified.
;
;	ZERO:	Set this keyword to place a dot at each point
;           with magnitude 0.
;
;	NOTAIL: Set this keyword to not draw any tail, just the arrowhead.
;           (note that HSIZE may then need to be adjusted)
;
;   NOASPECT: Do not compute the aspect ratio. Normally, OVECTOR takes
;             the ratio of plot window size to data coordinate 
;             window size and uses this to scale the V component,
;             so the scale in the V-direction is the same in the U.
;
; OUTPUTS:
;	None.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	Plotting on the selected device is performed.  System
;	variables concerning plotting are changed.
;
; RESTRICTIONS:
;	None.
;
; MODIFICATION HISTORY:
;	Created 2/16/93  C. Torrence
;	
;-
;
on_error,2      ;Return to caller IF an error occurs
s = size(u)
t = size(v)
x1 = N_ELEMENTS(xx)
y1 = N_ELEMENTS(yy)
mess0='U and V parameters must be 2D and same size.'
mess1='U('+StrCompress(s(1))+','+StrCompress(s(2))+')'
mess2=',   V('+StrCompress(t(1))+','+StrCompress(t(2))+')'
mess3=',   X('+StrCompress(x1)+')'
mess4=',   Y('+StrCompress(y1)+')'

IF s(0) ne 2 THEN begin 
	baduv:   MESSAGE,mess0,/inform
			 MESSAGE,mess1+mess2
ENDIF
IF total(abs(s(0:2)-t(0:2))) NE 0 THEN goto,baduv

IF N_PARAMS(0) lt 3 THEN x = FINDGEN(s(1)) $
	ELSE IF x1 NE s(1) THEN $
		  MESSAGE,'X array has incorrect size. '+mess1+mess3 $
		ELSE x=Float(xx)
IF N_PARAMS(1) lt 4 THEN y = FINDGEN(s(2)) $
	ELSE IF y1 NE s(2) THEN $
		  MESSAGE,'Y array has incorrect size. '+mess1+mess4 $
		ELSE y=Float(yy)
IF KEYWORD_SET(noAspect) EQ 0 THEN noAspect=0
IF N_ELEMENTS(missing) LE 0 THEN missing=1E30
IF N_ELEMENTS(length) LE 0 THEN length = 1.
IF N_ELEMENTS(hsize) LE 0 THEN hsize = 0.2
IF KEYWORD_SET(notail) LE 0 THEN notail=0
IF KEYWORD_SET(zero) THEN zero=1 ELSE zero=0
IF N_ELEMENTS(color) eq 0 THEN color = !P.COLOR

x0 = min(x) & x1 = max(x)		;get scaling
y0 = min(y) & y1 = max(y)

vv=v
IF NOT noAspect THEN BEGIN
	xw=!X.WINDOW
	yw=!Y.WINDOW
	aspectratio=(x1-x0)/(y1-y0)*(yw(1)-yw(0))/(xw(1)-xw(0))
	vv=vv/aspectratio
ENDIF

mag = SQRT(u^2+vv^2)     	;magnitude.
nbad = 0	;# of missing points
IF zero THEN zeroes=WHERE(mag EQ 0)
good = WHERE((mag NE 0) AND (mag LT missing))
IF KEYWORD_SET(dots) THEN bad = WHERE(mag GE missing, nbad)

mag = mag(good) 	;Discard missing values
ugood = u(good)
vgood = vv(good)
IF N_ELEMENTS(standard) LE 0 THEN standard = MAX(mag)

dx = (x1-x0)/s(1)/standard*length*ugood/2. 	;sin & cos comps.
dy = (y1-y0)/s(2)/standard*length*vgood/2.


angle=20*!Dtor
theta = ATAN(dy,dx)
cosa = COS(theta)
sina = SIN(theta)
dd=MAX(SQRT(dx^2+dy^2))
st = hsize*dd*SIN(angle)		;sin angle * length of head
ct = hsize*dd*COS(angle)

IF zero THEN PLOTS,x(zeroes MOD s(1)),y(zeroes/s(1)),psym=3,color=color

CASE notail OF

0: BEGIN
	x0 = x(good MOD s(1)) - dx		;get coords of start & end
	x1 = x0 + 2*dx
	y0 = y(good/s(1)) - dy
	y1 = y0 + 2*dy
	x2 = x1-ct*cosa+st*sina
	x3 = x1-ct*cosa-st*sina
	y2 = y1-ct*sina-st*cosa
	y3 = y1-ct*sina+st*cosa

	FOR i=0,N_ELEMENTS(good)-1L DO BEGIN	;Each point
		PLOTS,[x0(i),x1(i),x2(i),x1(i),x3(i)], $
			[y0(i),y1(i),y2(i),y1(i),y3(i)],color=color
; ** uncomment next line for filled arrow heads
;		POLYFILL,[x1(i),x2(i),x3(i)],[y1(i),y2(i),y3(i)],color=color
	ENDFOR
	IF nbad GT 0 THEN $     ;Dots for missing?
		OPLOT,x(bad mod s(1)),y(bad/s(1)),psym=3,color=color
   END

1: BEGIN
	x0 = x(good MOD s(1))		;get coords of start & end
	y0 = y(good/s(1))
	x2 = x0-ct*cosa+st*sina
	x3 = x0-ct*cosa-st*sina
	y2 = y0-ct*sina-st*cosa
	y3 = y0-ct*sina+st*cosa

	FOR i=0,N_ELEMENTS(good)-1L DO $     ;Each point
		PLOTS,[x2(i),x0(i),x3(i)], $
			[y2(i),y0(i),y3(i)],color=color
	IF nbad GT 0 THEN $     ;Dots for missing?
		OPLOT,x(bad mod s(1)),y(bad/s(1)),psym=3,color=color
   END

ENDCASE

END
