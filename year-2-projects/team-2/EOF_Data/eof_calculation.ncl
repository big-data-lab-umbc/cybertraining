;vname="sst"
neof=30
optEOF=False
optETS=False
optEVX=-1

fils=systemfunc("ls *.nc")
f=addfiles(fils,"r")
ListSetType(f,"cat")
var=short2flt(f[:]->$vname$)
var1=rmMonAnnCycTLL(var)
var2=dtrend_msg_n(var&time,var1,True,False,0)
printVarSummary(var2)
rad=4.0*atan(1.0)/180.0
clat=cos(rad*var&latitude)
clat=where(clat.lt.0,0.0,clat)
clat=sqrt(clat)
printVarSummary(clat)
wvar=var2*conform(var2,clat,1)
copy_VarCoords(var,wvar)
printVarSummary(wvar)

eof=eofunc_n(wvar,neof,optEOF,0)
eof_rot=eofunc_varimax(eof,optEVX)
eofunc_varimax_reorder(eof_rot)
eof_ts=eofunc_ts_n(wvar,eof_rot,optETS,0)

do gg=0,neof-1
eof_ts(gg,:)=dim_standardize(eof_ts(gg,:),0)
asciiwrite("eof_"+vname+"-"+sprinti("%0.2i",gg),eof_rot(gg,:,:))
asciiwrite("eof_ts_"+vname+"-"+sprinti("%0.2i",gg),eof_ts(gg,:))
end do
asciiwrite("eof_pcvar_"+vname,eof_rot@pcvar_varimax)
