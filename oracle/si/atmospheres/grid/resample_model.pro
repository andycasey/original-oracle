
readcol, 'LIST', atm, format='a'
natm = n_elements(atm)
ttemp=intarr(natm)
gtemp=fltarr(natm)
mtemp=fltarr(natm)

for i = 0, natm-1 do atm(i) = streplace(atm(i),'bin_','')
for i = 0, natm-1 do begin
 txt = strsplit(atm(i),'_',/extract)
 ttemp(i) = fix(strmid(txt[0],1)) 
 gtemp(i) =      float(txt[1])
 mtemp(i) =      float(txt[2])
endfor
for i = 0, natm-1 do atm(i) = streplace(atm(i), '_80','')
ttemp = ttemp(sort(ttemp)) & ind = uniq(ttemp) & tnew = ttemp(ind)
gtemp = gtemp(sort(gtemp)) & ind = uniq(gtemp) & gnew = gtemp(ind)
mtemp = mtemp(sort(mtemp)) & ind = uniq(mtemp) & mnew = mtemp(ind)

path = '~/siu/lineform_fug/code_pc/'
openr, lun, path + 'atmgrid.inc',/get_lun

header1 = strarr(13)
header2 = strarr(2)
header3 = strarr(2)
readf, lun, header1, format='(a80)'

txt    = extr_words(header1[4])
 txt[3] = streplace(txt[3],')','') & nteff  = fix(txt[3])
txt    = extr_words(header1[5])
 txt[3] = streplace(txt[3],')','') & nlogg  = fix(txt[3]) 
txt    = extr_words(header1[6])
 txt[3] = streplace(txt[3],')','') & nmeta  = fix(txt[3]) 
ntot = nteff*nlogg*nmeta

teff   = fltarr(nteff) & tind = teff
logg   = fltarr(nlogg) & gind = logg
meta   = fltarr(nmeta) & mind = meta

ft = '(6x,14x,i3,4x,f7.1,a1)'
fg = '(6x,14x,i3,4x,f6.3,a1)'
fm = '(6x,14x,i3,4x,f6.3,a1)'
fa = '(6x,18x,3(i3,1x),4x,a25)'

eps = 1.d-4

a = 0 & b = 0. & c = 0 &  tfile = ''
for i = 0, nteff-1 do begin 
  readf, lun, a, b, format=ft
  tind(i) = a &  teff(i) = b
endfor
readf, lun, header2, format='(a80)'
for i = 0, nlogg-1 do begin
  readf, lun, a, b, format=fg
  gind(i) = a &  logg(i) = b
endfor
readf, lun, header2, format='(a80)'
for i = 0, nmeta-1 do begin
  readf, lun, a, b, format=fm
  mind(i) = a &  meta(i) = b
endfor
readf, lun, header3, format='(a80)'

i1 = intarr(ntot)  & t1=intarr(ntot)
i2 = intarr(ntot)  & g1=fltarr(ntot)
i3 = intarr(ntot)  & m1=fltarr(ntot)
model = strarr(ntot)
for i = 0, ntot-1 do begin
  readf, lun, a, b, c, tfile, format=fa
  i1(i) = a &  i2(i) = b & i3(i) = c & model(i) = tfile
  t1(i) = teff(a-1)
  g1(i) = logg(b-1)
  m1(i) = meta(c-1)
;  print, t1(i), g1(i), m1(i),' ', model(i) 
endfor

t = [teff, tnew]
g = [logg, gnew]
m = [meta, mnew]
t=t(sort(t)) & g=g(sort(g)) & m=m(sort(m))
ind = uniq(t) & t   = t(ind)
ind = uniq(g) & g   = g(ind)
ind = uniq(m) & m   = m(ind)

nt = fix(n_elements(t))
ng = fix(n_elements(g))
nm = fix(n_elements(m))
print, t,g,m
nnew = nt*ng*nm

modarr = strarr(nt,ng,nm)

str1='C --- INCLUDE FILE ATMGRID.INC FOR LINEFORM.F ---------'
str2='C --- This Tables provides all accessible atmosphere grid points'
str3='C'
str4='INTEGER*4 NTEFFS,NLOGGS,NZS'
str5a='PARAMETER (NTEFFS ='
str5b=') ! Number of Teff-gridpoints'
str6a='PARAMETER (NLOGGS ='
str6b=') ! Number of log(g)-gridpoints'
str7a='      PARAMETER (NZS    ='
str7b=') ! Number of metallicity-gridpoints'
str8 ='REAL*4 atm_teff(NTEFFS)'
str9 ='REAL*4 atm_logg(NLOGGS)'
str10='REAL*4 atm_z   (NZS)'
str11='CHARACTER*25 atm_filename (NTEFFS,NLOGGS,NZS)'
str12='C --- AVAILABLE EFFECTIVE TEMPERATURES ----------------'
str13='C --- AVAILABLE GRAVITIES -----------------------------'
str14='C --- AVAILABLE METALLICITIES -------------------------'
str15='C --- AVAILABLE ATMOSPHERES ---------------------------'
s16='DATA atm_teff('
s17='DATA atm_logg('
s18='DATA atm_z   ('
s19='DATA atm_filename('
s1=')  /'
s2='/'
s3=','
s4="'"

listnm = '/afs/mpa/data/mbergema/archive/atmos/grid/LIST'
 readcol, listnm, newmodels, form='a'
nn = n_elements(newmodels)
t2=fltarr(nn)
g2=fltarr(nn)
m2=fltarr(nn)
newname=strarr(nn)
for i=0, nn-1 do begin
  getnum, newmodels(i), pars, ct
  t2(i) = pars(0)
  g2(i) = pars(1)
  m2(i) = pars(2)
  if m2(i) le 0. then charm='-' else charm = 'p'
  if g2(i) lt 0. then charg='-' else charg = ''

  if abs(g2(i)) lt 1. then ch1 = '0' else ch1 = ''
  if abs(m2(i)) lt 1. then ch2 = '0' else ch2 = ''
       if g2(i) eq 0. then ch3 = '0' else ch3  = ''
       if m2(i) eq 0. then ch4 = '0' else ch4  = ''

  temp  = strtrim(string(fix((t2(i)/10.))),2)
  grav  = strtrim(string(fix(abs(g2(i))*100)),2)
  meta  = strtrim(string(fix(abs(m2(i))*100)),2)
  if m2(i) eq 0. then meta = '0'

  newname(i) = 't'+ temp + charg + ch1 + ch3 + grav + charm + ch2 + ch4 + meta + '.dat'
  newname(i) = strtrim(newname(i),2)
  print, temp, ' ', grav, ' ', meta, '   ', newname(i)
  spawn, 'mv '+ newmodels(i) + ' ' + newname(i)
endfor

openw, lun,'/afs/mpa/home/mbergema/siu/lineform_fug/code_pc/atmgrid_new.inc', /get_lun
;for i = 0, nnew -1 do begin
    printf, lun, str1
    printf, lun, str2
    printf, lun, str3
    printf, lun, str4, form='(6x,a-30)'
;
    printf, lun, str5a, strtrim(string(nt),2), str5b, form='(6x,a-18,a3,a30)'
    printf, lun, str6a, strtrim(string(ng),2), str6b, form='(6x,a-18,a3,a30)'
    printf, lun, str7a, strtrim(string(nm),2), str7b, form='(6x,a-18,a3,a30)'

    printf, lun, str8,  form='(6x,a-29)'
    printf, lun, str9,  form='(6x,a-29)'
    printf, lun, str10, form='(6x,a-29)'
    printf, lun, str11, form='(6x,a-50)'
    printf, lun, str3

    printf, lun, str12, form='(a-50)'
    for i=0, nt-1 do printf, lun, s16, i+1, s1, t(i), s2, $
               form='(6x,a-14,i3,a4,f7.1,a1)'
    printf, lun, str3

    printf, lun, str13, form='(a-50)'
    for i=0, ng-1 do printf, lun, s17, i+1, s1, g(i), s2, $
               form='(6x,a-14,i3,a4,f7.3,a1)'
    printf, lun, str3

    printf, lun, str14, form='(a-50)'
    for i=0, nm-1 do printf, lun, s18, i+1, s1, string(m(i)), s2, $
               form='(6x,a-14,i3,a4,f7.3,a1)'
    printf, lun, str3

    printf, lun, str15, form='(a-50)'
    for i=0, nt-1 do begin
     for j=0, ng-1 do begin
      for k=0, nm-1 do begin

        stmod = ''
        ind = where(t(i) eq t1 and g(j) eq g1 and m(k) eq m1, nind)    ; old models
        if nind eq 1 then stmod=model(ind)
        ind = where(t(i) eq t2 and g(j) eq g2 and m(k) eq m2, nind)    ; new models
        if nind eq 1 then stmod=newname(ind) 

        printf, lun, s19, i+1, s3, j+1,s3, k+1,s1,s4,stmod,s4,s2, $
               form='(6x,a-18,i3,a1,i3,a1,i3,a4,a1,a-25,a1,a1)'
      endfor
     endfor
    endfor

close, lun
; endfor


end
