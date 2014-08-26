
            
C --- Voigt profile function ---------------------------------------------------
      REAL*8 FUNCTION VOIGT(a,v)        

      REAL*8    a           
      REAL*8    aa        
      REAL*8    v      
      REAL*8    vv,h0,h1,h2,hh1,hh2,hh3,hh4,div_vv,u1,uu

      vv= v * v
      IF ((a.LT.0.2D0) .AND. (v.GT.5.0D0))                     GOTO 70
      IF ((a.GE.0.2D0) .AND. (a.GT.1.4D0 .OR. (a+v).GT.3.2D0)) GOTO 80
      h0=EXP(-vv)
      h2=(1.D0-2.D0*vv)*h0
      IF (v.GT.2.4D0) GOTO 30
      IF (v.LE.1.3D0) GOTO 20
      h1=(-.220416D0*vv+1.989196D0*v-6.61487D0)*vv+9.39456D0*v-4.4848D0
      GOTO 40
20    h1=(.42139D0*vv-2.34358D0*v+3.28868D0)*vv-.15517D0*v-1.1247D0
      GOTO 40
30    h1=((-.0032783D0*vv+.0429913D0*v-.188326D0)*vv+.278712D0*v+
     *   .55415D0)/(vv-1.5D0)
40    CONTINUE
      IF (a.GE..2D0) GOTO 60
50    voigt=(h2*a+h1)*a+h0
      RETURN
60    hh1=h1+h0*1.12838D0
      hh2=h2+hh1*1.12838D0-h0
      hh3=(1.D0-h2)*0.37613D0-hh1*0.66667D0*vv+hh2*1.12838D0
      hh4=(3.D0*hh3-hh1)*0.37613D0+h0*0.66667D0*vv*vv
      voigt=((((hh4*a+hh3)*a+hh2)*a+hh1)*a+h0)*
     *      (((-.122727278D0*a+.532770573D0)*a-.962846325D0)*a+
     *         .979895023D0)
      RETURN
70    div_vv = 1.D0 / vv
      voigt=((2.12D0*div_vv+.8463D0)*div_vv+.5642D0)*a*div_vv       
      RETURN
80    aa=a*a
      u1=(aa+vv)*1.4142D0
      uu=u1*u1
      voigt=((((aa-10.D0*vv)*aa*3.D0+15.D0*vv*vv)/uu+3.D0*vv-aa)
     *      /uu+1.D0) * .79788D0*a/u1
      RETURN
      END
