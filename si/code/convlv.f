C  (C) Copr. 1986-92 Numerical Recipes Software ,2:-5K#R..
      SUBROUTINE convlv(data,respns,n,isign,ans,fft)
CU    USES realft,twofft
      INTEGER*4 isign,n
      INTEGER*4 i,no2
      REAL*4    data(n),respns(n)
      COMPLEX*8 ans(n)
      COMPLEX*8 fft(n)
      call twofft(data,respns,fft,ans,n)

      no2=n/2
      if (isign.eq.1) then               ! convolution
        do i=1,no2+1
          ans(i)=fft(i)*ans(i)/no2
        end do      
      else if (isign.eq.-1) then         ! deconvolution
        do i=1,no2+1
           if (abs(ans(i)).eq.0.0) pause
     *        'deconvolving at response zero in convlv'
           ans(i)=fft(i)/ans(i)/no2
        end do
      else
          pause 'no meaning for isign in convlv'
      endif
      ans(1)=cmplx(real(ans(1)),real(ans(no2+1)))
      call realft(ans,n,-1)
      return
      END
