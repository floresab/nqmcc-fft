      subroutine costi(n,wsave)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (pi=3.14159265358979d0)
c***begin prologue  costi
c***date written   790601   (yymmdd)
c***revision date  830401   (yymmdd)
c***category no.  j1a3
c***keywords  fourier transform
c***author  swarztrauber, p. n., (ncar)
c***purpose  initialize for cost.
c***description
c
c  subroutine costi initializes the array wsave which is used in
c  subroutine cost.  the prime factorization of n together with
c  a tabulation of the trigonometric functions are computed and
c  stored in wsave.
c
c  input parameter
c
c  n       the length of the sequence to be transformed.  the method
c          is most efficient when n-1 is a product of small primes.
c
c  output parameter
c
c  wsave   a work array which must be dimensioned at least 3*n+15.
c          different wsave arrays are required for different values
c          of n.  the contents of wsave must not be changed between
c          calls of cost.
c***references  (none)
c***routines called  rffti
c***end prologue  costi
      dimension       wsave(*)
c***first executable statement  costi
      if (n .le. 3) return
      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      dt = pi/nm1
      fk = zero
      do 101 k=2,ns2
         kc = np1-k
         fk = fk+one
         wsave(k) = two*sin(fk*dt)
         wsave(kc) = two*cos(fk*dt)
  101 continue
      call rffti (nm1,wsave(n+1))
      return
      end
      subroutine cost(n,x,wsave)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c***begin prologue  cost
c***date written   790601   (yymmdd)
c***revision date  851219   (yymmdd)
c***category no.  j1a3
c***keywords  fourier transform
c***author  swarztrauber, p. n., (ncar)
c***purpose  cosine transform of a real, even sequence.
c***description
c
c  subroutine cost computes the discrete fourier cosine transform
c  of an even sequence x(i).  the transform is defined below at output
c  parameter x.
c
c  cost is the unnormalized inverse of itself since a call of cost
c  followed by another call of cost will multiply the input sequence
c  x by 2*(n-1).  the transform is defined below at output parameter x.
c
c  the array wsave which is used by subroutine cost must be
c  initialized by calling subroutine costi(n,wsave).
c
c  input parameters
c
c  n       the length of the sequence x.  n must be greater than 1.
c          the method is most efficient when n-1 is a product of
c          small primes.
c
c  x       an array which contains the sequence to be transformed
c
c  wsave   a work array which must be dimensioned at least 3*n+15
c          in the program that calls cost.  the wsave array must be
c          initialized by calling subroutine costi(n,wsave), and a
c          different wsave array must be used for each different
c          value of n.  this initialization does not have to be
c          repeated so long as n remains unchanged.  thus subsequent
c          transforms can be obtained faster than the first.
c
c  output parameters
c
c  x       for i=1,...,n
c
c             x(i) = x(1)+(-1)**(i-1)*x(n)
c
c               + the sum from k=2 to k=n-1
c
c                 2*x(k)*cos((k-1)*(i-1)*pi/(n-1))
c
c               a call of cost followed by another call of
c               cost will multiply the sequence x by 2*(n-1).
c               hence cost is the unnormalized inverse
c               of itself.
c
c  wsave   contains initialization calculations which must not be
c          destroyed between calls of cost.
c***references  (none)
c***routines called  rfftf
c***end prologue  cost
      dimension       x(*)       ,wsave(*)
c***first executable statement  cost
      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      if (n-2) 106,101,102
  101 x1h = x(1)+x(2)
      x(2) = x(1)-x(2)
      x(1) = x1h
      return
  102 if (n .gt. 3) go to 103
      x1p3 = x(1)+x(3)
      tx2 = x(2)+x(2)
      x(2) = x(1)-x(3)
      x(1) = x1p3+tx2
      x(3) = x1p3-tx2
      return
  103 c1 = x(1)-x(n)
      x(1) = x(1)+x(n)
      do 104 k=2,ns2
         kc = np1-k
         t1 = x(k)+x(kc)
         t2 = x(k)-x(kc)
         c1 = c1+wsave(kc)*t2
         t2 = wsave(k)*t2
         x(k) = t1-t2
         x(kc) = t1+t2
  104 continue
      modn = mod(n,2)
      if (modn .ne. 0) x(ns2+1) = x(ns2+1)+x(ns2+1)
      call rfftf (nm1,x,wsave(n+1))
      xim2 = x(2)
      x(2) = c1
      do 105 i=4,n,2
         xi = x(i)
         x(i) = x(i-2)-x(i-1)
         x(i-1) = xim2
         xim2 = xi
  105 continue
      if (modn .ne. 0) x(n) = xim2
  106 return
      end
      subroutine sinti(n,wsave)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (pi=3.14159265358979d0)
c***begin prologue  sinti
c***date written   790601   (yymmdd)
c***revision date  830401   (yymmdd)
c***category no.  j1a3
c***keywords  fourier transform
c***author  swarztrauber, p. n., (ncar)
c***purpose  initialize for sint.
c***description
c
c  subroutine sinti initializes the array wsave which is used in
c  subroutine sint.  the prime factorization of n together with
c  a tabulation of the trigonometric functions are computed and
c  stored in wsave.
c
c  input parameter
c
c  n       the length of the sequence to be transformed.  the method
c          is most efficient when n+1 is a product of small primes.
c
c  output parameter
c
c  wsave   a work array with at least int(3.5*n+16) locations.
c          different wsave arrays are required for different values
c          of n.  the contents of wsave must not be changed between
c          calls of sint.
c***references  (none)
c***routines called  rffti
c***end prologue  sinti
      dimension       wsave(*)
c***first executable statement  sinti
      if (n .le. 1) return
      np1 = n+1
      ns2 = n/2
      dt = pi/np1
      ks = n+2
      kf = ks+ns2-1
      fk = 0.
      do 101 k=ks,kf
         fk = fk+one
         wsave(k) = two*sin(fk*dt)
  101 continue
      call rffti (np1,wsave(kf+1))
      return
      end
      subroutine sint(n,x,wsave)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (sqrt3=1.73205080756888d0)
c***begin prologue  sint
c***date written   790601   (yymmdd)
c***revision date  830401   (yymmdd)
c***category no.  j1a3
c***keywords  fourier transform
c***author  swarztrauber, p. n., (ncar)
c***purpose  sine transform of a real, odd sequence.
c***description
c
c  subroutine sint computes the discrete fourier sine transform
c  of an odd sequence x(i).  the transform is defined below at
c  output parameter x.
c
c  sint is the unnormalized inverse of itself since a call of sint
c  followed by another call of sint will multiply the input sequence
c  x by 2*(n+1).
c
c  the array wsave which is used by subroutine sint must be
c  initialized by calling subroutine sinti(n,wsave).
c
c  input parameters
c
c  n       the length of the sequence to be transformed.  the method
c          is most efficient when n+1 is the product of small primes.
c
c  x       an array which contains the sequence to be transformed
c
c
c  wsave   a work array with dimension at least int(3.5*n+16)
c          in the program that calls sint.  the wsave array must be
c          initialized by calling subroutine sinti(n,wsave), and a
c          different wsave array must be used for each different
c          value of n.  this initialization does not have to be
c          repeated so long as n remains unchanged.  thus subsequent
c          transforms can be obtained faster than the first.
c
c  output parameters
c
c  x       for i=1,...,n
c
c               x(i)= the sum from k=1 to k=n
c
c                    2*x(k)*sin(k*i*pi/(n+1))
c
c               a call of sint followed by another call of
c               sint will multiply the sequence x by 2*(n+1).
c               hence sint is the unnormalized inverse
c               of itself.
c
c  wsave   contains initialization calculations which must not be
c          destroyed between calls of sint.
c***references  (none)
c***routines called  rfftf
c***end prologue  sint
      dimension       x(*)       ,wsave(*)
c***first executable statement  sint
      if (n-2) 101,102,103
  101 x(1) = x(1)+x(1)
      return
  102 xh = sqrt3*(x(1)+x(2))
      x(2) = sqrt3*(x(1)-x(2))
      x(1) = xh
      return
  103 np1 = n+1
      ns2 = n/2
      wsave(1) = zero
      kw = np1
      do 104 k=1,ns2
1        kw = kw+1
         kc = np1-k
         t1 = x(k)-x(kc)
         t2 = wsave(kw)*(x(k)+x(kc))
         wsave(k+1) = t1+t2
         wsave(kc+1) = t2-t1
  104 continue
      modn = mod(n,2)
      if (modn .ne. 0) wsave(ns2+2) = four*x(ns2+1)
      nf = np1+ns2+1
      call rfftf (np1,wsave,wsave(nf))
      x(1) = half*wsave(1)
      do 105 i=3,n,2
         x(i-1) = -wsave(i)
         x(i) = x(i-2)+wsave(i-1)
  105 continue
      if (modn .ne. 0) return
      x(n) = -wsave(n+1)
      return
      end
      subroutine sinqf(n,x,wsave)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c***begin prologue  sinqf
c***date written   790601   (yymmdd)
c***revision date  830401   (yymmdd)
c***category no.  j1a3
c***keywords  fourier transform
c***author  swarztrauber, p. n., (ncar)
c***purpose  forward sine transform with odd wave numbers.
c***description
c
c  subroutine sinqf computes the fast fourier transform of quarter
c  wave data.  that is, sinqf computes the coefficients in a sine
c  series representation with only odd wave numbers.  the transform
c  is defined below at output parameter x.
c
c  sinqb is the unnormalized inverse of sinqf since a call of sinqf
c  followed by a call of sinqb will multiply the input sequence x
c  by 4*n.
c
c  the array wsave which is used by subroutine sinqf must be
c  initialized by calling subroutine sinqi(n,wsave).
c
c
c  input parameters
c
c  n       the length of the array x to be transformed.  the method
c          is most efficient when n is a product of small primes.
c
c  x       an array which contains the sequence to be transformed
c
c  wsave   a work array which must be dimensioned at least 3*n+15
c          in the program that calls sinqf.  the wsave array must be
c          initialized by calling subroutine sinqi(n,wsave), and a
c          different wsave array must be used for each different
c          value of n.  this initialization does not have to be
c          repeated so long as n remains unchanged.  thus subsequent
c          transforms can be obtained faster than the first.
c
c  output parameters
c
c  x       for i=1,...,n
c
c               x(i) = (-1)**(i-1)*x(n)
c
c                  + the sum from k=1 to k=n-1 of
c
c                  2*x(k)*sin((2*i-1)*k*pi/(2*n))
c
c               a call of sinqf followed by a call of
c               sinqb will multiply the sequence x by 4*n.
c               therefore sinqb is the unnormalized inverse
c               of sinqf.
c
c  wsave   contains initialization calculations which must not
c          be destroyed between calls of sinqf or sinqb.
c***references  (none)
c***routines called  cosqf
c***end prologue  sinqf
      dimension       x(*)       ,wsave(*)
c***first executable statement  sinqf
      if (n .eq. 1) return
      ns2 = n/2
      do 101 k=1,ns2
         kc = n-k
         xhold = x(k)
         x(k) = x(kc+1)
         x(kc+1) = xhold
  101 continue
      call cosqf (n,x,wsave)
      do 102 k=2,n,2
         x(k) = -x(k)
  102 continue
      return
      end
      subroutine sinqb(n,x,wsave)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c***begin prologue  sinqb
c***date written   790601   (yymmdd)
c***revision date  830401   (yymmdd)
c***category no.  j1a3
c***keywords  fourier transform
c***author  swarztrauber, p. n., (ncar)
c***purpose  unnormalized inverse of sinqf.
c***description
c
c  subroutine sinqb computes the fast fourier transform of quarter
c  wave data.  that is, sinqb computes a sequence from its
c  representation in terms of a sine series with odd wave numbers.
c  the transform is defined below at output parameter x.
c
c  sinqf is the unnormalized inverse of sinqb since a call of sinqb
c  followed by a call of sinqf will multiply the input sequence x
c  by 4*n.
c
c  the array wsave which is used by subroutine sinqb must be
c  initialized by calling subroutine sinqi(n,wsave).
c
c
c  input parameters
c
c  n       the length of the array x to be transformed.  the method
c          is most efficient when n is a product of small primes.
c
c  x       an array which contains the sequence to be transformed
c
c  wsave   a work array which must be dimensioned at least 3*n+15
c          in the program that calls sinqb.  the wsave array must be
c          initialized by calling subroutine sinqi(n,wsave), and a
c          different wsave array must be used for each different
c          value of n.  this initialization does not have to be
c          repeated so long as n remains unchanged.  thus subsequent
c          transforms can be obtained faster than the first.
c
c  output parameters
c
c  x       for i=1,...,n
c
c               x(i)= the sum from k=1 to k=n of
c
c                 4*x(k)*sin((2k-1)*i*pi/(2*n))
c
c               a call of sinqb followed by a call of
c               sinqf will multiply the sequence x by 4*n.
c               therefore sinqf is the unnormalized inverse
c               of sinqb.
c
c  wsave   contains initialization calculations which must not
c          be destroyed between calls of sinqb or sinqf.
c***references  (none)
c***routines called  cosqb
c***end prologue  sinqb
      dimension       x(*)       ,wsave(*)
c***first executable statement  sinqb
      if (n .gt. 1) go to 101
      x(1) = four*x(1)
      return
  101 ns2 = n/2
      do 102 k=2,n,2
         x(k) = -x(k)
  102 continue
      call cosqb (n,x,wsave)
      do 103 k=1,ns2
         kc = n-k
         xhold = x(k)
         x(k) = x(kc+1)
         x(kc+1) = xhold
  103 continue
      return
      end
      subroutine sinqi(n,wsave)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c***begin prologue  sinqi
c***date written   790601   (yymmdd)
c***revision date  830401   (yymmdd)
c***category no.  j1a3
c***keywords  fourier transform
c***author  swarztrauber, p. n., (ncar)
c***purpose  initialize for sinqf and sinqb.
c***description
c
c  subroutine sinqi initializes the array wsave which is used in
c  both sinqf and sinqb.  the prime factorization of n together with
c  a tabulation of the trigonometric functions are computed and
c  stored in wsave.
c
c  input parameter
c
c  n       the length of the sequence to be transformed.  the method
c          is most efficient when n is a product of small primes.
c
c  output parameter
c
c  wsave   a work array which must be dimensioned at least 3*n+15.
c          the same work array can be used for both sinqf and sinqb
c          as long as n remains unchanged.  different wsave arrays
c          are required for different values of n.  the contents of
c          wsave must not be changed between calls of sinqf or sinqb.
c***references  (none)
c***routines called  cosqi
c***end prologue  sinqi
      dimension       wsave(*)
c***first executable statement  sinqi
      call cosqi (n,wsave)
      return
      end
      subroutine cosqf(n,x,wsave)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (sqrt2=1.4142135623731d0)
c***begin prologue  cosqf
c***date written   790601   (yymmdd)
c***revision date  830401   (yymmdd)
c***category no.  j1a3
c***keywords  fourier transform
c***author  swarztrauber, p. n., (ncar)
c***purpose  forward cosine transform with odd wave numbers.
c***description
c
c  subroutine cosqf computes the fast fourier transform of quarter
c  wave data. that is, cosqf computes the coefficients in a cosine
c  series representation with only odd wave numbers.  the transform
c  is defined below at output parameter x
c
c  cosqf is the unnormalized inverse of cosqb since a call of cosqf
c  followed by a call of cosqb will multiply the input sequence x
c  by 4*n.
c
c  the array wsave which is used by subroutine cosqf must be
c  initialized by calling subroutine cosqi(n,wsave).
c
c
c  input parameters
c
c  n       the length of the array x to be transformed.  the method
c          is most efficient when n is a product of small primes.
c
c  x       an array which contains the sequence to be transformed
c
c  wsave   a work array which must be dimensioned at least 3*n+15
c          in the program that calls cosqf.  the wsave array must be
c          initialized by calling subroutine cosqi(n,wsave), and a
c          different wsave array must be used for each different
c          value of n.  this initialization does not have to be
c          repeated so long as n remains unchanged.  thus subsequent
c          transforms can be obtained faster than the first.
c
c  output parameters
c
c  x       for i=1,...,n
c
c               x(i) = x(1) plus the sum from k=2 to k=n of
c
c                  2*x(k)*cos((2*i-1)*(k-1)*pi/(2*n))
c
c               a call of cosqf followed by a call of
c               cosqb will multiply the sequence x by 4*n.
c               therefore cosqb is the unnormalized inverse
c               of cosqf.
c
c  wsave   contains initialization calculations which must not
c          be destroyed between calls of cosqf or cosqb.
c***references  (none)
c***routines called  cosqf1
c***end prologue  cosqf
      dimension       x(*)       ,wsave(*)
c***first executable statement  cosqf
      if (n-2) 102,101,103
  101 tsqx = sqrt2*x(2)
      x(2) = x(1)-tsqx
      x(1) = x(1)+tsqx
  102 return
  103 call cosqf1 (n,x,wsave,wsave(n+1))
      return
      end
      subroutine cosqb(n,x,wsave)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (tsqrt2=2.82842712474619d0)
c***begin prologue  cosqb
c***date written   790601   (yymmdd)
c***revision date  830401   (yymmdd)
c***category no.  j1a3
c***keywords  fourier transform
c***author  swarztrauber, p. n., (ncar)
c***purpose  unnormalized inverse of cosqf.
c***description
c
c  subroutine cosqb computes the fast fourier transform of quarter
c  wave data. that is, cosqb computes a sequence from its
c  representation in terms of a cosine series with odd wave numbers.
c  the transform is defined below at output parameter x.
c
c  cosqb is the unnormalized inverse of cosqf since a call of cosqb
c  followed by a call of cosqf will multiply the input sequence x
c  by 4*n.
c
c  the array wsave which is used by subroutine cosqb must be
c  initialized by calling subroutine cosqi(n,wsave).
c
c
c  input parameters
c
c  n       the length of the array x to be transformed.  the method
c          is most efficient when n is a product of small primes.
c
c  x       an array which contains the sequence to be transformed
c
c  wsave   a work array that must be dimensioned at least 3*n+15
c          in the program that calls cosqb.  the wsave array must be
c          initialized by calling subroutine cosqi(n,wsave), and a
c          different wsave array must be used for each different
c          value of n.  this initialization does not have to be
c          repeated so long as n remains unchanged.  thus subsequent
c          transforms can be obtained faster than the first.
c
c  output parameters
c
c  x       for i=1,...,n
c
c               x(i)= the sum from k=1 to k=n of
c
c                 4*x(k)*cos((2*k-1)*(i-1)*pi/(2*n))
c
c               a call of cosqb followed by a call of
c               cosqf will multiply the sequence x by 4*n.
c               therefore cosqf is the unnormalized inverse
c               of cosqb.
c
c  wsave   contains initialization calculations which must not
c          be destroyed between calls of cosqb or cosqf.
c***references  (none)
c***routines called  cosqb1
c***end prologue  cosqb
      dimension       x(*)       ,wsave(*)
c***first executable statement  cosqb
      if (n-2) 101,102,103
  101 x(1) = four*x(1)
      return
  102 x1 = four*(x(1)+x(2))
      x(2) = tsqrt2*(x(1)-x(2))
      x(1) = x1
      return
  103 call cosqb1 (n,x,wsave,wsave(n+1))
      return
      end
      subroutine cosqi(n,wsave)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (pih=1.57079632679491d0)
c***begin prologue  cosqi
c***date written   790601   (yymmdd)
c***revision date  830401   (yymmdd)
c***category no.  j1a3
c***keywords  fourier transform
c***author  swarztrauber, p. n., (ncar)
c***purpose  initialize for cosqf and cosqb.
c***description
c
c  subroutine cosqi initializes the array wsave which is used in
c  both cosqf and cosqb.  the prime factorization of n together with
c  a tabulation of the trigonometric functions are computed and
c  stored in wsave.
c
c  input parameter
c
c  n       the length of the array to be transformed.  the method
c          is most efficient when n is a product of small primes.
c
c  output parameter
c
c  wsave   a work array which must be dimensioned at least 3*n+15.
c          the same work array can be used for both cosqf and cosqb
c          as long as n remains unchanged.  different wsave arrays
c          are required for different values of n.  the contents of
c          wsave must not be changed between calls of cosqf or cosqb.
c***references  (none)
c***routines called  rffti
c***end prologue  cosqi
      dimension       wsave(*)
c***first executable statement  cosqi
      dt = pih/dfloat(n)
      fk = zero
      do 101 k=1,n
         fk = fk+one
         wsave(k) = cos(fk*dt)
  101 continue
      call rffti (n,wsave(n+1))
      return
      end
      subroutine cosqb1(n,x,w,xh)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c***begin prologue  cosqb1
c***refer to  cosqb
c***routines called  rfftb
c***end prologue  cosqb1
      dimension       x(*)       ,w(*)       ,xh(*)
c***first executable statement  cosqb1
      ns2 = (n+1)/2
      np2 = n+2
      do 101 i=3,n,2
         xim1 = x(i-1)+x(i)
         x(i) = x(i)-x(i-1)
         x(i-1) = xim1
  101 continue
      x(1) = x(1)+x(1)
      modn = mod(n,2)
      if (modn .eq. 0) x(n) = x(n)+x(n)
      call rfftb (n,x,xh)
      do 102 k=2,ns2
         kc = np2-k
         xh(k) = w(k-1)*x(kc)+w(kc-1)*x(k)
         xh(kc) = w(k-1)*x(k)-w(kc-1)*x(kc)
  102 continue
      if (modn .eq. 0) x(ns2+1) = w(ns2)*(x(ns2+1)+x(ns2+1))
      do 103 k=2,ns2
         kc = np2-k
         x(k) = xh(k)+xh(kc)
         x(kc) = xh(k)-xh(kc)
  103 continue
      x(1) = x(1)+x(1)
      return
      end
      subroutine rfftb(n,r,wsave)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c***begin prologue  rfftb
c***date written   790601   (yymmdd)
c***revision date  830401   (yymmdd)
c***category no.  j1a1
c***keywords  fourier transform
c***author  swarztrauber, p. n., (ncar)
c***purpose  backward transform of a real coefficient array.
c***description
c
c  subroutine rfftb computes the real perodic sequence from its
c  fourier coefficients (fourier synthesis).  the transform is defined
c  below at output parameter r.
c
c  input parameters
c
c  n       the length of the array r to be transformed.  the method
c          is most efficient when n is a product of small primes.
c          n may change so long as different work arrays are provided.
c
c  r       a real array of length n which contains the sequence
c          to be transformed
c
c  wsave   a work array which must be dimensioned at least 2*n+15
c          in the program that calls rfftb.  the wsave array must be
c          initialized by calling subroutine rffti(n,wsave), and a
c          different wsave array must be used for each different
c          value of n.  this initialization does not have to be
c          repeated so long as n remains unchanged.  thus subsequent
c          transforms can be obtained faster than the first.
c          the same wsave array can be used by rfftf and rfftb.
c
c
c  output parameters
c
c  r       for n even and for i = 1,...,n
c
c               r(i) = r(1)+(-1)**(i-1)*r(n)
c
c                    plus the sum from k=2 to k=n/2 of
c
c                     2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)
c
c                    -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
c
c          for n odd and for i = 1,...,n
c
c               r(i) = r(1) plus the sum from k=2 to k=(n+1)/2 of
c
c                    2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)
c
c                   -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
c
c   *****  note:
c               this transform is unnormalized since a call of rfftf
c               followed by a call of rfftb will multiply the input
c               sequence by n.
c
c  wsave   contains results which must not be destroyed between
c          calls of rfftb or rfftf.
c***references  (none)
c***routines called  rfftb1
c***end prologue  rfftb
      dimension       r(*)       ,wsave(*)
c***first executable statement  rfftb
      if (n .eq. 1) return
      call rfftb1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
      return
      end
      subroutine rfftb1(n,c,ch,wa,ifac)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c***begin prologue  rfftb1
c***refer to  rfftb
c***routines called  radb2,radb3,radb4,radb5,radbg
c***end prologue  rfftb1
      dimension       ch(*)      ,c(*)       ,wa(*)      ,ifac(*)
c***first executable statement  rfftb1
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call radb4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call radb4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call radb2 (ido,l1,c,ch,wa(iw))
         go to 105
  104    call radb2 (ido,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip .ne. 3) go to 109
         ix2 = iw+ido
         if (na .ne. 0) go to 107
         call radb3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call radb3 (ido,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip .ne. 5) go to 112
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 110
         call radb5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call radb5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na .ne. 0) go to 113
         call radbg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call radbg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (ido .eq. 1) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*ido
  116 continue
      if (na .eq. 0) return
      do 117 i=1,n
         c(i) = ch(i)
  117 continue
      return
      end
      subroutine cosqf1(n,x,w,xh)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c***begin prologue  cosqf1
c***refer to  cosqf
c***routines called  rfftf
c***end prologue  cosqf1
      dimension       x(*)       ,w(*)       ,xh(*)
c***first executable statement  cosqf1
      ns2 = (n+1)/2
      np2 = n+2
      do 101 k=2,ns2
         kc = np2-k
         xh(k) = x(k)+x(kc)
         xh(kc) = x(k)-x(kc)
  101 continue
      modn = mod(n,2)
      if (modn .eq. 0) xh(ns2+1) = x(ns2+1)+x(ns2+1)
      do 102 k=2,ns2
         kc = np2-k
         x(k) = w(k-1)*xh(kc)+w(kc-1)*xh(k)
         x(kc) = w(k-1)*xh(k)-w(kc-1)*xh(kc)
  102 continue
      if (modn .eq. 0) x(ns2+1) = w(ns2)*xh(ns2+1)
      call rfftf (n,x,xh)
      do 103 i=3,n,2
         xim1 = x(i-1)-x(i)
         x(i) = x(i-1)+x(i)
         x(i-1) = xim1
  103 continue
      return
      end
      subroutine rfftf(n,r,wsave)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c***begin prologue  rfftf
c***date written   790601   (yymmdd)
c***revision date  830401   (yymmdd)
c***category no.  j1a1
c***keywords  fourier transform
c***author  swarztrauber, p. n., (ncar)
c***purpose  forward transform of a real, periodic sequence.
c***description
c
c  subroutine rfftf computes the fourier coefficients of a real
c  perodic sequence (fourier analysis).  the transform is defined
c  below at output parameter r.
c
c  input parameters
c
c  n       the length of the array r to be transformed.  the method
c          is most efficient when n is a product of small primes.
c          n may change so long as different work arrays are provided
c
c  r       a real array of length n which contains the sequence
c          to be transformed
c
c  wsave   a work array which must be dimensioned at least 2*n+15
c          in the program that calls rfftf.  the wsave array must be
c          initialized by calling subroutine rffti(n,wsave), and a
c          different wsave array must be used for each different
c          value of n.  this initialization does not have to be
c          repeated so long as n remains unchanged.  thus subsequent
c          transforms can be obtained faster than the first.
c          the same wsave array can be used by rfftf and rfftb.
c
c
c  output parameters
c
c  r       r(1) = the sum from i=1 to i=n of r(i)
c
c          if n is even set l = n/2; if n is odd set l = (n+1)/2
c
c            then for k = 2,...,l
c
c               r(2*k-2) = the sum from i = 1 to i = n of
c
c                    r(i)*cos((k-1)*(i-1)*2*pi/n)
c
c               r(2*k-1) = the sum from i = 1 to i = n of
c
c                   -r(i)*sin((k-1)*(i-1)*2*pi/n)
c
c          if n is even
c
c               r(n) = the sum from i = 1 to i = n of
c
c                    (-1)**(i-1)*r(i)
c
c   *****  note:
c               this transform is unnormalized since a call of rfftf
c               followed by a call of rfftb will multiply the input
c               sequence by n.
c
c  wsave   contains results which must not be destroyed between
c          calls of rfftf or rfftb.
c***references  (none)
c***routines called  rfftf1
c***end prologue  rfftf
      dimension       r(*)       ,wsave(*)
c***first executable statement  rfftf
      if (n .eq. 1) return
      call rfftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
      return
      end
      subroutine rfftf1(n,c,ch,wa,ifac)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c***begin prologue  rfftf1
c***refer to  rfftf
c***routines called  radf2,radf3,radf4,radf5,radfg
c***end prologue  rfftf1
      dimension       ch(*)      ,c(*)       ,wa(*)      ,ifac(*)
c***first executable statement  rfftf1
      nf = ifac(2)
      na = 1
      l2 = n
      iw = n
      do 111 k1=1,nf
         kh = nf-k1
         ip = ifac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
         if (ip .ne. 4) go to 102
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call radf4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 110
  101    call radf4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
         go to 110
  102    if (ip .ne. 2) go to 104
         if (na .ne. 0) go to 103
         call radf2 (ido,l1,c,ch,wa(iw))
         go to 110
  103    call radf2 (ido,l1,ch,c,wa(iw))
         go to 110
  104    if (ip .ne. 3) go to 106
         ix2 = iw+ido
         if (na .ne. 0) go to 105
         call radf3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 110
  105    call radf3 (ido,l1,ch,c,wa(iw),wa(ix2))
         go to 110
  106    if (ip .ne. 5) go to 108
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 107
         call radf5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  107    call radf5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  108    if (ido .eq. 1) na = 1-na
         if (na .ne. 0) go to 109
         call radfg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         na = 1
         go to 110
  109    call radfg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
         na = 0
  110    l2 = l1
  111 continue
      if (na .eq. 1) return
      do 112 i=1,n
         c(i) = ch(i)
  112 continue
      return
      end
      subroutine rffti(n,wsave)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c***begin prologue  rffti
c***date written   790601   (yymmdd)
c***revision date  830401   (yymmdd)
c***category no.  j1a1
c***keywords  fourier transform
c***author  swarztrauber, p. n., (ncar)
c***purpose  initialize for rfftf and rfftb.
c***description
c
c  subroutine rffti initializes the array wsave which is used in
c  both rfftf and rfftb.  the prime factorization of n together with
c  a tabulation of the trigonometric functions are computed and
c  stored in wsave.
c
c  input parameter
c
c  n       the length of the sequence to be transformed.
c
c  output parameter
c
c  wsave   a work array which must be dimensioned at least 2*n+15.
c          the same work array can be used for both rfftf and rfftb
c          as long as n remains unchanged.  different wsave arrays
c          are required for different values of n.  the contents of
c          wsave must not be changed between calls of rfftf or rfftb.
c***references  (none)
c***routines called  rffti1
c***end prologue  rffti
      dimension       wsave(*)
c***first executable statement  rffti
      if (n .eq. 1) return
      call rffti1 (n,wsave(n+1),wsave(2*n+1))
      return
      end
      subroutine rffti1(n,wa,ifac)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (tpi = 6.28318530717959d0)
c***begin prologue  rffti1
c***refer to  rffti
c***routines called  (none)
c***end prologue  rffti1
      dimension       wa(*)      ,ifac(*)    ,ntryh(4)
      data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/4,2,3,5/
c***first executable statement  rffti1
      nl = n
      nf = 0
      j = 0
  101 j = j+1
      if (j-4) 102,102,103
  102 ntry = ntryh(j)
      go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr) 101,105,101
  105 nf = nf+1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         ifac(ib+2) = ifac(ib+1)
  106 continue
      ifac(3) = 2
  107 if (nl .ne. 1) go to 104
      ifac(1) = n
      ifac(2) = nf
      argh = tpi/dfloat(n)
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1 .eq. 0) return
      do 110 k1=1,nfm1
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do 109 j=1,ipm
            ld = ld+l1
            i = is
            argld = dfloat(ld)*argh
            fi = zero
            do 108 ii=3,ido,2
               i = i+2
               fi = fi+one
               arg = fi*argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
  108       continue
            is = is+ido
  109    continue
         l1 = l2
  110 continue
      return
      end
      subroutine radb2(ido,l1,cc,ch,wa1)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c***begin prologue  radb2
c***refer to  rfftb
c***routines called  (none)
c***end prologue  radb2
      dimension       cc(ido,2,l1)           ,ch(ido,l1,2)           ,
     1                wa1(*)
c***first executable statement  radb2
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(ido,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(ido,2,k)
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      if((ido-1)/2.lt.l1) go to 108
      do 104 k=1,l1
cdir$ ivdep
         do 103 i=3,ido,2
            ic = idp2-i
            ch(i-1,k,1) = cc(i-1,1,k)+cc(ic-1,2,k)
            tr2 = cc(i-1,1,k)-cc(ic-1,2,k)
            ch(i,k,1) = cc(i,1,k)-cc(ic,2,k)
            ti2 = cc(i,1,k)+cc(ic,2,k)
            ch(i-1,k,2) = wa1(i-2)*tr2-wa1(i-1)*ti2
            ch(i,k,2) = wa1(i-2)*ti2+wa1(i-1)*tr2
  103    continue
  104 continue
      go to 111
  108 do 110 i=3,ido,2
         ic = idp2-i
cdir$ ivdep
         do 109 k=1,l1
            ch(i-1,k,1) = cc(i-1,1,k)+cc(ic-1,2,k)
            tr2 = cc(i-1,1,k)-cc(ic-1,2,k)
            ch(i,k,1) = cc(i,1,k)-cc(ic,2,k)
            ti2 = cc(i,1,k)+cc(ic,2,k)
            ch(i-1,k,2) = wa1(i-2)*tr2-wa1(i-1)*ti2
            ch(i,k,2) = wa1(i-2)*ti2+wa1(i-1)*tr2
  109    continue
  110 continue
  111 if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(ido,k,1) = cc(ido,1,k)+cc(ido,1,k)
         ch(ido,k,2) = -(cc(1,2,k)+cc(1,2,k))
  106 continue
  107 return
      end
      subroutine radb3(ido,l1,cc,ch,wa1,wa2)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (taur=-.5d0,taui=.866025403784439d0)
c***begin prologue  radb3
c***refer to  rfftb
c***routines called  (none)
c***end prologue  radb3
      dimension       cc(ido,3,l1)           ,ch(ido,l1,3)           ,
     1                wa1(*)     ,wa2(*)
c***first executable statement  radb3
      do 101 k=1,l1
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         cr2 = cc(1,1,k)+taur*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ci3 = taui*(cc(1,3,k)+cc(1,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      if((ido-1)/2.lt.l1) go to 104
      do 103 k=1,l1
cdir$ ivdep
         do 102 i=3,ido,2
            ic = idp2-i
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,3,k)-cc(ic-1,2,k))
            ci3 = taui*(cc(i,3,k)+cc(ic,2,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
  102    continue
  103 continue
      return
  104 do 106 i=3,ido,2
         ic = idp2-i
cdir$ ivdep
         do 105 k=1,l1
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,3,k)-cc(ic-1,2,k))
            ci3 = taui*(cc(i,3,k)+cc(ic,2,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
  105    continue
  106 continue
      return
      end
      subroutine radb4(ido,l1,cc,ch,wa1,wa2,wa3)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (sqrt2=1.414213562373095d0)
c***begin prologue  radb4
c***refer to  rfftb
c***routines called  (none)
c***end prologue  radb4
      dimension       cc(ido,4,l1)           ,ch(ido,l1,4)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)
c***first executable statement  radb4
      do 101 k=1,l1
         tr1 = cc(1,1,k)-cc(ido,4,k)
         tr2 = cc(1,1,k)+cc(ido,4,k)
         tr3 = cc(ido,2,k)+cc(ido,2,k)
         tr4 = cc(1,3,k)+cc(1,3,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,2) = tr1-tr4
         ch(1,k,3) = tr2-tr3
         ch(1,k,4) = tr1+tr4
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      if((ido-1)/2.lt.l1) go to 108
      do 104 k=1,l1
cdir$ ivdep
         do 103 i=3,ido,2
            ic = idp2-i
            ti1 = cc(i,1,k)+cc(ic,4,k)
            ti2 = cc(i,1,k)-cc(ic,4,k)
            ti3 = cc(i,3,k)-cc(ic,2,k)
            tr4 = cc(i,3,k)+cc(ic,2,k)
            tr1 = cc(i-1,1,k)-cc(ic-1,4,k)
            tr2 = cc(i-1,1,k)+cc(ic-1,4,k)
            ti4 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr3 = cc(i-1,3,k)+cc(ic-1,2,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1-tr4
            cr4 = tr1+tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-2)*cr2-wa1(i-1)*ci2
            ch(i,k,2) = wa1(i-2)*ci2+wa1(i-1)*cr2
            ch(i-1,k,3) = wa2(i-2)*cr3-wa2(i-1)*ci3
            ch(i,k,3) = wa2(i-2)*ci3+wa2(i-1)*cr3
            ch(i-1,k,4) = wa3(i-2)*cr4-wa3(i-1)*ci4
            ch(i,k,4) = wa3(i-2)*ci4+wa3(i-1)*cr4
  103    continue
  104 continue
      go to 111
  108 do 110 i=3,ido,2
         ic = idp2-i
cdir$ ivdep
         do 109 k=1,l1
            ti1 = cc(i,1,k)+cc(ic,4,k)
            ti2 = cc(i,1,k)-cc(ic,4,k)
            ti3 = cc(i,3,k)-cc(ic,2,k)
            tr4 = cc(i,3,k)+cc(ic,2,k)
            tr1 = cc(i-1,1,k)-cc(ic-1,4,k)
            tr2 = cc(i-1,1,k)+cc(ic-1,4,k)
            ti4 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr3 = cc(i-1,3,k)+cc(ic-1,2,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1-tr4
            cr4 = tr1+tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-2)*cr2-wa1(i-1)*ci2
            ch(i,k,2) = wa1(i-2)*ci2+wa1(i-1)*cr2
            ch(i-1,k,3) = wa2(i-2)*cr3-wa2(i-1)*ci3
            ch(i,k,3) = wa2(i-2)*ci3+wa2(i-1)*cr3
            ch(i-1,k,4) = wa3(i-2)*cr4-wa3(i-1)*ci4
            ch(i,k,4) = wa3(i-2)*ci4+wa3(i-1)*cr4
  109    continue
  110 continue
  111 if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ti1 = cc(1,2,k)+cc(1,4,k)
         ti2 = cc(1,4,k)-cc(1,2,k)
         tr1 = cc(ido,1,k)-cc(ido,3,k)
         tr2 = cc(ido,1,k)+cc(ido,3,k)
         ch(ido,k,1) = tr2+tr2
         ch(ido,k,2) = sqrt2*(tr1-ti1)
         ch(ido,k,3) = ti2+ti2
         ch(ido,k,4) = -sqrt2*(tr1+ti1)
  106 continue
  107 return
      end
      subroutine radb5(ido,l1,cc,ch,wa1,wa2,wa3,wa4)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (tr11=.309016994374947d0,ti11=.951056516295154d0)
      parameter (tr12=-.809016994374947d0,ti12=.587785252292473d0)
c***begin prologue  radb5
c***refer to  rfftb
c***routines called  (none)
c***end prologue  radb5
      dimension       cc(ido,5,l1)           ,ch(ido,l1,5)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
c***first executable statement  radb5
      do 101 k=1,l1
         ti5 = cc(1,3,k)+cc(1,3,k)
         ti4 = cc(1,5,k)+cc(1,5,k)
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         tr3 = cc(ido,4,k)+cc(ido,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
         cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
         ci5 = ti11*ti5+ti12*ti4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(1,k,5) = cr2+ci5
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      if((ido-1)/2.lt.l1) go to 104
      do 103 k=1,l1
cdir$ ivdep
         do 102 i=3,ido,2
            ic = idp2-i
            ti5 = cc(i,3,k)+cc(ic,2,k)
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ti4 = cc(i,5,k)+cc(ic,4,k)
            ti3 = cc(i,5,k)-cc(ic,4,k)
            tr5 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            tr4 = cc(i-1,5,k)-cc(ic-1,4,k)
            tr3 = cc(i-1,5,k)+cc(ic-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
            ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
            ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
            ch(i-1,k,4) = wa3(i-2)*dr4-wa3(i-1)*di4
            ch(i,k,4) = wa3(i-2)*di4+wa3(i-1)*dr4
            ch(i-1,k,5) = wa4(i-2)*dr5-wa4(i-1)*di5
            ch(i,k,5) = wa4(i-2)*di5+wa4(i-1)*dr5
  102    continue
  103 continue
      return
  104 do 106 i=3,ido,2
         ic = idp2-i
cdir$ ivdep
         do 105 k=1,l1
            ti5 = cc(i,3,k)+cc(ic,2,k)
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ti4 = cc(i,5,k)+cc(ic,4,k)
            ti3 = cc(i,5,k)-cc(ic,4,k)
            tr5 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            tr4 = cc(i-1,5,k)-cc(ic-1,4,k)
            tr3 = cc(i-1,5,k)+cc(ic-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
            ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
            ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
            ch(i-1,k,4) = wa3(i-2)*dr4-wa3(i-1)*di4
            ch(i,k,4) = wa3(i-2)*di4+wa3(i-1)*dr4
            ch(i-1,k,5) = wa4(i-2)*dr5-wa4(i-1)*di5
            ch(i,k,5) = wa4(i-2)*di5+wa4(i-1)*dr5
  105    continue
  106 continue
      return
      end
      subroutine radbg(ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (tpi=6.28318530717959d0)
c***begin prologue  radbg
c***refer to  rfftb
c***routines called  (none)
c***end prologue  radbg
      dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          ,
     1                c1(ido,l1,ip)          ,c2(idl1,ip),
     2                ch2(idl1,ip)           ,wa(*)
c***first executable statement  radbg
      arg = tpi/dfloat(ip)
      dcp = cos(arg)
      dsp = sin(arg)
      idp2 = ido+2
      nbd = (ido-1)/2
      ipp2 = ip+2
      ipph = (ip+1)/2
      if (ido .lt. l1) go to 103
      do 102 k=1,l1
         do 101 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  101    continue
  102 continue
      go to 106
  103 do 105 i=1,ido
         do 104 k=1,l1
            ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
  106 do 108 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 107 k=1,l1
            ch(1,k,j) = cc(ido,j2-2,k)+cc(ido,j2-2,k)
            ch(1,k,jc) = cc(1,j2-1,k)+cc(1,j2-1,k)
  107    continue
  108 continue
      if (ido .eq. 1) go to 116
      if (nbd .lt. l1) go to 112
      do 111 j=2,ipph
         jc = ipp2-j
         do 110 k=1,l1
cdir$ ivdep
            do 109 i=3,ido,2
               ic = idp2-i
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  109       continue
  110    continue
  111 continue
      go to 116
  112 do 115 j=2,ipph
         jc = ipp2-j
cdir$ ivdep
         do 114 i=3,ido,2
            ic = idp2-i
            do 113 k=1,l1
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  113       continue
  114    continue
  115 continue
  116 ar1 = one
      ai1 = zero
      do 120 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 117 ik=1,idl1
            c2(ik,l) = ch2(ik,1)+ar1*ch2(ik,2)
            c2(ik,lc) = ai1*ch2(ik,ip)
  117    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 119 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 118 ik=1,idl1
               c2(ik,l) = c2(ik,l)+ar2*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)+ai2*ch2(ik,jc)
  118       continue
  119    continue
  120 continue
      do 122 j=2,ipph
         do 121 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  121    continue
  122 continue
      do 124 j=2,ipph
         jc = ipp2-j
         do 123 k=1,l1
            ch(1,k,j) = c1(1,k,j)-c1(1,k,jc)
            ch(1,k,jc) = c1(1,k,j)+c1(1,k,jc)
  123    continue
  124 continue
      if (ido .eq. 1) go to 132
      if (nbd .lt. l1) go to 128
      do 127 j=2,ipph
         jc = ipp2-j
         do 126 k=1,l1
cdir$ ivdep
            do 125 i=3,ido,2
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  125       continue
  126    continue
  127 continue
      go to 132
  128 do 131 j=2,ipph
         jc = ipp2-j
         do 130 i=3,ido,2
            do 129 k=1,l1
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  129       continue
  130    continue
  131 continue
  132 continue
      if (ido .eq. 1) return
      do 133 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  133 continue
      do 135 j=2,ip
         do 134 k=1,l1
            c1(1,k,j) = ch(1,k,j)
  134    continue
  135 continue
      if (nbd .gt. l1) go to 139
      is = -ido
      do 138 j=2,ip
         is = is+ido
         idij = is
         do 137 i=3,ido,2
            idij = idij+2
            do 136 k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  136       continue
  137    continue
  138 continue
      go to 143
  139 is = -ido
      do 142 j=2,ip
         is = is+ido
         do 141 k=1,l1
            idij = is
cdir$ ivdep
            do 140 i=3,ido,2
               idij = idij+2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  140       continue
  141    continue
  142 continue
  143 return
      end
      subroutine radf2(ido,l1,cc,ch,wa1)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c***begin prologue  radf2
c***refer to  rfftf
c***routines called  (none)
c***end prologue  radf2
      dimension       ch(ido,2,l1)           ,cc(ido,l1,2)           ,
     1                wa1(*)
c***first executable statement  radf2
      do 101 k=1,l1
         ch(1,1,k) = cc(1,k,1)+cc(1,k,2)
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,2)
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      if((ido-1)/2.lt.l1) go to 108
      do 104 k=1,l1
cdir$ ivdep
         do 103 i=3,ido,2
            ic = idp2-i
            tr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ti2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            ch(i,1,k) = cc(i,k,1)+ti2
            ch(ic,2,k) = ti2-cc(i,k,1)
            ch(i-1,1,k) = cc(i-1,k,1)+tr2
            ch(ic-1,2,k) = cc(i-1,k,1)-tr2
  103    continue
  104 continue
      go to 111
  108 do 110 i=3,ido,2
         ic = idp2-i
cdir$ ivdep
         do 109 k=1,l1
            tr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ti2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            ch(i,1,k) = cc(i,k,1)+ti2
            ch(ic,2,k) = ti2-cc(i,k,1)
            ch(i-1,1,k) = cc(i-1,k,1)+tr2
            ch(ic-1,2,k) = cc(i-1,k,1)-tr2
  109    continue
  110 continue
  111 if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(1,2,k) = -cc(ido,k,2)
         ch(ido,1,k) = cc(ido,k,1)
  106 continue
  107 return
      end
      subroutine radf3(ido,l1,cc,ch,wa1,wa2)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (taur=-.5d0,taui=.866025403784439d0)
c***begin prologue  radf3
c***refer to  rfftf
c***routines called  (none)
c***end prologue  radf3
      dimension       ch(ido,3,l1)           ,cc(ido,l1,3)           ,
     1                wa1(*)     ,wa2(*)
c***first executable statement  radf3
      do 101 k=1,l1
         cr2 = cc(1,k,2)+cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2
         ch(1,3,k) = taui*(cc(1,k,3)-cc(1,k,2))
         ch(ido,2,k) = cc(1,k,1)+taur*cr2
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      if((ido-1)/2.lt.l1) go to 104
      do 103 k=1,l1
cdir$ ivdep
         do 102 i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr2 = dr2+dr3
            ci2 = di2+di3
            ch(i-1,1,k) = cc(i-1,k,1)+cr2
            ch(i,1,k) = cc(i,k,1)+ci2
            tr2 = cc(i-1,k,1)+taur*cr2
            ti2 = cc(i,k,1)+taur*ci2
            tr3 = taui*(di2-di3)
            ti3 = taui*(dr3-dr2)
            ch(i-1,3,k) = tr2+tr3
            ch(ic-1,2,k) = tr2-tr3
            ch(i,3,k) = ti2+ti3
            ch(ic,2,k) = ti3-ti2
  102    continue
  103 continue
      return
  104 do 106 i=3,ido,2
         ic = idp2-i
cdir$ ivdep
         do 105 k=1,l1
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr2 = dr2+dr3
            ci2 = di2+di3
            ch(i-1,1,k) = cc(i-1,k,1)+cr2
            ch(i,1,k) = cc(i,k,1)+ci2
            tr2 = cc(i-1,k,1)+taur*cr2
            ti2 = cc(i,k,1)+taur*ci2
            tr3 = taui*(di2-di3)
            ti3 = taui*(dr3-dr2)
            ch(i-1,3,k) = tr2+tr3
            ch(ic-1,2,k) = tr2-tr3
            ch(i,3,k) = ti2+ti3
            ch(ic,2,k) = ti3-ti2
  105    continue
  106 continue
      return
      end
      subroutine radf4(ido,l1,cc,ch,wa1,wa2,wa3)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (hsqt2=.7071067811865475d0)
c***begin prologue  radf4
c***refer to  rfftf
c***routines called  (none)
c***end prologue  radf4
      dimension       cc(ido,l1,4)           ,ch(ido,4,l1)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)
c***first executable statement  radf4
      do 101 k=1,l1
         tr1 = cc(1,k,2)+cc(1,k,4)
         tr2 = cc(1,k,1)+cc(1,k,3)
         ch(1,1,k) = tr1+tr2
         ch(ido,4,k) = tr2-tr1
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,3)
         ch(1,3,k) = cc(1,k,4)-cc(1,k,2)
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      if((ido-1)/2.lt.l1) go to 111
      do 104 k=1,l1
cdir$ ivdep
         do 103 i=3,ido,2
            ic = idp2-i
            cr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ci2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            cr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            ci3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            ci4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            tr1 = cr2+cr4
            tr4 = cr4-cr2
            ti1 = ci2+ci4
            ti4 = ci2-ci4
            ti2 = cc(i,k,1)+ci3
            ti3 = cc(i,k,1)-ci3
            tr2 = cc(i-1,k,1)+cr3
            tr3 = cc(i-1,k,1)-cr3
            ch(i-1,1,k) = tr1+tr2
            ch(ic-1,4,k) = tr2-tr1
            ch(i,1,k) = ti1+ti2
            ch(ic,4,k) = ti1-ti2
            ch(i-1,3,k) = ti4+tr3
            ch(ic-1,2,k) = tr3-ti4
            ch(i,3,k) = tr4+ti3
            ch(ic,2,k) = tr4-ti3
  103    continue
  104 continue
      go to 110
  111 do 109 i=3,ido,2
         ic = idp2-i
cdir$ ivdep
         do 108 k=1,l1
            cr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ci2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            cr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            ci3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            ci4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            tr1 = cr2+cr4
            tr4 = cr4-cr2
            ti1 = ci2+ci4
            ti4 = ci2-ci4
            ti2 = cc(i,k,1)+ci3
            ti3 = cc(i,k,1)-ci3
            tr2 = cc(i-1,k,1)+cr3
            tr3 = cc(i-1,k,1)-cr3
            ch(i-1,1,k) = tr1+tr2
            ch(ic-1,4,k) = tr2-tr1
            ch(i,1,k) = ti1+ti2
            ch(ic,4,k) = ti1-ti2
            ch(i-1,3,k) = ti4+tr3
            ch(ic-1,2,k) = tr3-ti4
            ch(i,3,k) = tr4+ti3
            ch(ic,2,k) = tr4-ti3
  108    continue
  109 continue
  110 if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ti1 = -hsqt2*(cc(ido,k,2)+cc(ido,k,4))
         tr1 = hsqt2*(cc(ido,k,2)-cc(ido,k,4))
         ch(ido,1,k) = tr1+cc(ido,k,1)
         ch(ido,3,k) = cc(ido,k,1)-tr1
         ch(1,2,k) = ti1-cc(ido,k,3)
         ch(1,4,k) = ti1+cc(ido,k,3)
  106 continue
  107 return
      end
      subroutine radf5(ido,l1,cc,ch,wa1,wa2,wa3,wa4)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (tr11=.309016994374947d0,ti11=.951056516295154d0)
      parameter (tr12=-.809016994374947d0,ti12=.587785252292473d0)
c***begin prologue  radf5
c***refer to  rfftf
c***routines called  (none)
c***end prologue  radf5
      dimension       cc(ido,l1,5)           ,ch(ido,5,l1)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
c***first executable statement  radf5
      do 101 k=1,l1
         cr2 = cc(1,k,5)+cc(1,k,2)
         ci5 = cc(1,k,5)-cc(1,k,2)
         cr3 = cc(1,k,4)+cc(1,k,3)
         ci4 = cc(1,k,4)-cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2+cr3
         ch(ido,2,k) = cc(1,k,1)+tr11*cr2+tr12*cr3
         ch(1,3,k) = ti11*ci5+ti12*ci4
         ch(ido,4,k) = cc(1,k,1)+tr12*cr2+tr11*cr3
         ch(1,5,k) = ti12*ci5-ti11*ci4
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      if((ido-1)/2.lt.l1) go to 104
      do 103 k=1,l1
cdir$ ivdep
         do 102 i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            dr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            di4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            dr5 = wa4(i-2)*cc(i-1,k,5)+wa4(i-1)*cc(i,k,5)
            di5 = wa4(i-2)*cc(i,k,5)-wa4(i-1)*cc(i-1,k,5)
            cr2 = dr2+dr5
            ci5 = dr5-dr2
            cr5 = di2-di5
            ci2 = di2+di5
            cr3 = dr3+dr4
            ci4 = dr4-dr3
            cr4 = di3-di4
            ci3 = di3+di4
            ch(i-1,1,k) = cc(i-1,k,1)+cr2+cr3
            ch(i,1,k) = cc(i,k,1)+ci2+ci3
            tr2 = cc(i-1,k,1)+tr11*cr2+tr12*cr3
            ti2 = cc(i,k,1)+tr11*ci2+tr12*ci3
            tr3 = cc(i-1,k,1)+tr12*cr2+tr11*cr3
            ti3 = cc(i,k,1)+tr12*ci2+tr11*ci3
            tr5 = ti11*cr5+ti12*cr4
            ti5 = ti11*ci5+ti12*ci4
            tr4 = ti12*cr5-ti11*cr4
            ti4 = ti12*ci5-ti11*ci4
            ch(i-1,3,k) = tr2+tr5
            ch(ic-1,2,k) = tr2-tr5
            ch(i,3,k) = ti2+ti5
            ch(ic,2,k) = ti5-ti2
            ch(i-1,5,k) = tr3+tr4
            ch(ic-1,4,k) = tr3-tr4
            ch(i,5,k) = ti3+ti4
            ch(ic,4,k) = ti4-ti3
  102    continue
  103 continue
      return
  104 do 106 i=3,ido,2
         ic = idp2-i
cdir$ ivdep
         do 105 k=1,l1
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            dr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            di4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            dr5 = wa4(i-2)*cc(i-1,k,5)+wa4(i-1)*cc(i,k,5)
            di5 = wa4(i-2)*cc(i,k,5)-wa4(i-1)*cc(i-1,k,5)
            cr2 = dr2+dr5
            ci5 = dr5-dr2
            cr5 = di2-di5
            ci2 = di2+di5
            cr3 = dr3+dr4
            ci4 = dr4-dr3
            cr4 = di3-di4
            ci3 = di3+di4
            ch(i-1,1,k) = cc(i-1,k,1)+cr2+cr3
            ch(i,1,k) = cc(i,k,1)+ci2+ci3
            tr2 = cc(i-1,k,1)+tr11*cr2+tr12*cr3
            ti2 = cc(i,k,1)+tr11*ci2+tr12*ci3
            tr3 = cc(i-1,k,1)+tr12*cr2+tr11*cr3
            ti3 = cc(i,k,1)+tr12*ci2+tr11*ci3
            tr5 = ti11*cr5+ti12*cr4
            ti5 = ti11*ci5+ti12*ci4
            tr4 = ti12*cr5-ti11*cr4
            ti4 = ti12*ci5-ti11*ci4
            ch(i-1,3,k) = tr2+tr5
            ch(ic-1,2,k) = tr2-tr5
            ch(i,3,k) = ti2+ti5
            ch(ic,2,k) = ti5-ti2
            ch(i-1,5,k) = tr3+tr4
            ch(ic-1,4,k) = tr3-tr4
            ch(i,5,k) = ti3+ti4
            ch(ic,4,k) = ti4-ti3
  105    continue
  106 continue
      return
      end
      subroutine radfg(ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (tpi=6.28318530717959d0)
c***begin prologue  radfg
c***refer to  rfftf
c***routines called  (none)
c***end prologue  radfg
      dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          ,
     1                c1(ido,l1,ip)          ,c2(idl1,ip),
     2                ch2(idl1,ip)           ,wa(*)
c***first executable statement  radfg
      arg = tpi/dfloat(ip)
      dcp = cos(arg)
      dsp = sin(arg)
      ipph = (ip+1)/2
      ipp2 = ip+2
      idp2 = ido+2
      nbd = (ido-1)/2
      if (ido .eq. 1) go to 119
      do 101 ik=1,idl1
         ch2(ik,1) = c2(ik,1)
  101 continue
      do 103 j=2,ip
         do 102 k=1,l1
            ch(1,k,j) = c1(1,k,j)
  102    continue
  103 continue
      if (nbd .gt. l1) go to 107
      is = -ido
      do 106 j=2,ip
         is = is+ido
         idij = is
         do 105 i=3,ido,2
            idij = idij+2
            do 104 k=1,l1
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  104       continue
  105    continue
  106 continue
      go to 111
  107 is = -ido
      do 110 j=2,ip
         is = is+ido
         do 109 k=1,l1
            idij = is
cdir$ ivdep
            do 108 i=3,ido,2
               idij = idij+2
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  108       continue
  109    continue
  110 continue
  111 if (nbd .lt. l1) go to 115
      do 114 j=2,ipph
         jc = ipp2-j
         do 113 k=1,l1
cdir$ ivdep
            do 112 i=3,ido,2
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  112       continue
  113    continue
  114 continue
      go to 121
  115 do 118 j=2,ipph
         jc = ipp2-j
         do 117 i=3,ido,2
            do 116 k=1,l1
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  116       continue
  117    continue
  118 continue
      go to 121
  119 do 120 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  120 continue
  121 do 123 j=2,ipph
         jc = ipp2-j
         do 122 k=1,l1
            c1(1,k,j) = ch(1,k,j)+ch(1,k,jc)
            c1(1,k,jc) = ch(1,k,jc)-ch(1,k,j)
  122    continue
  123 continue
c
      ar1 = one
      ai1 = zero
      do 127 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 124 ik=1,idl1
            ch2(ik,l) = c2(ik,1)+ar1*c2(ik,2)
            ch2(ik,lc) = ai1*c2(ik,ip)
  124    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 126 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 125 ik=1,idl1
               ch2(ik,l) = ch2(ik,l)+ar2*c2(ik,j)
               ch2(ik,lc) = ch2(ik,lc)+ai2*c2(ik,jc)
  125       continue
  126    continue
  127 continue
      do 129 j=2,ipph
         do 128 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+c2(ik,j)
  128    continue
  129 continue
c
      if (ido .lt. l1) go to 132
      do 131 k=1,l1
         do 130 i=1,ido
            cc(i,1,k) = ch(i,k,1)
  130    continue
  131 continue
      go to 135
  132 do 134 i=1,ido
         do 133 k=1,l1
            cc(i,1,k) = ch(i,k,1)
  133    continue
  134 continue
  135 do 137 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 136 k=1,l1
            cc(ido,j2-2,k) = ch(1,k,j)
            cc(1,j2-1,k) = ch(1,k,jc)
  136    continue
  137 continue
      if (ido .eq. 1) return
      if (nbd .lt. l1) go to 141
      do 140 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 139 k=1,l1
cdir$ ivdep
            do 138 i=3,ido,2
               ic = idp2-i
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  138       continue
  139    continue
  140 continue
      return
  141 do 144 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 143 i=3,ido,2
            ic = idp2-i
            do 142 k=1,l1
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  142       continue
  143    continue
  144 continue
      return
      end
