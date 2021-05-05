      program findheadw
c	
c	Made from tracerays2.f to separate those time peacks that are from
c	head waves from that from direct waves
c	Back-propagates ray till it hits basement. If so writes the time peack
c	to head waves peacks file, otherwise - to direct waves peack file
c	
c	Slightly modified from cover.f by J.Hole to provide information about
c       all rays paths
c	S.Tikhotski, Strasbourg,06.2004.
c
c     -----------COMMENT TO COVER.F by J.HOLE------------------------------
c     Given shotpoint location within a sampled 3d traveltime grid,
c        find the ray from the shotpoint to the receiver(source),
c        and find how the traveltime residual affects each cell.
c        Input:   3d travel time field:   nx,ny,nz,t(i,j,k),h,x0,y0,z0,
c                                         xs,ys,zs
c                 shotpoints:   line,spn,xsh,ysh,zsh
c                 traveltime picks:   line,spn,tpick
c                 optional:  ray coverage from previous receivers(sources)
c                               to be added to this rec.:   nrays(i,j,k)
c                            slowness purturbation from previous:   du(i,j,k)
c        Output:  length of each ray
c                 ray coverage array, nrays(i,j,k), defined as the number of 
c                     rays through the model cell between
c                     t(i,j,k) and t(i+1,j+1,k+1)
c                 slowness purturbation array,du(i,j,k), defined as the sum
c                     of all rays' purturbations within each cell
c                     For each ray, assign a slowness purturbation along the
c                         entire ray, of value
c                                du = (traveltime residual)/(raylength)
c
c     J.Hole  March, 1991
c             April, 1991    changed source cube coverage definitions
c             April, 1991    changed for rays caught on faces and edges
c             October, 1991  changed to allow time statics (from location file)
c             February, 1994 changed to properly handle 2d modelling
c                             note: dimensions must be big enough to handle
c                              ny=2 when using 2d ny=1 (same for nx,nz)
c             February, 1994 changed to 1d indexing for convenience
c             late, 1995     changed to consistent units (km)
c                            allowed max offset in timefield (maxoff in punch.c)
c                            changed which statistics are computed
c             April, 1996    added weight factors
c
c     ****** replace all recl=4* (on Sun) with recl=1* (on DEC)
c
      implicit logical(a-z)
      integer      nxyzmx,nsmax,nxymax
      parameter    (nxyzmx=12000000,nsmax=50000,nxymax=140000)
      integer      nx,ny,nz,i,j,k,is,js,ks,ish,jsh,ksh,iii,iiiii,
     *             nshot,ngshot,line,spn,spnp,nseg,cr1,
     *             md,addcov,iscell,jscell,kscell,ist,nstat,
     *             nrays(nxyzmx),iseg(0:nsmax),nk,nkj,nk2,nj,nj2,
     *             pwflag,r1flag,qq,qqq,nshot_start,
     *             bflag,crossed_flag,nbshot,ix2d,iy2d,m,ip,jp,kp
      real         t(nxyzmx),du(nxyzmx),zb(nxymax),v1(nxyzmx),
     *             dum,iwght,ipwght,weight
      double precision x0,y0,z0,x,y,z,h,xi,yj,zk,xs,ys,zs,
     *             gradt(3),dd(3),ddb(3),d,length,vred,dist,tstat(5),
     *             dlen,tpick,dt,dtsum,dt2sum,dusum,du2sum,fx,fy,fz,
     *             duray,length_add,ray1(3),zbavg,xbas,ybas,zbas,tt,
     *             v,cos1,cos2,cosn,vel1,vel2,dtdz,nb(3),gradb(3),
     *             xip,yjp,zkp,ray2(3)
      character*80 sfile,tfile,rfile,ufile,pfile,zfile,
     *             crossfile,v1file,v2file,hwpeakfile,dirpeakfile
c
      write (*,*) 'findheadw (based on cover by J.Hole)'
      write (*,*) 'upgraded to be a part of GRAS algorithm'
c
c     NOW READING THE MODEL DIMENTIONS AND GEOMETRY FROM FILE 
      open (21,file='model-geometry.gras',status='old')
      read (21,*) nx,ny,nz
      if  ((nx*ny*nz.gt.nxyzmx).or.(nx*ny.gt.nxymax))  then
         write (*,*) '***** ERROR:  dimensions are too big'
         stop
      end if
      read (21,*) x0,y0,z0
      read (21,*) h
      close (21)
c      
      write (*,*) 'input the shotpoint locations filename'
      read (*,1010) sfile
      write (*,*) 'input # of time statics (max 5) in the above file'
      read (*,*) nstat
      if  (nstat.gt.5)  then
         write (*,*) 'too many statics'
         stop
      endif
      write (*,*) 'input the traveltime picks filename'
      read (*,1010) pfile
      write (*,*) 'input the reducing velocity for the picks;'
      write (*,*) 'also, input integer weight for data'
      write (*,*) '  if iwght<0, use abs(iwght)*ipwght'
      write (*,*) '  where ipwght is read from fifth column of pickfile'
      read (*,*) vred, iwght
      if (iwght.eq.0) then
         write (*,*) '***** ERROR:  data weight = 0 is inappropriate'
         stop
      endif
      write (*,*) 'input the 3d traveltime filename'
      read (*,1010) tfile
      write (*,*) 'input source xs,ys,zs from forward modelling'
      read (*,*) xs,ys,zs
      write (*,*) xs,ys,zs
      write (*,*) 'input basement depth file'
      read (*,*) zfile
      write (*,*) 'input the filename for the head wave peaks'
      read (*,*) hwpeakfile
      write (*,*) 'input the filename for the direct waves peaks'
      read (*,*) dirpeakfile
      write (*,*) 'input sediments velocity model file name'
      read (*,*) v1file
      write (*,*) 'input crustal velocity model file name'
      read (*,*) v2file
      write (*,*) 'input start shot number'
      read (*,*) nshot_start
c
      if  (vred.lt.0.001) then
         write (*,*) 'Warning:  reducing velocity not used'
      endif
      if (iwght.lt.0) then
         iwght = iwght * -1
         write (*,*) 'data weight =',iwght,' * fifth column of pickfile'
         pwflag = 1
      else
         ipwght = 1
         write (*,*) 'data weight =',iwght
         pwflag = 0
      endif
c
c     OPEN SHOT POSITION FILE
      open (21,file=sfile,status='old')
c
C     OPEN AND READ 3D TRAVELTIME FILE
      open (22,file=tfile,form='unformatted',access='direct',
     *         recl=4*nx,status='old')
      do 10, k=1,nz
         nk = nx*ny*(k-1)
         nk2 = (k-1)*ny
         do 5, j=1,ny
            nkj = nk+nx*(j-1)
            read (22,rec=nk2+j) (t(nkj+i),i=1,nx)
 5       continue
 10   continue
      close (22)
c
c     TEST FOR 2D MODELS
      ix2d=0
      iy2d=0
      if  (nx.eq.1)  then
         x0 = x0 - h/2.
         nx = 2
	 ix2d=1
         write (*,*) '2d model encountered (nx=1)'
         if  (nx*ny*nz.gt.nxyzmx)  then
            write (*,*) '***** ERROR:  dimensions are too big (for 2d)'
            stop
         end if
         do 310, i=1,ny*nz
            du(i) = t(i)
 310     continue
         do 320, k=1,nz
            nk = nx*ny*(k-1)
            nk2 = 1*ny*(k-1)
            do 325, j=1,ny
               iiiii = nk+nx*(j-1)+1
               t(iiiii) = du(nk2+1*(j-1)+1)
               iiiii = iiiii+1
               t(iiiii) = du(nk2+1*(j-1)+1)
 325        continue
 320     continue
      endif
      if  (ny.eq.1)  then
         y0 = y0 - h/2.
         ny = 2
	 iy2d = 1
         write (*,*) '2d model encountered (ny=1)'
         if  (nx*ny*nz.gt.nxyzmx)  then
            write (*,*) '***** ERROR:  dimensions are too big (for 2d)'
            stop
         end if
         do 330, i=1,nx*nz
            du(i) = t(i)
 330     continue
         do 340, k=1,nz
            nk = nx*ny*(k-1)
            nk2 = nx*1*(k-1)
            do 345, i=1,nx
               iiiii = nk+i
               t(iiiii) = du(nk2+i)
               iiiii = iiiii+nx
               t(iiiii) = du(nk2+i)
 345        continue
 340     continue
      endif
      if  (nz.eq.1)  then
         z0 = z0 - h/2.
         nz = 2
         write (*,*) '2d model encountered (nz=1)'
         if  (nx*ny*nz.gt.nxyzmx)  then
            write (*,*) '***** ERROR:  dimensions are too big (for 2d)'
            stop
         end if
         nk = nx*ny
         do 350, i=1,nx
            do 360, j=1,ny
               iiiii = nk+nx*(j-1)+i
               t(iiiii) = t(nx*(j-1)+i)
 360        continue
 350     continue
      endif
c
c     OPEN THE TRAVELTIME PICKS FILE
      open (23,file=pfile,status='old')
c
c     FIND THE INTEGER SOURCE POINT
      is = nint((xs-x0)/h) + 1
      js = nint((ys-y0)/h) + 1
      ks = nint((zs-z0)/h) + 1
c     FIND THE CELL CONTAINING THE SOURCE POINT
      iscell = int((xs-x0)/h) + 1
      jscell = int((ys-y0)/h) + 1
      kscell = int((zs-z0)/h) + 1
c
      nshot = nshot_start
      ngshot = 0
      nbshot = 0
      dtsum = 0.
      dt2sum = 0.
      dusum = 0.
      du2sum = 0.

C     OPEN AND READ 2D BASEMENT DEPTH FILE
      open (22,file=zfile,form='unformatted',access='direct',
     *         recl=4*nx,status='old')
      do 45, j=1,ny
         nkj = nx*(j-1)
         read (22,rec=j) (zb(nkj+i),i=1,nx)
 45   continue
      close (22)
c
c     Open file for head wave time peaks
      open (30,file=hwpeakfile,status='unknown')
c     Open file for direct wave time peaks
      open (31,file=dirpeakfile,status='unknown')
c
c
c
c     ******* START A SHOT (FIND RAY AND BASEMENT INTERSECTION) *******
 200  continue
c
      nseg = 0
      iseg(0) = 0
c
 250  continue
C     READ AND FIND SHOT TRAVELTIME PICK
      if (pwflag.eq.1) then
         read (23,*,end=300) spnp,dist,dum,tpick,ipwght
      else
         read (23,*,end=300) spnp,dist,dum,tpick
      endif
C     READ AND FIND SHOT POSITION
      if  (nstat.gt.0)  then
         read (21,*,end=300) line,spn,x,y,z,(tstat(ist),ist=1,nstat)
      else
         read (21,*,end=300) line,spn,x,y,z
      endif
c
 240  continue
      if  (spnp.gt.spn)  then
C        READ AND FIND SHOT POSITION
         if  (nstat.gt.0)  then
            read (21,*,end=300) line,spn,x,y,z,(tstat(ist),ist=1,nstat)
         else
            read (21,*,end=300) line,spn,x,y,z
         endif
         goto240
      endif
      if  (spn.gt.spnp)  then
C        READ AND FIND SHOT TRAVELTIME PICK
         if (pwflag.eq.1) then
            read (23,*,end=300) spnp,dist,dum,tpick,ipwght
         else
            read (23,*,end=300) spnp,dist,dum,tpick
         endif
         goto240
      endif
      if  (abs(tpick).lt.1.e-6)  goto 250
      if  (abs(dist).lt.1.e-3)  goto 250
c
c
c     CHANGED 17.06.2004 FOR USE WITH GRAS FROM if (ipwght.lt.0) goto 250
c     i.e. IN GRAS VERSION ipwgt*abs(iwght) HAS THE MEANING OF PEACK'S STANDART
c     IF ipwght<0 then peak is missed
      if  (ipwght.lt.0) goto 250

c     LOCATE SHOT IN MODEL
c      write (29,*) '******************************************'
c      write (29,*) 'line',line,'  spn',spn
c      write (*,*) 'line',line,'  spn',spn
c      write (29,fmt='(''>'')')
c      write (31,fmt='(''>'')')
c      write (32,fmt='(''>'')')
c      x = x/1000.
c      y = y/1000.
c      z = z/1000.
      nshot = nshot+1
C     FIND CELL
      i = int((x-x0)/h) + 1
      if  (((x-x0).lt.0.) .or. (i.ge.nx))  goto 200
      j = int((y-y0)/h) + 1
      if  (((y-y0).lt.0.) .or. (j.ge.ny))  goto 200
      k = int((z-z0)/h) + 1
      if  (((z-z0).lt.0.) .or. (k.ge.nz))  goto 200
      ish = i
      jsh = j
      ksh = k
c      write (29,*) x,y,z,i,j,k,' shot'
c      write (29,fmt='(2f10.3)') z,x
c      write (31,fmt='(2f10.3)') x,-z
c      write (32,fmt='(2f10.3)') x,y
      xi = h*(i-1) + x0
      yj = h*(j-1) + y0
      zk = h*(k-1) + z0
c
c     CALCULATE THE TRAVELTIME AT THE SHOTPOINT (TRILINEAR INTERPOLATION)
      nk = nx*ny*(k-1)
      nj = nx*(j-1)
      nk2 = nx*ny*k
      nj2 = nx*j
c     DON'T USE SHOT IF IT'S IN REGION OF MODEL WHERE TIMES WERE NOT CALCULATED
C        (SEE PARAMETER maxoff IN punch.c)
      if (t(nk+nj+i).gt.1.e9   .or. t(nk+nj+i+1).gt.1.e9  .or.
     *    t(nk+nj2+i).gt.1.e9  .or. t(nk+nj2+i+1).gt.1.e9 .or.
     *    t(nk2+nj+i).gt.1.e9  .or. t(nk2+nj+i+1).gt.1.e9 .or.
     *    t(nk2+nj2+i).gt.1.e9 .or. t(nk2+nj2+i+1).gt.1.e9 ) goto200
      fx = (x-xi)/h
      fy = (y-yj)/h
      fz = (z-zk)/h
      dt = (1.-fx)*(1.-fy)*(1.-fz)*t(nk+nj+i)
      dt = dt + fx*(1.-fy)*(1.-fz)*t(nk+nj+i+1)
      dt = dt + (1.-fx)*fy*(1.-fz)*t(nk+nj2+i)
      dt = dt + (1.-fx)*(1.-fy)*fz*t(nk2+nj+i)
      dt = dt + fx*fy*(1.-fz)*t(nk+nj2+i+1)
      dt = dt + fx*(1.-fy)*fz*t(nk2+nj+i+1)
      dt = dt + (1.-fx)*fy*fz*t(nk2+nj2+i)
      dt = dt + fx*fy*fz*t(nk2+nj2+i+1)
c      dist = sqrt((xs-x)**2+(ys-y)**2+(zs-z)**2)
c      write (29,*) spn,dist,0.,dt
c
c      tt=dt
c     CALCULATE DT, THE TRAVELTIME RESIDUAL
      if  (vred.gt.0.001)  tpick = tpick + abs(dist)/vred
      if  (nstat.gt.0)  then
         do 50, ist=1,nstat
            tpick = tpick + tstat(ist)
 50      continue
      endif
c     MY DT IS THE NEGATIVE OF THE USUAL DEFINITION --- WHY?
c        BECAUSE T(CALC)-T(OBS) IMPLIES DU = U(CALC)-U(EARTH)
C           BUT THEN U(EARTH) = U(CALC) - DU   [I.E. DU HAS
C           OPPOSITE SIGN THAN INTUITION STATES]
C        THEREFORE: LEAVE LIKE THIS
      dt = tpick - dt
      tt = dt
c
      length = 0.
c
      r1flag=0
      bflag=0
      crossed_flag=0
      nb(3)=-1.
      cr1=0
c
C     *** BACK-PROPAGATE RAY ***
 100  continue
c
c     FIND LOCATION OF CELL
      xip = xi;
      yjp = yj;
      zkp = zk;
      xi = h*(i-1) + x0
      yj = h*(j-1) + y0
      zk = h*(k-1) + z0
c
c     COPIED FROM BACKRAY.F, 06.2004--------------------------------
c     IF LAST CELL ABOVE BASEMENT, SET RAY1=(GRADT IN LAST CELL)
      if  (r1flag.eq.0)  then
         ray2(1) = ray1(1)
	 ray2(2) = ray2(2)
	 ray2(3) = ray2(3)
         ray1(1) = gradt(1)
         ray1(2) = gradt(2)
         ray1(3) = gradt(3)
      end if
c     --------------------------------------------------------------
c
c     IF RAY IN SOURCE CELL, FINISH RAY BY ADDING TO COVERAGE
      if  ((i.eq.iscell .and. j.eq.jscell .and. k.eq.kscell) .or.
     *     (((x-xs)**2+(y-ys)**2+(z-zs)**2).lt.(h/1000.)))  then
c        *** LENGTH OF RAY IN SOURCE CUBE APPROXIMATED BY STRAIGHT LINE
         length_add = sqrt((x-xs)**2+(y-ys)**2+(z-zs)**2)
         length = length + length_add
         nk = (nx-1)*(ny-1)
         iiiii = nk*(ksh-1)+(nx-1)*(jsh-1)+ish
c
c
c         write (29,*) xs,ys,zs,' source'
c         write (29,fmt='(2f10.3)') zs,xs
c         write (31,fmt='(2f10.3)') xs,-zs
c         write (32,fmt='(2f10.3)') xs,ys
         weight = iwght * ipwght
         duray = dt/length
         do 151, iii=0,nseg
c           SINCE ISEG(0)=0, SHOTPOINT IS PROPERLY CONSIDERED
            if  (iseg(iii).eq.-1)  iiiii=iiiii-1
            if  (iseg(iii).eq.1)   iiiii=iiiii+1
            if  (iseg(iii).eq.-20)  iiiii=iiiii-(nx-1)
            if  (iseg(iii).eq.20)   iiiii=iiiii+(nx-1)
            if  (iseg(iii).eq.-300)  iiiii=iiiii-nk
            if  (iseg(iii).eq.300)   iiiii=iiiii+nk
            nrays(iiiii) = nrays(iiiii) + weight
            du(iiiii) = du(iiiii) + duray*weight
 151     continue
c         write (29,1020) line,spn,dt,length
         dtsum = dtsum + abs(dt)
         dt2sum = dt2sum + dt**2
         dusum = dusum + duray
         du2sum = du2sum + (duray)**2
         ngshot = ngshot + 1
c
c	 Write peak to direct waves peaks file
	if (crossed_flag.eq.0) then 
	 if (pwflag.eq.1) then
	    write (31,*) spnp,dist,dum,tpick,ipwght
	 else
	    write (31,*) spnp,dist,dum,tpick
	 end if
	end if
c
c
         goto 200
      end if
c
c     FIND THE RAY GRAD(T)
      nk = nx*ny*(k-1)
      nj = nx*(j-1)
      nk2 = nx*ny*k
      nj2 = nx*j
      gradt(1) = (t(nk+nj+i+1)+t(nk+nj2+i+1)
     *     +t(nk2+nj+i+1)+t(nk2+nj2+i+1)
     *     -t(nk+nj+i)-t(nk+nj2+i)
     *     -t(nk2+nj+i)-t(nk2+nj2+i)) / (4.*h)
      gradt(2) = (t(nk+nj2+i)+t(nk+nj2+i+1)
     *     +t(nk2+nj2+i)+t(nk2+nj2+i+1)
     *     -t(nk+nj+i)-t(nk+nj+i+1)
     *     -t(nk2+nj+i)-t(nk2+nj+i+1)) / (4.*h)
      gradt(3) = (t(nk2+nj+i)+t(nk2+nj+i+1)
     *     +t(nk2+nj2+i)+t(nk2+nj2+i+1)
     *     -t(nk+nj+i)-t(nk+nj+i+1)
     *     -t(nk+nj2+i)-t(nk+nj2+i+1)) / (4.*h)
c     IF RAY IN SOURCE CUBE, USE STRAIGHT RAY FROM SOURCE
      if  (i.ge.(is-2) .and. i.lt.(is+2) .and.
     *     j.ge.(js-2) .and. j.lt.(js+2) .and.
     *     k.ge.(ks-2) .and. k.lt.(ks+2))  then
         gradt(1) = x - xs
         gradt(2) = y - ys
         gradt(3) = z - zs
      endif
c      write (29,*) 'gradt',gradt
c
c
c
c     FIND THE POSSIBLE STEP SCALARS
      if  (gradt(1).gt.0.)  then
         dd(1) = (xi-x)/gradt(1)
         dlen = dd(1)*sqrt(gradt(1)**2+gradt(2)**2+gradt(3)**2)
         if  (dlen.gt.-h/1000.)  then
            if  ((iseg(nseg).eq.1) .or.
     *           ((iseg(nseg)+iseg(nseg-1)+iseg(nseg-2).eq.1).and.
     *            (nseg.ge.3)))  then
c               write (*,*) 'dd=0 fix spn',spn,' ijk',i,j,k,' -x'
               gradt(1) = 0.
               dd(1) = -1e20
            endif
	  end if
      else if  (gradt(1).lt.0.)  then
         dd(1) = (xi+h-x)/gradt(1)
         dlen = dd(1)*sqrt(gradt(1)**2+gradt(2)**2+gradt(3)**2)
         if  (dlen.gt.-h/1000.)  then
            if  ((iseg(nseg).eq.-1) .or.
     *           ((iseg(nseg)+iseg(nseg-1)+iseg(nseg-2).eq.-1).and.
     *            (nseg.ge.3)))  then
c               write (*,*) 'dd=0 fix spn',spn,' ijk',i,j,k,' +x'
               gradt(1) = 0.
               dd(1) = -1e20
            endif
	  end if
      else
         dd(1) = -1e20
      end if
      if  (gradt(2).gt.0.)  then
         dd(2) = (yj-y)/gradt(2)
         dlen = dd(2)*sqrt(gradt(1)**2+gradt(2)**2+gradt(3)**2)
         if  (dlen.gt.-h/1000.)  then
            if  ((iseg(nseg).eq.20) .or.
     *           ((iseg(nseg)+iseg(nseg-1)+iseg(nseg-2).eq.20).and.
     *            (nseg.ge.3)))  then
c               write (*,*) 'dd=0 fix spn',spn,' ijk',i,j,k,' -y'
               gradt(2) = 0.
               dd(2) = -1e20
            endif
	  end if
      else if  (gradt(2).lt.0.)  then
         dd(2) = (yj+h-y)/gradt(2)
         dlen = dd(2)*sqrt(gradt(1)**2+gradt(2)**2+gradt(3)**2)
         if  (dlen.gt.-h/1000.)  then
            if  ((iseg(nseg).eq.-20) .or.
     *           ((iseg(nseg)+iseg(nseg-1)+iseg(nseg-2).eq.-20).and.
     *            (nseg.ge.3)))  then
c               write (*,*) 'dd=0 fix spn',spn,' ijk',i,j,k,' +y'
               gradt(2) = 0.
               dd(2) = -1e20
            endif
          end if
      else
         dd(2) = -1e20
      end if
      if  (gradt(3).gt.0.)  then
         dd(3) = (zk-z)/gradt(3)
         dlen = dd(3)*sqrt(gradt(1)**2+gradt(2)**2+gradt(3)**2)
         if  (dlen.gt.-h/1000.)  then
            if  ((iseg(nseg).eq.300) .or.
     *           ((iseg(nseg)+iseg(nseg-1)+iseg(nseg-2).eq.300).and.
     *            (nseg.ge.3)))  then
c               write (*,*) 'dd=0 fix spn',spn,' ijk',i,j,k,' -z'
               gradt(3) = 0.
               dd(3) = -1e20
            endif
	  end if
      else if  (gradt(3).lt.0.)  then
         dd(3) = (zk+h-z)/gradt(3)
         dlen = dd(3)*sqrt(gradt(1)**2+gradt(2)**2+gradt(3)**2)
         if  (dlen.gt.-h/1000.)  then
            if  ((iseg(nseg).eq.-300) .or.
     *           ((iseg(nseg)+iseg(nseg-1)+iseg(nseg-2).eq.-300).and.
     *            (nseg.ge.3)))  then
c               write (*,*) 'dd=0 fix spn',spn,' ijk',i,j,k,' +z'
               gradt(3) = 0.
               dd(3) = -1e20
            endif
          end if
      else
         dd(3) = -1e20
      end if
c      write (29,*) 'dd',dd
c     DETERMINE WHICH IS DESIRED (SMALLEST ABSOLUTE VALUE; ALL SHOULD
c        BE NEGATIVE, BUT NEAR-ZERO MIGHT BE POSITIVE)
      v = sqrt(gradt(1)**2+gradt(2)**2+gradt(3)**2);
      if  (abs(dd(1)).le.abs(dd(2)).and.abs(dd(1)).le.abs(dd(3)))  then
         md = 1
         d = dd(1)
      else if (abs(dd(2)).le.abs(dd(3)))  then
         md = 2
         d = dd(2)
      else
         md = 3
         d = dd(3)
      end if
c      write (29,*) 'md,d',md,d
c
c
c     ------------------NOW LOOKING FOR CROSSING POINTS--------------------------
c     CROSSPOINT NEAR RECEIVER
      if (crossed_flag.eq.0) then
c      IF BASEMENT IS IN CELL, BUT NOT YET CROSSED...
       if  (((zk+h).ge.
     *     max(zb(nx*(j-1)+i),zb(nx*(j-1)+i+1),zb(nx*j+i),zb(nx*j+i+1)))
     *     .and. (bflag.eq.0))  then
C        IF FIRST TIME, SAVE PREVIOUS GRADT AS RAY ABOVE BASEMENT
         r1flag = 1
c        FIND NORMAL TO BASEMENT
         zbavg = (zb(nx*(j-1)+i)+zb(nx*(j-1)+i+1)
     *            +zb(nx*j+i)+zb(nx*j+i+1))/4.
         nb(1) = (-zb(nx*(j-1)+i)+zb(nx*(j-1)+i+1)
     *            -zb(nx*j+i)+zb(nx*j+i+1))/(2.*h)
         nb(2) = (-zb(nx*(j-1)+i)-zb(nx*(j-1)+i+1)
     *            +zb(nx*j+i)+zb(nx*j+i+1))/(2.*h)
c         write (29,*) 'zbavg,nb',zbavg,nb
C        FIND STEP SCALAR TO BASEMENT
          ddb(1) = ( nb(1)*(x-xi-h/2.) + nb(2)*(y-yj-h/2.)
     *	             - z + zbavg )
     *	          /( gradt(3) - gradt(1)*nb(1) - gradt(2)*nb(2) )
c         write (29,*) 'basement d',dd(1)
         if  (ddb(1).ge.d)  then
c           RAY HIT BASEMENT!
            if  (ddb(1).lt.0.)  then
               xbas = x + ddb(1)*gradt(1)
               ybas = y + ddb(1)*gradt(2)
               zbas = z + ddb(1)*gradt(3)
            else
               xbas = x
               ybas = y
               zbas = z
            end if
c            write (29,*) xbas,ybas,zbas,' basement'
            bflag = 1
         end if
c
c      ELSE IF BASEMENT ABOVE CELL, SAVING THE CROSSPOINT
       else
        if (((zk).ge.
     *     max(zb(nx*(j-1)+i),zb(nx*(j-1)+i+1),zb(nx*j+i),zb(nx*j+i+1)))
     *     .and. (bflag.eq.1))  then
         v = sqrt(ray1(1)**2+ray1(2)**2+ray1(3)**2)
         ray1(1) = ray1(1)/v
         ray1(2) = ray1(2)/v
         ray1(3) = ray1(3)/v
         v = sqrt(gradt(1)**2+gradt(2)**2+gradt(3)**2)
         gradb(1) = gradt(1)/v
         gradb(2) = gradt(2)/v
         gradb(3) = gradt(3)/v
         v = sqrt(nb(1)**2+nb(2)**2+nb(3)**2)
         nb(1) = nb(1)/v
         nb(2) = nb(2)/v
         nb(3) = nb(3)/v
         cosn = -nb(3)
         cos1 = ray1(1)*nb(1) + ray1(2)*nb(2) + ray1(3)*nb(3)
         cos2 = gradb(1)*nb(1) + gradb(2)*nb(2) + gradb(3)*nb(3)
c         write (29,*) line,spn,xbas,ybas,zbas,cos1,cos2,cosn
c         write (29,*)
c        resetting flags to use with second crosspoint
	 r1flag = 0
	 bflag = 0
c        Designate that first crosspoint is done
	 crossed_flag = 1
         nbshot = nbshot + 1
         cr1 = 1;
c
c	 Write peak to headwaves peaks file
	 if (pwflag.eq.1) then
	    write (30,*) spnp,dist,dum,tpick,ipwght
	 else
	    write (30,*) spnp,dist,dum,tpick
	 end if
c
c	Job done, go to next peak
	goto 200
	end if
c      END IF BASEMENT IS IN CELL, BUT NOT YET CROSSED...
       end if
c     ELSE IF CROSSPOINT NEAR RECEIVER
      else
c      ******************************************************************
c      CROSSPOINT NEAR SOURCE
c      ******************************************************************
c
       if (crossed_flag.eq.1) then
c       FIRST TIME IN THE CELL ABOVE BASEMENT
        if  (((zk+h).le.
     *    min(zb(nx*(j-1)+i),zb(nx*(j-1)+i+1),zb(nx*j+i),zb(nx*j+i+1)))
     *     .and.(bflag.eq.0)) then
c        THE PREVIOUS CELL AND GRADIENT WILL BE USED
C        FIND STEP SCALAR TO BASEMENT (BAKWARD DIRECTION)
c        FIND NORMAL TO BASEMENT
         zbavg = (zb(nx*(jp-1)+ip)+zb(nx*(jp-1)+ip+1)
     *            +zb(nx*jp+ip)+zb(nx*jp+ip+1))/4.
         nb(1) = (-zb(nx*(jp-1)+ip)+zb(nx*(jp-1)+ip+1)
     *            -zb(nx*jp+ip)+zb(nx*jp+ip+1))/(2.*h)
         nb(2) = (-zb(nx*(jp-1)+ip)-zb(nx*(jp-1)+ip+1)
     *            +zb(nx*jp+ip)+zb(nx*jp+ip+1))/(2.*h)
          ddb(1) = (nb(1)*(x-xip-h/2.) + nb(2)*(y-yjp-h/2.)
     *	               - z + zbavg )
     *	          /( - gradt(3) - gradt(1)*nb(1) - gradt(2)*nb(2) )
          if  (ddb(1).lt.0.)  then
               xbas = x - ddb(1)*ray2(1)
               ybas = y - ddb(1)*ray2(2)
c              zbas = z - ddb(1)*ray2(3)
	       zbas = min(zb(nx*(j-1)+i),zb(nx*(j-1)+i+1),zb(nx*j+i),zb(nx*j+i+1))
	       write (*,*) '?'
          else
               xbas = x
               ybas = y
               zbas = z
          end if
          bflag = 1
          write (*,*) 'rec ',spn
          v = sqrt(ray2(1)**2+ray2(2)**2+ray2(3)**2)
          ray1(1) = ray2(1)/v
          ray1(2) = ray2(2)/v
          ray1(3) = ray2(3)/v
          v = sqrt(gradt(1)**2+gradt(2)**2+gradt(3)**2)
          gradb(1) = gradt(1)/v
          gradb(2) = gradt(2)/v
          gradb(3) = gradt(3)/v
          v = sqrt(nb(1)**2+nb(2)**2+nb(3)**2)
          nb(1) = nb(1)/v
          nb(2) = nb(2)/v
          nb(3) = nb(3)/v
c         cosn = -nb(3)
c         cos1 = ray2(1)*nb(1) + ray2(2)*nb(2) + ray2(3)*nb(3)
c         cos2 = gradb(1)*nb(1) + gradb(2)*nb(2) + gradb(3)*nb(3)
c         write (29,*) line,spn,xbas,ybas,zbas,cos1,cos2,cosn
c         write (29,*)
	  crossed_flag = 2
c         nbshot = nbshot + 1
	  bflag = 1
          cr1 = 3
        end if
       end if
c     *****************************************************************************
c     END IF CROSSPOINT NEAR SOURCE
c     *****************************************************************************
      end if
c     ----------------------END LOOKING FOR CROSSING POINTS-------------------------
c
c     CONTINUE FINDING RAY
      length_add = -d*sqrt(gradt(1)**2+gradt(2)**2+gradt(3)**2)
      length = length + length_add
      x = x + d*gradt(1)
      y = y + d*gradt(2)
      z = z + d*gradt(3)
c
      if ((cr1.eq.1)) then
        cr1=2
      end if
      if ((cr1.eq.3)) then
        cr1=4
      end if
c
c      write (29,*) x,y,z,i,j,k
c      write (29,fmt='(2f10.3)') z,x
c      write (31,fmt='(2f10.3)') x,-z
c      write (32,fmt='(2f10.3)') x,y
      nseg = nseg + 1
      if  (nseg.gt.nsmax)  then
         write (*,*) 'line',line,' spn',spn,'  ray too long'
         goto 200
      endif
      ip=i
      jp=j
      kp=k
      if  (md.eq.1)  then
         if  (gradt(1).ge.0) then
            i=i-1
            if  (i.lt.1)  goto 200
            iseg(nseg) = -1
         else
            i=i+1
            if  (i.ge.nx)  goto 200
            iseg(nseg) = 1
         end if
      else if  (md.eq.2)  then
         if  (gradt(2).ge.0) then
            j=j-1
            if  (j.lt.1)  goto 200
            iseg(nseg) = -20
         else
            j=j+1
            if  (j.ge.ny)  goto 200
            iseg(nseg) = 20
         end if
      else
         if  (gradt(3).ge.0) then
            k=k-1
            if  (k.lt.1)  goto 200
            iseg(nseg) = -300
         else
            k=k+1
            if  (k.ge.nz)  goto 200
            iseg(nseg) = 300
         end if
      end if
c
c
      goto 100
c 
 300  continue
c
      close (22)
      close (31)
c
      write (*,*) 'finished source:',xs,ys,zs
      write (*,*) '         line:',line
      write (*,*) '         number of shots:',nshot-nshot_start
      write (*,*) '         goodshots:',ngshot,' could be raytraced'
      write (*,*) '   avgabs dt =',real(dtsum/real(ngshot))
      write (*,*) '   rms dt =',real(sqrt(dt2sum/real(ngshot)))
      write (*,*) '   avg du =',real(dusum/real(ngshot))
      write (*,*) '   rms du =',real(sqrt(du2sum/real(ngshot)))
      if (iwght.gt.1 .or. pwflag.eq.1) then
         write (*,*) '   *** statistics do NOT include pick weights'
      endif
      write (*,*) nbshot, 'cross the basement'
c
 1010 format (a80)
 1020 format (1x,i2,i5,f7.3,f8.3)
 1030 format (1x,i3,i6,3f10.4,f15.10,f15.10)
 1040 format (1x,i5,i5,i5,i5,f10.5,f10.5,f10.5,f15.10,f15.10)
 1050 format (1x,i5,i5,i5,i5,f10.5,f10.5,f10.5,f10.5,f10.5,i5)
 1060 format (1x,i5,i5,f10.5,f10.5,f10.5,f15.10,f15.10,f15.10,f15.10)
c
      stop
      end
