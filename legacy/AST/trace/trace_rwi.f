      program tracerays
c
c
c     ---------------------------------------------------------------------
c
c     Made from tracerays2.f to trace reflected rays. Apr., 2005
c     S. Tikhotski, Strasbourg. 
c     Improved for use in the iterative inversion process with multiple
c     sources. June, 2005.
c
c     ---------------------------------------------------------------------	
c       Slightly modified from cover.f by J.Hole to provide information about
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
      parameter    (nxyzmx=6000000,nsmax=50000,nxymax=50000)
      integer      nx,ny,nz,i,j,k,is,js,ks,ish,jsh,ksh,iii,iiiii,
     *             nshot,ngshot,line,spn,spnp,nseg,cr1,
     *             md,addcov,iscell,jscell,kscell,ist,nstat,
     *             nrays(nxyzmx),iseg(0:nsmax),nk,nkj,nk2,nj,nj2,
     *             pwflag,r1flag,qq,qqq,start_sn,
     *             bflag,crossed_flag,nbshot,ix2d,iy2d,m,ip,jp,kp,
     *             kc,ic,jc,just_reflected
      real         t(nxyzmx),du(nxyzmx),zb(nxymax),v1(nxyzmx),
     *             dum,iwght,ipwght,weight
      double precision x0,y0,z0,x,y,z,h,xi,yj,zk,xs,ys,zs,
     *             gradt(3),dd(3),ddb(3),d,length,vred,dist,tstat(5),
     *             dlen,tpick,dt,dtsum,dt2sum,dusum,du2sum,fx,fy,fz,
     *             duray,length_add,ray1(3),zbavg,xbas,ybas,zbas,tt,
     *             v,cos_r,cosn,vel1,dtdz,nb(3),gradb(3),
     *             xip,yjp,zkp,ray2(3)
      character*80 sfile,tinc_file,trlw_file,rfile,ufile,pfile,
     *             tracefile,zfile,crossfile,v1file
c
      write (*,*) 'tracerrlw (based on cover by J.Hole)'
      write (*,*) 'trace reflected rays. S.Tikhotski, 2005'
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
      write (*,*) 'input the 3d INCIDENT traveltime filename'
      read (*,1010) tinc_file
      write (*,*) 'input the 3d REFLECTED traveltime filename'
      read (*,1010) trlw_file
      write (*,*) 'input source xs,ys,zs from forward modelling'
      read (*,*) xs,ys,zs
      write (*,*) xs,ys,zs
c     write (*,*) 'input the filename for the ray coverage'
c     read (*,1010) rfile
c     write (*,*) 'input the filename for the du coverage'
c     read (*,1010) ufile
c     write (*,*) 'enter 1 if coverage files contain data to which'
c     write (*,*) '  the current coverage should be added; else the'
c     write (*,*) '  coverage files will be overwritten'
c     read (*,*) addcov
      write (*,*) 'input reflector depth file'
      read (*,*) zfile
      write (*,*) 'input the filename for the ray traces'
      read (*,*) tracefile
      write (*,*) 'input the filename for the reflection points'
      read (*,*) crossfile
      write (*,*) 'input upper velocity model file name'
      read (*,*) v1file
      write (*,*) 'input the number of rays in preceeding fans'
      read (*,*) start_sn
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
c     OPEN AND READ (IF REQUIRED) 3D COVERAGE FILES
c     if  (addcov .eq. 1)  then
c        open (22,file=rfile,form='unformatted',access='direct',
c    *        recl=4*(nx-1),status='old')
c        open (23,file=ufile,form='unformatted',access='direct',
c    *        recl=4*(nx-1),status='old')
c        do 20, k=1,nz-1
c           nk = (nx-1)*(ny-1)*(k-1)
c           nk2 = (k-1)*(ny-1)
c           do 15, j=1,ny-1
c              nkj = nk+(nx-1)*(j-1)
c              read (22,rec=nk2+j) (nrays(nkj+i),i=1,nx-1)
c              read (23,rec=nk2+j) (du(nkj+i),i=1,nx-1)
c15        continue
c20      continue
c        close (22)
c        close (23)
c     else
c        do 40, k=1,nz-1
c           nk = (nx-1)*(ny-1)*(k-1)
c           do 35, j=1,ny-1
c              nkj = nk+(nx-1)*(j-1)
c              do 30, i=1,nx-1
c                 iiiii = nkj+i
c                 nrays(iiiii) = 0
c                 du(iiiii) = 0.
c30            continue
c35         continue
c40      continue
c     endif
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
      nshot = 0
      ngshot = 0
      nbshot = 0
      dtsum = 0.
      dt2sum = 0.
      dusum = 0.
      du2sum = 0.
c
c     OPEN AND READ 2D BASEMENT DEPTH FILE
      open (22,file=zfile,form='unformatted',access='direct',
     *         recl=4*nx,status='old')
      do 45, j=1,ny
         nkj = nx*(j-1)
         read (22,rec=j) (zb(nkj+i),i=1,nx)
 45   continue
      close (22)
c
c     Open file for ray traces
      open (30,file='rays.temp',status='unknown')
      open (24,status='scratch')
c     Open file for reflection points
      open (31,file=crossfile,status='unknown')
c
c
c
c     ******* START A SHOT (FIND RAY AND REFLECTION POINT) *******
c     ******* HERE STARTS LOOP OVER RAYS *******
 200  continue
c
      nseg = 0
      iseg(0) = 0
c
c     ******* RETURN POINT FOR FINDING PROPER RAYS (SHOTPOINTS) TO CONSIDER *******
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
c     if  (abs(tpick).lt.1.e-6)  goto 250
c     if  (abs(dist).lt.1.e-3)  goto 250
c
c     CHANGED 17.06.2004 FOR USE WITH GRAS FROM if (ipwght.lt.0) goto 250
c     i.e. IN GRAS VERSION ipwgt*abs(iwght) HAS THE MEANING OF PEACK'S STANDART
c     IF ipwght<0 then peak is missed
      if  (ipwght.lt.0) goto 250
c
c     *********** AT THIS POINT PROPER SHOTPOINT SHOULD BE FOUND ***********
c
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
c     FIND CELL THAT CONTAINS THE SHOTPOINT
c     GO TO NEXT SHOTPOINT IF THIS IS OUTSIDE THE VOLUME
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
c     FIND COORDINATES OF THE TOP LEFT FRONT ANGLE OF THE CELL THAT CONTAINS THE SHOTPOINT
      xi = h*(i-1) + x0
      yj = h*(j-1) + y0
      zk = h*(k-1) + z0
c
c     WE WILL START EACH RAY FROM THE SHOT AND FIRST 
c     PROPAGATE IT THROUGH REFLECTED TIMEFIELD TO THE REFLECTOR
c     THUS FIRST OPEN AND READ 3D REFLECTED TRAVELTIME FILE INTO MEMORY
      open (22,file=trlw_file,form='unformatted',access='direct',
     *         recl=4*nx,status='old')
      do 241, kc=1,nz
         nk = nx*ny*(kc-1)
         nk2 = (kc-1)*ny
         do 242, jc=1,ny
            nkj = nk+nx*(jc-1)
            read (22,rec=nk2+jc) (t(nkj+ic),ic=1,nx)
 242       continue
 241   continue
      close (22)
c     ***********************************************************
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
      write (*,*) 'shotpoint',spn
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
c     SAVE INFORMATION ON THE SHOTPOINT IN THE TRACE FILE
      write (30,*) nshot
      if (pwflag.eq.1) then
        write (30,1040) spn,i,j,k,x,y,z,dt,ipwght*iwght
      else
        write (30,1040) spn,i,j,k,x,y,z,dt,iwght
      end if
c
c     BASEMENT IS NOT YET HITED
      bflag=0
      nb(3)=-1.
c     WE START AFTER THE REFLECTION POINT, TRACING UPGOING RAY: cr1=1
      cr1=1
      just_reflected=0
c
c     *** BACK-PROPAGATE RAY ***
c     *** HERE STARTS THE LOOP OVER CELLS HITED BY RAY ***
 100  continue
c
c     FIND LOCATION OF CELL
c
c     SAVE PREVIOUS VALUES HERE
      xip = xi;
      yjp = yj;
      zkp = zk;
c     AND FIND THE NEW ONES
      xi = h*(i-1) + x0
      yj = h*(j-1) + y0
      zk = h*(k-1) + z0
c
c     IF RAY IN SOURCE CELL AND BASEMENT WAS ALREADY HITED, FINISH RAY BY ADDING TO COVERAGE
      if  ((bflag.eq.1) .and. 
     *     ((i.eq.iscell .and. j.eq.jscell .and. k.eq.kscell) .or.
     *     (((x-xs)**2+(y-ys)**2+(z-zs)**2).lt.(h/1000.))))  then
c        *** LENGTH OF RAY IN SOURCE CUBE APPROXIMATED BY STRAIGHT LINE
         length_add = sqrt((x-xs)**2+(y-ys)**2+(z-zs)**2)
         length = length + length_add
         nk = (nx-1)*(ny-1)
         iiiii = nk*(ksh-1)+(nx-1)*(jsh-1)+ish
c
         write (30,1050) nseg,i,j,k,xs,ys,zs,length_add,length,cr1
         write (30,1050) -1,nseg+1,0,0,0.,0.,0.,0.,0.,0
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
         goto 200
      endif
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
      if   ((bflag.eq.1) .and. 
     *     i.ge.(is-2) .and. i.lt.(is+2) .and.
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
           if ((just_reflected.eq.1).and.(-d*gradt(3).ge.1e-6*h)) then
c           NOW FIND THE REFLECTION ANGLE
              v = sqrt(gradt(1)**2+gradt(2)**2+gradt(3)**2)
              gradb(1) = -gradt(1)/v
              gradb(2) = -gradt(2)/v
              gradb(3) = -gradt(3)/v
              v = sqrt(nb(1)**2+nb(2)**2+nb(3)**2)
              nb(1) = nb(1)/v
              nb(2) = nb(2)/v
              nb(3) = nb(3)/v
              cosn = -nb(3)
              cos_r = gradb(1)*nb(1) + gradb(2)*nb(2) + gradb(3)*nb(3)
c             SAVE FOUND VALUES TO FILE
              nbshot = nbshot + 1
	      write (24,1060) nbshot,tt,iwght*ipwght,line,spn,xbas,ybas,
     *	                      zbas,cos_r,cosn
              just_reflected=0
            end if
c
c     ------------------NOW LOOKING FOR REFLECTION POINT--------------------------
c
c     IF BASEMENT IS IN CELL, BUT NOT YET CROSSED...
       if  (((zk+h).ge.
     *     min(zb(nx*(j-1)+i),zb(nx*(j-1)+i+1),zb(nx*j+i),zb(nx*j+i+1)))
     *     .and. (bflag.eq.0))  then
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
c        IF REFLECTION POINT FOUND...
         if  (ddb(1).ge.d)  then
c           AHA! RAY HIT BASEMENT!
            if  (ddb(1).lt.0.)  then
               xbas = x + ddb(1)*gradt(1)
               ybas = y + ddb(1)*gradt(2)
               zbas = z + ddb(1)*gradt(3)
            else
               xbas = x
               ybas = y
               zbas = z
            end if
            write (*,*) xbas,ybas,zbas,' basement'
c           BASEMENT HITED, SET FLAG
            bflag = 1
	    just_reflected=1
c           
c 
c	    NOW WE WILL PROPAGATE RAY BACK TO THE SOURCE THROUGH THE INCIDENCE TIMEFIELD	    
c           THUS OPEN AND READ 3D INCIDENCE TRAVELTIME FILE INTO MEMORY
            open (22,file=tinc_file,form='unformatted',access='direct',
     *               recl=4*nx,status='old')
           do 341, kc=1,nz
             nk = nx*ny*(kc-1)
             nk2 = (kc-1)*ny
             do 342, jc=1,ny
               nkj = nk+nx*(jc-1)
               read (22,rec=nk2+jc) (t(nkj+ic),ic=1,nx)
 342         continue
 341       continue
           close (22)
c          ***********************************************************
c
           length_add = sqrt((x-xbas)**2+(y-ybas)**2+(z-zbas)**2)
           length = length + length_add
           x = xbas
           y = ybas
           z = zbas
c     FIND CELL THAT CONTAINS THE REFLECTION POINT
c     GO TO NEXT SHOTPOINT IF THIS IS OUTSIDE THE VOLUME
           i = int((x-x0)/h) + 1
           j = int((y-y0)/h) + 1
           k = int((z-z0)/h) + 1
c     FIND COORDINATES OF THE TOP LEFT FRONT ANGLE OF THE CELL THAT CONTAINS THE SHOTPOINT
           xi = h*(i-1) + x0
           yj = h*(j-1) + y0
           zk = h*(k-1) + z0	   
c          Saving the cell indices, coordinates of ray's exit point, length of ray in the
c          current cell and total length of ray before the exit point
           write (30,1050) nseg,i,j,k,x,y,z,length_add,length,cr1
c          AND NOW WE WILL TRACE THE DOWNGOUING PART OF THE RAY: cr1=0
           cr1=0
c
c          AND CONTINUE FINDING RAY, NOW IN THE INCIDENCE TIMEFIELD
           goto 100
	 end if
       end if
c      ----------- END OF LOOKING FOR REFLECTION POINT ---------------
c
c     CONTINUE FINDING RAY
      length_add = -d*sqrt(gradt(1)**2+gradt(2)**2+gradt(3)**2)
      length = length + length_add
      x = x + d*gradt(1)
      y = y + d*gradt(2)
      z = z + d*gradt(3)
c
c     Saving the cell indices, coordinates of ray's exit point, length of ray in the
c     current cell and total length of ray before the exit point
      write (30,1050) nseg,i,j,k,x,y,z,length_add,length,cr1
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
c
 300  continue
c     ******* END OF SHOTS *******
c
      close (21)
      close (23)
      rewind (30)
      open (21,file=tracefile,status='unknown')
c     write (21,*) nshot
      do 505, qqq=1,nshot
        write (*,*) 'loop ',qqq
        read (30,*) qq
	qq=qq+start_sn
	write (21,*) qq
	write (*,*) qq,' of ',nshot
	read (30,1040) spn,i,j,k,x,y,z,dt,iwght
	write (21,1040) spn,i,j,k,x,y,z,dt,iwght
 501    continue
          read (30,1050) nseg,i,j,k,x,y,z,length_add,length,cr1
	  write (21,1050) nseg,i,j,k,x,y,z,length_add,length,cr1
	if (nseg.ne.-1) goto 501
 505  continue
c     write (21,*) nshot
      close (21)
      close (30)
c
c     WRITE OUT THE COVERAGE
c     open (22,file=rfile,form='unformatted',access='direct',
c    *recl=4*(nx-1),status='unknown')
c     open (23,file=ufile,form='unformatted',access='direct',
c    *     recl=4*(nx-1),status='unknown')
c     do 520, k=1,nz-1
c        nk = (nx-1)*(ny-1)*(k-1)
c    nk2 = (k-1)*(ny-1)
c        do 510, j=1,ny-1
c           nkj = nk+(nx-1)*(j-1)
c           write (22,rec=nk2+j) (nrays(nkj+i),i=1,nx-1)
c           write (23,rec=nk2+j) (du(nkj+i),i=1,nx-1)
c510     continue
c520  continue
c     close (22)
c     close (23)

c     ******* NOW FINDING DT/DZ *******
      write (*,*) 'finished finding rays, now find dt/dz'
c
      rewind (24)
c
      if (ix2d.eq.1) nx=1
      if (iy2d.eq.1) ny=1
c
c     READ UPPER VELOCITY FILE
      open (21,file=v1file,form='unformatted',access='direct',
     *     recl=4*nx,status='old')
      do 610, k=1,nz
         nk = nx*ny*(k-1)
         nk2 = (k-1)*ny
         do 605, j=1,ny
            nkj = nk+nx*(j-1)
            read (21,rec=nk2+j) (v1(nkj+i),i=1,nx)
 605     continue
 610  continue
      close (21)
      write (*,*) 'Velocities read...'
c
c
c     TEST FOR 2D MODELS
      if  (nx.eq.1)  then
         nx = 2
         do 620, k=nz,1,-1
            nk = nx*ny*(k-1)
            nk2 = 1*ny*(k-1)
            do 625, j=ny,1,-1
               iiiii = nk+nx*(j-1)+1
               v1(iiiii) = v1(nk2+1*(j-1)+1)
               iiiii = iiiii+1
               v1(iiiii) = v1(nk2+1*(j-1)+1)
 625        continue
 620     continue
      endif
      if  (ny.eq.1)  then
         ny = 2
         do 640, k=nz,1,-1
            nk = nx*ny*(k-1)
            nk2 = nx*1*(k-1)
            do 645, i=nx,1,-1
               iiiii = nk+i
               v1(iiiii) = v1(nk2+i)
               iiiii = iiiii+nx
               v1(iiiii) = v1(nk2+i)
 645        continue
 640     continue
      endif
c
c     ******* FIND DT/DZ FOR BASEMENT FOR EACH SHOT *******
c
      write (31,*) nbshot
      
      do 700, m=1,nbshot
c
c        READ RAYS AND BASEMENT INTERSECTION
         read (24,1060) nshot,tt,iwght,line,spn,xbas,ybas,zbas,
     *                  cos_r,cosn
c
c        CALCULATE DT/DZ AND WRITE
         i = int((xbas-x0)/h) + 1
         j = int((ybas-y0)/h) + 1
         k = int((zbas-z0)/h) + 1
         iiiii = nx*ny*(k-1)+nx*(j-1)+i
         vel1 = (v1(iiiii)+v1(iiiii+1)+v1(iiiii+nx)+v1(iiiii+nx*ny)+
     *        v1(iiiii+nx+1)+v1(iiiii+nx*ny+1)+v1(iiiii+nx*ny+nx)+
     *        v1(iiiii+nx*ny+nx+1))/8.
         dtdz = 2.* cosn * ( cos_r/vel1 )
         nshot = nshot + start_sn
         write (31,1070) nshot,tt,iwght,line,spn,xbas,ybas,zbas,
     *                   cos_r,cosn,dtdz
c
 700  continue
c
      close (22)
      close (31)
c
      write (*,*) 'finished source:',xs,ys,zs
      write (*,*) '         line:',line
      write (*,*) '         number of shots:',nshot
      write (*,*) '         goodshots:',ngshot,' could be raytraced'
      write (*,*) '   avgabs dt =',real(dtsum/real(ngshot))
      write (*,*) '   rms dt =',real(sqrt(dt2sum/real(ngshot)))
      write (*,*) '   avg du =',real(dusum/real(ngshot))
      write (*,*) '   rms du =',real(sqrt(du2sum/real(ngshot)))
      if (iwght.gt.1 .or. pwflag.eq.1) then
         write (*,*) '   *** statistics do NOT include pick weights'
      endif
      write (*,*) nbshot, 'rays hit the reflector'
c
 1010 format (a80)
 1020 format (1x,i2,i5,f7.3,f8.3)
 1030 format (1x,i3,i6,3f10.4,f15.10,f15.10,f15.10,f15.10)
 1040 format (1x,i5,i5,i5,i5,f10.5,f10.5,f10.5,f15.10,f15.10)
 1050 format (1x,i5,i5,i5,i5,f10.5,f10.5,f10.5,f10.5,f10.5,i5)
 1060 format (1x,i5,f10.5,f10.5,i5,i5,f10.5,f10.5,f10.5,f15.10,f15.10)
 1070 format (1x,i5,f10.5,f10.5,i5,i5,f10.5,f10.5,f10.5,f15.10,f15.10,
     *           f15.10)
c
      stop
      end
