
!! (wanhui 2003.07.07)
!! (wanhui 2003.11.04)
!! (b.wang 2004.02.15)


       real*8  du  (nx,ny,nl)  !
       real*8  dv  (nx,ny,nl)  ! tendency at present step
       real*8  dtt (nx,ny,nl)  !

       real*8  du0 (nx,ny,nl)  !
       real*8  dv0 (nx,ny,nl)  !
       real*8  dtt0(nx,ny,nl)  ! tendency at step n-1
       real*8  dps0(nx,ny)     !

       real*8  du1 (nx,ny,nl)  !
       real*8  dv1 (nx,ny,nl)  !
       real*8  dtt1(nx,ny,nl)  ! tendency at step n
       real*8  dps1(nx,ny)     !

       real*8  uu  (nx,ny,nl)  !
       real*8  vv  (nx,ny,nl)  !
       real*8  tt  (nx,ny,nl)  ! variables at present step
       real*8  p   (nx,ny)     !
       real*8  ply (nx,ny,nl)  !

       real*8  up  (nx,ny,nl)  !
       real*8  vp  (nx,ny,nl)  !
       real*8  ttp (nx,ny,nl)  ! variables at step n-1
       real*8  pps (nx,ny)     !
       real*8  dps (nx,ny)     !

       real*8  uk  (nx,ny,nl)  !
       real*8  vk  (nx,ny,nl)  !
       real*8  ttk (nx,ny,nl)  ! variables at step n
       real*8  psk (nx,ny)     !

       real*8  dlt1
       real*8  dlt2

       real*8  tb  (nx,ny,nl)
       real*8  cb  (nx,ny,nl)
       real*8  dcb (nx,ny,nl)

       real*8  hps (nx,ny)
       real*8  c0  (nx,ny)
       real*8  cs0 (nx,ny)
       real*8  cb0 (NX,NY,NL)
       real*8  cbs (NX,NY,NL)
       integer nigw(ny),nzad(ny)

       common/comfm2/ du,dv,dtt, du0,dv0,dtt0,dps0, du1,dv1,dtt1,dps1
       common/comfm2/ uu,vv,tt,p,ply, up,vp,ttp,pps,dps, uk,vk,ttk,psk
       common/comfm2/ dlt1,dlt2, tb,cb,dcb
       common/comfm2/ hps,c0,cb0,nigw
