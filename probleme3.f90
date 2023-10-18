! Copyritht 2009 by Richard Marchand,
! Department of Physics, University of Alberta.
! All rights reserved. This program or any part thereof may not be
! copied, modified or transferred to other parties without the explicit
! written consent of the author.
! The program may be used as is, and the author may not be held
! responsible for any loss or damage caused by its use.


! A ameliorer:
!
!  Procedure a suivre pour definir un probleme:
!  1) Definir les indices d'equation dans la routine inipro
!  2) Definir les equations dans "coefeq"
!  3) Definir les coefficients de conditions aux limites dans "coefcf"
!
!=======================================================================
      include 'mod_proion.foc90'
      include 'mod_promag.foc90'
!     include 'mod_pronav.foc90'
      include 'mod_prorat.foc90'
      include 'mod_pronmeq.foc90'
      include 'mod_propro.foc90'
      include 'mod_proequil.foc90'
!     include 'mod_prophy.foc90'
!=======================================================================
      subroutine coefel(v,e,cotali,vxy,nv,volum,gradxy,nequat,um,nun)
      use mod_proion
      use mod_promag
!     use mod_pronav
      use mod_prorat
      use mod_propro
      use mod_prophy
      use mod_comctr
      use mod_proequil
      implicit none
!
!  calcul des coefficients numeriques et physiques propres a un element
!
!  arguments
      integer v(DimP1),e(DimP1),cotali,nv,nequat,nun
      real vxy(DIM,nv),volum,gradxy(DIM,DimP1),um(nun)
!
!  variables en common
!     include 'prophy.foc90'
      include 'npamax.foc90'
      include 'compar.foc90'
      include 'com_pronav.foc90'
!
!  variables locales
      integer icall,ifail,i,j,k,iposdn,iposps,ii,iposdx,iposdy
      real boris2
! ad hoc alignment variables
      integer iv1(3),iv2(3),imax,v1,v2
      real :: prodmax,vx,vy,vz,uv,ux,uy,uz,vxMax,vyMax,vzMax
! end of ad hoc alignment variables
      real :: teev,dn,ps,dx,dy,rEarth,rh0
      save icall,teev,rEarth,rh0,boris2
!
!  procedures
      real valpar,tcolei,prod
      external valpar,tcolei
!
!  data
      data icall/0/
      data iv1/1,2,3/
      data iv2/2,3,1/
!
!  calcul
!
!
!  1.   Coefficients numeriques specifiques a un element
!
      if(icall.eq.0) then
        icall=1
        teev=valpar('tPlasma','$parametres physiqu',nompa,physpa,npa, &
          ifail)
        rEarth=valpar('rEarth','$parametres physiq',nompa,physpa,npa, &
          ifail)
        rh0=valpar('rh0','$parametres physiq',nompa,physpa,npa,ifail)
        ifail=0
        boris2=valpar('boris','$parametres physiq',nompa,physpa,npa,ifail)
        if(ifail == 1) then
          boris2=1.
        else
          boris2=boris2*boris2
        endif
      endif
!
!  2.   Coefficients physiques specifiques a un element
!  2.1  Champ magnetique
      if(cotali > 0 .and. 1 == 1) then
        v1=v(iv1(cotali))
        v2=v(iv2(cotali))
        vx=vxy(1,v1)-vxy(1,v2)
        vy=vxy(2,v1)-vxy(2,v2)
        uv=sqrt(vx*vx+vy*vy)
        vx=vx/uv
        vy=vy/uv
        prod=(magxyz(1,v1)+magxyz(1,v2))*vx &
            +(magxyz(2,v1)+magxyz(2,v2))*vy
        do i=1,DimP1
          ii=v(i)
          bMag(i)=magxyz(1,ii)*magxyz(1,ii)+magxyz(2,ii)*magxyz(2,ii)
          bMag(i)=sqrt(bMag(i))
          uMag(1,i)=vx*sign(1.,prod)
          uMag(2,i)=vy*sign(1.,prod)
          vAlfv2(i)=bMag(i)*bMag(i)/rh0T(v(i)) *vnorm**2 !physical units
        enddo
      else
        do i=1,DimP1
          ii=v(i)
          bMag(i)=magxyz(1,ii)*magxyz(1,ii)+magxyz(2,ii)*magxyz(2,ii)
          bMag(i)=sqrt(bMag(i))
          uMag(1,i)=magxyz(1,ii)/bMag(i)
          uMag(2,i)=magxyz(2,ii)/bMag(i)
          vAlfv2(i)=bMag(i)*bMag(i)/rh0T(v(i)) *vnorm**2
        enddo
      endif
! Ad hoc check of alignment
      if(1 == 2) then
        imax=-1
        prodmax=0.
        do i=1,3
          v1=v(iv1(i))
          v2=v(iv2(i))
          vx=vxy(1,v1)-vxy(1,v2)
          vy=vxy(2,v1)-vxy(2,v2)
          uv=sqrt(vx*vx+vy*vy)
!         write(6,*)'uv=',uv
          vx=vx/uv
          vy=vy/uv
          ux=magxyz(1,v1)+magxyz(1,v2)
          uy=magxyz(2,v1)+magxyz(2,v2)
          uv=sqrt(ux*ux+uy*uy)
!         write(6,*)'uv=',uv
          ux=ux/uv
          uy=uy/uv
          prod=vx*ux+vy*uy
          write(6,102)i,prod
 102      format('  i=',i2,'  prod=',s,e14.6)
          if(abs(prod) > prodmax)then
            prodmax=abs(prod)
            imax=i
            vxMax=vx
            vyMax=vy
            vzMax=vz
          endif
        enddo
        write(6,101)cotali,imax,prodmax
 101    format('cotaly=',i2,' imax=',i2,'prodmax=',s,e15.6)
        if(prodmax > 0.99) then
          prod=0.
          do i=1,DimP1
            ii=v(i)
            prod=-prod+vxMax*magxyz(1,ii)+vyMax*magxyz(2,ii) &
                 +vzMax*magxyz(DIM,ii)
            prod=sign(1.,prod)
          enddo
          do i=1,DimP1
            uMag(1,i)=vxMax*prod
            uMag(2,i)=vxMax*prod
            uMag(DIM,i)=vxMax*prod
          enddo
        endif
      endif
! end of ad hoc check

!  2.  compute gradient of equilibrium quantities
!  2.1 equilibrium density
      grdRh0=0. !vector
      do j=1,DimP1
      do i=1,DIM
        grdRh0(i)=grdRh0(i)+rh0T(v(j))*gradxy(i,j)
      enddo
      enddo
 
!  2.2 equilibrium pressure
      grdPr0=0. !vector
      do j=1,DimP1
! Hossna add: DIM is dimonsion which is 2 here (x,y)
      do i=1,DIM 
        grdPr0(i)=grdPr0(i)+pr0T(v(j))*gradxy(i,j)
      enddo
      enddo
 
!  2.3 equilibrium magnetic field
      grdB0=0. !vector
      do j=1,DimP1
      do i=1,DIM
        grdB0(i)=grdB0(i)+bmag(j)*gradxy(i,j)
      enddo
      enddo

!  2.4 unit magnetic field
      grdumag=0.
      do k=1,Dim   !summation over components of the unit vector
      do j=1,DimP1 !summation over nodes
      do i=1,DIM   !loop over the components of the gradient
        grdumag(i,k)=grdumag(i,k)+umag(k,j)*gradxy(i,j)
      enddo
      enddo
      enddo
 

      return
      end
!=======================================================================
! Hossna add: this subroutine calculate the coeefficients of equation
      subroutine coefeq(iequat,nequat,temps,vxy,v,um,gradxy,a,b, &
        c,d,z,ifail)
      use mod_proion
      use mod_promag
!     use mod_pronav
      use mod_propro
      use mod_prorat
      use mod_prophy
      use mod_comctr
      use mod_proequil
      implicit none
!
!  calcul des coefficients qui definissent chaque equation du systeme
!  d'equation aux derivees partielles. Pour chaque equation k, ce
!  systeme est de la forme:
!       dU_l
!  z_kl ---- + div[A_kl*U_l + B_kl.grad(U_l)] + D_kl.grad(U_l)
!       dt
!
!       + C_kl*U_l + C_k0 = 0
!
!  ou il y a sommation sur les indices repetes. On note que le terme
!  B_kl est un tenseur d'ordre 2. Les coefficients sont definis en
!  chaque noeud du triangle t
!
!  on considere deux types de conditions frontieres: essentielles
!                                                    naturelles
!  1) conditions essentielles: de la forme q0 = valeur
!  2) conditions naturelles: Le flux sortant a la frontiere est egal a
!     q0 + q1_k*u_k + q2_kl*grad_l(u_k)
!     ou k = indice de l'inconnue
!        l = indice de la coordonnee spatiale

!
!  arguments
      integer iequat,nequat,v(DimP1),ifail
      real temps,vxy(DIM,*),um(*),gradxy(DIM,DimP1), &
        a(DIM,DimP1,nequat),b(DIM,DIM,DimP1,nequat),c(DimP1,0:nequat), &
        d(DIM,DimP1,nequat),z(DimP1,nequat)
!  iequat: indice de l'equation
!  nequat: nombre total d'inconnues (d'equations).
!  vxy(i,j): coordonnees x,y (i=1,2) du noeud j dans un element donne.
!  a(i,j,ieq): vecteur vitesses
!  b(i1,i2,j,ieq): tenseur diffusivite
!  c(j,ieq): vecteur sources
!  d(i,j,ieq): vecteur projection
!
!  variables common
!     include 'prophy.foc90'
      include 'npamax.foc90'
      include 'compar.foc90'
      include 'com_pronav.foc90'
!
!  variables locales
      integer icall,iv,i,j
      real, parameter :: epsi = 0.00
      real geo,rh,gr,gt,gz,vr,vt,vz,br,bt,bz,rPepsi,r,psi,coef,delta
      REAL dt,ps,fact,tPlasma,artvis,difLpArt,artLap,fudLap,damping, &
           dpdpsi,FdFdPsi,psiUp,psiDn,jt,x,y,zz,rEarth,rh0,dipmom,boris2, &
           deltaAlf
      real, parameter :: gamma=5./3.
!     real, parameter :: del=0.020,x0=0.620613,y0=0.846155,period=2.5
!     real, parameter :: del=0.020,x0=0.5835,y0=0.8099,period=2.5
!     real, parameter :: del= 0.20,x0=3.0000,y0=0.0000,period=380.0
      real, parameter :: del= 0.20,x0=3.0000,y0=0.0000,period=13.74
!     real, parameter :: del=0.20,x0=7.5000,y0=0.0000,period=3701.
      REAL temp1,phi,rhmin
      real, parameter :: Mach=6.0
! tempo
      parameter(difLpArt=0.1)
      save icall,dt,geo,tPlasma,artvis,artlap,damping,dpdpsi,FdFdPsi, &
           psiUp,psiDn,rEarth,rh0,dipmom,boris2,deltaAlf,fudLap
!
!  procedures
      real valrpa,valpar
      external valrpa,valpar
!
!  data
      data icall/0/
!
!  calcul
!
!  0.   initialisation de constantes
!
      if(icall.eq.0) then
        icall=1
        dt=valrpa('deltat',nompa,physpa,npa,ifail)
        geo=valrpa('geo',nompa,physpa,npa,ifail)
        tPlasma=valpar('tPlasma','$parametres physiq',nompa,physpa,npa,ifail)
        artvis=valpar('artvis','$parametres physiq',nompa,physpa,npa,ifail)
        artLap=valpar('artLap','$parametres physiq',nompa,physpa,npa,ifail)
        fudLap=valpar('fudLap','$parametres physiq',nompa,physpa,npa,ifail)
        damping=valpar('damping','$parametres physiq',nompa,physpa,npa,ifail)
        dpdpsi=valpar('dpdpsi','$parametres physiq',nompa,physpa,npa,ifail)
        FdFdPsi=valpar('FdFdPsi','$parametres physiq',nompa,physpa,npa,ifail)
        psiUp=valpar('psiUp','$parametres physiq',nompa,physpa,npa,ifail)
        psiDn=valpar('psiDn','$parametres physiq',nompa,physpa,npa,ifail)
        rEarth=valpar('rEarth','$parametres physiq',nompa,physpa,npa,ifail)
        ifail=0
        boris2=valpar('boris','$parametres physiq',nompa,physpa,npa,ifail)
        deltaAlf=valpar('deltaAlf','$parametres physiq',nompa,physpa,npa,ifail)
        if(ifail == 1) then
          boris2=1.
        else
          boris2=boris2*boris2
        endif
        rhmin=valpar('rhmin','$parametres physiqu',nompa,physpa,npa,ifail)

      endif
!
!  0.   Initialisation des coefficients a zero
!
      a=0. !(vecteur)
      b=0. !(vecteur)
      c=0. !(vecteur)
      d=0. !(vecteur)
      z=0. !(vecteur)
      do iv=1,DimP1
        z(iv,iequat)=1.
      enddo
      

!Hossna add:
! iequat means the index of equation
! this below means if the euation is ieqn1(the equation of n1)
      if(iequat.eq.ieqn1) then 
!
!  1.   EQUATION ieqn1 : a-parallel
!
        do i=1,DimP1
! Hossna add
! vxy has two elements: the first one indicates the dimonsion(x or y)
! the second one indicates the index of each point of triangle
          x=vxy(1,v(i))
          y=vxy(2,v(i))
          
!  time derivative = 0
!  Hossna Add: z is the coefficient of time derivative which is one in the 
!  equation for n1
          z(i,iequat)=1.

!  parallel derivative
!  Hossna add:
!  the first element of each coefficient means dimonsion:1 is x, 2 is y.
!  the second element correspond to each points of triangle
!  the third one indicates the corresponding equation
!  ieqm1 is equation corresponding to Mparallel 1
          !!!a(1,i,iequat)=a(1,i,iequat)+Mach
         !! d(2,i,iequat)=d(2,i,iequat)+um((v(i)-1)*nequat+ieqm1)
         !! d(2,i,ieqm1)=d(2,i,ieqm1)+um((v(i)-1)*nequat+iequat)
         !! d(1,i,iequat)=d(1,i,iequat)+um((v(i)-1)*nequat+ieqm2)
         !! d(1,i,ieqm2)=d(1,i,ieqm2)+um((v(i)-1)*nequat+iequat)

          a(1,i,iequat)=a(1,i,iequat)+um((v(i)-1)*nequat+ieqm2)
          a(2,i,iequat)=a(2,i,iequat)+um((v(i)-1)*nequat+ieqm1)

!  artificial parallel diffusion
          b(1,1,i,iequat)=b(1,1,i,iequat)-artlap
          b(2,2,i,iequat)=b(2,2,i,iequat)-artlap

        enddo

      else if(iequat.eq.ieqm1) then
!
!  2.   EQUATION ieqm1 : M1
!
        do i=1,DimP1

!  time derivative = 0
          z(i,iequat)=max(um((v(i)-1)*nequat+ieqn1),rhmin)

!  spatial derivatives
          d(1,i,iequat)=d(1,i,iequat)+um((v(i)-1)*nequat+ieqn1) &
                        *um((v(i)-1)*nequat+ieqm2)
          d(2,i,iequat)=d(2,i,iequat)+um((v(i)-1)*nequat+ieqn1) &
                        *um((v(i)-1)*nequat+iequat)

!  pressure gradient
!  Hossna add: if one of the below is commented out means the interaction between
!  fluid 1 and fluid 2 is ignored
!  If both is in effect means the interaction between two fluids are considered
          d(2,i,ieqn1)=d(2,i,ieqn1)+1.
 !        d(2,i,ieqn2)=d(2,i,ieqn2)+1.

!  artificial parallel diffusion
          b(1,1,i,iequat)=b(1,1,i,iequat)-artvis
          b(2,2,i,iequat)=b(2,2,i,iequat)-artvis

        enddo
      else if(iequat.eq.ieqn2) then
!
!  3.   EQUATION ieqn2: n2
!
        do i=1,DimP1
          
!  time derivative = 0
          z(i,iequat)=1.

!  parallel derivative
          a(1,i,iequat)=a(1,i,iequat)+Mach
 !        a(2,i,iequat)=a(2,i,iequat)+0.5*um((v(i)-1)*nequat+ieqm2)
 !        a(2,i,ieqm2)=a(2,i,ieqm2)+0.5*um((v(i)-1)*nequat+iequat)
          d(2,i,iequat)=d(2,i,iequat)+um((v(i)-1)*nequat+ieqm2)
          d(2,i,ieqm2)=d(2,i,ieqm2)+um((v(i)-1)*nequat+iequat)

!  artificial parallel diffusion
          b(1,1,i,iequat)=b(1,1,i,iequat)-artlap
          b(2,2,i,iequat)=b(2,2,i,iequat)-artlap

        enddo

      else if(iequat.eq.ieqm2) then
!
!  4.   EQUATION ieqm2: M2
!
        do i=1,DimP1

!  time derivative = 0
          z(i,iequat)=max(um((v(i)-1)*nequat+ieqn1),rhmin)

!  spatial derivatives
          d(1,i,iequat)=d(1,i,iequat)+um((v(i)-1)*nequat+ieqn1) &
                        *um((v(i)-1)*nequat+iequat)
          d(2,i,iequat)=d(2,i,iequat)+um((v(i)-1)*nequat+ieqn1) &
                        *um((v(i)-1)*nequat+ieqm1)
 !        d(1,i,iequat)=d(1,i,iequat)+Mach*max(um((v(i)-1)*nequat+ieqn2),rhmin)
 !        d(2,i,iequat)=d(2,i,iequat)+max(um((v(i)-1)*nequat+ieqn2),rhmin) &
 !                      *um((v(i)-1)*nequat+iequat)

!  pressure gradient
!  Hossna add: if one of the below is commented out means the interaction between
!  fluid 1 and fluid 2 is ignored
!  If both is in effect means the interaction between two fluids are considered

!          d(2,i,ieqn1)=d(2,i,ieqn1)+1.
          d(1,i,ieqn1)=d(1,i,ieqn1)+1.

!  artificial parallel diffusion
          b(1,1,i,iequat)=b(1,1,i,iequat)-artvis
          b(2,2,i,iequat)=b(2,2,i,iequat)-artvis

        enddo

      else
        write(6,*)'indice d''equation non prevu dans coefeq: iequat=', &
          iequat
        ifail=1
      endif
!***
!     print*,'sortie de coefeq'
!***
      return

      CONTAINS

        real function step(x)
        real x
        if(x<0) then
          step=0.
        else
          step=1.
        endif
        end function step

      END
!=======================================================================
! to do: modify for v1, v2, i1, vv, ur, uz, ...
      subroutine coefcf(ifr,iequat,nequat,temps,vxy,vv,i1,um,q0,q1, &
        q2,ifail)
      use mod_proion
      use mod_promag
!     use mod_pronav
      use mod_comctr
      use mod_propro
      use mod_prorat
      use mod_prophy
      use mod_proequil
      implicit none
!  coefficients de condition frontiere essentielle et naturelle
!  1) conditions essentielles: de la forme q0 = valeur
!  2) conditions naturelles: Le flux non convectif sortant a la
!     frontiere est egal a
!     q0 + q1_k*u_k + q2_kl*grad_l(u_k)
!     ou k = indice de l'inconnue
!        l = indice de la coordonnee spatiale
!
!  arguments
      integer ifr,iequat,nequat,vv(DIM),i1,ifail
      real temps,vxy(DIM,*),um(*),q0(DIM),q1(DIM,nequat), &
        q2(DIM,DIM,nequat)
!  vv: array of absolute indices of the nodes on the boundary
!  i1: indice (de 1 a 3) du premier noeud: v(i1,t)=v1.
!      cette position est requise pour utiliser les bons indices des
!      quantites calculees dans coefel (comme br, bz, etc.)

!
!  variables en common
!     include 'prophy.foc90'
      include 'npamax.foc90'
      include 'compar.foc90'
      include 'com_pronav.foc90'
!
!  variables locales
      integer i,ieq,lni(DIM),icall
      real, parameter :: epsi = 0.01,del=0.015,x0=0.5835,y0=0.8099,&
                         mach=2.0
      real ux,uy,uz,uu,x,y,z,phi,temp1,vx,vy,vz,fact
      real artvis,rh0,R0,bt0,bp0,psiUp,psiDn
      real bM,vAlfvNorm,boris2,period
      save icall,rh0,artvis,R0,bt0,bp0,psiUp,psiDn,boris2,period
!
!  procedures
      real valpar
      external valpar
      intrinsic abs,sqrt
!
!  data
      data icall/0/,period/2.5/

!
!  calcul
!
!  0.   initialisation de coefficients numeriques
!
      if(icall.eq.0) then
        rh0=valpar('rh0','$parametres physiq',nompa,physpa,npa,ifail)
        bt0=valpar('bt0','$parametres physiq',nompa,physpa,npa,ifail)
        bp0=valpar('bp0','$parametres physiq',nompa,physpa,npa,ifail)
        R0=valpar('R0','$parametres physiq',nompa,physpa,npa,ifail)
        artvis=valpar('artvis','$parametres physiq',nompa,physpa,npa,ifail)
        psiUp=valpar('psiUp','$parametres physiq',nompa,physpa,npa,ifail)
        psiDn=valpar('psiDn','$parametres physiq',nompa,physpa,npa,ifail)
        boris2=valpar('boris','$parametres physiq',nompa,physpa,npa,ifail)
        if(ifail == 1) then
          boris2=1.
        else
          boris2=boris2*boris2
        endif
        icall=1
      endif


!  1.   initialisation des coefficients a zero
!
      q0=0.  !(vecteur)
      q1=0.  !(vecteur)
      q2=0.  !(vecteur)

      do i=1,DIM
        lni(i)=vvoi(i,i1)
      enddo

!
!  1.1  calculs preliminaires
!
      ux=vxy(2,vv(2))-vxy(2,vv(1))
      uz=-(vxy(1,vv(2))-vxy(1,vv(1)))
      uu=sqrt(ux*ux+uz*uz)
      ux=ux/uu
      uz=uz/uu
!
!  2.   calcul des coefficients
!
      if(iequat.eq.ieqn1) then
!
!  2.1  EQUATION ieqn1 : parallel equation
!
      do i=1,DIM
        if(ifr == 1) then! Hossna :this is moon
          q0(i)=1.
        endif
      enddo

      elseif(iequat.eq.ieqm1) then
!
!  2.1  EQUATION ieqm1 : electrostatic potential
!
      do i=1,DIM
        if (ifr == 1) then
           q0(i)=0.
        elseif (ifr == 3) then
               q1(i,iequat)=uz
               !print*,'uz',uz,q1(i,iequat)
               q1(i,ieqm2)=ux
               !print*,'ux',ux,q1(i,ieqm2)
               !stop
        endif
      enddo

      elseif(iequat.eq.ieqn2) then
!
!  2.1  EQUATION ieqn2 : psi1 = perp-laplacian of A-parallel
!
      do i=1,DIM
        if(ifr == 5) then! Hossna :this is moon
          q0(i)=0.
        elseif(ifr == 4) then
          q0(i)=0.
        elseif(ifr == 6) then
          q0(i)=0.
        endif
      enddo

      elseif(iequat.eq.ieqm2) then
!
!  2.1  EQUATION ieqm2 : psi2 = div((1+va^2/c^2)grad(phi))
!
      do i=1,DIM
        if (ifr == 1) then
           q0(i)=Mach
           !print*,q0(i),mach
           !stop
        endif
      enddo


      endif

!
!  Pour garantir la compatibilite avec les versions anterieures du
!  code, on met la diagonale a 1 dans le cas d'une condition
!  essentielle, si tous les q1 sont nuls.
!
      if(typcfr(ifr,iequat) == 1) then
        x=0.
        do ieq=1,nequat
        do i=1,2
          x=x+abs(q1(i,ieq))
        enddo
        enddo
        if(x == 0.) then
          q1(1:DIM,iequat)=1.
        endif
      endif


!     if(ifr.ge.1) then
!       q0(1)=q01(1,ifr,iequat)
!       q0(2)=q01(1,ifr,iequat)
!       q1(1,iequat)=q01(2,ifr,iequat)
!       q1(2,iequat)=q01(2,ifr,iequat)
!     else
!       write(6,*)'Indice de frontiere non prevu dans coefcf: ifr=',
!    .    ifr
!       ifail=1
!     endif




!***
!     print*,'sortie de coefcf'
!***
      return
      end
!=======================================================================
      subroutine ecrpro0(isor,gradxy,vxy,v,e,cotali,bndrid,nv,nt,um, &
        volum)
      use mod_comctr
      use mod_prorat
!     use mod_pronav
      use mod_propro
      use mod_proion
      use mod_promag
      use mod_prophy
      use mod_proequil
      implicit none
!  ecriture des resultats specifiques au probleme
!
!  arguments
      integer nt,nv,isor,v(DimP1,nt),e(DimP1,nt),cotali(nt),bndrid(DimP1,nt)
      real vxy(DIM,nv),um(*),gradxy(DIM,DimP1,nt),volum(nt)

!  variables common
!     include 'prophy.foc90'
      include 'npamax.foc90'
      include 'compar.foc90'
      include 'com_pronav.foc90'
!
!  variables locales
      integer i,nequat,ifail,nvoisi
      real zz
      integer nflx,nun,geo, &
        nfrnts
      real temps
      real, allocatable :: grad(:,:),a(:,:,:),b(:,:,:,:),c(:,:), &
        d(:,:,:),z(:,:),q1(:,:),q2(:,:,:),flux(:),fluint(:,:)
      real dpdpsi,psiUp,rh0
      real xm,ym,zm,vyvy
!
!  procedures
      real valrpa,valpar
      external coefeq,coefcf,coefel,valrpa,valpar
!
!  calcul
!
!  1.   preliminaires
      nequat=valrpa('nequat',nompa,physpa,npa,ifail)
      nvoisi=valrpa('nvoisi',nompa,physpa,npa,ifail)
      temps=valrpa('temps',nompa,physpa,npa,ifail)
      geo=valrpa('geo',nompa,physpa,npa,ifail)
      nfrnts=valrpa('nfrnts',nompa,physpa,npa,ifail)
      dpdpsi=valpar('dpdpsi','$parametres physiq',nompa,physpa,npa,ifail)
      psiUp=valpar('psiUp','$parametres physiq',nompa,physpa,npa,ifail)
      rh0=valpar('rh0','$parametres physiq',nompa,physpa,npa,ifail)
      nun=nv*nequat
!
!  2.   allocation

!     Les flux: qe, et (gp, qi) pour chaque espece ionique.

      nflx=2+2*nsmax
      allocate(grad(DIM,nequat))
      allocate(a(DIM,DimP1,nequat))
      allocate(b(DIM,DIM,DimP1,nequat))
      allocate(c(0:DimP1,nequat))
      allocate(d(DIM,DimP1,nequat))
      allocate(z(DimP1,nequat))
      allocate(q1(DIM,nequat))
      allocate(q2(DIM,DIM,nequat))
      allocate(flux(nflx))
      allocate(fluint(nflx,nfrnts))

 104  format(i5,2x,es15.6)
 106  format(a,es9.2,a)
 107  format(es15.6)
!
!  3.   ecriture des solutions
!
!
!  3.00 ecriture des etiquettes de solution pour Vu
      if(outputFormat == 'Vu') then
        write(isor,*)'SOLUTION MaSolution( ) ='
        write(isor,*)'{'
!       write(isor,*)'VARIABLE Va( LagrTetra04, Va+0%1, Connec1, Zone1 );'
 !      write(isor,*)'VARIABLE B( LagrTetra04, B+0%1, Connec1, Zone1 );'
        write(isor,*)'VARIABLE vx( LagrTetra04, vx+0%1, Connec1, Zone1 );'
        write(isor,*)'VARIABLE vy( LagrTetra04, vy+0%1, Connec1, Zone1 );'
        write(isor,*)'VARIABLE vz( LagrTetra04, vz+0%1, Connec1, Zone1 );'
        write(isor,*)'VARIABLE bx( LagrTetra04, bx+0%1, Connec1, Zone1 );'
        write(isor,*)'VARIABLE by( LagrTetra04, by+0%1, Connec1, Zone1 );'
        write(isor,*)'VARIABLE bz( LagrTetra04, bz+0%1, Connec1, Zone1 );'
        write(isor,*)'VARIABLE pr( LagrTetra04, pr+0%1, Connec1, Zone1 );'
        write(isor,*)'VARIABLE rh( LagrTetra04, rh+0%1, Connec1, Zone1 );'
        write(isor,*)'};'
      endif
!
!  3.01 Temperature Te
!
 !    write(isor,106)'$dependent variable:Te T=',temps,'$'
 !    do i=1,nv
 !      zz=um((i-1)*nequat+ieqTe)
 !      write(isor,104)i,zz
 !    enddo
!
!  3.01 Alfven velocity
!
!     if(outputFormat == 'Vu') then
!       write(isor,*)' '
!       write(isor,*)'CHAMP Va( ) ='
!       write(isor,*)'{'
!       do i=1,nv
!!        zz=(magxyz(1,i)*magxyz(1,i)+magxyz(2,i)*magxyz(2,i) &
!!           +magxyz(3,i)*magxyz(3,i))/(rh0T(i)*mu0)
!!        zz=zz*bNorm*bNorm/rNorm
!         zz=(magxyz(1,i)*magxyz(1,i)+magxyz(2,i)*magxyz(2,i) &
!            +magxyz(3,i)*magxyz(3,i))/rh0T(i)
!         zz=sqrt(zz)
!         if(abs(zz) < 1.e-99) zz=0.
!         write(isor,107)zz
!       enddo
!       write(isor,*)'};'
!     endif
!
!  3.01 normalized perturbed magnetic field
!
 !    if(outputFormat == 'Vu') then
 !      write(isor,*)' '
 !      write(isor,*)'CHAMP B( ) ='
 !      write(isor,*)'{'
 !      do i=1,nv
 !        zz=(magxyz(1,i)*magxyz(1,i)+magxyz(2,i)*magxyz(2,i) &
 !           +magxyz(DIM,i)*magxyz(DIM,i))
 !        zz=sqrt(zz)
 !        zz=2.*(um((i-1)*nequat+ieqbx)*magxyz(1,i) &
 !              +um((i-1)*nequat+ieqby)*magxyz(2,i) &
 !              +um((i-1)*nequat+ieqbz)*magxyz(DIM,i))/zz
 !        if(abs(zz) < 1.e-99) zz=0.
 !        write(isor,107)zz
 !      enddo
 !      write(isor,*)'};'
 !    endif
!
!  3.01 x component of the velocity
!
!     if(outputFormat == 'Vu') then
!       write(isor,*)' '
!       write(isor,*)'CHAMP vx( ) ='
!       write(isor,*)'{'
!       do i=1,nv
!         zz=um((i-1)*nequat+ieqvx)*vNorm
!         if(abs(zz) < 1.e-99) zz=0.
!         write(isor,107)zz
!       enddo
!       write(isor,*)'};'
!     else
!       write(isor,106)'$dependent variable:vx T=',temps,'$'
!       do i=1,nv
!         zz=um((i-1)*nequat+ieqvx)*vNorm
!         if(abs(zz) < 1.e-99) zz=0.
!         write(isor,104)i,zz
!       enddo
!     endif
!
!  3.01 y component of the velocity
!

 !    if(outputFormat == 'Vu') then
 !      write(isor,*)' '
 !      write(isor,*)'CHAMP vy( ) ='
 !      write(isor,*)'{'
 !      do i=1,nv
 !        zz=um((i-1)*nequat+ieqvy)*vNorm
 !        if(abs(zz) < 1.e-99) zz=0.
 !        write(isor,107)zz
 !      enddo
 !      write(isor,*)'};'
 !    else
 !      write(isor,106)'$dependent variable:vy T=',temps,'$'
 !      do i=1,nv
 !        zz=um((i-1)*nequat+ieqvy)*vNorm
 !        if(abs(zz) < 1.e-99) zz=0.
 !        write(isor,104)i,zz
 !      enddo
 !    endif
!
!  3.01 z component of the velocity
!
 !    if(outputFormat == 'Vu') then
 !      write(isor,*)' '
 !      write(isor,*)'CHAMP vz( ) ='
 !      write(isor,*)'{'
 !      do i=1,nv
 !        zz=um((i-1)*nequat+ieqvz)*vNorm
 !        if(abs(zz) < 1.e-99) zz=0.
 !        write(isor,107)zz
 !      enddo
 !      write(isor,*)'};'
 !    else
 !      write(isor,106)'$dependent variable:vz T=',temps,'$'
 !      do i=1,nv
 !        zz=um((i-1)*nequat+ieqvz)*vNorm
 !        if(abs(zz) < 1.e-99) zz=0.
 !        write(isor,104)i,zz
 !      enddo
 !    endif
!
!  3.01 x component of the magnetic field
!
      if(outputFormat == 'Vu') then
        write(isor,*)' '
        write(isor,*)'CHAMP bx( ) ='
        write(isor,*)'{'
        do i=1,nv
          zz=um((i-1)*nequat+ieqbx)*bNorm
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,107)zz
        enddo
        write(isor,*)'};'
      else
        write(isor,106)'$dependent variable:bx T=',temps,'$'
        do i=1,nv
          zz=um((i-1)*nequat+ieqbx)*bNorm
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,104)i,zz
        enddo
      endif
!
!  3.01 y component of the magnetic field
!
      if(outputFormat == 'Vu') then
        write(isor,*)' '
        write(isor,*)'CHAMP by( ) ='
        write(isor,*)'{'
        do i=1,nv
          zz=um((i-1)*nequat+ieqby)*bNorm
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,107)zz
        enddo
        write(isor,*)'};'
      else
        write(isor,106)'$dependent variable:by T=',temps,'$'
        do i=1,nv
          zz=um((i-1)*nequat+ieqby)*bNorm
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,104)i,zz
        enddo
      endif
!
!  3.01 z component of the magnetic field
!
      if(outputFormat == 'Vu') then
        write(isor,*)' '
        write(isor,*)'CHAMP bz( ) ='
        write(isor,*)'{'
        do i=1,nv
          zz=um((i-1)*nequat+ieqvt)*bNorm
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,107)zz
        enddo
        write(isor,*)'};'
      else
        write(isor,106)'$dependent variable:bz T=',temps,'$'
        do i=1,nv
          zz=um((i-1)*nequat+ieqvt)*bNorm
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,104)i,zz
        enddo
      endif
!
!  3.01 pressure
!
!     if(outputFormat == 'Vu') then
!       write(isor,*)' '
!       write(isor,*)'CHAMP pr( ) ='
!       write(isor,*)'{'
!       do i=1,nv
!         zz=um((i-1)*nequat+ieqpr)*rNorm
!!zz=pr0T(i)*pNorm
!         if(abs(zz) < 1.e-99) zz=0.
!         write(isor,107)zz
!       enddo
!       write(isor,*)'};'
!     else
!       write(isor,106)'$dependent variable:pr T=',temps,'$'
!       do i=1,nv
!         zz=um((i-1)*nequat+ieqpr)*rNorm
!         if(abs(zz) < 1.e-99) zz=0.
!         write(isor,104)i,zz
!       enddo
!     endif
!
!  3.01 perturbed B toroidal
!
      if(outputFormat == 'Vu') then
        write(isor,*)' '
        write(isor,*)'CHAMP rh( ) ='
        write(isor,*)'{'
        do i=1,nv
          zz=um((i-1)*nequat+ieqbt)*rNorm
!zz=rh0T(i)*rNorm
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,107)zz
        enddo
        write(isor,*)'};'
      else
        write(isor,106)'$dependent variable:rh T=',temps,'$'
        do i=1,nv
          zz=um((i-1)*nequat+ieqbt)*rNorm
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,104)i,zz
        enddo
      endif
!
!  7.   Retour
!
      return
      end
!=======================================================================
      subroutine ecrpro(isor,gradxy,vxy,v,e,cotali,bndrid,nv,nt,um, &
        volum)
      use mod_comctr
      use mod_prorat
!     use mod_pronav
      use mod_propro
      use mod_proion
      use mod_promag
      use mod_prophy
      use mod_proequil
      implicit none
!  ecriture des resultats specifiques au probleme
!
!  arguments
      integer nt,nv,isor,v(DimP1,nt),e(DimP1,nt),cotali(nt),bndrid(DimP1,nt)
      real vxy(DIM,nv),um(*),gradxy(DIM,DimP1,nt),volum(nt)

!  variables common
!     include 'prophy.foc90'
      include 'npamax.foc90'
      include 'compar.foc90'
      include 'com_pronav.foc90'
!
!  variables locales
      integer i,nequat,ifail,nvoisi
      real zz
      integer nflx,nun,geo, &
        nfrnts
      real temps
      real, allocatable :: grad(:,:),a(:,:,:),b(:,:,:,:),c(:,:), &
        d(:,:,:),z(:,:),q1(:,:),q2(:,:,:),flux(:),fluint(:,:)
      real dpdpsi,psiUp,rh0,tPlasma
      real xm,ym,zm,vyvy
!
!  procedures
      real valrpa,valpar
      external coefeq,coefcf,coefel,valrpa,valpar
!
!  calcul
!
!  1.   preliminaires
      nequat=valrpa('nequat',nompa,physpa,npa,ifail)
      nvoisi=valrpa('nvoisi',nompa,physpa,npa,ifail)
      temps=valrpa('temps',nompa,physpa,npa,ifail)
      geo=valrpa('geo',nompa,physpa,npa,ifail)
      nfrnts=valrpa('nfrnts',nompa,physpa,npa,ifail)
      dpdpsi=valpar('dpdpsi','$parametres physiq',nompa,physpa,npa,ifail)
      psiUp=valpar('psiUp','$parametres physiq',nompa,physpa,npa,ifail)
      rh0=valpar('rh0','$parametres physiq',nompa,physpa,npa,ifail)
      tPlasma=valpar('tPlasma','$parametres physiq',nompa,physpa,npa,ifail)
      nun=nv*nequat
!
!  2.   allocation

!     Les flux: qe, et (gp, qi) pour chaque espece ionique.

      nflx=2+2*nsmax
      allocate(grad(DIM,nequat))
      allocate(a(DIM,DimP1,nequat))
      allocate(b(DIM,DIM,DimP1,nequat))
      allocate(c(0:DimP1,nequat))
      allocate(d(DIM,DimP1,nequat))
      allocate(z(DimP1,nequat))
      allocate(q1(DIM,nequat))
      allocate(q2(DIM,DIM,nequat))
      allocate(flux(nflx))
      allocate(fluint(nflx,nfrnts))

 101  format(a)
 102  format(a,i6,a)
 103  format(a,2i6)
 104  format(99e15.7)
 105  format(a,es8.2,a)
!
!  3.   ecriture des solutions
!  n1
      write(isor,105)'SCALARS n1_',temps,' float 1'
      write(isor,101)'LOOKUP_TABLE default'
      do i=1,nv
        zz=um((i-1)*nequat+ieqn1) !&
        if(abs(zz) < 1.e-30) zz=0.
        write(isor,104)zz
      enddo

!  m1
      write(isor,101)
      write(isor,105)'SCALARS m1_',temps,' float 1'
      write(isor,101)'LOOKUP_TABLE default'
      do i=1,nv
        zz=um((i-1)*nequat+ieqm1) !&
        if(abs(zz) < 1.e-30) zz=0.
        write(isor,104)zz
      enddo

!  n2
      write(isor,105)'SCALARS n2_',temps,' float 1'
      write(isor,101)'LOOKUP_TABLE default'
      do i=1,nv
        zz=um((i-1)*nequat+ieqn2)
        if(abs(zz) < 1.e-30) zz=0.
        write(isor,104)zz
      enddo

!  m2
      write(isor,105)'SCALARS m2_',temps,' float 1'
      write(isor,101)'LOOKUP_TABLE default'
      do i=1,nv
        zz=um((i-1)*nequat+ieqm2)
        if(abs(zz) < 1.e-30) zz=0.
        write(isor,104)zz
      enddo


!  7.   Retour
!
      return
      end
!=======================================================================
      subroutine inipro(vxy,v,e,cotali,bndrid,nv,nt,u0,um,up,option)
      use mod_proion
      use mod_promag
!     use mod_pronav
      use mod_pronmeq
      use mod_prorat
      use mod_propro
      use mod_comctr
      use mod_prophy
      use mod_proequil
      implicit none
!  initialisation propre au probleme physique
!
!  arguments
      integer, intent(in) :: nv,nt,v(DimP1,nt),e(DimP1,nt),option
      integer, intent(inout) :: cotali(nt),bndrid(DimP1,nt)
      real, intent(in) :: vxy(DIM,nv)
      real, intent(inout) :: u0(*),um(*),up(*)
!
!  variables en common
      include 'npamax.foc90'
      include 'compar.foc90'
!     include 'prophy.foc90'
      include 'com_pronav.foc90'
!
!  variables locales
      integer ngesti,sor2
      parameter(ngesti=10,sor2=8)
      integer :: i,iv,ifail,nequat,nin,iflag,lngpop,iequat,ipos,is
      real :: redmar,zz,tPlasma,dipMom,rEarth,rh0,pr0,boris2
      real :: x,y,fac,L,Lmp,r,rhmin
!     real, parameter :: del=0.020,x0=0.620613,y0=0.846155
      real, parameter :: del=0.20,x0=3.0,y0=0.
      character (LEN=132) :: ligne
      INTEGER ngroup
!
!  procedures
      integer numreq,indtch
      real valrpa,valpar
      external entete,inttr2,numreq,popin,valrpa,valpar,indtch
 
!     Initialisation des parametres
 
      redmar=valrpa('redmar',nompa,physpa,npa,ifail)
      nequat=valrpa('nequat',nompa,physpa,npa,ifail)

      ngroup=valpar('ngroup','$valeurs numeriques',nompa,physpa, &
        npa,ifail)
      dipMom=valpar('dipMom','$parametres physiqu',nompa,physpa, &
        npa,ifail)
      rEarth=valpar('rEarth','$parametres physiqu',nompa,physpa, &
        npa,ifail)
      rh0=valpar('rh0','$parametres physiq',nompa,physpa,npa,ifail)
      tPlasma=valpar('tPlasma','$parametres physiq',nompa,physpa, &
        npa,ifail)
      ifail=0
      boris2=valpar('boris','$parametres physiq',nompa,physpa,npa,ifail)
      if(ifail == 1) then
        boris2=1.
      else
        boris2=boris2*boris2
      endif
      rhmin=valpar('rhmin','$parametres physiqu',nompa,physpa,npa,ifail)
!  calcul

      if(option.eq.1) then  !-------------------------------------------


!  0.1  allocation de la memoire pour proion
        nsmax=0
        nin=1
        call entete(nin,'$ion species',iflag)
        if(iflag == 0) then
          npa=npa+1
          nompa(npa)=' $ion species'
          read(nin,101)ligne
 101      format(a)
          do while(index(ligne,'$fin') == 0)
            if(index(ligne,'atomic_no') /= 0) then
              i=index(ligne,'=') + index(ligne,':')
              if(i == 0) exit
              call rdfrin(99,ligne(i+1:80),iv,ifail)
              npa=npa+1
              nompa(npa)=' atomic_no'
              physpa(npa)=iv
              read(nin,101)ligne
            elseif(index(ligne,'mass') /= 0) then
              i=index(ligne,'=') + index(ligne,':')
              if(i == 0) exit
              call rdfrre(99,ligne(i+1:80),zz,ifail)
              npa=npa+1
              nompa(npa)=' mass'
              physpa(npa)=zz
              read(nin,101)ligne
            elseif(index(ligne,'charge') /= 0) then
              i=index(ligne,'=') + index(ligne,':')
              if(i == 0) exit
              lngpop=132-i
              call popin(ligne(i+1:132),lngpop,iv,ifail)
              do while(ifail == 0)
                nsmax=nsmax+1
                npa=npa+1
                nompa(npa)=' charge'
                physpa(npa)=iv
                call popin(ligne(i+1:132),lngpop,iv,ifail)
              enddo
              read(nin,101)ligne
            else
              write(6,*)'Variable non prevue dans le bloc impurete:'
              write(6,*)ligne
              read(nin,101)ligne
            endif
          enddo
          npa=npa+1
          nompa(npa)=' $fin'
        endif

        allocate(mi(nsmax))
        allocate(zi(nsmax))
        allocate(za(nsmax))


!  0.2  allocation de la memoire pour promag
        allocate(magxyz(DIM,nv))
        allocate(rh0T(nv))
        allocate(pr0T(nv))


!  0.3  allocation de la memoire pour pronav


!  0.4  allocation pour nmequa
        allocate(nmequa(nequat))

!
!  1.   definition des noms de variables pour chaque equation
        nmequa(1)='n1'
        ieqn1=numreq('n1')
        nmequa(2)='m1'
        ieqm1=numreq('m1')
        nmequa(3)='n2'
        ieqn2=numreq('n2')
        nmequa(4)='m2'
        ieqm2=numreq('m2')
!
!  1.1  definition des masses et des charges des ions
        ipos=indtch(nompa,npa,'$ion species')
        is=0
        do while(ipos > 0)
          is=is+1
          ipos=ipos+1
          za(is)=valrpa('atomic_no',nompa(ipos),physpa(ipos), &
            npa-ipos+1,ifail)
          ipos=ipos+1
          mi(is)=valrpa('mass',nompa(ipos),physpa(ipos),npa-ipos+1, &
            ifail)
          ipos=ipos+1
          zi(is)=valrpa('charge',nompa(ipos),physpa(ipos),npa-ipos+1, &
            ifail)
          ipos=ipos+1
          do while(index(nompa(ipos),'charge') /= 0) 
            is=is+1
            za(is)=za(is-1)
            mi(is)=mi(is-1)
            zi(is)=physpa(ipos)
            ipos=ipos+1
          enddo
          i=ipos
          ipos=indtch(nompa(ipos),npa-ipos+1,'no_atomique')-1
          if(ipos >= 0) ipos=ipos+i-1

        enddo

        pr0=rh0*tPlasma*qelec/(mprot*mi(1))
        rNorm=rh0
        pNorm=(rNorm/mprot)*tPlasma*qelec
        lNorm=1.0*rEarth
        bNorm=1.0*3.*mu0*dipMom/(4.*picst*lNorm**3)
        vNorm=sqrt(bNorm*bNorm/(mu0*rNorm))
!pr0=1.
!rNorm=1.
!pNorm=1.
!lNorm=1.
!bNorm=1.
!vNorm=1.
        tNorm=lNorm/vNorm
        print*,'rNorm=',rNorm
        print*,'lNorm=',lNorm
        print*,'bNorm=',bNorm
        print*,'vNorm=',vNorm
        print*,'tNorm=',tNorm

        do iv=1,nv
          x=vxy(1,iv)
          y=vxy(2,iv)
          Lmp=10.
          r=sqrt(x*x+y*y)
          L=r*r*r/(x*x+y*y)
          rh0T(iv)=rh0*(Lmp/L)**3*(L/r)**6 &
            *(1.+0.e2*exp(-(((r-1.)*lnorm-350.e3)/200.e3)**2))  /rNorm
          pr0T(iv)=pr0*(Lmp/L)**3*(L/r)**6 &
            *(1.+0.e2*exp(-(((r-1.)*lnorm-350.e3)/200.e3)**2))  /pNorm
!rh0T(iv)=rh0
!pr0T(iv)=pr0
        enddo

!  7.1  Assign a dipole magnetic field
!       Assume a dipole moment dipMom oriented along the z axis,
!       centered at the origin of the coordinates
        do iv=1,nv
          x=vxy(1,iv)
          y=vxy(2,iv)
          fac=mu0*dipMom/(4.*picst*(x*x+y*y)**2.5)/rEarth**3
          magxyz(1,iv)=fac*3.*x*y /bnorm
          magxyz(DIM,iv)=fac*(2.*y*y-x*x) /bnorm
        enddo


      elseif(option.eq.2) then  !---------------------------------------
!
!  Initialisation specifique des variables dependantes
!

        if(redmar /= 1) then

print*,'rhmin=',rhmin
          do iv=1,nv

!  comment out the following lines when driving from coefep (within volume)
!  or coefcf (from the boundary)
!           x=vxy(1,iv)
!           y=vxy(2,iv)
!           fac=exp(-((x-x0)**2+(y-y0)**2)/del**2)
            fac=0.
            um((iv-1)*nequat+ieqn1)=rhmin
            up((iv-1)*nequat+ieqn1)=rhmin
            u0((iv-1)*nequat+ieqn1)=rhmin
            um((iv-1)*nequat+ieqm1)=0.
            up((iv-1)*nequat+ieqm1)=0.
            u0((iv-1)*nequat+ieqm1)=0.
            um((iv-1)*nequat+ieqn2)=rhmin
            up((iv-1)*nequat+ieqn2)=rhmin
            u0((iv-1)*nequat+ieqn2)=rhmin
            um((iv-1)*nequat+ieqm2)=0.
            up((iv-1)*nequat+ieqm2)=0.
            u0((iv-1)*nequat+ieqm2)=0.
          enddo
        endif

!
!  7.   Initialisation des profils
!


      endif  !----------------------------------------------------------


      end subroutine inipro
!=======================================================================
      subroutine spcpro(vxy,v,e,cotali,bndrid,nv,nt,nvoi,indeq,neqs, &
        nequat,indlin,indvoi,u0,up,nun)
      use mod_propro
      use mod_proion
      use mod_prorat
      use mod_prophy
      use mod_comctr
      implicit none
!  traitement specifique au probleme a chaque pas de temps, dans chaque
!  groupe de solution
!
!  arguments
      integer nv,nvoi,neqs,nequat,nun,indeq(neqs),indlin(nv+1), &
        indvoi(nvoi)
      real u0(nun),up(nun)
      integer, intent(in) :: nt,v(DimP1,nt),e(DimP1,nt),cotali(nt), &
                             bndrid(DimP1,nt)
      real, intent(in) :: vxy(DIM,nv)
!
!  variables en common
!     include 'prophy.foc90'
      include 'npamax.foc90'
      include 'compar.foc90'
!
!  variables locales
      integer iequat,ipos,iposn,iv,icall,ifail,iposdn
      real rhmin,vmax,dt,ps0,dpdpsi,psiUp
      save icall,rhmin,vmax,dt,dpdpsi,psiUp
!  ad hoc
      integer i
      real xm,ym,zm,vyvy
!  end ad hoc
!
!  procedures
      real valpar,valrpa
      external valpar,valrpa
!
!  data
      data icall/0/
!
!  calcul
!
!  0.   Initialisation
      if(icall.eq.0) then
        icall=1
        rhmin=valpar('rhmin','$parametres physiqu',nompa,physpa,npa,ifail)
        vmax=valpar('vmax','$parametres physiqu',nompa,physpa,npa,ifail)
        dt=valrpa('deltat',nompa,physpa,npa,ifail)
        dpdpsi=valpar('dpdpsi','$parametres physiq',nompa,physpa,npa,ifail)
        psiUp=valpar('psiUp','$parametres physiq',nompa,physpa,npa,ifail)
        icall=1
      endif

!  1.   Balayage sur toutes les equations
      do iequat=1,nequat
        
 !      if(iequat == ieqn1 .or. iequat == ieqn2) then
!  1.1   plasma pressure
 !        do iv=1,nv
 !          ipos=(iv-1)*nequat+iequat
 !          up(ipos)=max(rhmin,up(ipos))
 !          u0(ipos)=max(rhmin,up(ipos))
 !        enddo
 !      endif
 !      if(iequat == ieqm1) then
!  1.1   M1
 !        do iv=1,nv
 !          ipos=(iv-1)*nequat+iequat
 !          iposn=(iv-1)*nequat+ieqn1
 !          if(up(iposn) < rhmin) then
 !            up(ipos)=max(0.,up(ipos))
 !            u0(ipos)=max(0.,up(ipos))
 !          endif
 !        enddo
 !      endif
 !      if(iequat == ieqm2) then
!  1.1   M2
 !        do iv=1,nv
 !          ipos=(iv-1)*nequat+iequat
 !          iposn=(iv-1)*nequat+ieqn2
 !          if(up(iposn) < rhmin) then
 !            up(ipos)=min(0.,up(ipos))
 !            u0(ipos)=min(0.,up(ipos))
 !          endif
 !        enddo
 !      endif
      enddo

      return
      end
!=======================================================================
!=======================================================================
!     real function alphacol(coef,zeff)
!     implicit none
!     if(coef == 0.71) then
!       alphacol=0.71*(0.473+zeff)/(1.+zeff*0.473)
!     elseif(coef == 1.96) then
!       alphacol=(1.+1.76*zeff)/(0.51*(1.76+zeff))
!     elseif(coef == 3.16) then
!       alphacol=3.16*(0.24+zeff)/(1.+0.24*zeff)
!     elseif(coef == 4.66) then
!       alphacol=4.66*(1.46+zeff)/(1.+1.46*zeff)
!     else
!       write(6,*)'cas non prevu dans alphacol: coef=',coef
!       alphacol=coef
!     return
!=======================================================================
      real function coulog(ne,teev)
      implicit none
!  logarithme de coulomb d'apres brajinskii
!
!  arguments
      real ne,teev
!
!  ne: densite electronique (particules / m^3)
!  teev: temperature electronique (ev)
!
!  variables locales
      real t0,clgmin
      parameter(t0=44.7,clgmin=1.)
!
!  procedures
      intrinsic log,max
!
!  calcul
      if(teev.ge.t0) then
        coulog=30.31-0.5*log(ne/(teev*teev*teev))
      else
        coulog=32.21-0.5*log(ne/(teev*teev))
      endif
      coulog=max(coulog,clgmin)
      return
      end
!=======================================================================
      integer function numreq(ch)
      use mod_pronmeq
      implicit none
!  calcul de l'indice d'equation d'equation correspondant a la variable
!  dependante de nom ch
!
!  arguments
      character ch*(*)
!
!  variables en common
      include 'npamax.foc90'
      include 'compar.foc90'
!
!  variables locales
      integer i,nequat,ifail
!
!  procedures
      INTERFACE

        real function valrpa(chaine,nompa,physpa,npa,ifail)
        INTEGER npa,ifail
        REAL physpa(npa)
        CHARACTER chaine*(*),nompa(npa)*20
        end function valrpa

      ENDINTERFACE
!
!  calcul
      nequat=valrpa('nequat',nompa,physpa,npa,ifail)
      ifail=0
      do i=1,nequat
        if(index(nmequa(i),ch).gt.0) then
          numreq=i
          return
        endif
      enddo
      write(6,*)'attention: chaine non trouvee dans numreq: ch=',ch
      do i=1,nequat
        write(6,*)'i, nmequa=',i,nmequa(i)
      enddo
      stop
      end
!=======================================================================
      subroutine grdnod(iv,t0,tx,ty,nv,v,e,nt,z,dzdxy,ifail)
      use mod_comctr
      implicit none
!  calcul de du gradient d'une fonction z definie sur une maille
!  triangulaire non structuree.
!
!  arguments
      integer iv,t0,nv,nt,ifail,v(DIM,nt),e(DIM,nt)
      real tx(nv),ty(nv),z(nv),dzdxy(2)
!  z: tableau de la variable dependante
!  dzdxy(j): tableau de la derivee de la variable dependante au point iv
!            j=1,2 pour la derivee par rapport a x ou y.
!
!  variables locales
      integer t,tdep,i,ii,j,idep,ntri,chnge,inouv,tnouv,is
      real denom,grdxy(2,3),x1,y1,x2,y2,x3,y3
!
!  procedures
      integer edg,triloc
      external edg,triloc
!
!  calcul
!
!  1.   On trouve d'abord un triangle dont iv est un sommet. On part
!       du triangle t0
!rm   do t=t0,nt
!rm     do i=1,3
!rm       if(v(i,t).eq.iv) go to 5
!rm     enddo
!rm   enddo
!rm   do t=1,t0-1
!rm     do i=1,3
!rm       if(v(i,t).eq.iv) go to 5
!rm     enddo
!rm   enddo
      t=triloc(tx(iv),ty(iv),tx,ty,v,e,t0)
      if(t.le.0) then
        write(6,*)'retour force dans grdnod: t non trouve'
        write(6,*)'iv, t0=',iv,t0
        ifail=1
        return
      endif
      do i=1,3
        if(v(i,t).eq.iv) go to 5
      enddo
!
      write(6,*)'retour force dans grdnod: t non trouve'
      write(6,*)'iv, t0=',iv,t0
      ifail=1
      return
!
 5    continue
      t0=t
      tdep=t
      idep=i
      ntri=0
      do j=1,2
        dzdxy(j)=0.
      enddo


      chnge=0
 10   continue
      ntri=ntri+1
      x1=tx(v(1,t))
      y1=ty(v(1,t))
      x2=tx(v(2,t))
      y2=ty(v(2,t))
      x3=tx(v(DIM,t))
      y3=ty(v(DIM,t))
      denom=(x3-x2)*(y1-y2)-(y3-y2)*(x1-x2)
      if(denom.le.0.) then
        write(6,*)'element singulier ou mal oriente: t=',t
        ifail=1
        return
      endif
      grdxy(1,1)=(y2-y3)/denom
      grdxy(2,1)=(x3-x2)/denom
      grdxy(1,2)=(y3-y1)/denom
      grdxy(2,2)=(x1-x3)/denom
      grdxy(1,3)=(y1-y2)/denom
      grdxy(2,3)=(x2-x1)/denom
      do ii=1,3
        is=v(ii,t)
        do j=1,2
          dzdxy(j)=dzdxy(j)+z(is)*grdxy(j,ii)
        enddo
      enddo
      if(chnge.eq.0) then
        tnouv=e(i,t)
        if(tnouv.eq.tdep) go to 15
        if(tnouv.gt.0) then
          inouv=edg(tnouv,t,e)
          i=mod(inouv,3)+1
          t=tnouv
          go to 10
        else
          chnge=1
          i=mod(idep+1,3)+1
          tnouv=e(i,tdep)
          if(tnouv.gt.0) then
            inouv=edg(tnouv,tdep,e)
            i=inouv
            t=tnouv
            go to 10
          endif
        endif
      else
        tnouv=e(mod(i+1,3)+1,t)
        if(tnouv.gt.0) then
          inouv=edg(tnouv,t,e)
          i=inouv
          t=tnouv
          go to 10
        endif
      endif
!
 15   continue
      do j=1,2
        dzdxy(j)=dzdxy(j)/ntri
      enddo


      ifail=0
      return
      end
!=======================================================================
      INTEGER FUNCTION ifind(x,t,n,inc)
      IMPLICIT NONE


!..  Cette fonction retourne l'indice i tel que t(i)<=x<t(i+1) sur un
!  tableau t d'indice 1, ..., n par bonds de inc.
!  On suppose ici que mod(n-1,inc)=0.


!  arguments
      INTEGER n,inc
      REAL x, t(n)


!  variables locales
      INTEGER j, jp


!=========================
!.. x: valeur pour laquelle on veut trouver l'indice de t
!.. t: tableau de coordonnees sur lequel on veut interpoler
!.. n: nombre de points dans le tableau
!.. inc: increment
!=========================


!---
!     print*,'x, t(1)=',x,t(1)
!---
      IF (x.LE.t(1)) THEN
         ifind=1
         RETURN


      ELSE IF (x.GE.t(1+n-inc)) THEN
         ifind=n-(2*inc)+1
         RETURN


      ELSE
         ifind=0
         j=n/inc+1


    2    IF (j-ifind.EQ.1) THEN
            ifind=1+inc*(ifind-1)
            RETURN


         ENDIF


         jp=(j+ifind)/2
         IF (x.GT.t(1+inc*(jp-1))) THEN
            ifind=jp


         ELSE
            j=jp


         ENDIF


         GO TO 2
      ENDIF


      END
!=======================================================================
      real function inlirc(x1,x2,ttaux,tx1,nx1,nx1max,tx2,nx2)
      implicit none
!  Interpolation lineaire sur un tableau rectangulaire
!
!  arguments
      integer nx1,nx2,nx1max
      real x1,x2,tx1(nx1),tx2(nx2),ttaux(nx1max,nx2)
!
!  variables locales
      integer ix1,ix2
      real zx1,zx2,w1,w2
!
!  precedures
      integer ifind
      external ifind
!
!  calcul
      ix1=ifind(x1,tx1,nx1,1)
      ix2=ifind(x2,tx2,nx2,1)
      zx1=min(tx1(ix1+1),max(x1,tx1(ix1)))
      zx2=min(tx2(ix2+1),max(x2,tx2(ix2)))
      w1=(zx1-tx1(ix1))/(tx1(ix1+1)-tx1(ix1))
      w2=(zx2-tx2(ix2))/(tx2(ix2+1)-tx2(ix2))
      inlirc=ttaux(ix1,ix2)*(1.-w1)*(1.-w2) &
        +ttaux(ix1+1,ix2)*w1*(1.-w2) &
        +ttaux(ix1,ix2+1)*(1.-w1)*w2 &
        +ttaux(ix1+1,ix2+1)*w1*w2
      return
      end
!=======================================================================
      real function inttr2(x,y,xm,ym,zm,vm,em,t)
      use mod_comctr
      implicit none
!  Interpolation d'une variable dependante zm sur une maille non
!  structuree xm, zm, en un point (x,y).
!
!  arguments
      integer vm(DIM,*),em(DIM,*),t
      real x,y,xm(*),ym(*),zm(*)
!
!  procedures
      integer triloc
      real inttri
      external inttri,triloc
!
!  calculs
!
!  1.   On trouve l'idice du triangle dans lequel se trouve le point
      t=triloc(x,y,xm,ym,vm,em,t)
!
!  2.   On calcule l'interpolation sur ce triangle.
      inttr2=inttri(x,y,xm,ym,zm,vm,t)
!
      return
      end
!=======================================================================
      real function inttri(xp,yp,xt,yt,delta,v,t)
      use mod_comctr
      implicit none


!  Cette fonction renvoit la valeur de la fonction delta au point xp,yp
!  demande. Pour ce faire, on fait passer un plan par les 3 points
!  formant le triangle qui contient xp,yp, et on trouve la valeur de la
!  fonction en determinant la valeur de z (fonction) pour que ce point
!  soit situe sur le plan. Soit v et w, 2 vecteurs directeurs d'un plan;
!  A, un point du plan; et P un point quelconque. Si (v x w).AP=0, alors
!  le point P est situe sur le plan.


!  arguments
      integer v(DIM,*),t
      real xp,yp,xt(*),yt(*),delta(*)


!  variables locales
      real vx,vy,vz,wx,wy,wz,apx,apy


!  procedures: aucune


!..xt,yt: coords des points de triangles.
!..xp,yp: coords du point dont on cherche la valeur de la fonction.
!..delta: valeurs du de la fonction en chaque point.
!..v : vecteur qui contient les sommets de chaque triangle,(ind. du
!..    sommet(1 a 3), ind. du triangle).
!..e : vecteur qui contient les triangles voisins de chaque triangle,
!..    (ind. du triangle voisin(1 a 3), ind. du triangle).
!..t : indice du triangle qui contient le point xp, yp.
!..vx,vy,vz: vecteur allant du sommet 1 au sommet 2 pour le triangle t.
!..wx,wy,wz: vecteur allant du sommet 1 au sommet 3 pour le triangle t.
!..apx,apy: vecteur allant du sommet 1 au point xp,yp.


!  definition des vecteurs v, w et AP.


      vx=xt(v(2,t))-xt(v(1,t))
      vy=yt(v(2,t))-yt(v(1,t))
      vz=delta(v(2,t))-delta(v(1,t))
      wx=xt(v(DIM,t))-xt(v(1,t))
      wy=yt(v(DIM,t))-yt(v(1,t))
      wz=delta(v(DIM,t))-delta(v(1,t))
      apx=xp-xt(v(1,t))
      apy=yp-yt(v(1,t))


!  en effectuant le produit mixte:(v x w).AP=0, et en isolant la
!  composante z (inttri):


      inttri=((-apx*(vy*wz-vz*wy)+apy*(vx*wz-vz*wx))/(vx*wy-vy*wx)) &
               +delta(v(1,t))


      return
      end
!=======================================================================
      subroutine lecpsi(nmapsi,ntrpsi,vpsi,epsi,xmapsi,ymapsi,psi)
      use mod_comctr
      implicit none


!    Cette sous-routine effectue la lecture des points ou psi est
!  defini.
!
!  arguments
      integer, intent(in) :: nmapsi,ntrpsi
      integer, intent(out) :: vpsi(DIM,ntrpsi),epsi(DIM,ntrpsi)
      real, intent(out) :: xmapsi(nmapsi),ymapsi(nmapsi),psi(nmapsi)
!
!  variables locales
      integer ifail,i,j,bidon
      character (LEN=80) :: vari


!  procedures
      external entete
      intrinsic index


!..nmapsi: nb de points de la maille de psi.
!..xmapsi,ymapsi: coords des points de la maille de psi.
!..psi: valeur de la fonction en chaque point de la maille de psi.
!..vari: chaine de caractere qui est lue.


      open(unit=14,file='mailpsi.out',status='old')
      ifail=0
      call entete(14,'$coord',ifail)
      read(14,2)vari
 2    format(A)
!     i=index(vari,'=')
!     read(vari(i+1:80),*)nmapsi
      do i=1,nmapsi
        read(14,*)bidon,xmapsi(i),ymapsi(i)
      enddo
!
      call entete(14,'$elements',ifail)
      read(14,2)vari
!     i=index(vari,'=')
!     read(vari(i+1:80),*)ntrpsi
      do i=1,ntrpsi
        read(14,*)bidon,(vpsi(j,i),j=1,3),(epsi(j,i),j=1,3)
      enddo
!
      call entete(14,'$dependent variable: psi',ifail)
      do i=1,nmapsi
        read(14,*)bidon,psi(i)
      enddo
!
      close(unit=14)
      return
      end
!=======================================================================
      real function tcolei(ne,ni,teev,zi)
      use mod_prophy
      implicit none
!  temps de collision des electrons sur des ions
!
!  arguments
      real ne,ni,teev,zi
!
!  variables en common
!     include 'prophy.foc90'
!
!  variables locales
      integer icall
      real coef
      save icall,coef
!
!  procedures
      real coulog
      external coulog
!
!  data
      data icall/0/
!
!  calcul
      if(icall.eq.0) then
        icall=1
!       coef=3.5e11
        coef=3.*(4.*picst*eps0)**2/(4.*sqrt(2.*picst/melec)*qelec**2.5)
      endif
      tcolei=coef*teev**1.5/(zi*zi*ni*coulog(ne,teev))
      return
      end
!=======================================================================
      real function tcoli0(ne,teev,tiev,mi,ai,zi,nj,tjev,mj,aj,zj)
      use mod_prophy
      implicit none
!  temps de collision des ions sur des ions
!
!  arguments
      real ne,teev,tiev,mi,ai,zi,nj,tjev,mj,aj,zj
!  ne: densite (si)
!  nj: densite des ions sur lesquels on calcule le temps de collision
!  teev: temperature electronique (ev)
!  tiev: temperature ionique (ev)
!  mi: masse ionique (en unites de masse protonique)
!  zi: charge protonique
!
!  variables en common
!     include 'prophy.foc90'
!
!  variables locales
      integer icall
      real coef,sigma
      save icall,coef,sigma
!
!  data
      data icall/0/
!
!  calcul
      if(icall.eq.0) then
        icall=1
        sigma=2.2e-19  !temporaire
        coef=1.018e-4  !sqrt(mprot/1eV)
      endif
      tcoli0=coef*sqrt((mi*tjev+mj*tiev)/(tiev*tiev))/(sigma*nj) &
        *0.25*(1.+max(0.,ai-zi-2.)+1.+max(0.,aj-zj-2))**2
      return
      end
!=======================================================================
      real function tcolii(ne,teev,tiev,mi,zi,nj,tjev,mj,zj)
      use mod_prophy
      implicit none
!  temps de collision des ions sur des ions
!
!  arguments
      real ne,teev,tiev,mi,zi,nj,tjev,mj,zj
!  ne: densite (si)
!  nj: densite des ions sur lesquels on calcule le temps de collision
!  teev: temperature electronique (ev)
!  tiev: temperature ionique (ev)
!  mi: masse ionique (en unites de masse protonique)
!  zi: charge protonique
!
!  variables en common
!     include 'prophy.foc90'
!
!  variables locales
      integer icall
      real coef
      save icall,coef
!
!  procedures
      real coulog
      external coulog
!
!  data
      data icall/0/
!
!  calcul
      if(icall.eq.0) then
        icall=1
!       coef=0.5*sqrt(0.5)* &
!         3.*(4.*picst*eps0)**2/(4.*sqrt(picst/mprot)*qelec**2.5)
        coef=1./1.39e-13
      endif
!     tcolii=coef*tiev**1.5*sqrt(mi)/(zi*zi*nj*coulog(ne,teev))
      tcolii=coef*(mi*tjev+mj*tiev)**1.5 &
        /(sqrt(mi*mj)*zi*zi*zj*zj*nj*coulog(ne,teev))
      return
      end
!=======================================================================
      FUNCTION TRILOC(XP,YP,X,Y,V,E,NUMTRI)
      use mod_comctr
      implicit none
!
!     FUNCTION TRILOC
!
!     PURPOSE:
!     --------
!
!     LOCATE TRIANGLE WHICH ENCLOSES POINT WITH COORDS (XP,YP) USING
!     LAWSON'S SEARCH
!
!     INPUT:
!     -------
!
!     'XP,YP'  - X-Y COORDINATES OF POINT
!
!     'X,Y'    - X-Y COORDINATES OF POINTS AND SUPERTRIANGLE VERTICES
!              - LISTS OF LENGTH NUMPTS+3
!              - LAST THREE LOCATIONS USED TO STORE COORDS OF
!                SUPERTRIANGLE
!
!     'V'      - VERTEX ARRAY FOR TRIANGULATION
!              - VERTICES LISTED IN ANTICLOCKWISE SEQUENCE
!              - VERTICES FOR TRIANGLE J ARE FOUND IN V(I,J) FOR I=1,2,3
!              - FIRST VERTEX IS AT POINT OF CONTACT OF FIRST AND THIRD
!                ADJACENT TRIANGLES
!              - V HAS DIMENSIONS V(3,2*N+1), WHERE N IS THE NUMBER OF
!                POINTS TO BE TRIANGULATED
!
!     'E'      - ADJACENCY ARRAY FOR TRIANGULATION
!              - TRIANGLES ADJACENT TO J ARE FOUND IN E(I,J) FOR I=1,2,3
!              - ADJACENT TRIANGLES LISTED IN ANTICLOCKWISE SEQUENCE
!              - ZERO DENOTES NO ADJACENT TRIANGLE
!              - E HAS DIMENSIONS E(3,2*N+1), WHERE N IS THE NUMBER OF
!                POINTS TO BE TRIANGULATED
!
!     'NUMTRI' - NUMBER OF TRIANGLES IN TRIANGULATION
!
!     'TRILOC' - NOT DEFINED
!
!     OUTPUT:
!     --------
!
!     'XP'YP'  - UNCHANGED
!
!     'X,Y'    - UNCHANGED
!
!     'V'      - UNCHANGED
!
!     'E'      - UNCHANGED
!
!     'NUMTRI' - UNCHANGED
!
!     'TRILOC' - NUMBER OF TRIANGLE CONTAINING POINT WITH COORDS (XP,YP)
!
!     PROGRAMMER:
!     -----------
!
!     S W SLOAN
!
!     LAST MODIFIED:
!
!
!     30 JAN 1986     S W SLOAN
!
!
      INTEGER V(DIM,*),E(DIM,*),NUMTRI,V1,V2,I,T,TRILOC
!
      REAL X(*),Y(*),XP,YP
!
      T=NUMTRI
   10 CONTINUE
      DO 20 I=1,3
        V1=V(I,T)
        V2=V(MOD(I,3)+1,T)
        IF((Y(V1)-YP)*(X(V2)-XP).GT.(X(V1)-XP)*(Y(V2)-YP))THEN
          T=E(I,T)
!rmarchand
          if(t.le.0) then
            triloc=t
            return
          endif
!rmarchand
          GOTO 10
        END IF
   20 CONTINUE
!
!     TRIANGLE HAS BEEN FOUND
!
      TRILOC=T
!
      END
!=======================================================================
      REAL FUNCTION vparcf(cs,x,VExBr,VExBz,ur,uz,bn0,CODE_CL_plaques)
      IMPLICIT NONE


!     Arguments
      INTEGER, INTENT(IN) :: CODE_CL_plaques
      REAL, INTENT(IN) :: cs,x,VExBr,VExBz,ur,uz,bn0


!     Variables locales
      REAL GCterm


      SELECT CASE (CODE_CL_plaques)
 
         CASE(1)
!           Condition de Bohm standard
            vparcf=cs*tanh(x)

         CASE(2)
            GCterm = 0.0
            IF (VExBr*uz-VExBz*ur >= 0.) THEN
!              Condition de Gerhauser-Claasen
               IF (x > -0.1 .AND. x < 0.1) THEN
                  GCterm=-(VExBr*uz-VExBz*ur)/bn0*(x-2./3.*x**3)
               ELSE
                  GCterm=-(VExBr*uz-VExBz*ur)/bn0*tanh(x)**2/x
               ENDIF
            ENDIF
            vparcf=cs*tanh(x) + GCterm
            !write(*,999),x,cs*tanh(x),VExBr,VExBz,GCterm
            !999 format(5(2x,g10.3))

      END SELECT

      RETURN
      END                                                                       
