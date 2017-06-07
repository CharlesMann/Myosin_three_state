c
c $Id: dyn21.F 86848 2014-01-29 17:11:37Z tsay $
c
      subroutine umat43v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsv,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,cma,qmat,elsizv,idelev)
c******************************************************************
c|  Livermore Software Technology Corporation (LSTC)              |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2007 LSTC                                      |
c|  All Rights Reserved                                           |
c******************************************************************
      include 'nlqparm'
      include 'iounits.inc'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension cm(*),epsps(*),hsv(nlq,*),dt1siz(*)
      dimension temps(*),crv(101,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      dimension idelev(*)
      logical failels(*)
      character*5 etype
c
c     eps(1)=local x  strain increment
c     eps(2)=local y  strain increment
c     eps(3)=local z  strain increment
c     eps(4)=local xy strain increment
c     eps(5)=local yz strain increment
c     eps(6)=local zx strain increment
c
c     sig(1)=local x  stress
c     sig(2)=local y  stress
c     sig(3)=local z  stress
c     sig(4)=local xy stress
c     sig(5)=local yz stress
c     sig(6)=local zx stress
c
c     material properties for orthotropic elastic material 
c     If transverse elastic properties are not provided, the transverse
c     directions take on the same properties as the longitudinal direction,
c     i.e. behaves as an isotropic material
c
      E1 = cm(1)
      E2 = cm(12)
      E3 = cm(13)
      v21 = cm(14)
      v31=cm(15)
      v32=cm(16)
      G12 = cm(17)
      G23 = cm(18)
      G31 = cm(19)            
c      
      v12 = v21*E1/E2
      v23 = v32*E2/E3
      v13 = v31*E1/E3
c
c     stiffness matrix 
c
      S = (1.-v12*v21-v23*v32-v31*v13-2*v21*v32*v13)/(E1*E2*E3)
c
      C11=(1. - v23*v32)/(E2*E3*S)
      C12=(v21 + v31*v23)/(E2*E3*S)
      C13=(v31 + v21*v32)/(E2*E3*S)
      C22=(1. - v31*v13)/(E1*E3*S)
      C23=(v32 + v31*v12)/(E1*E3*S)
      C33=(1. - v12*v21)/(E1*E2*S)
      C44=G12
      C55=G23
      C66=G31
c
c     write(59,*) 'Stiffness Matrix'
c     write(59,*) C11, '  ',C12,  '  ',  C13      
c     write(59,*) C22, '  ',C23,  '  ',  C33
c     write(59,*) C44, '  ',C55,  '  ',  C66
c     write(59,*) '      '
c
      do  10 i=lft,llt
      if (etype.eq.'solid'.or.etype.eq.'shl_t'.or.
     1     etype.eq.'sld2d'.or.etype.eq.'tshel') then
c
c       integrate 3-D strains (six components)
c       +++++++++++++++++++++++++++++++++++++++++++++++++
c       hsv(i,1)=hsv(i,1)+d1(i)     
c       hsv(i,6)=hsv(i,6)+d2(i)
c       hsv(i,7)=hsv(i,7)+d3(i)
c       hsv(i,8)=hsv(i,8)+d4(i)
c       hsv(i,9)=hsv(i,9)+d5(i)
c       hsv(i,10)=hsv(i,10)+d6(i)
c       +++++++++++++++++++++++++++++++++++++++++++++++++
c
c       Computes stresses 
c      
c       sig1(i) = C11*hsv(i,1)+C12*hsv(i,6)+C13*hsv(i,7)
c       sig2(i) = C12*hsv(i,1)+C22*hsv(i,6)+C23*hsv(i,7)
c       sig3(i) = C13*hsv(i,1)+C23*hsv(i,6)+C33*hsv(i,7)
c       sig4(i) = C44*hsv(i,8)
c       sig5(i) = C55*hsv(i,9)
c       sig6(i) = C66*hsv(i,10)
c
c       sig1(i) = sig1(i)+(C11*d1(i)+C12*d2(i)+C13*d3(i))
c       sig2(i) = sig2(i)+(C12*d1(i)+C22*d2(i)+C23*d3(i))
c       sig3(i) = sig3(i)+(C13*d1(i)+C23*d2(i)+C33*d3(i))
c       sig4(i) = sig4(i)+C44*d4(i)
c       sig5(i) = sig5(i)+C55*d5(i)
c       sig6(i) = sig6(i)+C66*d6(i)
c
c       S1 = C11*d1(i)+C12*d2(i)+C13*d3(i)
c       S2 = C12*d1(i)+C22*d2(i)+C23*d3(i)
c       S3 = C13*d1(i)+C23*d2(i)+C33*d3(i)
c       S4 = C44*d4(i)
c       S5 = C55*d5(i)
c       S6 = C66*d6(i)
c
c       hsv(i,1) = 2.0
c       hsv(i,2) = 0.9
c       hsv(i,3) = 1.2
c       hsv(i,4) = 0.5
c       hsv(i,5) = 3.0
c       hsv(i,6) = 5.0
c       hsv(i,7) = 0.7
c       hsv(i,8) = 3.4
c       hsv(i,9) = 8.0
c
c       Right Green-St. Venant Strain Tensor
        G1 = 0.5*(hsv(i,1)*hsv(i,1)+hsv(i,2)*hsv(i,2)+
     .       hsv(i,3)*hsv(i,3)-1.)
        G2 = 0.5*(hsv(i,4)*hsv(i,4)+hsv(i,5)*hsv(i,5)+
     .       hsv(i,6)*hsv(i,6)-1.)
        G3 = 0.5*(hsv(i,7)*hsv(i,7)+hsv(i,8)*hsv(i,8)+
     .       hsv(i,9)*hsv(i,9)-1.)
        G4 = 0.5*(hsv(i,1)*hsv(i,4)+hsv(i,2)*hsv(i,5)+
     .       hsv(i,3)*hsv(i,6))
        G5 = 0.5*(hsv(i,4)*hsv(i,7)+hsv(i,5)*hsv(i,8)+
     .       hsv(i,6)*hsv(i,9))
        G6 = 0.5*(hsv(i,1)*hsv(i,7)+hsv(i,2)*hsv(i,8)+
     .       hsv(i,3)*hsv(i,9))
c
c       Left Green-St. Venant Strain Tensor
c       G1 = 0.5*(hsv(i,1)*hsv(i,1)+hsv(i,4)*hsv(i,4)+
c     .      hsv(i,7)*hsv(i,7)-1.)
c       G2 = 0.5*(hsv(i,2)*hsv(i,2)+hsv(i,5)*hsv(i,5)+
c     .      hsv(i,8)*hsv(i,8)-1.)
c       G3 = 0.5*(hsv(i,3)*hsv(i,3)+hsv(i,6)*hsv(i,6)+
c     .      hsv(i,9)*hsv(i,9)-1.)
c       G4 = 0.5*(hsv(i,1)*hsv(i,2)+hsv(i,4)*hsv(i,5)+
c     .      hsv(i,7)*hsv(i,8))
c       G5 = 0.5*(hsv(i,2)*hsv(i,3)+hsv(i,5)*hsv(i,6)+
c     .      hsv(i,8)*hsv(i,9))
c       G6 = 0.5*(hsv(i,1)*hsv(i,3)+hsv(i,4)*hsv(i,6)+
c     .      hsv(i,7)*hsv(i,9))                     
c
c       write(59,*) 'Green-St. Venant Strain Tensor'
c       write(59,*) G1, '  ' ,G4,  '  ' ,G6      
c       write(59,*) G4, '  ' ,G2,  '  ' ,G5
c       write(59,*) G6, '  ' ,G5,  '  ' ,G3
c       write(59,*) '      '
c
c       Second Piola-Kirchhoff Stress Tensor
        S1 = C11*G1 + C12*G2 + C13*G3
        S2 = C12*G1 + C22*G2 + C23*G3
        S3 = C13*G1 + C23*G2 + C33*G3
        S4 = 2.*C44*G4
        S5 = 2.*C55*G5
        S6 = 2.*C66*G6
c
c       write(59,*) 'Second Piola-Kirchhoff Stress Tensor'
c       write(59,*) S1, '  ' ,S4,  '  ' ,S6      
c       write(59,*) S4, '  ' ,S2,  '  ' ,S5
c       write(59,*) S6, '  ' ,S5,  '  ' ,S3
c       write(59,*) '      '
c
c       Jacobian
        detF=hsv(i,1)*(hsv(i,5)*hsv(i,9)-hsv(i,6)*hsv(i,8))-
     .       hsv(i,2)*(hsv(i,4)*hsv(i,9)-hsv(i,6)*hsv(i,7))+
     .       hsv(i,3)*(hsv(i,4)*hsv(i,8)-hsv(i,5)*hsv(i,7))
c
c       write(59,*) 'Jacobian'
c       write(59,*) detF
c       write(59,*) '      '
c
c       Cauchy Stresses (F*S*F^T)
        FS11 = hsv(i,1)*S1 + hsv(i,4)*S4 + hsv(i,7)*S6
        FS12 = hsv(i,1)*S4 + hsv(i,4)*S2 + hsv(i,7)*S5
        FS13 = hsv(i,1)*S6 + hsv(i,4)*S5 + hsv(i,7)*S3
        FS21 = hsv(i,2)*S1 + hsv(i,5)*S4 + hsv(i,8)*S6
        FS22 = hsv(i,2)*S4 + hsv(i,5)*S2 + hsv(i,8)*S5
        FS23 = hsv(i,2)*S6 + hsv(i,5)*S5 + hsv(i,8)*S3
        FS31 = hsv(i,3)*S1 + hsv(i,6)*S4 + hsv(i,9)*S6
        FS32 = hsv(i,3)*S4 + hsv(i,6)*S2 + hsv(i,9)*S5
        FS33 = hsv(i,3)*S6 + hsv(i,6)*S5 + hsv(i,9)*S3
c
c       write(59,*) 'F*S'
c       write(59,*) FS11, '  ' ,FS12,  '  ' ,FS13      
c       write(59,*) FS21, '  ' ,FS22,  '  ' ,FS23
c       write(59,*) FS31, '  ' ,FS32,  '  ' ,FS33
c       write(59,*) '      '
c
        sig1(i) = 1./detF*(FS11*hsv(i,1)+FS12*hsv(i,4)+FS13*hsv(i,7))
        sig2(i) = 1./detF*(FS21*hsv(i,2)+FS22*hsv(i,5)+FS23*hsv(i,8))
        sig3(i) = 1./detF*(FS31*hsv(i,3)+FS32*hsv(i,6)+FS33*hsv(i,9))
        sig4(i) = 1./detF*(FS11*hsv(i,2)+FS12*hsv(i,5)+FS13*hsv(i,8))
        sig5(i) = 1./detF*(FS21*hsv(i,3)+FS22*hsv(i,6)+FS23*hsv(i,9))
        sig6(i) = 1./detF*(FS11*hsv(i,3)+FS12*hsv(i,6)+FS13*hsv(i,9))
c
c       write(59,*) 'Cauchy Stresses '
c       write(59,*) sig1(i), '  ' ,sig4(i),  '  ' ,sig6(i)    
c       write(59,*) sig4(i), '  ' ,sig2(i),  '  ' ,sig5(i)
c       write(59,*) sig6(i), '  ' ,sig5(i),  '  ' ,sig3(i)
c       write(59,*) '      '
c
c       write(59,*) '-----------------------------'
c
      endif
   10 continue
c       
      return
      end
      subroutine umat42v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,cma,qmat,elsizv,idelev)
c
c******************************************************************
c|  Livermore Software Technology Corporation (LSTC)              |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2007 LSTC                                      |
c|  All Rights Reserved                                           |
c******************************************************************
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*)
      dimension temps(*),crv(lq1,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      dimension idelev(*)
      logical failels(*)
      character*5 etype
c     
c     eps(1)=local x  strain increment
c     eps(2)=local y  strain increment
c     eps(3)=local z  strain increment
c     eps(4)=local xy strain increment
c     eps(5)=local yz strain increment
c     eps(6)=local zx strain increment
c
c     sig(1)=local x  stress
c     sig(2)=local y  stress
c     sig(3)=local z  stress
c     sig(4)=local xy stress
c     sig(5)=local yz stress
c     sig(6)=local zx stress
c
c     material properties for orthotropic elastic material 
c     If transverse elastic properties are not provided, the transverse
c     directions take on the same properties as the longitudinal direction,
c     i.e. behaves as an isotropic material
c     
      E1 = cm(1)
      E2 = cm(12)
      E3 = cm(13)
      v21 = cm(14)
      v31=cm(15)
      v32=cm(16)
      G12 = cm(17)
      G23 = cm(18)
      G31 = cm(19)            
c      
      v12 = v21*E1/E2
      v23 = v32*E2/E3
      v13 = v31*E1/E3
c
c     stiffness matrix 
c
      S = 1.-v12*v21-v23*v32-v31*v13-2.*v21*v32*v13
c
      C11=(1. - v23*v32)*E1/S
      C12=(v21 + v31*v23)*E1/S
      C13=(v31 + v21*v32)*E1/S
c     C13=(v31 + v21*v23)/(E2*E3*S)
      C22=(1. - v31*v13)*E2/S
      C23=(v32 + v31*v12)*E2/S
      C33=(1. - v12*v21)*E3/S
      C44=G12
      C55=G23
      C66=G31
c
c     write(59,*) 'Stiffness Matrix'
c     write(59,*) C11, '  ',C12,  '  ',  C13      
c     write(59,*) C22, '  ',C23,  '  ',  C33
c     write(59,*) C44, '  ',C55,  '  ',  C66
c     write(59,*) '      '
c
      do  10 i=lft,llt
      if (etype.eq.'solid'.or.etype.eq.'shl_t'.or.
     1     etype.eq.'sld2d'.or.etype.eq.'tshel') then
c
c       Computes stresses 
c
        sig1(i) = sig1(i)+(C11*d1(i)+C12*d2(i)+C13*d3(i))
        sig2(i) = sig2(i)+(C12*d1(i)+C22*d2(i)+C23*d3(i))
        sig3(i) = sig3(i)+(C13*d1(i)+C23*d2(i)+C33*d3(i))
        sig4(i) = sig4(i)+C44*d4(i)
        sig5(i) = sig5(i)+C55*d5(i)
        sig6(i) = sig6(i)+C66*d6(i) 
c
c       sig1(i) = sig1(i)+(C11*d1(i)+C12*d2(i)+C13*d3(i))*dt1siz(i)
c       sig2(i) = sig2(i)+(C12*d1(i)+C22*d2(i)+C23*d3(i))*dt1siz(i)
c       sig3(i) = sig3(i)+(C13*d1(i)+C23*d2(i)+C33*d3(i))*dt1siz(i)
c       sig4(i) = sig4(i)+C44*d4(i)*dt1siz(i)
c       sig5(i) = sig5(i)+C55*d5(i)*dt1siz(i)
c       sig6(i) = sig6(i)+C66*d6(i)*dt1siz(i)
c
      else if (etype.eq.'shell') then
c       hsv(i,7)=-(hsv(i,4)*C13*hsv(i,1) + C23*hsv(i,6))/C33
        d3(i) = -(C13*d1(i) + C23*d2(i))/C33
        sig1(i) = sig1(i)+(C11*d1(i)+C12*d2(i)+C13*d3(i))
        sig2(i) = sig2(i)+(C12*d1(i)+C22*d2(i)+C23*d3(i))               
        sig3(i) = 0.0
        sig4(i) = sig4(i)+C44*d4(i)
        sig5(i) = sig5(i)+capa*C55*d5(i)
        sig6(i) = sig6(i)+capa*C66*d6(i)
         endif
   10 continue
c
      return
      end
      subroutine urmats(lft,llt,cm,capa,mt,crv,ipt,rcoor,scoor,tcoor,
     . nnm1,nconstp,nip,ipt_thk,npc,plc,eltype,nshbwp)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
      integer
     1 nh01,nh02,nh03,nh04,nh05,nh06,nh07,nh08,nh09,nh10,
     2 nb01,nb02,nb03,nb04,nb05,nb06,nb07,nb08,nb09,nb10,
     3 ns01,ns02,ns03,ns04,ns05,ns06,ns07,ns08,ns09,ns10,
     4 nt01,nt02,nt03,nt04,nt05,nt06,nt07,nt08,nt09,nt10,
     5 nu01,nu02,nu03,nu04,nu05,nu06,nu07,nu08,nu09,nu10,
     6 nd01,nd02,nd03,nd04,nd05,nd06,nd07,nd08,nd09,nd10,
     7 ntbqs
      common/bk05/
     1 nh01,nh02,nh03,nh04,nh05,nh06,nh07,nh08,nh09,nh10,
     2 nb01,nb02,nb03,nb04,nb05,nb06,nb07,nb08,nb09,nb10,
     3 ns01,ns02,ns03,ns04,ns05,ns06,ns07,ns08,ns09,ns10,
     4 nt01,nt02,nt03,nt04,nt05,nt06,nt07,nt08,nt09,nt10,
     5 nu01,nu02,nu03,nu04,nu05,nu06,nu07,nu08,nu09,nu10,
     6 nd01,nd02,nd03,nd04,nd05,nd06,nd07,nd08,nd09,nd10,
     7 ntbqs
c
      integer*8 dm_hmtnum,dm_hnfegp
      integer*8 dm_bmtnum,dm_bnfegp
      integer*8 dm_smtnum,dm_snfegp
      integer*8 dm_tmtnum,dm_tnfegp
c
      common/dmbk05/
     1 dm_hmtnum,dm_hnfegp,
     2 dm_bmtnum,dm_bnfegp,
     3 dm_smtnum,dm_snfegp,
     4 dm_tmtnum,dm_tnfegp
      include 'bk06.inc'
      integer n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      common/bk07/n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      integer lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      common/bk13/lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      include 'bk19.inc'
      include 'iounits.inc'
      include 'memaia.inc'
      include 'umatss.inc'
      include 'shlioc.inc'
c
c     dynamic memory allocation stuff.  This has to be a C style 
c     include because of the use of integer*8
c
      integer i_mem(0:1)
      real r_mem(0:1)
      integer*8 i8_mem(0:1)
      real*8 r8_mem(0:1)
      real*8 real8_mem
      integer*8 m_mem(0:1)
      common/dynmem/real8_mem(0:1)
      equivalence (i_mem,r_mem,i8_mem,r8_mem,real8_mem,m_mem)
      integer*8 mem_alloc, mems_alloc,mem_realloc, memh_ptr
      integer mem_length,memh_alloc,memh_realloc,memh_length
      integer memh_enum,memh_enums,memh_attach
c      integer mem_newmark,mem_getmark,memh_getmark
      integer*8 memh_detach,memh_ralloc
      integer*8 m_to_m8
      integer*8 m_to_mm
      external mem_alloc, mems_alloc, mem_realloc, mem_length
      external memh_alloc, memh_realloc, memh_length, memh_ptr
      external memh_enum,memh_attach,memh_detach,memh_ralloc
      external memh_enums
      external m_to_m8
      external m_to_mm
      integer*8 dm_x,dm_xr,dm_v,dm_vr,dm_a,dm_ar,dm_rots,dm_rotr
      integer*8 dm_xms,dm_xmr,dm_xm2,dm_xmz,dm_xma
      integer*8 dm_xtpz,dm_disp
      integer*8 dm_me1,dm_me2,dm_me3
      integer*8 dm_masp,dm_masr
      integer*8 dm_xrb,dm_yrb,dm_zrb
      integer*8 dm_axrb,dm_ayrb,dm_azrb
      integer*8 dm_rbfx,dm_rbfy,dm_rbfz
      integer*8 dm_n45,dm_n46,dm_n47
      integer len_n45,len_n46,len_n47
      common/dynmem1/ 
     . dm_x,dm_xr,dm_v,dm_vr,dm_a,dm_ar,dm_rots,dm_rotr,
     . dm_xms,dm_xmr,dm_xm2,dm_xmz,dm_xma,
     . dm_xtpz,dm_disp,
     . dm_me1,dm_me2,dm_me3,
     . dm_masp,dm_masr,
     . dm_xrb,dm_yrb,dm_zrb,
     . dm_axrb,dm_ayrb,dm_azrb,
     . dm_rbfx,dm_rbfy,dm_rbfz,
     . dm_n45,dm_n46,dm_n47
      common/dynmem1l/len_n45,len_n46,len_n47
c
      integer*8 mem_allocm,mems_allocm,mem_reallocm
      integer memh_allocm,memh_reallocm
      external mem_allocm,mems_allocm,mem_reallocm
      external memh_allocm,memh_reallocm
      integer*8 mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      integer*8 mem_reallocmchk,mem_reallocchk
      integer memh_allocmchk,memh_allocchk
      external mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      external mem_reallocmchk,mem_reallocchk
      external memh_allocmchk,memh_allocchk
      common/frozloc/yhat(nlq,3,4),qloc(nlq,3,3),qmat(nlq,3,3)
      common/bux13loc/ eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      common/aux2loc/e1(nlq),e2(nlq),d3(nlq),e4(nlq),e5(nlq),e6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common /aux8cloc/
     & d1(nlq),d2(nlq),d4(nlq),d5(nlq),d6(nlq),
     & dfg11(nlq),dfg21(nlq),dfg31(nlq),dfg12(nlq),dfg22(nlq),
     & dfg32(nlq),dfg13(nlq),dfg23(nlq),dfg33(nlq),stg5(nlq),
     & stg6(nlq),a11(nlq),a12(nlq),a21(nlq),a22(nlq)
      common/aux14loc/
     & sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),
     & sig5(nlq),sig6(nlq),epsps(nlq),hsvs(nlq,NHISVAR)
      common/aux19loc/
     1 sign0(nlq),sign1(nlq),sign2(nlq),sign3(nlq),sign4(nlq),
     2 sign5(nlq),sign6(nlq)
      common/aux33loc/
     1 ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ixs(nlq,4),mxt(nlq)
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,numelh10,numels8
      common/bk36loc/index
      common/failcmloc/ifail(nlq)
      common/failuloc/sieu(nlq),fail(nlq),ifaili(nlq),nfipts(nlq)
      common/hourgloc/ym(nlq),gm(nlq),ifsv(nlq)
      common/shloptloc/ibelyt
c
      common/bk01/itherm,itemp,ntmp0,ntmp1,itempan,itempdr,itmpe(44)
      common/bk26/begtim,nintcy
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie,selie,selke,
     . erodehg,eintcnt,engjnt
      logical failur,failgl
      common/failcm/failur,failgl
      common/numcpu/ncpu,ncpua,ncpub,lenvec(8)
      integer istrn,istupd,ibelyts,miter,ipstpd,intsts,nodsts,intstn,
     1 nodstn,jstrn,nddet,nsout01,nsout02,nsout03,nsout04,
     2 nsout05,nsout06
      real wrpang
      common/shlopt/istrn,istupd,ibelyts,miter,wrpang,ipstpd,intsts,
     1 nodsts,intstn,nodstn,jstrn,nddet,nsout01,nsout02,nsout03,
     2 nsout04,nsout05,nsout06
      common/soundloc/ss(nlq),sndsp(nlq),diagm(nlq),sarea(nlq),
     . dxl(nlq)
      common/thropt/itopaz(101)
      common/subtssloc/dt1siz(nlq)
      common/sidesloc/sidmn(nlq)
c
      dimension cm(*),crv(lq1,2,*),nconstp(*),npc(*),plc(*)
c     dimension eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      dimension thhsv(nlq,100),thhsvi(100)
      dimension elsizv(nlq),idelev(nlq)
      logical failel,failels(nlq)
      dimension qmats(3,3)
      character*(*)eltype
      dimension nshbwp(*)
c
c     data temps/nlq*0.0/
c
c     processing elements
c
c$omp threadprivate (/aux2loc/)
c$omp threadprivate (/aux8cloc/)
c$omp threadprivate (/aux14loc/)
c$omp threadprivate (/aux19loc/)
c$omp threadprivate (/aux33loc/)
c$omp threadprivate (/bk36loc/)
c$omp threadprivate (/failuloc/)
c$omp threadprivate (/failcmloc/)
c$omp threadprivate (/hourgloc/)
c$omp threadprivate (/shloptloc/)
c$omp threadprivate (/soundloc/)
c$omp threadprivate (/subtssloc/)
c$omp threadprivate (/bux13loc/)
c$omp threadprivate (/frozloc/)
c$omp threadprivate (/sidesloc/)
      mx=48*(mxt(lft)-1)
      gm(lft)=cm(mx+abs(ishrmp(mt)))
      bk=cm(mx+ibulkp(mt))
      pr=(3.*bk-2.*gm(lft))/(2.*(3.*bk+gm(lft)))
      ym(lft)=2.*(1.+pr)*gm(lft)
      ss(lft)=ym(lft)/(1.-pr*pr)
c
c     Determining number of requested history variables
c
      no_hsvs=nconstp(mxt(lft))
      if (ivumat(mt).ne.0) no_hsvs=no_hsvs-9-6
      if (lenvec(5) .eq.1) no_hsvs=no_hsvs-6
      if (iorien(mt).ne.0) no_hsvs=no_hsvs-2
      if (ihyper(mt).ne.0) no_hsvs=no_hsvs-9
      if (iabs(ihyper(mt)).eq.4) no_hsvs=no_hsvs-18
      if (ifipts(mt).lt.0) no_hsvs=no_hsvs-1
      if (ioshl(67) .eq.1) no_hsvs=no_hsvs-1
c
c     initializing indeces
c
      if(ivumat(mt).ne.0) then
        if (iorien(mt).ne.0) then
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_ortho=ind_defgrad+9
          else
            ind_ortho=1+no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_ortho
            ind_ortho=ind_ortho+1
          endif
          ind_q1=ind_ortho
          ind_q2=ind_ortho+1
          ind_last=ind_q2
        else
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_last=ind_defgrad+8
          else
            ind_last=no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_last+1
            ind_last=ind_last+1
          endif             
        endif
        ind_vumat=ind_last+1
        ind_last=ind_vumat+9+6-1
c        ind_last=ind_last+1
      else
        if (iorien(mt).ne.0) then
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_ortho=ind_defgrad+9
            if (iabs(ihyper(mt)).eq.4) ind_ortho=ind_ortho+18
          else
            ind_ortho=1+no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_ortho
            ind_ortho=ind_ortho+1
          endif
          ind_q1=ind_ortho
          ind_q2=ind_ortho+1
          ind_last=ind_q2
        else
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_last=ind_defgrad+8
            if (iabs(ihyper(mt)).eq.4) ind_last=ind_last+18
          else
            ind_last=no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_last+1
            ind_last=ind_last+1
          endif        
        endif 
      endif
c
      if (ind_last.gt.NHISVAR) then
c       write(iotty,900) lqfmiv(mxt(lft))
c       write(iohsp,900) lqfmiv(mxt(lft))
c       write(iomsg,900) lqfmiv(mxt(lft))
c       call adios(2)
        ierdat(1)=lqfmiv(mxt(lft))
        call lsmsg(3,MSG_SOL+1148,ioall,ierdat,rerdat,cerdat,0)
      endif
c$$$     if (tt.eq.0.0.and.nintcy.eq.0) then
c$$$c
c$$$c      move transformation matrix in history variables array
c$$$c
c$$$       if (iorien(mt).ne.0) then
c$$$         do i=lft,llt
c$$$           hsvs(i,ind_q2)=hsvs(i,2)
c$$$           hsvs(i,ind_q1)=hsvs(i,1)
c$$$         enddo
c$$$         if (no_hsvs.gt.0) then
c$$$           do j=1,no_hsvs
c$$$             do i=lft,llt
c$$$               hsvs(i,j)=0.0
c$$$             enddo
c$$$           enddo
c$$$         endif
c$$$       endif
c$$$     endif
c
      if (ipt.eq.1.and.ifipts(mt).ne.0) then
        do i=lft,llt
          nfipts(i)=0
        enddo
      endif
c
      num_nods=4
      if (ibelyt.eq.0) then
        num_nods=8
        do 10 i=lft,llt
        sign0(i)=sign3(i)+ym(lft)*d3(i)
        sign3(i)=sign0(i)
 10     continue
      endif
c
c     transform stresses and strains to material system
c
      if (iorien(mt).ne.0.and.ihyper(mt).ge.0) then
        do 20 i=lft,llt
        stg5(i)=sig5(i)
        stg6(i)=sig6(i)
        if (ncycle.le.1) then
          rr = 1./sqrt(hsvs(i,ind_q1)**2+hsvs(i,ind_q2)**2)
          hsvs(i,ind_q1)=hsvs(i,ind_q1)*rr
          hsvs(i,ind_q2)=hsvs(i,ind_q2)*rr
        end if
        a11(i) =hsvs(i,ind_q1)*sig1(i)-hsvs(i,ind_q2)*sig4(i)
        a12(i) =hsvs(i,ind_q2)*sig1(i)+hsvs(i,ind_q1)*sig4(i)
        a21(i) =hsvs(i,ind_q1)*sig4(i)-hsvs(i,ind_q2)*sig2(i)
        a22(i) =hsvs(i,ind_q2)*sig4(i)+hsvs(i,ind_q1)*sig2(i)
        sig1(i)=hsvs(i,ind_q1)*a11(i)-hsvs(i,ind_q2)*a21(i)
        sig2(i)=hsvs(i,ind_q2)*a12(i)+hsvs(i,ind_q1)*a22(i)
        sig4(i)=hsvs(i,ind_q1)*a12(i)-hsvs(i,ind_q2)*a22(i)
        sig5(i)=hsvs(i,ind_q2)*stg6(i)+hsvs(i,ind_q1)*stg5(i)
        sig6(i)=hsvs(i,ind_q1)*stg6(i)-hsvs(i,ind_q2)*stg5(i)
 20     continue
c
        do 30 i=lft,llt
        d4(i)  =.5*e4(i)
        a11(i) =hsvs(i,ind_q1)*e1(i)-hsvs(i,ind_q2)*d4(i)
        a12(i) =hsvs(i,ind_q2)*e1(i)+hsvs(i,ind_q1)*d4(i)
        a21(i) =hsvs(i,ind_q1)*d4(i)-hsvs(i,ind_q2)*e2(i)
        a22(i) =hsvs(i,ind_q2)*d4(i)+hsvs(i,ind_q1)*e2(i)
        d1(i)  =hsvs(i,ind_q1)*a11(i)-hsvs(i,ind_q2)*a21(i)
        d2(i)  =hsvs(i,ind_q2)*a12(i)+hsvs(i,ind_q1)*a22(i)
        d4(i)  =2.*(hsvs(i,ind_q1)*a12(i)-hsvs(i,ind_q2)*a22(i))
        d5(i)  =hsvs(i,ind_q2)*e6(i)+hsvs(i,ind_q1)*e5(i)
        d6(i)  =hsvs(i,ind_q1)*e6(i)-hsvs(i,ind_q2)*e5(i)
 30     continue
c
c     element local strains
c
      else
        do i=lft,llt
          d1(i)  =e1(i)
          d2(i)  =e2(i)
          d4(i)  =e4(i)
          d5(i)  =e5(i)
          d6(i)  =e6(i)
        enddo
      endif
c
c     storage of deformation gradient
c
      if (ivumat(mt).ne.0) then
        ifo11=ind_vumat
        ifo21=ind_vumat+1
        ifo31=ind_vumat+2
        ifo12=ind_vumat+3
        ifo22=ind_vumat+4
        ifo32=ind_vumat+5
        ifo13=ind_vumat+6
        ifo23=ind_vumat+7
        ifo33=ind_vumat+8
      endif
      if (ihyper(mt).ne.0) then
        if11=ind_defgrad
        if21=ind_defgrad+1
        if31=ind_defgrad+2
        if12=ind_defgrad+3
        if22=ind_defgrad+4
        if32=ind_defgrad+5
        if13=ind_defgrad+6
        if23=ind_defgrad+7
        if33=ind_defgrad+8
c
c       deformation gradient in element system at time = 0
c
        if (tt.eq.begtim.and.nintcy.eq.0) then
          do i=lft,llt
            hsvs(i,if11)=1.0
            hsvs(i,if21)=0.0
            hsvs(i,if31)=0.0
            hsvs(i,if12)=0.0
            hsvs(i,if22)=1.0
            hsvs(i,if32)=0.0
            hsvs(i,if13)=0.0
            hsvs(i,if23)=0.0
            hsvs(i,if33)=1.0
          enddo
          if (ivumat(mt).ne.0) then
            do i=lft,llt
              hsvs(i,ifo11)=1.0
              hsvs(i,ifo21)=0.0
              hsvs(i,ifo31)=0.0
              hsvs(i,ifo12)=0.0
              hsvs(i,ifo22)=1.0
              hsvs(i,ifo32)=0.0
              hsvs(i,ifo13)=0.0
              hsvs(i,ifo23)=0.0
              hsvs(i,ifo33)=1.0
            enddo
          endif
        endif
c
c       store deformation gradient in element system
c
        if (iorien(mt).eq.0) then
          do i=lft,llt
            dfg11(i)=hsvs(i,if11)
            dfg21(i)=hsvs(i,if21)
            dfg31(i)=hsvs(i,if31)
            dfg12(i)=hsvs(i,if12)
            dfg22(i)=hsvs(i,if22)
            dfg32(i)=hsvs(i,if32)
            dfg13(i)=hsvs(i,if13)
            dfg23(i)=hsvs(i,if23)
            dfg33(i)=hsvs(i,if33)
           enddo
        endif
      endif
c
      if (iorien(mt).ne.0) then
c
c       store deformation gradient in material system
c       as F_bar=Q^T*F
c
      if (ihyper(mt).gt.0) then
c        if (ihyper(mt).eq.1) then
          do i=lft,llt
            dfg11(i)=hsvs(i,ind_q1)*hsvs(i,if11)-
     1               hsvs(i,ind_q2)*hsvs(i,if21)
            dfg21(i)=hsvs(i,ind_q2)*hsvs(i,if11)+
     1               hsvs(i,ind_q1)*hsvs(i,if21)
            dfg31(i)=hsvs(i,if31)
            dfg12(i)=hsvs(i,ind_q1)*hsvs(i,if12)-
     1               hsvs(i,ind_q2)*hsvs(i,if22)
            dfg22(i)=hsvs(i,ind_q2)*hsvs(i,if12)+
     1               hsvs(i,ind_q1)*hsvs(i,if22)
            dfg32(i)=hsvs(i,if32)
            dfg13(i)=hsvs(i,ind_q1)*hsvs(i,if13)-
     1               hsvs(i,ind_q2)*hsvs(i,if23)
            dfg23(i)=hsvs(i,ind_q2)*hsvs(i,if13)+
     1               hsvs(i,ind_q1)*hsvs(i,if23)
            dfg33(i)=hsvs(i,if33)
c
            hsvs(i,if11)=dfg11(i)
            hsvs(i,if21)=dfg21(i)
            hsvs(i,if12)=dfg12(i)
            hsvs(i,if22)=dfg22(i)
            hsvs(i,if13)=dfg13(i)
            hsvs(i,if23)=dfg23(i)
          enddo
c
c       store deformation gradient in material system
c       as F_bar=Q^T*F*Q
c
        else if (ihyper(mt).eq.-2) then
          do i=lft,llt
            dfg11(i)=hsvs(i,ind_q1)*hsvs(i,if11)-
     1               hsvs(i,ind_q2)*hsvs(i,if21)
            dfg21(i)=hsvs(i,ind_q2)*hsvs(i,if11)+
     1               hsvs(i,ind_q1)*hsvs(i,if21)
            dfg31(i)=hsvs(i,if31)
            dfg12(i)=hsvs(i,ind_q1)*hsvs(i,if12)-
     1               hsvs(i,ind_q2)*hsvs(i,if22)
            dfg22(i)=hsvs(i,ind_q2)*hsvs(i,if12)+
     1               hsvs(i,ind_q1)*hsvs(i,if22)
            dfg32(i)=hsvs(i,if32)
            dfg13(i)=hsvs(i,ind_q1)*hsvs(i,if13)-
     1               hsvs(i,ind_q2)*hsvs(i,if23)
            dfg23(i)=hsvs(i,ind_q2)*hsvs(i,if13)+
     1               hsvs(i,ind_q1)*hsvs(i,if23)
            dfg33(i)=hsvs(i,if33)
c
            hsvs(i,if11)=dfg11(i)*hsvs(i,ind_q1)-dfg12(i)*hsvs(i,ind_q2)
            hsvs(i,if12)=dfg11(i)*hsvs(i,ind_q2)+dfg12(i)*hsvs(i,ind_q1)
            hsvs(i,if13)=dfg13(i)
            hsvs(i,if21)=dfg21(i)*hsvs(i,ind_q1)-dfg22(i)*hsvs(i,ind_q2)
            hsvs(i,if22)=dfg21(i)*hsvs(i,ind_q2)+dfg22(i)*hsvs(i,ind_q1)
            hsvs(i,if23)=dfg23(i)
            hsvs(i,if31)=dfg31(i)*hsvs(i,ind_q1)-dfg32(i)*hsvs(i,ind_q2)
            hsvs(i,if32)=dfg31(i)*hsvs(i,ind_q2)+dfg32(i)*hsvs(i,ind_q1)
            hsvs(i,if33)=dfg33(i)
            dfg11(i)=hsvs(i,if11)
            dfg22(i)=hsvs(i,if22)
            dfg33(i)=hsvs(i,if33)
            dfg12(i)=hsvs(i,if12)
            dfg21(i)=hsvs(i,if21)
            dfg13(i)=hsvs(i,if13)
            dfg31(i)=hsvs(i,if31)
            dfg23(i)=hsvs(i,if23)
            dfg32(i)=hsvs(i,if32)
          enddo
c
c       store deformation gradient in material system
c       as F_bar=F*Q
c
c        else if (ihyper(mt).lt.0) then
        else if (ihyper(mt).eq.-1) then
          do i=lft,llt
            dfg11(i)=hsvs(i,ind_q1)*hsvs(i,if11)-
     1               hsvs(i,ind_q2)*hsvs(i,if12)
            dfg12(i)=hsvs(i,ind_q2)*hsvs(i,if11)+
     1               hsvs(i,ind_q1)*hsvs(i,if12)
            dfg13(i)=hsvs(i,if13)
            dfg21(i)=hsvs(i,ind_q1)*hsvs(i,if21)-
     1               hsvs(i,ind_q2)*hsvs(i,if22)
            dfg22(i)=hsvs(i,ind_q2)*hsvs(i,if21)+
     1               hsvs(i,ind_q1)*hsvs(i,if22)
            dfg23(i)=hsvs(i,if23)
            dfg31(i)=hsvs(i,ind_q1)*hsvs(i,if31)-
     1               hsvs(i,ind_q2)*hsvs(i,if32)
            dfg32(i)=hsvs(i,ind_q2)*hsvs(i,if31)+
     1               hsvs(i,ind_q1)*hsvs(i,if32)
            dfg33(i)=hsvs(i,if33)
c
            hsvs(i,if11)=dfg11(i)
            hsvs(i,if12)=dfg12(i)
            hsvs(i,if21)=dfg21(i)
            hsvs(i,if22)=dfg22(i)
            hsvs(i,if31)=dfg31(i)
            hsvs(i,if32)=dfg32(i)
          enddo
        endif
      endif
c
c     update first two rows of deformation gradient at time = n+1
c
      if (ihyper(mt).ne.0) then
        call integf_s(hsvs(1,if11),hsvs(1,if21),hsvs(1,if31),
     &   hsvs(1,if12),hsvs(1,if22),hsvs(1,if32),
     &   hsvs(1,if13),hsvs(1,if23),hsvs(1,if33),lft,llt)
      endif
c
      if (iabs(ihyper(mt)).eq.3.or.iabs(ihyper(mt)).eq.4) then
         idxshl=5
         if (numels8.gt.0) idxshl=9
         call usrshl_defgrad(hsvs(1,if11),ihyper(mt),
     1        rcoor,scoor,tcoor,
     2        r_mem(dm_x),r_mem(dm_rots),
     3        a(ns03+3*(idxshl-1)*numels+12*nnm1),
     4        a(ns03+3*(idxshl+3)*numels+ 4*nnm1),
     5        a(ns05+9*nnm1),
     6        lft,llt)         
         do k=1,3
            do j=1,3               
               do i=lft,llt
                  qmat(i,j,k)=qloc(i,j,k)
               enddo
            enddo
         enddo
         if (iorien(mt).ne.0) then
            do j=1,3
               do i=lft,llt
                  q1=hsvs(i,ind_q1)*qmat(i,j,1)-
     1                 hsvs(i,ind_q2)*qmat(i,j,2)
                  q2=hsvs(i,ind_q2)*qmat(i,j,1)+
     1                 hsvs(i,ind_q1)*qmat(i,j,2)
                  qmat(i,j,1)=q1
                  qmat(i,j,2)=q2
               enddo
            enddo
         endif
      endif
c
c     if ishrmp(mt).lt.0 then the flag for temperature is active
c     temperatures are stored in the temps array.
c
      if (ishrmp(mt).lt.0) then
        call usr_temps (a(ntmp0+1),a(n19),itemp,temps,lft,llt,ix1,ix2,
     .   ix3,ix4,ixs,ixs(1,2),ixs(1,3),ixs(1,4),num_nods,itopaz,
     .   nnm1,rcoor,scoor,tcoor)
      else
        do i=lft,llt
          temps(i)=0.0
        enddo
      endif
c
c     Get thermal history variables from thermal
c     user material
c
      call get_thhsv(nthhsv,thhsv,itemp,lft,llt,
     1 num_nods,itopaz,nnm1,rcoor,scoor,tcoor)
c
c     energy calculations
c
      do 33 i=lft,llt
      einc(i)=d1(i)*sig1(i)+d2(i)*sig2(i)+d4(i)*sig4(i)+d5(i)*sig5(i)
     1       +d6(i)*sig6(i)
 33   continue
c
c     characteristic element size and element id
      do i=lft,llt
        elsizv(i)=sqrt(sidmn(i))
        ielm=nshbwp(nnm1+i)
        if (ielm.gt.0) idelev(i)=lqfinv(ielm,4)
      enddo
c
c     memory pointer for extra memory data
      lcma=ncma+nmmat+ia(ncma+mxt(lft)-1)-1
c
c     vectorization flag
c
      ivect=ivectr(mt)
c
      mte40=mt-40
c
c     scalar umat
c
      if (ivect.eq.0) then
        do 120 i=lft,llt
c
c       gathering of variables
c
        index=i
        dt1   =dt1siz(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
        epsp=epsps(i)
        temper=temps(i)
        failel=.false.
        if (ifipts(mt).lt.0) failel=hsvs(i,ind_fipts).ne.0.
        if (no_hsvs.gt.0) then
          do j=1,no_hsvs
            hsv(j)=hsvs(i,j)
          enddo
        endif
        if (ihyper(mt).ne.0) then
          do j=1,9
            hsv(no_hsvs+j)=hsvs(i,ind_defgrad-1+j)
          enddo
          if (iabs(ihyper(mt)).eq.4) then
          do j=1,18
            hsv(no_hsvs+j+9)=hsvs(i,ind_defgrad-1+j+9)
          enddo
          endif
          if (iabs(ihyper(mt)).eq.3.or.iabs(ihyper(mt)).eq.4) then
             do k=1,3
                do j=1,3
                   qmats(j,k)=qmat(i,j,k)
                enddo
             enddo
          endif
        endif
        if (nthhsv.gt.0) then
          do j=1,nthhsv
            thhsvi(j)=thhsv(i,j)
          enddo
        endif        
        elsiz=elsizv(i)
        idele=idelev(i)
c
c       call user developed subroutines here
c
        go to (41,42,43,44,45,46,47,48,49,50), mte40
   41   call umat41 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,a(lcma),qmats,elsiz,idele)
        go to 60
   42   call umat42 (cm(mx+1),eps,sig,hsv,dt1,capa,'shell',tt,crv,
     .   a(lcma),temper,qmats,elsiz,idele)
        go to 60
   43   call umat43 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,a(lcma),qmats,elsiz,idele)
        go to 60
   44   call umat44 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,a(lcma),qmats,elsiz,idele)
        go to 60
   45   call umat45 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,a(lcma),qmats,elsiz,idele)
        go to 60
   46   call umat46 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,a(lcma),thhsvi,nthhsv,qmats,elsiz,idele)
        go to 60
   47   call umat47 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,a(lcma),qmats,elsiz,idele)
        go to 60
   48   call umat48 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,a(lcma),qmats,elsiz,idele)
        go to 60
   49   call umat49 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,a(lcma),qmats,elsiz,idele)
        go to 60
   50   call umat50 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,a(lcma),qmats,elsiz,idele)
c
c       scattering of variables
c
 60     if (no_hsvs.gt.0) then
          do j=1,no_hsvs
            hsvs(i,j)=hsv(j)
          enddo
        endif
c
        if (ihyper(mt).ne.0) then
          do j=1,9
            hsvs(i,ind_defgrad-1+j)=hsv(no_hsvs+j)
          enddo
          if (iabs(ihyper(mt)).eq.4) then
          do j=1,18
            hsvs(i,ind_defgrad-1+j+9)=hsv(no_hsvs+j+9)
          enddo
          endif
        endif
c
        d3(i)  =eps(3)
c
        if (ifipts(mt).eq.0) failel=.false.
c
        if (failel) then
          sig1(i)=0.
          sig2(i)=0.
          sig3(i)=0.
          sig4(i)=0.
          sig5(i)=0.
          sig6(i)=0.
          if (ifipts(mt).lt.0) hsvs(i,ind_fipts)=1.
          nfipts(i)=nfipts(i)+1
          if (ifipts(mt).ge.0) then
          failur=.true.
          ifail(i)=1
          fail(i)=0.
          elseif (nfipts(i).ge.nint(cm(mx-ifipts(mt)))) then
          failur=.true.
          ifail(i)=1
          fail(i)=0.            
          endif          
        else
          sig1(i)=sig(1)
          sig2(i)=sig(2)
          sig3(i)=sig(3)
          sig4(i)=sig(4)
          sig5(i)=sig(5)
          sig6(i)=sig(6)
        endif
        epsps(i)=epsp
 120    continue
c
c     vector umat
c
      else
c
c       assume no failed elements
c
        if (ifipts(mt).ge.0) then
        do i=lft,llt
          failels(i)=.false.
        enddo
        else
          do i=lft,llt
            failels(i)=hsvs(i,ind_fipts).ne.0.
          enddo
        endif
c
c       call user developed subroutines here
c       first history variable is the one requested
c
        if(ivumat(mt).eq.0) then
c
        if (mt.eq.41) then
          call umat41v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,a(lcma),qmat,
     .     elsizv,idelev)
        elseif (mt.eq.42) then
          call umat42v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,a(lcma),qmat,
     .     elsizv,idelev)
        elseif (mt.eq.43) then
          call umat43v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,a(lcma),qmat,
     .     elsizv,idelev)
        elseif (mt.eq.44) then
          call umat44v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,a(lcma),qmat,
     .     elsizv,idelev)
        elseif (mt.eq.45) then
          call umat45v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,a(lcma),qmat,
     .     elsizv,idelev)
        elseif (mt.eq.46) then
          call umat46v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,a(lcma),
     .     thhsv,nthhsv,qmat,elsizv,idelev)
        elseif (mt.eq.47) then
          call umat47v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,a(lcma),qmat,
     .     elsizv,idelev)
        elseif (mt.eq.48) then
          call umat48v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,ipt,a(n8),a(n9),a(lcma),qmat,elsizv,
     .     idelev)
        elseif (mt.eq.49) then
          call umat49v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,ipt,temps,a(lcma),qmat,elsizv,idelev)
        elseif (mt.eq.50) then
          call umat50v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,a(lcma),qmat,
     .     elsizv,idelev)
c
        elseif (mt.eq.281) then
          call mat281s(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
        elseif (mt.eq.282) then
          call mat282s(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
c       elseif (mt.eq.283) then
c         call umat283s(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
c    .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
c    .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,
c    .     a(lcma))
        elseif (mt.eq.284) then
          call mat284s(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
        elseif (mt.eq.285) then
          call mat285s(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
        elseif (mt.eq.286) then
          call mat286s(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
        elseif (mt.eq.287) then
          call mat287s(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,
     .    ipt, npc, plc, nnm1,a(lcma))
        endif
c
        endif
c
        if (ivumat(mt).ne.0) then 
          ndi=3
          nsh=1
          nblock=llt-lft+1
          curveid=ivumat(mt)
          icid=int(curveid)
          if (icid.lt.0) then
            iicid=lcidi(-icid)
            k2=npc(iicid)
            k3=npc(iicid+1)
            nprops=nint((k3-k2)/2.0)
          else
            nprops=40+icmadd(mt)
          endif
          if (eltype.eq.'shell') if13=ind_q1
          call ext_vecumat(lft,llt,cm,a(lcma),bqs,mt,crv,npc,plc,nnm1,
     $         curveid,rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,
     $         nprops,ndi,nsh,
     $         if11,if22,if33,if12,if21,if23,if32,if13,if31,no_hsvs,
     $         ifo11,ifo22,ifo33,ifo12,ifo21,ifo23,ifo32,ifo13,ifo31)
c          call ext_vecumat(lft,llt,cm,bqs,mt,crv,npc,plc,nnm1,curveid,
c     $         rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,nprops,
c     $         ndi,nsh,if11,if22,if33,if12,if21,if23,if32,if13,if31)
        endif
c
c       check for failure
c
        if (ifipts(mt).ne.0) then
        do i=lft,llt
          if (failels(i)) then
            sig1(i)=0.
            sig2(i)=0.
            sig3(i)=0.
            sig4(i)=0.
            sig5(i)=0.
            sig6(i)=0.
            if (ifipts(mt).lt.0) hsvs(i,ind_fipts)=1.
            nfipts(i)=nfipts(i)+1
            if (ifipts(mt).ge.0) then
            failur=.true.
            ifail(i)=1
            fail(i)=0.
            elseif (nfipts(i).ge.nint(cm(mx-ifipts(mt)))) then
            failur=.true.
            ifail(i)=1
            fail(i)=0.
            endif
           endif
        enddo
      endif
      endif
c
      if (iabs(ihyper(mt)).eq.3.or.iabs(ihyper(mt)).eq.4) then
         do i=lft,llt
            d3(i)=0.
         enddo
      endif
c
c     energy calculations
c
      do 130 i=lft,llt
      einc(i)=d1(i)*sig1(i)+d2(i)*sig2(i)+d4(i)*sig4(i)+d5(i)*sig5(i)
     1       +d6(i)*sig6(i)+einc(i)
 130  continue
c
      if (iorien(mt).ne.0) then
c
c       transform stresses to element system
c
        if (ihyper(mt).ge.0) then
          do 140 i=lft,llt
          stg5(i)=sig5(i)
          stg6(i)=sig6(i)
          a11(i) = sig1(i)*hsvs(i,ind_q1)+sig4(i)*hsvs(i,ind_q2)
          a12(i) =-sig1(i)*hsvs(i,ind_q2)+sig4(i)*hsvs(i,ind_q1)
          a21(i) = sig4(i)*hsvs(i,ind_q1)+sig2(i)*hsvs(i,ind_q2)
          a22(i) =-sig4(i)*hsvs(i,ind_q2)+sig2(i)*hsvs(i,ind_q1)
          sig1(i)= hsvs(i,ind_q1)*a11(i)+hsvs(i,ind_q2)*a21(i)
          sig2(i)=-hsvs(i,ind_q2)*a12(i)+hsvs(i,ind_q1)*a22(i)
          sig4(i)= hsvs(i,ind_q1)*a12(i)+hsvs(i,ind_q2)*a22(i)
          sig5(i)=-hsvs(i,ind_q2)*stg6(i)+hsvs(i,ind_q1)*stg5(i)
          sig6(i)= hsvs(i,ind_q1)*stg6(i)+hsvs(i,ind_q2)*stg5(i)
 140      continue
        endif
c
c       transform deformation gradient back to element system
c
      if (ihyper(mt).gt.0) then
c        if (ihyper(mt).eq.1) then
          do i=lft,llt
            ftemp=hsvs(i,ind_q1)*hsvs(i,if11)+
     1            hsvs(i,ind_q2)*hsvs(i,if21)
            hsvs(i,if21)=-hsvs(i,ind_q2)*hsvs(i,if11)+
     1                    hsvs(i,ind_q1)*hsvs(i,if21)
            hsvs(i,if11)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if12)+
     1            hsvs(i,ind_q2)*hsvs(i,if22)
            hsvs(i,if22)=-hsvs(i,ind_q2)*hsvs(i,if12)+
     1                    hsvs(i,ind_q1)*hsvs(i,if22)
            hsvs(i,if12)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if13)+
     1            hsvs(i,ind_q2)*hsvs(i,if23)
            hsvs(i,if23)=-hsvs(i,ind_q2)*hsvs(i,if13)+
     1                    hsvs(i,ind_q1)*hsvs(i,if23)
            hsvs(i,if13)=ftemp
          enddo
c
c
c        transform deformation gradient back to element system
c
        else if (ihyper(mt).eq.-2) then
          do i=lft,llt
            f11=hsvs(i,ind_q1)*hsvs(i,if11)+
     1            hsvs(i,ind_q2)*hsvs(i,if21)
            f21=-hsvs(i,ind_q2)*hsvs(i,if11)+
     1                    hsvs(i,ind_q1)*hsvs(i,if21)
            f31=hsvs(i,if31)
c            hsvs(i,if11)=ftemp
            f12=hsvs(i,ind_q1)*hsvs(i,if12)+
     1            hsvs(i,ind_q2)*hsvs(i,if22)
            f22=-hsvs(i,ind_q2)*hsvs(i,if12)+
     1                    hsvs(i,ind_q1)*hsvs(i,if22)
            f32=hsvs(i,if32)
c            hsvs(i,if12)=ftemp
            f13=hsvs(i,ind_q1)*hsvs(i,if13)+
     1            hsvs(i,ind_q2)*hsvs(i,if23)
            f23=-hsvs(i,ind_q2)*hsvs(i,if13)+
     1                    hsvs(i,ind_q1)*hsvs(i,if23)
            f33=hsvs(i,if33)
c            hsvs(i,if13)=ftemp
            hsvs(i,if11)=f11*hsvs(i,ind_q1)+f12*hsvs(i,ind_q2)
            hsvs(i,if12)=-f11*hsvs(i,ind_q2)+f12*hsvs(i,ind_q1)
            hsvs(i,if13)=f13
            hsvs(i,if21)=f21*hsvs(i,ind_q1)+f22*hsvs(i,ind_q2)
            hsvs(i,if22)=-f21*hsvs(i,ind_q2)+f22*hsvs(i,ind_q1)
            hsvs(i,if23)=f23
            hsvs(i,if31)=f31*hsvs(i,ind_q1)+f32*hsvs(i,ind_q2)
            hsvs(i,if32)=-f31*hsvs(i,ind_q2)+f32*hsvs(i,ind_q1)
            hsvs(i,if33)=f33
          enddo
c
c
c        transform deformation gradient back to element system
c
c         elseif (ihyper(mt).lt.0) then
         elseif (ihyper(mt).eq.-1) then
           do i=lft,llt
            ftemp=hsvs(i,ind_q1)*hsvs(i,if11)+
     1            hsvs(i,ind_q2)*hsvs(i,if12)
            hsvs(i,if12)=-hsvs(i,ind_q2)*hsvs(i,if11)+
     1                    hsvs(i,ind_q1)*hsvs(i,if12)
            hsvs(i,if11)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if21)+
     1            hsvs(i,ind_q2)*hsvs(i,if22)
            hsvs(i,if22)=-hsvs(i,ind_q2)*hsvs(i,if21)+
     1                    hsvs(i,ind_q1)*hsvs(i,if22)
            hsvs(i,if21)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if31)+
     1            hsvs(i,ind_q2)*hsvs(i,if32)
            hsvs(i,if32)=-hsvs(i,ind_q2)*hsvs(i,if31)+
     1                    hsvs(i,ind_q1)*hsvs(i,if32)
            hsvs(i,if31)=ftemp
c
            enddo
         endif
      endif
c
      return
c900  format(/
c    1 ' *** Error User defined material in part',i12,
c    2 ' exceeds the current limit'/
c    3 10x,'of history variables for shells.')
      end
      subroutine urmatb (lft,llt,cm,capa,mt,crv,nnm1,nconstp,npc,plc,
     $     eltype,nbmbwp)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
      integer n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      common/bk07/n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      integer lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      common/bk13/lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      include 'bk19.inc'
      include 'memaia.inc'
      include 'umatss.inc'
c
c     dynamic memory allocation stuff.  This has to be a C style 
c     include because of the use of integer*8
c
      integer i_mem(0:1)
      real r_mem(0:1)
      integer*8 i8_mem(0:1)
      real*8 r8_mem(0:1)
      real*8 real8_mem
      integer*8 m_mem(0:1)
      common/dynmem/real8_mem(0:1)
      equivalence (i_mem,r_mem,i8_mem,r8_mem,real8_mem,m_mem)
      integer*8 mem_alloc, mems_alloc,mem_realloc, memh_ptr
      integer mem_length,memh_alloc,memh_realloc,memh_length
      integer memh_enum,memh_enums,memh_attach
c      integer mem_newmark,mem_getmark,memh_getmark
      integer*8 memh_detach,memh_ralloc
      integer*8 m_to_m8
      integer*8 m_to_mm
      external mem_alloc, mems_alloc, mem_realloc, mem_length
      external memh_alloc, memh_realloc, memh_length, memh_ptr
      external memh_enum,memh_attach,memh_detach,memh_ralloc
      external memh_enums
      external m_to_m8
      external m_to_mm
      integer*8 dm_x,dm_xr,dm_v,dm_vr,dm_a,dm_ar,dm_rots,dm_rotr
      integer*8 dm_xms,dm_xmr,dm_xm2,dm_xmz,dm_xma
      integer*8 dm_xtpz,dm_disp
      integer*8 dm_me1,dm_me2,dm_me3
      integer*8 dm_masp,dm_masr
      integer*8 dm_xrb,dm_yrb,dm_zrb
      integer*8 dm_axrb,dm_ayrb,dm_azrb
      integer*8 dm_rbfx,dm_rbfy,dm_rbfz
      integer*8 dm_n45,dm_n46,dm_n47
      integer len_n45,len_n46,len_n47
      common/dynmem1/ 
     . dm_x,dm_xr,dm_v,dm_vr,dm_a,dm_ar,dm_rots,dm_rotr,
     . dm_xms,dm_xmr,dm_xm2,dm_xmz,dm_xma,
     . dm_xtpz,dm_disp,
     . dm_me1,dm_me2,dm_me3,
     . dm_masp,dm_masr,
     . dm_xrb,dm_yrb,dm_zrb,
     . dm_axrb,dm_ayrb,dm_azrb,
     . dm_rbfx,dm_rbfy,dm_rbfz,
     . dm_n45,dm_n46,dm_n47
      common/dynmem1l/len_n45,len_n46,len_n47
c
      integer*8 mem_allocm,mems_allocm,mem_reallocm
      integer memh_allocm,memh_reallocm
      external mem_allocm,mems_allocm,mem_reallocm
      external memh_allocm,memh_reallocm
      integer*8 mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      integer*8 mem_reallocmchk,mem_reallocchk
      integer memh_allocmchk,memh_allocchk
      external mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      external mem_reallocmchk,mem_reallocchk
      external memh_allocmchk,memh_allocchk
      common/bux13loc/ eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      common/aux2loc/d1(nlq),d2(nlq),d3(nlq),d4(nlq),d5(nlq),d6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common/aux14loc/
     1 sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),
     2 sig5(nlq),sig6(nlq),epsps(nlq),hsvs(nlq,NHISVAR)
      common/aux33loc/
     1 ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ix5(nlq),mxt(nlq)
      common/aux35loc/rhoa(nlq),cb(nlq),davg(nlq),p(nlq)
      common/soundloc/ss(nlq),sndsp(nlq),diagm(nlq),sarea(nlq),
     . dxl(nlq)
      common/subtssloc/dt1siz(nlq)
      logical failur
      common/failcm/failur
      common/failcmloc/ifail(nlq)
      common/failuloc/sieu(nlq),fail(nlq)
c
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
      common/bk01/itherm,itemp,ntmp0,ntmp1,itempan,itempdr,itmpe(44)
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie
      common/thropt/itopaz(101)
c
      dimension cm(*),crv(lq1,2,*),nconstp(*),npc(*),plc(*)
c     dimension eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      dimension elsizv(nlq),idelev(nlq)
      logical failel,failels(nlq)
      character*(*)eltype
      dimension nbmbwp(*)
c
c     data temps/nlq*0.0/
      data failels/nlq*.false./
c
c$omp threadprivate (/aux2loc/)
c$omp threadprivate (/aux14loc/)
c$omp threadprivate (/aux33loc/)
c$omp threadprivate (/subtssloc/)
c$omp threadprivate (/aux35loc/)
c$omp threadprivate (/soundloc/)
c$omp threadprivate (/failcmloc/)
c$omp threadprivate (/failuloc/)
c$omp threadprivate (/bux13loc/)
      failel=.false.
c
      ivect=ivectr(mt)
      mx=48*(mxt(lft)-1)
      gm=cm(mx+abs(ishrmp(mt)))
      bk=cm(mx+ibulkp(mt))
      pr=(3.*bk-2.*gm)/(2.*(3.*bk+gm))
      ym=2.*(1.+pr)*gm
      ss(lft)=ym
      do 10 i=lft,llt
      einc(i)=d1(i)*sig1(i)+d4(i)*sig4(i)+d6(i)*sig6(i)
      cb(i)=ss(lft)
   10 continue
c
c     Determining number of requested history variables
c
      no_hsvs=nconstp(mxt(lft))
c
c     if ishrmp(mt).lt.0 then the flag for temperature is active
c     temperatures are stored in the temps array.
c
      if (ishrmp(mt).lt.0) then
        num_nods=2
        call usr_temps (a(ntmp0+1),a(n19),itemp,temps,lft,llt,ix1,
     .   ix2,ix3,ix4,ix4,ix4,ix4,ix4,num_nods,itopaz,0,0.,0.,0.)
      else
        do i=lft,llt
          temps(i)=0.0
        enddo
      endif
c
c     characteristic element size and element id
      do i=lft,llt
        elsizv(i)=diagm(i)
        ielm=nbmbwp(nnm1+i)
        if (ielm.gt.0) idelev(i)=lqfinv(ielm,3)
      enddo
c
c     memory pointer for extra memory data
      lcma=ncma+nmmat+ia(ncma+mxt(lft)-1)-1
c
      mte40=mt-40
c
      if (ivect.eq.0) then
        do 120 i=lft,llt
        dt1   =dt1siz(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
        epsp=epsps(i)
        temper=temps(i)
        if (no_hsvs.gt.0) then
          do 20 j=1,no_hsvs
          hsv(j)=hsvs(i,j)
 20       continue
        endif
        elsiz=elsizv(i)
        idele=idelev(i)
c
c       call user developed subroutines here
c
        go to (41,42,43,44,45,46,47,48,49,50), mte40
   41   call umat41 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   42   call umat42 (cm(mx+1),eps,sig,hsv,dt1,capa,'beam ',tt,crv,
     1   a(lcma),temper,0,elsiz,idele)
        go to 60
   43   call umat43 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   44   call umat44 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   45   call umat45 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   46   call umat46 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,a(lcma),0.0,0,0,elsiz,idele)
        go to 60
   47   call umat47 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   48   call umat48 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   49   call umat49 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   50   call umat50 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
 60     continue
c
        if (no_hsvs.gt.0) then
          do 100 j=1,no_hsvs
          hsvs(i,j)=hsv(j)
 100      continue
        endif
        d2(i)  =eps(2)
        d3(i)  =eps(3)
        epsps(i)=epsp
        if (failel) then
          sig1(i)=0.
          sig2(i)=0.
          sig3(i)=0.
          sig4(i)=0.
          sig5(i)=0.
          sig6(i)=0.
          failur=.true.
          ifail(i)=1
          fail(i)=0.
        else
          sig1(i)=sig(1)
          sig2(i)=sig(2)
          sig3(i)=sig(3)
          sig4(i)=sig(4)
          sig5(i)=sig(5)
          sig6(i)=sig(6)
        endif
 120    continue
c
c     vector umat
c
      else
c
        if (ivumat(mt).eq.0) then
        if (mt.eq.41) then
          call umat41v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.42) then
          call umat42v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.43) then
          call umat43v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.44) then
          call umat44v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .   sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .   'beam ',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.45) then
          call umat45v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.46) then
          call umat46v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,temps,failels,nlq,crv,a(lcma),0.0,0,0,elsizv,
     .     idelev)
        elseif (mt.eq.47) then
          call umat47v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.48) then
          call umat48v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,ipt,a(n8),a(n9),a(lcma),0,elsizv,idelev)
        elseif (mt.eq.49) then
          call umat49v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,ipt,temps,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.50) then
          call umat50v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        endif
        endif
        if (ivumat(mt).ne.0) then
          ndi=1
          nsh=1
          nblock=llt-lft+1
          curveid=ivumat(mt)
          icid=int(curveid)
          if (icid.lt.0) then
            iicid=lcidi(-icid)
            k2=npc(iicid)
            k3=npc(iicid+1)
            nprops=nint((k3-k2)/2.0)
          else
            nprops=40+icmadd(mt)
          endif
          if11=0.0
          if22=0.0
          if33=0.0
          if12=0.0
          if21=0.0
          if23=0.0
          if32=0.0
          if13=0.0
          if31=0.0
          ifo11=0.0
          ifo22=0.0
          ifo33=0.0
          ifo12=0.0
          ifo21=0.0
          ifo23=0.0
          ifo32=0.0
          ifo13=0.0
          ifo31=0.0
          call ext_vecumat(lft,llt,cm,a(lcma),bqs,mt,crv,npc,plc,nnm1,
     $         curveid,rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,
     $         nprops,ndi,nsh,
     $         if11,if22,if33,if12,if21,if23,if32,if13,if31,no_hsvs,
     $         ifo11,ifo22,ifo33,ifo12,ifo21,ifo23,ifo32,ifo13,ifo31)
c          call ext_vecumat(lft,llt,cm,bqs,mt,crv,npc,plc,nnm1,curveid,
c     $         rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,nprops,
c     $         ndi,nsh,if11,if22,if33,if12,if21,if23,if32,if13,if31)
        endif
      endif
c
      do 130 i=lft,llt
      einc(i)=d1(i)*sig1(i)+d4(i)*sig4(i)+d6(i)*sig6(i)+einc(i)
 130  continue
c
      return
      end
      subroutine urmatd (lft,llt,cm,capa,mt,crv,nnm1,nconstp,npc,plc,
     $     eltype,nbmbwp)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     For now, call only vectorized umat routines for discrte beam.
c     umat47v includes curves
c
      include 'nlqparm'
      include 'nhisparm.inc'
      integer n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      common/bk07/n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      integer lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      common/bk13/lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      include 'bk19.inc'
      include 'memaia.inc'
      include 'umatss.inc'
c
c     dynamic memory allocation stuff.  This has to be a C style 
c     include because of the use of integer*8
c
      integer i_mem(0:1)
      real r_mem(0:1)
      integer*8 i8_mem(0:1)
      real*8 r8_mem(0:1)
      real*8 real8_mem
      integer*8 m_mem(0:1)
      common/dynmem/real8_mem(0:1)
      equivalence (i_mem,r_mem,i8_mem,r8_mem,real8_mem,m_mem)
      integer*8 mem_alloc, mems_alloc,mem_realloc, memh_ptr
      integer mem_length,memh_alloc,memh_realloc,memh_length
      integer memh_enum,memh_enums,memh_attach
c      integer mem_newmark,mem_getmark,memh_getmark
      integer*8 memh_detach,memh_ralloc
      integer*8 m_to_m8
      integer*8 m_to_mm
      external mem_alloc, mems_alloc, mem_realloc, mem_length
      external memh_alloc, memh_realloc, memh_length, memh_ptr
      external memh_enum,memh_attach,memh_detach,memh_ralloc
      external memh_enums
      external m_to_m8
      external m_to_mm
      integer*8 dm_x,dm_xr,dm_v,dm_vr,dm_a,dm_ar,dm_rots,dm_rotr
      integer*8 dm_xms,dm_xmr,dm_xm2,dm_xmz,dm_xma
      integer*8 dm_xtpz,dm_disp
      integer*8 dm_me1,dm_me2,dm_me3
      integer*8 dm_masp,dm_masr
      integer*8 dm_xrb,dm_yrb,dm_zrb
      integer*8 dm_axrb,dm_ayrb,dm_azrb
      integer*8 dm_rbfx,dm_rbfy,dm_rbfz
      integer*8 dm_n45,dm_n46,dm_n47
      integer len_n45,len_n46,len_n47
      common/dynmem1/ 
     . dm_x,dm_xr,dm_v,dm_vr,dm_a,dm_ar,dm_rots,dm_rotr,
     . dm_xms,dm_xmr,dm_xm2,dm_xmz,dm_xma,
     . dm_xtpz,dm_disp,
     . dm_me1,dm_me2,dm_me3,
     . dm_masp,dm_masr,
     . dm_xrb,dm_yrb,dm_zrb,
     . dm_axrb,dm_ayrb,dm_azrb,
     . dm_rbfx,dm_rbfy,dm_rbfz,
     . dm_n45,dm_n46,dm_n47
      common/dynmem1l/len_n45,len_n46,len_n47
c
      integer*8 mem_allocm,mems_allocm,mem_reallocm
      integer memh_allocm,memh_reallocm
      external mem_allocm,mems_allocm,mem_reallocm
      external memh_allocm,memh_reallocm
      integer*8 mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      integer*8 mem_reallocmchk,mem_reallocchk
      integer memh_allocmchk,memh_allocchk
      external mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      external mem_reallocmchk,mem_reallocchk
      external memh_allocmchk,memh_allocchk
      common/bux13loc/ eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      common/aux2loc/d1(nlq),d2(nlq),d3(nlq),d4(nlq),d5(nlq),d6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common/aux14loc/
     1 sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),
     2 sig5(nlq),sig6(nlq),epsps(nlq),hsvs(nlq,NHISVAR)
      common/aux33loc/
     1 ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ix5(nlq),mxt(nlq)
      common/aux35loc/rhoa(nlq),cb(nlq),davg(nlq),p(nlq)
      common/soundloc/ss(nlq),sndsp(nlq),diagm(nlq),sarea(nlq),
     . dxl(nlq)
      common/subtssloc/dt1siz(nlq)
c
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
      common/bk01/itherm,itemp,ntmp0,ntmp1,itempan,itempdr,itmpe(44)
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie
      common/slcntr/islcnt(21)
      common/thropt/itopaz(101)
c
c     capa here is ies in conmd6 and bqs in usrmat (all are internal energy)
c     real capa
      dimension capa(*),cm(*),crv(lq1,2,*),nconstp(*),npc(*),plc(*)
c     dimension eps(6),sig(6),hsv(6),temps(nlq)
      dimension elsizv(nlq),idelev(nlq)
      logical failel,failels(nlq)
      character*(*)eltype
      dimension nbmbwp(*)
c
c     data temps/nlq*0.0/
      data failels/nlq*.false./
c
c$omp threadprivate (/aux2loc/)
c$omp threadprivate (/aux14loc/)
c$omp threadprivate (/aux33loc/)
c$omp threadprivate (/aux35loc/)
c$omp threadprivate (/soundloc/)
c$omp threadprivate (/subtssloc/)
c$omp threadprivate (/bux13loc/)
      failel=.false.
      ivect=ivectr(mt)
      mx=48*(mxt(lft)-1)
c
c     Determining number of requested history variables
c
      no_hsvs=nconstp(mxt(lft))
c
c     if ishrmp(mt).lt.0 then the flag for temperature is active
c     temperatures are stored in the temps array.
c
      if (ishrmp(mt).lt.0) then
        num_nods=2
        call usr_temps (a(ntmp0+1),a(n19),itemp,temps,lft,llt,
     .   ix1,ix2,ix3,ix4,ix4,ix4,ix4,ix4,num_nods,itopaz,0,
     .   0.,0.,0.)
      else
        do i=lft,llt
          temps(i)=0.0
        enddo
      endif
c
c     characteristic element size and element id
      do i=lft,llt
        elsizv(i)=diagm(i)
        ielm=nbmbwp(nnm1+i)
        if (ielm.gt.0) idelev(i)=lqfinv(ielm,3)
      enddo
c
c     memory pointer for extra memory data
      lcma=ncma+nmmat+ia(ncma+mxt(lft)-1)-1
c
      mte40=mt-40
c
      if (ivect.eq.0) then
        do 120 i=lft,llt
        dt1   =dt1siz(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=wxxdt(i)
        eps(5)=wyydt(i)
        eps(6)=wzzdt(i)
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
        temper=temps(i)
        epsp=epsps(i)
        if (no_hsvs.gt.0) then
          do 20 j=1,no_hsvs
          hsv(j)=hsvs(i,j)
 20       continue
        endif
        elsiz=elsizv(i)
        idele=idelev(i)
c
c       call user developed subroutines here
c
        go to (41,42,43,44,45,46,47,48,49,50), mte40
   41   call umat41 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   42   call umat42 (cm(mx+1),eps,sig,hsv,dt1,capa,'dbeam',tt,crv,
     1   a(lcma),temper,0,elsiz,idele)
        go to 60
   43   call umat43 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   44   call umat44 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   45   call umat45 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   46   call umat46 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,0.0,0,a(lcma),0,elsiz,idele)
        go to 60
   47   call umat47 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   48   call umat48 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   49   call umat49 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   50   call umat50 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
 60     continue
c
        if (no_hsvs.gt.0) then
          do 100 j=1,no_hsvs
          hsvs(i,j)=hsv(j)
 100      continue
        endif
        sig1(i)=sig(1)
        sig2(i)=sig(2)
        sig3(i)=sig(3)
        sig4(i)=sig(4)
        sig5(i)=sig(5)
        sig6(i)=sig(6)
        epsps(i)=epsp
 120  continue
c
c     vector umat
c
      else
        if (ivumat(mt).eq.0) then
        if (mt.eq.41) then
          call umat41v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.42) then
          call umat42v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.43) then
          call umat43v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.44) then
          call umat44v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.45) then
          call umat45v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.46) then
          call umat46v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,a(lcma),0.0,0,0,elsizv,
     .     idelev)
        elseif (mt.eq.47) then
          call umat47v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.48) then
          call umat48v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,ipt,a(n8),a(n9),a(lcma),0,elsizv,idelev)
        elseif (mt.eq.49) then
          call umat49v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,ipt,temps,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.50) then
          call umat50v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        endif
        endif
        if (ivumat(mt).ne.0) then
          ndi=1
          nsh=1
          nblock=llt-lft+1
          curveid=ivumat(mt)
          icid=int(curveid)
          if (icid.lt.0) then
            iicid=lcidi(-icid)
            k2=npc(iicid)
            k3=npc(iicid+1)
            nprops=nint((k3-k2)/2.0)
          else
            nprops=40+icmadd(mt)
          endif
          call ext_vecumat(lft,llt,cm,a(lcma),bqs,mt,crv,npc,plc,nnm1,
     $         curveid,rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,
     $         nprops,ndi,nsh,
     $         if11,if22,if33,if12,if21,if23,if32,if13,if31,no_hsvs,
     $         ifo11,ifo22,ifo33,ifo12,ifo21,ifo23,ifo32,ifo13,ifo31)
c          call ext_vecumat(lft,llt,cm,bqs,mt,crv,npc,plc,nnm1,curveid,
c     $         rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,nprops,
c     $         ndi,nsh,if11,if22,if33,if12,if21,if23,if32,if13,if31)
        endif
      endif
c
      do 130 i=lft,llt
      einc(i)=d1(i)*sig1(i)+d4(i)*sig4(i)+d6(i)*sig6(i)+einc(i)
  130 continue
c
      return
      end
      subroutine urmatt (lft,llt,cm,capa,mt,crv,nnm1,nconstp,npc,plc,
     $     eltype,nbmbwp)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
      integer n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      common/bk07/n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      integer lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      common/bk13/lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      include 'bk19.inc'
      include 'memaia.inc'
      include 'umatss.inc'
c
c     dynamic memory allocation stuff.  This has to be a C style 
c     include because of the use of integer*8
c
      integer i_mem(0:1)
      real r_mem(0:1)
      integer*8 i8_mem(0:1)
      real*8 r8_mem(0:1)
      real*8 real8_mem
      integer*8 m_mem(0:1)
      common/dynmem/real8_mem(0:1)
      equivalence (i_mem,r_mem,i8_mem,r8_mem,real8_mem,m_mem)
      integer*8 mem_alloc, mems_alloc,mem_realloc, memh_ptr
      integer mem_length,memh_alloc,memh_realloc,memh_length
      integer memh_enum,memh_enums,memh_attach
c      integer mem_newmark,mem_getmark,memh_getmark
      integer*8 memh_detach,memh_ralloc
      integer*8 m_to_m8
      integer*8 m_to_mm
      external mem_alloc, mems_alloc, mem_realloc, mem_length
      external memh_alloc, memh_realloc, memh_length, memh_ptr
      external memh_enum,memh_attach,memh_detach,memh_ralloc
      external memh_enums
      external m_to_m8
      external m_to_mm
      integer*8 dm_x,dm_xr,dm_v,dm_vr,dm_a,dm_ar,dm_rots,dm_rotr
      integer*8 dm_xms,dm_xmr,dm_xm2,dm_xmz,dm_xma
      integer*8 dm_xtpz,dm_disp
      integer*8 dm_me1,dm_me2,dm_me3
      integer*8 dm_masp,dm_masr
      integer*8 dm_xrb,dm_yrb,dm_zrb
      integer*8 dm_axrb,dm_ayrb,dm_azrb
      integer*8 dm_rbfx,dm_rbfy,dm_rbfz
      integer*8 dm_n45,dm_n46,dm_n47
      integer len_n45,len_n46,len_n47
      common/dynmem1/ 
     . dm_x,dm_xr,dm_v,dm_vr,dm_a,dm_ar,dm_rots,dm_rotr,
     . dm_xms,dm_xmr,dm_xm2,dm_xmz,dm_xma,
     . dm_xtpz,dm_disp,
     . dm_me1,dm_me2,dm_me3,
     . dm_masp,dm_masr,
     . dm_xrb,dm_yrb,dm_zrb,
     . dm_axrb,dm_ayrb,dm_azrb,
     . dm_rbfx,dm_rbfy,dm_rbfz,
     . dm_n45,dm_n46,dm_n47
      common/dynmem1l/len_n45,len_n46,len_n47
c
      integer*8 mem_allocm,mems_allocm,mem_reallocm
      integer memh_allocm,memh_reallocm
      external mem_allocm,mems_allocm,mem_reallocm
      external memh_allocm,memh_reallocm
      integer*8 mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      integer*8 mem_reallocmchk,mem_reallocchk
      integer memh_allocmchk,memh_allocchk
      external mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      external mem_reallocmchk,mem_reallocchk
      external memh_allocmchk,memh_allocchk
      common/bux13loc/ eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      common/aux2loc/d1(nlq),d2(nlq),d3(nlq),d4(nlq),d5(nlq),d6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common/aux14loc/
     1 sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),
     2 sig5(nlq),sig6(nlq),epsps(nlq),hsvs(nlq,NHISVAR)
      common/aux33loc/
     1 ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ix5(nlq),mxt(nlq)
      common/aux35loc/rhoa(nlq),cb(nlq),davg(nlq),p(nlq)
      common/failcmloc/ifail(nlq)
      common/failuloc/sieu(nlq),fail(nlq)
      common/soundloc/ss(nlq),sndsp(nlq),diagm(nlq),sarea(nlq),
     . dxl(nlq)
      common/subtssloc/dt1siz(nlq)
c
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
      common/bk01/itherm,itemp,ntmp0,ntmp1,itempan,itempdr,itmpe(44)
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie
      logical failur
      common/failcm/failur
      common/thropt/itopaz(101)
      common/slcntr/islcnt(21)
c
c     real capa
      dimension capa(*),cm(*),crv(lq1,2,*),nconstp(*),npc(*),plc(*)
c     dimension eps(6),sig(6),hsv(6),temps(nlq)
      dimension elsizv(nlq),idelev(nlq)
      logical failel,failels(nlq)
      character*(*)eltype
      dimension nbmbwp(*)
c
c     data temps/nlq*0.0/
      data failels/nlq*.false./
c
c$omp threadprivate (/aux2loc/)
c$omp threadprivate (/aux14loc/)
c$omp threadprivate (/aux33loc/)
c$omp threadprivate (/aux35loc/)
c$omp threadprivate (/failcmloc/)
c$omp threadprivate (/failuloc/)
c$omp threadprivate (/soundloc/)
c$omp threadprivate (/subtssloc/)
c$omp threadprivate (/bux13loc/)
      ivect=ivectr(mt)
      mx=48*(mxt(lft)-1)
      gm=cm(mx+abs(ishrmp(mt)))
      bk=cm(mx+ibulkp(mt))
      ym=9.*bk*gm/(3.*bk+gm)
      ss(lft)=ym
      do 10 i=lft,llt
      einc(i)=d1(i)*sig1(i)
      cb(i)=ss(lft)
   10 continue
c
c     Determining number of requested history variables
c
      no_hsvs=nconstp(mxt(lft))
c
c     if ishrmp(mt).lt.0 then the flag for temperature is active
c     temperatures are stored in the temps array.
c
      if (ishrmp(mt).lt.0) then
        num_nods=2
        call usr_temps (a(ntmp0+1),a(n19),itemp,temps,lft,llt,
     .   ix1,ix2,ix3,ix4,ix4,ix4,ix4,ix4,num_nods,itopaz,
     .   0,0.,0.,0.)
      else
        do i=lft,llt
          temps(i)=0.0
        enddo
      endif
c
c     characteristic element size
      do i=lft,llt
        elsizv(i)=diagm(i)
        ielm=nbmbwp(nnm1+i)
        if (ielm.gt.0) idelev(i)=lqfinv(ielm,3)
      enddo
c
c     memory pointer for extra memory data
      lcma=ncma+nmmat+ia(ncma+mxt(lft)-1)-1
c
      mte40=mt-40
c
      if (ivect.eq.0) then
        do 120 i=lft,llt
        dt1   =dt1siz(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        epsp=epsps(i)
        temper=temps(i)
        failel=.false.
        if (no_hsvs.gt.0) then
          do 20 j=1,no_hsvs
          hsv(j)=hsvs(i,j)
 20       continue
        endif
        elsiz=elsizv(i)
        idele=idelev(i)
c
c       call user developed subroutines here
c
        go to (41,42,43,44,45,46,47,48,49,50), mte40
   41   call umat41 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   42   call umat42 (cm(mx+1),eps,sig,hsv,dt1,capa,'tbeam',tt,crv,
     1   a(lcma),temper,0,elsiz,idele)
        go to 60
   43   call umat43 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   44   call umat44 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   45   call umat45 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   46   call umat46 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,a(lcma),0.0,0,0,elsiz,idele)
        go to 60
   47   call umat47 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   48   call umat48 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   49   call umat49 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
   50   call umat50 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,a(lcma),0,elsiz,idele)
 60   continue
c
      if (no_hsvs.gt.0) then
        do 100 j=1,no_hsvs
        hsvs(i,j)=hsv(j)
 100   continue
      endif
c
      if (failel) then
        einc(i)=0.
        sig1(i)=0.
        sig2(i)=0.
        sig3(i)=0.
        failur=.true.
        ifail(i)=1
        fail(i)=0.
      else
        sig1(i)=sig(1)
        sig2(i)=sig(2)
        sig3(i)=sig(3)
      endif
      d2(i)  =eps(2)
      d3(i)  =eps(3)
      epsps(i)=epsp
 120  continue
c
c     vector umat
c
      else
        if (ivumat(mt).eq.0) then
        if (mt.eq.41) then
          call umat41v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.42) then
          call umat42v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.43) then
          call umat43v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.44) then
          call umat44v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.45) then
          call umat45v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.46) then
          call umat46v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,a(lcma),0.0,0,0,elsizv,
     .     idelev)
        elseif (mt.eq.47) then
          call umat47v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.48) then
          call umat48v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,ipt,a(n8),a(n9),a(lcma),0,elsizv,idelev)
        elseif (mt.eq.49) then
          call umat49v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,ipt,temps,a(lcma),0,elsizv,idelev)
        elseif (mt.eq.50) then
          call umat50v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,a(lcma),0,elsizv,idelev)
        endif
        endif
        if (ivumat(mt).ne.0) then
          ndi=1
          nsh=1
          nblock=llt-lft+1
          curveid=ivumat(mt)
          icid=int(curveid)
          if (icid.lt.0) then
            iicid=lcidi(-icid)
            k2=npc(iicid)
            k3=npc(iicid+1)
            nprops=nint((k3-k2)/2.0)
          else
            nprops=40+icmadd(mt)
          endif
          call ext_vecumat(lft,llt,cm,a(lcma),bqs,mt,crv,npc,plc,nnm1,
     $         curveid,rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,
     $         nprops,ndi,nsh,
     $         if11,if22,if33,if12,if21,if23,if32,if13,if31,no_hsvs,
     $         ifo11,ifo22,ifo33,ifo12,ifo21,ifo23,ifo32,ifo13,ifo31)
c          call ext_vecumat(lft,llt,cm,bqs,mt,crv,npc,plc,nnm1,curveid,
c     $         rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,nprops,
c     $         ndi,nsh,if11,if22,if33,if12,if21,if23,if32,if13,if31)
        endif
c
        do i=lft,llt
          if (failels(i)) then
            sig1(i)=0.
            sig2(i)=0.
            sig3(i)=0.
            fail(i)=0.
            ifail(i)=1
            failur=.true.
         endif
        enddo
      endif
c
      do 130 i=lft,llt
      einc(i)=d1(i)*sig1(i)+einc(i)
 130  continue
c
      return
      end
      subroutine umat41 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,cma,qmat,elsiz,idele)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     isotropic elastic material (sample user subroutine)
c
c     Variables
c
c     cm(1)=first material constant, here young's modulus
c     cm(2)=second material constant, here poisson's ratio
c        .
c        .
c        .
c     cm(n)=nth material constant
c
c     eps(1)=local x  strain increment
c     eps(2)=local y  strain increment
c     eps(3)=local z  strain increment
c     eps(4)=local xy strain increment
c     eps(5)=local yz strain increment
c     eps(6)=local zx strain increment
c
c     sig(1)=local x  stress
c     sig(2)=local y  stress
c     sig(3)=local z  stress
c     sig(4)=local xy stress
c     sig(5)=local yz stress
c     sig(6)=local zx stress
c
c     hsv(1)=1st history variable
c     hsv(2)=2nd history variable
c        .
c        .
c        .
c        .
c     hsv(n)=nth history variable
c
c     dt1=current time step size
c     capa=reduction factor for transverse shear
c     etype:
c        eq."solid" for solid elements
c        eq."sld2d" for shell forms 13, 14, and 15 (2D solids)
c        eq."shl_t" for shell forms 25, 26, and 27 (shells with thickness stretch)
c        eq."shell" for all other shell elements plus thick shell forms 1 and 2
c        eq."tshel" for thick shell forms 3 and 5
c        eq."hbeam" for beam element forms 1 and 11
c        eq."tbeam" for beam element form 3 (truss)
c        eq."dbeam" for beam element form 6 (discrete)
c        eq."beam " for all other beam elements
c
c     tt=current problem time.
c
c     temper=current temperature
c
c     failel=flag for failure, set to .true. to fail an integration point,
c            if .true. on input the integration point has failed earlier
c
c     crv=array representation of curves in keyword deck
c
c     cma=additional memory for material data defined by LMCA at 
c       6th field of 2nd crad of *DATA_USER_DEFINED
c
c     elsiz=characteristic element size
c
c     idele=element id
c
c     All transformations into the element local system are
c     performed prior to entering this subroutine.  Transformations
c     back to the global system are performed after exiting this
c     routine.
c
c     All history variables are initialized to zero in the input
c     phase. Initialization of history variables to nonzero values
c     may be done during the first call to this subroutine for each
c     element.
c
c     Energy calculations for the dyna3d energy balance are done
c     outside this subroutine.
c
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      logical failel
      character*5 etype
c
      if (ncycle.eq.1) then
        if (cm(16).ne.1234567) then
          call usermsg('mat41')
        endif
      endif
c
c     compute shear modulus, g
c
      g2 =abs(cm(1))/(1.+cm(2))
      g  =.5*g2
c
      if (etype.eq.'solid'.or.etype.eq.'shl_t'.or.
     1     etype.eq.'sld2d'.or.etype.eq.'tshel') then
        if (cm(16).eq.1234567) then
          call mitfail3d(cm,eps,sig,epsp,hsv,dt1,capa,failel,tt,crv)
        else
          if (.not.failel) then
          davg=(-eps(1)-eps(2)-eps(3))/3.
          p=-davg*abs(cm(1))/(1.-2.*cm(2))
          sig(1)=sig(1)+p+g2*(eps(1)+davg)
          sig(2)=sig(2)+p+g2*(eps(2)+davg)
          sig(3)=sig(3)+p+g2*(eps(3)+davg)
          sig(4)=sig(4)+g*eps(4)
          sig(5)=sig(5)+g*eps(5)
          sig(6)=sig(6)+g*eps(6)
          if (cm(1).lt.0.) then            
            if (sig(1).gt.cm(5)) failel=.true.
          endif
          endif
        end if
c
      else if (etype.eq.'shell') then
        if (cm(16).eq.1234567) then
          call mitfailure(cm,eps,sig,epsp,hsv,dt1,capa,failel,tt,crv)
        else
          if (.not.failel) then
          gc    =capa*g
          q1    =abs(cm(1))*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
          q3    =1./(q1+g2)
          eps(3)=-q1*(eps(1)+eps(2))*q3
          davg  =(-eps(1)-eps(2)-eps(3))/3.
          p     =-davg*abs(cm(1))/(1.-2.*cm(2))
          sig(1)=sig(1)+p+g2*(eps(1)+davg)
          sig(2)=sig(2)+p+g2*(eps(2)+davg)
          sig(3)=0.0
          sig(4)=sig(4)+g *eps(4)
          sig(5)=sig(5)+gc*eps(5)
          sig(6)=sig(6)+gc*eps(6)
          if (cm(1).lt.0.) then
            if (sig(1).gt.cm(5)) failel=.true.
          endif
          endif
        end if
      elseif (etype.eq.'beam ' ) then
          q1    =cm(1)*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
          q3    =q1+2.0*g
          gc    =capa*g
          deti  =1./(q3*q3-q1*q1)
          c22i  = q3*deti
          c23i  =-q1*deti
          fac   =(c22i+c23i)*q1
          eps(2)=-eps(1)*fac-sig(2)*c22i-sig(3)*c23i
          eps(3)=-eps(1)*fac-sig(2)*c23i-sig(3)*c22i
          davg  =(-eps(1)-eps(2)-eps(3))/3.
          p     =-davg*cm(1)/(1.-2.*cm(2))
          sig(1)=sig(1)+p+g2*(eps(1)+davg)
          sig(2)=0.0
          sig(3)=0.0
          sig(4)=sig(4)+gc*eps(4)
          sig(5)=0.0
          sig(6)=sig(6)+gc*eps(6)
c
      elseif (etype.eq.'tbeam') then
        q1    =cm(1)*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
        q3    =q1+2.0*g
        deti  =1./(q3*q3-q1*q1)
        c22i  = q3*deti
        c23i  =-q1*deti
        fac   =(c22i+c23i)*q1
        eps(2)=-eps(1)*fac
        eps(3)=-eps(1)*fac
        davg  =(-eps(1)-eps(2)-eps(3))/3.
        p     =-davg*cm(1)/(1.-2.*cm(2))
        sig(1)=sig(1)+p+g2*(eps(1)+davg)
        sig(2)=0.0
        sig(3)=0.0
c
      else
c       write(iotty,10) etype
c       write(iohsp,10) etype
c       write(iomsg,10) etype
c       call adios(2)
        cerdat(1)=etype
        call lsmsg(3,MSG_SOL+1150,ioall,ierdat,rerdat,cerdat,0)
      endif
c
c10   format(/
c    1 ' *** Error element type ',a,' can not be',
c    2 '           run with the current material model.')
      return
      end
c------------------------------------------------------------------
c      subroutine umat45 (cm,eps,sig,epsp,hsv,dt1,capa,
c     . etype,time,temp,failel,crv,cma,qmat,elsiz,idele)
      subroutine umat45 (cm,eps,sig,epsp,hsv,dt1,capa,
     . etype,time,temp,failel,crv,xM,xS,xN,vrat,Eb1)
c*******************************************************************
c     cardiac mechanics using Fung's exponential function          *
c                                                                  *
c     follows Guccione et al, "Mechanism Underlying Mechanical     *
c     Dysfunction in Border Zone of LV Aneurysm", 2001, Society of *
c     Thoracic Surgeons.                                           *
c                                                                  *
c     follows Gucione et al., Finite Element Stress Analysis of    *
c     LV mechanics in the Beating Dog Heart, 1995, J. Biomech.     *
c                                                                  *
c     created by Kay Sun                                           *
c     02/2007                                                      *
c                                                                  *
c*******************************************************************

c -------------------------------------------------------------- 
      include 'iounits.inc'
c      common/bk06/idmmy,iaddp,ifil,maxsiz,ncycle,ctime(2,30)
      character*(*) etype
      dimension  cm(*),sig(*),hsv(*),xM(*),xS(*),xN(*)
c
      real sp1,sp2,sp3,sp4,sp5,sp6
      real lr,lo
c
c     terms added by XZ for N_overlap
      real x_overlap, x_maxoverlap

c     terms added by LC Lee to include time dependency
c      real t0, tr
c     define t0 as an input 8/22/13
c      t0 = 0.3

c     terms added to include Ken's model (XZ 4/10/14)
c     A: number of attached heads; D: number of detached heads; N_bound:total of A at previous step
      real*8 A(21), N_overlap, N_bound,D
      real*8 N_total,N_total0

c     l(21):21 bins for cross-bridge lengths; k_a=k1; k_d=k-1
c     S_dk:sum of k_a(i)
      real*8 l(21), k_a(21),k_d(21),S_dk,l_mv(21)

c     beta: used to calculate calcium
      real Ca,aCa,pCa 
c     
c     matrix used to calculate ODEs    
      real*8 x_cb0(21),x_cb(21),x_cb_new(21)
      real*8 dxdt1(21),dxdt2(21),dxdt3(21),dxdt4(21)
      real*8 f_dxdt1(21),f_dxdt2(21),f_dxdt3(21)
      real*8 Ak1(21),Ak2(21),Ak3(21),Ak4(21)
      real*8 dNdt
      real*8 p_time,time1,time2,t_p
      real*8 ka1(21),ka2,ka3(21),ka4(21)
c
c     passive material parameters (1 to 4)
      ta=cm(1)
      b2=cm(2)
      b3=cm(3)
      b4=cm(4)
c      B=cm(5)
c     
c     active material parameters and bulk modulus (20)
      lo=cm(5)
      lr=cm(6)
c     
      t_act=cm(7)
c
c     parameters determining attachment rate functions April 2016, XZ
c     cross_bridge density per cm*2
      rho_c=cm(8)
c     Boltzmann's constant (kg*cm2/s2/K)
      k_b=cm(9)
c     Experimental temperature (K)
      T_k=cm(10)
c
c      Ca_max=cm(7)
c      m=cm(9)
c      bb=cm(10)
c      Ca_o=cm(11)
c      T_max=cm(12)
c
c     a_on, a_off: relate to change rate of N_total
      a_on=cm(11)
      a_off=cm(12)
c
c     compliance factor for distribution shift
      Cp=cm(13)    
c
c     lthin,lthick,lbare,l_falloff: relates to N_overlap
      lthin=cm(14)
      lthick=cm(15)
      lbare=cm(16)
      l_falloff=cm(17)
c     K---- stiffness of elastic spring
      K=cm(18)
c     other active terms
c     bulk modulus
      bp=cm(20)
c
c     coefficient for k1(x) calculation 
      Cka=cm(21)
      f_p2=cm(22)
c     coefficients for k-1(x) calculation
      kd1=cm(23)
      kdb1=cm(24)
      kd2=cm(25)
      kdb2=cm(26)
c
      x_ps=cm(27)
c
c     constants: strengths of two potential cooperative efffects
      kplus=cm(28)
      kminus=cm(29)

c     added two cm for another two stress components, March, 2016 by XZ
      s_2=cm(30)
      s_3=cm(31)
c
      Cka_2=cm(32)
c     factor scaling calcium transient
      fCa=cm(33)
      fCa_2=cm(34)
c
c     deformation gradient from hsv 13 to 21
c     IHYPER=1, F=R'*F, in global frames
c
      F1=hsv(37)
      F4=hsv(38)
      F7=hsv(39)
      F2=hsv(40)
      F5=hsv(41)
      F8=hsv(42)
      F3=hsv(43)
      F6=hsv(44)
      F9=hsv(45)
c
c     save fiber directions in hsv
c     xM is fiber direction, xN is cross-fiber
      if (time.eq.0.0) then
        hsv(1)=xM(1)
        hsv(2)=xM(2)
        hsv(3)=xM(3)
c        hsv(5)=xS(1)
c        hsv(6)=xS(2)
c        hsv(7)=xS(3)
        hsv(4)=xN(1)
        hsv(5)=xN(2)
        hsv(6)=xN(3)
c      save initial sarcomere length 11/18/13
        hsv(11)=lr
        hsv(12)=lr
c      save initial distribution of myosin heads 4/16/14
        N_total=0
        do 900 i=1, 21
           A(i)=0
900     continue 
        hsv(14)=A(1)
        hsv(15)=A(2)
        hsv(16)=A(3)
        hsv(17)=A(4)
        hsv(18)=A(5)
        hsv(19)=A(6)
        hsv(20)=A(7)
        hsv(21)=A(8)
        hsv(22)=A(9)
        hsv(23)=A(10)
        hsv(24)=A(11)
        hsv(25)=A(12)
        hsv(26)=A(13)
        hsv(27)=A(14)
        hsv(28)=A(15)
        hsv(29)=A(16)
        hsv(30)=A(17)
        hsv(31)=A(18)
        hsv(32)=A(19)
        hsv(33)=A(20)
        hsv(34)=A(21)
        hsv(35)=N_total
c
      endif
c
c     cross product of fiber and cross fiber to get orthogonal    
c     coordinate transformation from material to global
      s11=hsv(1)
      s12=hsv(2)
      s13=hsv(3)
      s21=hsv(5)*hsv(3)-hsv(2)*hsv(6)
      s22=hsv(6)*hsv(1)-hsv(3)*hsv(4)
      s23=hsv(4)*hsv(2)-hsv(1)*hsv(5)
      s31=hsv(4)
      s32=hsv(5)
      s33=hsv(6)
c
c     preparing material parameters
      tb2=2.*b2
      tb3=2.*b3
      tb4=2.*b4
c
cdjb  limit magnitude of exponential term to prevent overflow errors
c     not use in this case
c      qlim=0.1
c      aqmx=log(qlim*bp/tb2)
c
c     calculation of sound speed for time step determination
      tc=0.5*ta
      bmbill=0.6666*max(b2,b3)
      smbill=2.*max(b3,b4)
      scb=max(bmbill,bp) + 1.33333*smbill
c
c     Jacobian of F
      df  =F1*(F5*F9-F8*F6)-
     .        F2*(F4*F9-F7*F6)+
     .        F3*(F4*F8-F7*F5)
      df23=df**(-2.0/3.0)     
c
c     green-st.venant strain  e=(F^t*F-1)/2
c
      e1=.5*(df23*(F1*F1+F4*F4+F7*F7)-1.)
      e2=.5*(df23*(F2*F2+F5*F5+F8*F8)-1.)
      e3=.5*(df23*(F3*F3+F6*F6+F9*F9)-1.)
      e4=.5*(df23*(F1*F2+F4*F5+F7*F8))
      e5=.5*(df23*(F2*F3+F5*F6+F8*F9))
      e6=.5*(df23*(F1*F3+F4*F6+F7*F9))

c     calculate total green-st.venant strain by Xiaoyan Zhang (06/17/2013)
      te1=.5*(df*(F1*F1+F4*F4+F7*F7)-1.)
      te2=.5*(df*(F2*F2+F5*F5+F8*F8)-1.)
      te3=.5*(df*(F3*F3+F6*F6+F9*F9)-1.)
      te4=.5*(df*(F1*F2+F4*F5+F7*F8))
      te5=.5*(df*(F2*F3+F5*F6+F8*F9))
      te6=.5*(df*(F1*F3+F4*F6+F7*F9))
c
c     transform green's strain in global to material axes
      a11   =e1*s11+e4*s12+e6*s13
      a12   =e1*s21+e4*s22+e6*s23
      a13   =e1*s31+e4*s32+e6*s33
      a21   =e4*s11+e2*s12+e5*s13
      a22   =e4*s21+e2*s22+e5*s23
      a23   =e4*s31+e2*s32+e5*s33
      a31   =e6*s11+e5*s12+e3*s13
      a32   =e6*s21+e5*s22+e3*s23
      a33   =e6*s31+e5*s32+e3*s33
      e1=    s11*a11   +s12*a21   +s13*a31
      e2=    s21*a12   +s22*a22   +s23*a32
      e3=    s31*a13   +s32*a23   +s33*a33
      e4=    s11*a12   +s12*a22   +s13*a32
      e5=    s21*a13   +s22*a23   +s23*a33
      e6=    s11*a13   +s12*a23   +s13*a33

c     transformation of total green's strain added by Xiaoyan Zhang (06/17/2013)
      a11   =te1*s11+te4*s12+te6*s13
      a12   =te1*s21+te4*s22+te6*s23
      a13   =te1*s31+te4*s32+te6*s33
      a21   =te4*s11+te2*s12+te5*s13
      a22   =te4*s21+te2*s22+te5*s23
      a23   =te4*s31+te2*s32+te5*s33
      a31   =te6*s11+te5*s12+te3*s13
      a32   =te6*s21+te5*s22+te3*s23
      a33   =te6*s31+te5*s32+te3*s33
      te1=    s11*a11   +s12*a21   +s13*a31
      te2=    s21*a12   +s22*a22   +s23*a32
      te3=    s31*a13   +s32*a23   +s33*a33
      te4=    s11*a12   +s12*a22   +s13*a32
      te5=    s21*a13   +s22*a23   +s23*a33
      te6=    s11*a13   +s12*a23   +s13*a33
c
c     partial derivative of W with repect to q
c     in material coordinate system
      aq=b2* e1**2+
     1   b3*(e2**2+e3**2+2.*e5**2)+
     1  tb4*(e4**2+e6**2)
c      aq=min(aq,aqmx)
      dwdq=tc*exp(aq)

c    adding violating case by 2/11/16 by Xiaoyan Zhang 
      hsv(11)=hsv(12)
      if (e1.gt. 0.5*((lo/lr)**2.0-1.0)) then
         saclen=lr*sqrt(2.0*e1+1.0)
      else
         saclen=hsv(11)
      endif
c
c     saving sarcomere length 11/18/13 by Xiaoyan Zhang
c      hsv(11)=hsv(12)
      hsv(12)=saclen
c     calculate the changes in half sarcomere length 
      hsv(13)=0.5*(hsv(12)-hsv(11))
c      hsv(13)=hsv(12)-hsv(11)
c
c***********************************************************************************
c     active contraction incorporating Dr. Campbell's cellular contraction model   *
c     modified by Xiaoyan Zhang 2014                                               *
c***********************************************************************************
c
c     initialization: reset x_cb0(i,1) at values from previous time step 
      x_cb0(1)=hsv(14)
      x_cb0(2)=hsv(15)
      x_cb0(3)=hsv(16)
      x_cb0(4)=hsv(17)
      x_cb0(5)=hsv(18)
      x_cb0(6)=hsv(19)
      x_cb0(7)=hsv(20)
      x_cb0(8)=hsv(21)
      x_cb0(9)=hsv(22)
      x_cb0(10)=hsv(23)
      x_cb0(11)=hsv(24)
      x_cb0(12)=hsv(25)
      x_cb0(13)=hsv(26)
      x_cb0(14)=hsv(27)
      x_cb0(15)=hsv(28)
      x_cb0(16)=hsv(29)
      x_cb0(17)=hsv(30)
      x_cb0(18)=hsv(31)
      x_cb0(19)=hsv(32)
      x_cb0(20)=hsv(33)
      x_cb0(21)=hsv(34)
      N_total0=hsv(35)
c 
      if (time.ge.t_act) then
c     make active force zero when reaching l0 Sep, 2016 by XZ
       if (e1.gt. 0.5*((lo/lr)**2.0-1.0)) then
         p_time=time-dt1 
         time1=p_time+0.5*dt1
         time2=p_time+dt1
c
c     calculate calcium concentration in nM
c     modified calcium curve May, 2016 by XZ
         t_p=t_act+0.01
         if (time.ge.t_p) then
           pCa=0.5*exp(-((time-t_p)*cm(33))**cm(34))
         else
           pCa=(time-t_act)/0.02
         endif
         Ca=0.1+1000*sin(3.14*pCa)
c
c
         aca=a_on*Ca
c
c     calculate N_overlap based on hsl at current time step
         hsl=0.5*hsv(12)
         if (hsl.le.cm(14)) then 
            N_overlap=1-cm(17)*(cm(14)-hsl)
            if (N_overlap.lt.0) then
               N_overlap=0.0
            endif
         else
            x_overlap=cm(14)+cm(15)-hsl
            x_maxoverlap=cm(15)-cm(16)
            if (x_overlap.lt.0) then
               N_overlap=0.0
            elseif ((x_overlap.ge.0) .and. 
     .           (x_overlap.le.x_maxoverlap)) then
               N_overlap=x_overlap/x_maxoverlap
            elseif (x_overlap.gt.x_maxoverlap) then
               N_overlap=1.0
            endif
         endif
c    
c     calculate N_bound at previous time step
         N_bound=0.0
         do 100 i=1,21
           N_bound=N_bound+x_cb0(i)
100      continue
c     added April 2016 by XZ
         if (N_bound.gt.N_total0) then
            N_bound=N_total0
         endif
c
c     calculate N_total at the current time step
         N_total=N_total0+dt1*(aCa*(N_overlap-N_total0)
     .            -cm(12)*(N_total0-N_bound)
     .            +cm(28)*N_bound
     .            -cm(29)*(N_overlap-N_total0))
c
c     prevent negative N_total April 2016 by XZ
         if (N_total.lt.0) then
            N_total=0.0
         endif 
c
c     calculate D at the current time step
         D=N_total-N_bound

c     assign attachment and detachment rates ------ parameters from Ken's figure 6 and table S4
c     l(i) in nm; k_a(i) and k_d(i) are converted to cm.
         do 1000 i=1,21
            l(i)=-11+i

c     move l(i) to l_mv(i);hsv(13) is in cm;
c     cm(13) is compliance factor
         l_mv(i)=l(i)-hsv(13)*1e7*cm(13)
c
            ka1(i)=-cm(18)*l(i)*l(i)*1e-14
            ka2=cm(22)*cm(9)*cm(10)
            ka3(i)=ka1(i)/ka2
            ka4(i)=exp(ka3(i))
c     two rate functions for attachment April 2016, XZ
            if (l(i).gt.0) then            
               k_a(i)=cm(21)*ka4(i)
            else
               k_a(i)=cm(32)*ka4(i)     
            endif         
c
c    limit k_a to 5000 (may need to check the unit of maximum k_a)
            if (k_a(i).gt.5000) then
               k_a(i)=5000
            endif 
c
            if (l(i).gt.0) then
               k_d(i)=cm(23)+cm(24)*l(i)**4
            else
               k_d(i)=cm(25)+cm(26)*l(i)**4
            endif
c    limit k_d to 5000
            if (k_d(i).gt.5000) then
               k_d(i)=5000
            endif 

c    4th-order Runge-Kutta method for x_cb at current time step 
            call derivsA(p_time,x_cb0(i),k_a(i),D,k_d(i),dxdt1(i))
            Ak1(i)=dt1*dxdt1(i)
            f_dxdt1(i)=x_cb0(i)+0.5*Ak1(i)
        
            call derivsA(time1,f_dxdt1(i),k_a(i),D,k_d(i),dxdt2(i))
            Ak2(i)=dt1*dxdt2(i)
            f_dxdt2(i)=x_cb0(i)+0.5*Ak2(i)

            call derivsA(time1,f_dxdt2(i),k_a(i),D,k_d(i),dxdt3(i))
            Ak3(i)=dt1*dxdt3(i)
            f_dxdt3(i)=x_cb0(i)+Ak3(i)

            call derivsA(time2,f_dxdt3(i),k_a(i),D,k_d(i),dxdt4(i))
            Ak4(i)=dt1*dxdt4(i)

            x_cb(i)=x_cb0(i)+Ak1(i)/6+Ak2(i)/3+Ak3(i)/3+Ak4(i)/6
c
1000     continue 
c 
c      move cross-bridge distribution (linear interpolation)
         call pwl_interp_1d(21,l,x_cb,21,l_mv,x_cb_new)     

c      save solutions to A(i), and calculate each T_a1 
         T_a1=0.0
         do 2000 i=1,21
            A(i)=x_cb_new(i)
            T_a1=T_a1+cm(18)*A(i)*rho_c*(l(i)+cm(27))*1e-7
2000     continue 
c      zero force when reaching l0, Sep 2016 by XZ
       else
         T_a1=0.0
       endif
c 
      else
c 
         ca=0.1
c
         T_a1=0.0
      endif
c
c     save T_a1 in material system April 2016 by XZ
      hsv(10)=T_a1
c
c***********************************************************************************
c     active contraction                                                           *
c      if (time.ge.t_act) then                                                     *
c     added the if statement on the next line (6/24/13)                            *
c        if (e1.gt. 0.5*((lo/lr)**2.0-1.0)) then                                   *
c        saclen=lr*sqrt(2.0*e1+1.0)                                                *
c         ECa_50=Ca_max/sqrt(exp(B*(saclen-lo))-1.0)                               *
c         tr=(m*saclen+bb)                                                         *
c                                                                                  *
c	change by LC Lee to include time dependency                                *
c	omega=3.14159*(0.25+tr)/tr                                                 *
c       corrected the typo of tr by Xiaoyan Zhang                                  *
c	    if(abs(time - t_act) < t0) then                                        *
c	 	omega = 3.14159*abs(time - t_act)/t0                               *
c	   elseif(abs(time - t_act) < t0+tr .and. abs(time - t_act) >= t0) then    *
c		omega = 3.14159*(abs(time - t_act) - t0 + tr)/tr                   *
c	    else                                                                   *
c		omega = 0                                                          *
c	    end if                                                                 *
c                                                                                  *
c                                                                                  *
c                                                                                  *
c        Ct=0.5*(1.0-cos(omega))                                                   *
c         T_a1=T_max*(Ca_o*Ca_o)*Ct/(Ca_o*Ca_o + ECa_50*ECa_50)                    *
c       added the next 3 lines (6/24/13)                                           *
c        else                                                                      *
c          T_a1=0.0                                                                *
c        endif                                                                     *
c      else                                                                        *
c          T_a1=0.0                                                                *
c      endif                                                                       *
c***********************************************************************************

c     partial derivative of W with repect to e1 to e6
c     convert to wrt to C (dW/dC = 0.5*dW/dE)
c     where c is the right cauchy-green tensor
c     material coordinate system
      sp1=0.5*dwdq*(tb2*e1)
      sp2=0.5*dwdq*(tb3*e2)
      sp3=0.5*dwdq*(tb3*e3)
      sp4=0.5*dwdq*tb4*e4
      sp5=0.5*dwdq*tb3*e5
      sp6=0.5*dwdq*tb4*e6
c
c     transform dWdC back to gobal axes
c     a=dWdCm*R
      a11   =sp1*s11+sp4*s21+sp6*s31
      a12   =sp1*s12+sp4*s22+sp6*s32
      a13   =sp1*s13+sp4*s23+sp6*s33
      a21   =sp4*s11+sp2*s21+sp5*s31
      a22   =sp4*s12+sp2*s22+sp5*s32
      a23   =sp4*s13+sp2*s23+sp5*s33
      a31   =sp6*s11+sp5*s21+sp3*s31
      a32   =sp6*s12+sp5*s22+sp3*s32
      a33   =sp6*s13+sp5*s23+sp3*s33
c     dWdCg=R'*a
      sp1=s11*a11   +s21*a21   +s31*a31
      sp2=s12*a12   +s22*a22   +s32*a32
      sp3=s13*a13   +s23*a23   +s33*a33
      sp4=s11*a12   +s21*a22   +s31*a32
      sp5=s12*a13   +s22*a23   +s32*a33
      sp6=s11*a13   +s21*a23   +s31*a33
c
c     active stresses
c     40% of fiber stress is assigned to cross-fiber
      T_a2=cm(30)*T_a1
      T_a3=cm(31)*T_a1
      T_a4=0.0
      T_a5=0.0
      T_a6=0.0
c    
c     transform active stress back to global axes
      a11   =T_a1*s11+T_a4*s21+T_a6*s31
      a12   =T_a1*s12+T_a4*s22+T_a6*s32
      a13   =T_a1*s13+T_a4*s23+T_a6*s33
      a21   =T_a4*s11+T_a2*s21+T_a5*s31
      a22   =T_a4*s12+T_a2*s22+T_a5*s32
      a23   =T_a4*s13+T_a2*s23+T_a5*s33
      a31   =T_a6*s11+T_a5*s21+T_a3*s31
      a32   =T_a6*s12+T_a5*s22+T_a3*s32
      a33   =T_a6*s13+T_a5*s23+T_a3*s33
c     T_a=R'*a
      T_a1=s11*a11   +s21*a21   +s31*a31
      T_a2=s12*a12   +s22*a22   +s32*a32
      T_a3=s13*a13   +s23*a23   +s33*a33
      T_a4=s11*a12   +s21*a22   +s31*a32
      T_a5=s12*a13   +s22*a23   +s32*a33
      T_a6=s11*a13   +s21*a23   +s31*a33
c
c*********************************************************************
c     mulitplicative split
c     right cauchy deformation tensor
      c1=F1*F1+F4*F4+F7*F7
      c2=F2*F2+F5*F5+F8*F8
      c3=F3*F3+F6*F6+F9*F9
      c4=F1*F2+F4*F5+F7*F8
      c5=F2*F3+F5*F6+F8*F9
      c6=F1*F3+F4*F6+F7*F9
c
c     determinant of C
      dc  =df*df
c
c     inverse of right cauchy deformation tensor
      cinv1=(c2*c3-c5*c5)/dc
      cinv2=(c1*c3-c6*c6)/dc
      cinv3=(c1*c2-c4*c4)/dc
      cinv4=(c6*c5-c3*c4)/dc
      cinv5=(c6*c4-c5*c1)/dc
      cinv6=(c4*c5-c6*c2)/dc
c
c     trace of dWdC and C 
      trspc=sp1*c1+sp2*c2+sp3*c3+2.0*sp6*c6+
     &       2.0*sp4*c4+2.0*sp5*c5
c
c     calculating the deviator
      dev1=sp1-trspc*cinv1/3.0
      dev2=sp2-trspc*cinv2/3.0
      dev3=sp3-trspc*cinv3/3.0
      dev4=sp4-trspc*cinv4/3.0
      dev5=sp5-trspc*cinv5/3.0
      dev6=sp6-trspc*cinv6/3.0
c
c     calculating the incompressibility
c     from bulk modulus and pressure relations
      p=bp*(vrat-1.0)
c
c     2nd piola kirchoff stress of diastole and systole
      sp1=p*df*cinv1+2.0*df23*dev1+T_a1
      sp2=p*df*cinv2+2.0*df23*dev2+T_a2
      sp3=p*df*cinv3+2.0*df23*dev3+T_a3
      sp4=p*df*cinv4+2.0*df23*dev4+T_a4
      sp5=p*df*cinv5+2.0*df23*dev5+T_a5
      sp6=p*df*cinv6+2.0*df23*dev6+T_a6
c
c     convert passive stress back to material system
c     April 2016 by XZ
      T_passive1=p*df*cinv1+2.0*df23*dev1
      T_passive2=p*df*cinv2+2.0*df23*dev2
      T_passive3=p*df*cinv3+2.0*df23*dev3
      T_passive4=p*df*cinv4+2.0*df23*dev4
      T_passive5=p*df*cinv5+2.0*df23*dev5
      T_passive6=p*df*cinv6+2.0*df23*dev6
c
      a11   =T_passive1*s11+T_passive4*s12+T_passive6*s13
      a12   =T_passive1*s21+T_passive4*s22+T_passive6*s23
      a13   =T_passive1*s31+T_passive4*s32+T_passive6*s33
      a21   =T_passive4*s11+T_passive2*s12+T_passive5*s13
      a22   =T_passive4*s21+T_passive2*s22+T_passive5*s23
      a23   =T_passive4*s31+T_passive2*s32+T_passive5*s33
      a31   =T_passive6*s11+T_passive5*s12+T_passive3*s13
      a32   =T_passive6*s21+T_passive5*s22+T_passive3*s23
      a33   =T_passive6*s31+T_passive5*s32+T_passive3*s33
c
      T_passive1=    s11*a11   +s12*a21   +s13*a31
      T_passive2=    s21*a12   +s22*a22   +s23*a32
c      T_passive3=    s31*a13   +s32*a23   +s33*a33
c      T_passive4=    s11*a12   +s12*a22   +s13*a32
c      T_passive5=    s21*a13   +s22*a23   +s23*a33
c      T_passive6=    s11*a13   +s12*a23   +s13*a33
c  
c     save passive stress in material system
c      hsv(37)=T_passive1
c    
c     cauchy stress
c
      rx=1./df
c
      q1=F1*sp1+F2*sp4+F3*sp6
      q2=F2*sp2+F1*sp4+F3*sp5
      q3=F3*sp3+F1*sp6+F2*sp5
      q4=F4*sp1+F5*sp4+F6*sp6
      q5=F5*sp2+F4*sp4+F6*sp5
      q6=F6*sp3+F4*sp6+F5*sp5
      q7=F7*sp1+F8*sp4+F9*sp6
      q8=F8*sp2+F7*sp4+F9*sp5
      q9=F9*sp3+F7*sp6+F8*sp5
      ssig1=rx*(F1*q1+F2*q2+F3*q3)
      ssig2=rx*(F4*q4+F5*q5+F6*q6)
      ssig3=rx*(F7*q7+F8*q8+F9*q9)
      ssig4=rx*(F4*q1+F5*q2+F6*q3)
      ssig5=rx*(F7*q4+F8*q5+F9*q6)
      ssig6=rx*(F7*q1+F8*q2+F9*q3)
c
c     updating the stress
      sig(1)=ssig1
      sig(2)=ssig2
      sig(3)=ssig3
      sig(4)=ssig4
      sig(5)=ssig5
      sig(6)=ssig6
c
c*********************************************************************
c     stress in fiber direction
      a11   =ssig1*s11+ssig4*s12+ssig6*s13
      a12   =ssig1*s21+ssig4*s22+ssig6*s23
c      a13   =ssig1*s31+ssig4*s32+ssig6*s33
      a21   =ssig4*s11+ssig2*s12+ssig5*s13
      a22   =ssig4*s21+ssig2*s22+ssig5*s23
c      a23   =ssig4*s31+ssig2*s32+ssig5*s33
      a31   =ssig6*s11+ssig5*s12+ssig3*s13
      a32   =ssig6*s21+ssig5*s22+ssig3*s23
c      a33   =ssig6*s31+ssig5*s32+ssig3*s33
      sigf1=    s11*a11   +s12*a21   +s13*a31
      sigf2=    s21*a12   +s22*a22   +s23*a32
c      sigf3=    s31*a13   +s32*a23   +s33*a33
c      sigf4=    s11*a12   +s12*a22   +s13*a32
c      sigf5=    s21*a13   +s22*a23   +s23*a33
c      sigf6=    s11*a13   +s12*a23   +s13*a33
c
c     saving fiber and cross fiber stress and strain in hsv
c
      hsv(7)=sigf1
c      hsv(8)=sigf2
      hsv(8)=N_overlap
      hsv(9)=e1
c      hsv(10)=e2
c     save T_a1 in global system
c      hsv(10)=T_a1
c
c     saving current x_cb(i,1) in hsv by XZ
      hsv(14)=A(1)
      hsv(15)=A(2)
      hsv(16)=A(3)
      hsv(17)=A(4)
      hsv(18)=A(5)
      hsv(19)=A(6)
      hsv(20)=A(7)
      hsv(21)=A(8)
      hsv(22)=A(9)
      hsv(23)=A(10)
      hsv(24)=A(11)
      hsv(25)=A(12)
      hsv(26)=A(13)
      hsv(27)=A(14)
      hsv(28)=A(15)
      hsv(29)=A(16)
      hsv(30)=A(17)
      hsv(31)=A(18)
      hsv(32)=A(19)
      hsv(33)=A(20)
      hsv(34)=A(21)
      hsv(35)=N_total
c
      hsv(36)=Ca
c
c
c ************************************************************
c     comments by Daniel Einstein
c --------------------- SOUND SPEED --------------------------
C     NOTE THE SOUND SPEED IS CRITICAL FOR THE TIME STEP CALCULATION
C     THOUGH YOU DO NOT NEED THE TANGENT MODULS FOR ITERATIONS, THE
C     MAXIMUM OF THE 9X9 TANGENT IS THE SOUND SPEED. OFTEN A
C     CONSERVATIVE ESTIMATE WILL DO. NOTE AS THIS IS THE NONLINEAR WORLD
C     THE SOUND SPEED IS ALSO A FUNCTION OF DEFORMATION. THE BULK
C     MODULUS IS A MINIMUM. HERE'S MY BEST GUESS. THERE WILL BE A
C     PASSIVE AND AN ACTIVE COMPONENT
c
c      Eb1=bq + (b2*dwdq)
c ***************************************************************
c     sound speed to pass into urmathn.f for time step determination
      Eb1=scb
c
c     
      return
      end     
c--------------------------------------------------------------------
c     subroutine derivs calculate the derivative of a function
c     sep 14
c--------------------------------------------------------------------
      subroutine derivsA(var_x,var_y,const1,const2,const3,
     . f_dydx)
      real*8 var_x,var_y,f_dydx
      real*8 const1,const2,const3
      f_dydx=const1*const2-const3*var_y
      return
      end
c--------------------------------------------------------------------
c     subroutine pwl_interp_1d: piecewise linear interpolant 
c     sep 14
c--------------------------------------------------------------------
      subroutine pwl_interp_1d(nd,xd,yd,ni,xi,yi)
      integer nd,ni,i,k
      real*8 xd(nd),yd(nd),xi(ni),yi(ni)
      real*8 ff
      do 20 i=1,ni
        yi(i)=0.0
20    continue 
      if (nd.eq.1) then
        do 21 i=1,ni
          yi(i)=yd(1)
21      continue      
        return
      endif
      do 22 i=1,ni
        if (xi(i).lt.xd(1).or.xi(i).gt.xd(nd)) then
           yi(i)=0.0
        else
          do 23 k=2,nd
            if (xd(k-1).le.xi(i).and.xi(i).le.xd(k)) then
              ff=(xi(i)-xd(k-1))/(xd(k)-xd(k-1))
              yi(i)=(1.0-ff)*yd(k-1)+ff*yd(k)
            endif
23        continue
        endif
22    continue 
      return
      end 

c--------------------------------------------------------------------
      subroutine umat43(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,cma,qmat,elsiz,idele)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      character*5 etype
      logical failel
c
      if (ncycle.eq.1) then
        call usermsg('mat43')
      endif
c
      return
      end
      subroutine umat44 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,cma,qmat,elsiz,idele)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      character*5 etype
      logical failel
c
c
      if (ncycle.eq.1) then
        call usermsg('mat44')
      endif
c
      return
      end
c$$$      subroutine umat42 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
c$$$     1 temper,failel,crv,cma,temper)
c$$$      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
c$$$      character*5 etype
c$$$      logical failel
c$$$c
c$$$      return
c$$$      end

c    Commented out umat46, replacing with three-state code


c    subroutine umat46?
      subroutine umat46 (cm,eps,sig,epsp,hsv,dt1,capa,
     . etype,time,temp,failel,crv,xM,xS,xN,vrat,Eb1).
c*******************************************************************
c     cardiac mechanics using Fung's exponential function          *
c                                                                  *
c     follows Guccione et al, "Mechanism Underlying Mechanical     *
c     Dysfunction in Border Zone of LV Aneurysm", 2001, Society of *
c     Thoracic Surgeons.                                           *
c                                                                  *
c     follows Gucione et al., Finite Element Stress Analysis of    *
c     LV mechanics in the Beating Dog Heart, 1995, J. Biomech.     *
c                                                                  *
c     created by Kay Sun                                           *
c     02/2007                                                      *
c                                                                  *
c*******************************************************************

c -------------------------------------------------------------- 
      include 'iounits.inc'
c     common/bk06/idmmy,iaddp,ifil,maxsiz,ncycle,ctime(2,30)
      character*(*) etype
      dimension  cm(*),sig(*),hsv(*),xM(*),xS(*),xN(*)
c
      real sp1,sp2,sp3,sp4,sp5,sp6
      real lr,lo

c    *--------------------------------------------------------------*
c     terms for Dr. Campbell's three-state contraction model
c     Assuming dt1 passed in is the timestep

c     terms for N_overlap
      real x_overlap, x_maxoverlap, N_overlap
  
c     terms for the number in states D1 and D2, as well as N_bound
c     N_a_active is the number of actin active at any given step. 
      real N_D1, N_D2, N_bound, N_a_active, N_bound0, N_a_active0

c     terms for rates
      real a_on_rate, a_off_rate, sum_of_rates
  
c     arrays to hold the activation and de-activation rates.
c     kA_D2 is detachment, kD2_A is attachment.
c     bin_position(41) is the bin displacement from -10 to 10 nm
      real*8 kA_D2(41), kD2_A(41), bin_position(41), moved_bins(41)
  
c     x_cb0(43) is the population in each of the 41 bins, 
c     and in the D1 and D2 states (1-41 are bins, 42 and 43 are D1 
c     D2 respectively). x_cb0 are initial for RK method, 
c     which returns the updated populations in x_cb
      real*8 x_cb0(43), x_cb(43), x_cb_interp(41)
  
c     M is the matrix that contains the rates to calculate dy/dt, where 
c     y is the array of populations. Ca is calcium concentration, T_a1
c     is the force generated by the given activated half-sarcomere
      real*8 M(43,43)
      real*8 Ca
      real*8 T_a1, T_a10, delta_D2

c     p is set to one, used for calculation of off rate of actin,
c     H is 0 when myosin regulatory light chains are not phosphorylated
c     and 1 when they are.
      integer p=1, H=1
c    *-----------------------------------------------------------------*
c     passive material parameters (1 to 4)
      ta=cm(1)
      b2=cm(2)
      b3=cm(3)
      b4=cm(4)

c     active material parameters and bulk modulus (20)
      lo=cm(5)
      lr=cm(6)
  
c     The following are to be set in the instruction file as material constants?
      k_mlcp = cm(7)
      k_force = cm(8)
c     Experimental Temp
      T_k = cm(9)
c     Boltzmann's Constant
      k_b = cm(10)
c     constants describing the length of thick and thin filaments
      l_thin = cm(11)
      l_thick = cm(12)
      l_bare = cm(13)
c     time of activation
      t_act = cm(14)
c     terms for calculating rates
      k1 = cm(15)
      k3 = cm(16)
      k40 = cm(17)
      k41 = cm(18)
      a_on = cm(19)
      a_off = cm(20)
      kD2_D1 = cm(21)
c     stiffness of cross-bridge spring
      k_cb = cm(22)
c     length of powerstroke
      x_ps = cm(23)
c     How much of filament movement is transferred to cross-bridges
      compliance = cm(24)
c     density of cross-bridges in a half sarcomere
      cb_density = cm(25)
c     Cooperative activation of binding sites
      k_thin_coop = cm(26)
  
c     Other cm's from XZ code
c     bulk modulus
      bp=cm(27)
  
c     added two cm for another two stress components, March, 2016 by XZ
      s_2=cm(28)
      s_3=cm(29)
  
c     factor scaling calcium transient
      fCa=cm(30)
      fCa_2=cm(31)
      l_falloff = cm(32)
  
c     hsv stuff will go here, initialize things
c    ------------------------------------------
c     deformation gradient from hsv 13 to 21
c     IHYPER=1, F=R'*F, in global frames
c
      F1=hsv(60)
      F4=hsv(61)
      F7=hsv(62)
      F2=hsv(63)
      F5=hsv(64)
      F8=hsv(65)
      F3=hsv(66)
      F6=hsv(67)
      F9=hsv(68)
c     save fiber directions in hsv
c     xM is fiber direction, xN is cross-fiber
c     Initializing at time = 0.0
      if (time.eq.0.0) then
        hsv(1)=xM(1)
        hsv(2)=xM(2)
        hsv(3)=xM(3)
c        hsv(5)=xS(1)
c        hsv(6)=xS(2)
c        hsv(7)=xS(3)
        hsv(4)=xN(1)
        hsv(5)=xN(2)
        hsv(6)=xN(3)
c      save initial sarcomere length 11/18/13
        hsv(11)=lr
        hsv(12)=lr

c      Initial distribution of bin population
c      Set to zero
        do 900 i = 1, 41
         x_cb0(i) = 0.0
900     continue
c       Dr. Campbell's model sets initial proportions
c       in D1 to be 0.9, D2 to be 0.1
        x_cb0(42) = 0.9
        x_cb0(43) = 0.1 
c       x_cb0(1-41) are number of attached heads in bins 1-41
c       x_cb0(42) is the population in D1
c       x_cb0(43) is the population in D2
c       hsv(14-56) represent this x_cb0 array
c       The flow of the program is:
c             hsv -> x_cb0 -> x_cb via RK -> x_cb_interp via
c             interpolation -> hsv
        hsv(14) = x_cb0(1)
        hsv(15) = x_cb0(2)
        hsv(16) = x_cb0(3)
        hsv(17) = x_cb0(4)
        hsv(18) = x_cb0(5)
        hsv(19) = x_cb0(6)
        hsv(20) = x_cb0(7)
        hsv(21) = x_cb0(8)
        hsv(22) = x_cb0(9)
        hsv(23) = x_cb0(10)
        hsv(24) = x_cb0(11)
        hsv(25) = x_cb0(12)
        hsv(26) = x_cb0(13)
        hsv(27) = x_cb0(14)
        hsv(28) = x_cb0(15)
        hsv(29) = x_cb0(16)
        hsv(30) = x_cb0(17)
        hsv(31) = x_cb0(18)
        hsv(32) = x_cb0(19)
        hsv(33) = x_cb0(20)
        hsv(34) = x_cb0(21)
        hsv(35) = x_cb0(22)
        hsv(36) = x_cb0(23)
        hsv(37) = x_cb0(24)
        hsv(38) = x_cb0(25)
        hsv(39) = x_cb0(26)
        hsv(40) = x_cb0(27)
        hsv(41) = x_cb0(28)
        hsv(42) = x_cb0(29)
        hsv(43) = x_cb0(30)
        hsv(44) = x_cb0(31)
        hsv(45) = x_cb0(32)
        hsv(46) = x_cb0(33)
        hsv(47) = x_cb0(34)
        hsv(48) = x_cb0(35)
        hsv(49) = x_cb0(36)
        hsv(50) = x_cb0(37)
        hsv(51) = x_cb0(38)
        hsv(52) = x_cb0(39)
        hsv(53) = x_cb0(40)
        hsv(54) = x_cb0(41)
        hsv(55) = x_cb0(42)
        hsv(56) = x_cb0(43)

c       hsv(57) is N_bound, initialized to zero
        hsv(57) = 0.0

c       hsv(58) is N_a_active, initialized to zero
        hsv(58) = 0.0

        hsv(10) = 0.0

c       hsv(10) represents T_a1, saved at the end of timestep
       endif

c     cross product of fiber and cross fiber to get orthogonal    
c     coordinate transformation from material to global
      s11=hsv(1)
      s12=hsv(2)
      s13=hsv(3)
      s21=hsv(5)*hsv(3)-hsv(2)*hsv(6)
      s22=hsv(6)*hsv(1)-hsv(3)*hsv(4)
      s23=hsv(4)*hsv(2)-hsv(1)*hsv(5)
      s31=hsv(4)
      s32=hsv(5)
      s33=hsv(6)
c
c     preparing material parameters
      tb2=2.*b2
      tb3=2.*b3
      tb4=2.*b4
c
c     limit magnitude of exponential term to prevent overflow errors
c     not use in this case
c      qlim=0.1
c      aqmx=log(qlim*bp/tb2)
c
c     calculation of sound speed for time step determination
      tc=0.5*ta
      bmbill=0.6666*max(b2,b3)
      smbill=2.*max(b3,b4)
      scb=max(bmbill,bp) + 1.33333*smbill
c
c     Jacobian of F
      df  =F1*(F5*F9-F8*F6)-
     .        F2*(F4*F9-F7*F6)+
     .        F3*(F4*F8-F7*F5)
      df23=df**(-2.0/3.0)     
c
c     green-st.venant strain  e=(F^t*F-1)/2
c
      e1=.5*(df23*(F1*F1+F4*F4+F7*F7)-1.)
      e2=.5*(df23*(F2*F2+F5*F5+F8*F8)-1.)
      e3=.5*(df23*(F3*F3+F6*F6+F9*F9)-1.)
      e4=.5*(df23*(F1*F2+F4*F5+F7*F8))
      e5=.5*(df23*(F2*F3+F5*F6+F8*F9))
      e6=.5*(df23*(F1*F3+F4*F6+F7*F9))

c     calculate total green-st.venant strain by Xiaoyan Zhang (06/17/2013)
      te1=.5*(df*(F1*F1+F4*F4+F7*F7)-1.)
      te2=.5*(df*(F2*F2+F5*F5+F8*F8)-1.)
      te3=.5*(df*(F3*F3+F6*F6+F9*F9)-1.)
      te4=.5*(df*(F1*F2+F4*F5+F7*F8))
      te5=.5*(df*(F2*F3+F5*F6+F8*F9))
      te6=.5*(df*(F1*F3+F4*F6+F7*F9))
c
c     transform green's strain in global to material axes
      a11   =e1*s11+e4*s12+e6*s13
      a12   =e1*s21+e4*s22+e6*s23
      a13   =e1*s31+e4*s32+e6*s33
      a21   =e4*s11+e2*s12+e5*s13
      a22   =e4*s21+e2*s22+e5*s23
      a23   =e4*s31+e2*s32+e5*s33
      a31   =e6*s11+e5*s12+e3*s13
      a32   =e6*s21+e5*s22+e3*s23
      a33   =e6*s31+e5*s32+e3*s33
      e1=    s11*a11   +s12*a21   +s13*a31
      e2=    s21*a12   +s22*a22   +s23*a32
      e3=    s31*a13   +s32*a23   +s33*a33
      e4=    s11*a12   +s12*a22   +s13*a32
      e5=    s21*a13   +s22*a23   +s23*a33
      e6=    s11*a13   +s12*a23   +s13*a33

c     transformation of total green's strain added by Xiaoyan Zhang (06/17/2013)
      a11   =te1*s11+te4*s12+te6*s13
      a12   =te1*s21+te4*s22+te6*s23
      a13   =te1*s31+te4*s32+te6*s33
      a21   =te4*s11+te2*s12+te5*s13
      a22   =te4*s21+te2*s22+te5*s23
      a23   =te4*s31+te2*s32+te5*s33
      a31   =te6*s11+te5*s12+te3*s13
      a32   =te6*s21+te5*s22+te3*s23
      a33   =te6*s31+te5*s32+te3*s33
      te1=    s11*a11   +s12*a21   +s13*a31
      te2=    s21*a12   +s22*a22   +s23*a32
      te3=    s31*a13   +s32*a23   +s33*a33
      te4=    s11*a12   +s12*a22   +s13*a32
      te5=    s21*a13   +s22*a23   +s23*a33
      te6=    s11*a13   +s12*a23   +s13*a33
c
c     partial derivative of W with respect to q
c     in material coordinate system
      aq=b2* e1**2+
     1   b3*(e2**2+e3**2+2.*e5**2)+
     1  tb4*(e4**2+e6**2)
c      aq=min(aq,aqmx)
      dwdq=tc*exp(aq)

c    adding violating case by 2/11/16 by Xiaoyan Zhang 
      hsv(11)=hsv(12)
      if (e1.gt. 0.5*((lo/lr)**2.0-1.0)) then
         saclen=lr*sqrt(2.0*e1+1.0)
      else
         saclen=hsv(11)
      endif
c
c     saving sarcomere length 11/18/13 by Xiaoyan Zhang
c      hsv(11)=hsv(12)
      hsv(12)=saclen
c     calculate the changes in half sarcomere length 
      hsv(13)=0.5*(hsv(12)-hsv(11))
c      hsv(13)=hsv(12)-hsv(11)
c
c*******************************************************************************************
c     active contraction incorporating Dr. Campbell's threestate cellular contraction model*
c     modified by Kurtis Mann 2017                                                         *
c*******************************************************************************************
c
c     initialization: reset x_cb0 from previous timestep
        x_cb0(1) = hsv(14) 
        x_cb0(2) = hsv(15) 
        x_cb0(3) = hsv(16)
        x_cb0(4) = hsv(17)
        x_cb0(5) = hsv(18)
        x_cb0(6) = hsv(19)
        x_cb0(7) = hsv(20)
        x_cb0(8) = hsv(21)
        x_cb0(9) = hsv(22)
        x_cb0(10) = hsv(23)
        x_cb0(11) = hsv(24)
        x_cb0(12) = hsv(25)
        x_cb0(13) = hsv(26)
        x_cb0(14) = hsv(27)
        x_cb0(15) = hsv(28)
        x_cb0(16) = hsv(29)
        x_cb0(17) = hsv(30)
        x_cb0(18) = hsv(31)
        x_cb0(19) = hsv(32)
        x_cb0(20) = hsv(33)
        x_cb0(21) = hsv(34)
        x_cb0(22) = hsv(35)
        x_cb0(23) = hsv(36)
        x_cb0(24) = hsv(37)
        x_cb0(25) = hsv(38)
        x_cb0(26) = hsv(39)
        x_cb0(27) = hsv(40)
        x_cb0(28) = hsv(41)
        x_cb0(29) = hsv(42)
        x_cb0(30) = hsv(43)
        x_cb0(31) = hsv(44)
        x_cb0(32) = hsv(45)
        x_cb0(33) = hsv(46)
        x_cb0(34) = hsv(47)
        x_cb0(35) = hsv(48)
        x_cb0(36) = hsv(49)
        x_cb0(37) = hsv(50)
        x_cb0(38) = hsv(51)
        x_cb0(39) = hsv(52)
        x_cb0(40) = hsv(53)
        x_cb0(41) = hsv(54)
        x_cb0(42) = hsv(55)
        x_cb0(43) = hsv(56)

        N_bound0 = hsv(57)
        N_a_active0 = hsv(58)

c     Need the force from last time step to calculate rate of D1 -> D2        
        T_a10 = hsv(10)


c     Once the header is generated, use cm(") instead of for ex. k3 or T

      if (time.ge.t_act) then

c     make active force zero when reaching l0 Sep, 2016 by XZ
       if (e1.gt. 0.5*((lo/lr)**2.0-1.0)) then
         p_time=time-dt1 
         time1=p_time+0.5*dt1
         time2=p_time+dt1

c      Calculate Calcium
c      -----------------
c      This is using XZ's modified calcium curve. Check units. Ca in nM
       t_p=t_act+0.01
       if (time.ge.t_p) then
c      cm(33) and cm(34) are set by XZ. Need more information on this
         pCa=0.5*exp(-((time-t_p)*cm(33))**cm(34))
       else
c      time is greater than t_act, but not to t_p               
         pCa=(time-t_act)/0.02
       endif
       Ca=0.1+1000*sin(3.14*pCa)
      
c      Calculate N_overlap
c      -------------------
       hsl=0.5*hsv(12)
       if (hsl.le.cm(11)) then 
            N_overlap=1-cm(32)*(cm(11)-hsl)
            if (N_overlap.lt.0) then
               N_overlap=0.0
            endif
         else
       x_overlap = l_thick + l_thin -hsl
       x_maxoverlap = l_thick - l_bare
       if (x_overlap.le.0.0) then
               N_overlap = 0.0
       elseif (x_overlap.ge.0.0 .AND. x_overlap.le. x_maxoverlap) then
               N_overlap = x_overlap/x_maxoverlap
       elseif (x_overlap.ge.x_maxoverlap) then
              N_overlap = 1.0
       endif
      
c      Calculate the number of activated actin
c      ---------------------------------------
       a_on_rate = a_on*Ca*(N_overlap-N_a_active0)*
     & (1+(k_thin_coop*(N_a_active0/N_overlap)))
      
       a_off_rate = a_off*(N_a_active0 - N_bound0)*(1+k_thin_coop*(((N
     & _overlap-N_a_active0)/N_overlap)**p))
        
       N_a_active = N_a_active0 + (a_on_rate - a_off_rate)*dt1

c      Calculate Rate Constants kA_D2, kD2_A
c      For now, constant detachment rate
c      -------------------------------------
       do 10 i = 1,41
          bin_position(i) = ((i+1)*0.5-11)  
          kD2_A(i) = k3*exp((-1*k_cb*bin_position(i)**2)/(2*k_b*T))
     & *(N_a_active-N_bound0)
          kA_D2(i) = 10
10     continue
       kD1_D2 = (k1+H*k_mlcp)*(1+k_force*T_a10)

c     Create M matrix to be used in the Runge-Kutta Method
c     ----------------------------------------------------
       do 20 i = 1,43
          do 25 j = 1,43
             M(i,j) = 0.0
25        continue
20     continue

       do 30 i = 1,41
          M(i,43) = kD2_A(i)
          M(i,i)  = -1*kA_D2(i)
          M(43,i) = kA_D2(i)
30     continue
       M(42,42) = -1*kD1_D2
       M(43,43) = -1*(kD2_D1 + sum_of_rates)
       M(42,43) = kD2_D1
       M(43,42) = kD1_D2

c      Runge Kutta for updating the number in each state/bin
c      -----------------------------------------------------
       call RK(M, time, dt1, x_cb0, x_cb)
       N_D1 = yout(42)
       N_D2 = yout(43)
       
       do 40 i = 1,41
       N_bound = N_bound + x_cb(i)
40     continue


c      Interpolate the cross-bridges (same as XZ code)
c      -----------------------------------------------
c      hsv(13) is the change in hsl, cm(13) is compliance factor
       do 60 i = 1,41
          moved_bins(i) = bin_position(i) - hsv(13)*1e7*cm(13)
60     continue

       call pwl_interp_1d(41,bin_position,x_cb,41,moved_bins,x_cb_interp
     .)

c      Put ripped cross bridges back into D2
c      -------------------------------------
       do 61 i = 1,41
          delta_D2 = delta_D2 + (x_cb(i) - x_cb_interp(i))
61     continue
       x_cb(43) = x_cb(43) + delta_D2

c      Calculate the half-sarc force
c      -----------------------------
       do 50 i = 1,41
          T_a1 = T_a1+cb_density*k_cb*x_cb_interp(i)*(bin_po
     & sition(i)+x_ps)
50     continue
       
c      Save everything in the history variables for the next step
c      ---------------------------------------------------------
       else
        T_a1=0.0
       endif
       
      else
c     Prior to activation, low calcium and no force
       Ca=0.1
       T_a1=0.0
      endif

      hsv(10) = T_a1


c      Insert stuff from umat45 that comes after active contraction
c      ------------------------------------------------------------

c     partial derivative of W with repect to e1 to e6
c     convert to wrt to C (dW/dC = 0.5*dW/dE)
c     where c is the right cauchy-green tensor
c     material coordinate system
      sp1=0.5*dwdq*(tb2*e1)
      sp2=0.5*dwdq*(tb3*e2)
      sp3=0.5*dwdq*(tb3*e3)
      sp4=0.5*dwdq*tb4*e4
      sp5=0.5*dwdq*tb3*e5
      sp6=0.5*dwdq*tb4*e6
c
c     transform dWdC back to gobal axes
c     a=dWdCm*R
      a11   =sp1*s11+sp4*s21+sp6*s31
      a12   =sp1*s12+sp4*s22+sp6*s32
      a13   =sp1*s13+sp4*s23+sp6*s33
      a21   =sp4*s11+sp2*s21+sp5*s31
      a22   =sp4*s12+sp2*s22+sp5*s32
      a23   =sp4*s13+sp2*s23+sp5*s33
      a31   =sp6*s11+sp5*s21+sp3*s31
      a32   =sp6*s12+sp5*s22+sp3*s32
      a33   =sp6*s13+sp5*s23+sp3*s33
c     dWdCg=R'*a
      sp1=s11*a11   +s21*a21   +s31*a31
      sp2=s12*a12   +s22*a22   +s32*a32
      sp3=s13*a13   +s23*a23   +s33*a33
      sp4=s11*a12   +s21*a22   +s31*a32
      sp5=s12*a13   +s22*a23   +s32*a33
      sp6=s11*a13   +s21*a23   +s31*a33
c
c     active stresses
c     40% of fiber stress is assigned to cross-fiber
      T_a2=cm(30)*T_a1
      T_a3=cm(31)*T_a1
      T_a4=0.0
      T_a5=0.0
      T_a6=0.0
c    
c     transform active stress back to global axes
      a11   =T_a1*s11+T_a4*s21+T_a6*s31
      a12   =T_a1*s12+T_a4*s22+T_a6*s32
      a13   =T_a1*s13+T_a4*s23+T_a6*s33
      a21   =T_a4*s11+T_a2*s21+T_a5*s31
      a22   =T_a4*s12+T_a2*s22+T_a5*s32
      a23   =T_a4*s13+T_a2*s23+T_a5*s33
      a31   =T_a6*s11+T_a5*s21+T_a3*s31
      a32   =T_a6*s12+T_a5*s22+T_a3*s32
      a33   =T_a6*s13+T_a5*s23+T_a3*s33
c     T_a=R'*a
      T_a1=s11*a11   +s21*a21   +s31*a31
      T_a2=s12*a12   +s22*a22   +s32*a32
      T_a3=s13*a13   +s23*a23   +s33*a33
      T_a4=s11*a12   +s21*a22   +s31*a32
      T_a5=s12*a13   +s22*a23   +s32*a33
      T_a6=s11*a13   +s21*a23   +s31*a33
c
c*********************************************************************
c     mulitplicative split
c     right cauchy deformation tensor
      c1=F1*F1+F4*F4+F7*F7
      c2=F2*F2+F5*F5+F8*F8
      c3=F3*F3+F6*F6+F9*F9
      c4=F1*F2+F4*F5+F7*F8
      c5=F2*F3+F5*F6+F8*F9
      c6=F1*F3+F4*F6+F7*F9
c
c     determinant of C
      dc  =df*df
c
c     inverse of right cauchy deformation tensor
      cinv1=(c2*c3-c5*c5)/dc
      cinv2=(c1*c3-c6*c6)/dc
      cinv3=(c1*c2-c4*c4)/dc
      cinv4=(c6*c5-c3*c4)/dc
      cinv5=(c6*c4-c5*c1)/dc
      cinv6=(c4*c5-c6*c2)/dc
c
c     trace of dWdC and C 
      trspc=sp1*c1+sp2*c2+sp3*c3+2.0*sp6*c6+
     &       2.0*sp4*c4+2.0*sp5*c5
c
c     calculating the deviator
      dev1=sp1-trspc*cinv1/3.0
      dev2=sp2-trspc*cinv2/3.0
      dev3=sp3-trspc*cinv3/3.0
      dev4=sp4-trspc*cinv4/3.0
      dev5=sp5-trspc*cinv5/3.0
      dev6=sp6-trspc*cinv6/3.0
c
c     calculating the incompressibility
c     from bulk modulus and pressure relations
      p=bp*(vrat-1.0)
c
c     2nd piola kirchoff stress of diastole and systole
      sp1=p*df*cinv1+2.0*df23*dev1+T_a1
      sp2=p*df*cinv2+2.0*df23*dev2+T_a2
      sp3=p*df*cinv3+2.0*df23*dev3+T_a3
      sp4=p*df*cinv4+2.0*df23*dev4+T_a4
      sp5=p*df*cinv5+2.0*df23*dev5+T_a5
      sp6=p*df*cinv6+2.0*df23*dev6+T_a6
c
c     convert passive stress back to material system
c     April 2016 by XZ
      T_passive1=p*df*cinv1+2.0*df23*dev1
      T_passive2=p*df*cinv2+2.0*df23*dev2
      T_passive3=p*df*cinv3+2.0*df23*dev3
      T_passive4=p*df*cinv4+2.0*df23*dev4
      T_passive5=p*df*cinv5+2.0*df23*dev5
      T_passive6=p*df*cinv6+2.0*df23*dev6
c
      a11   =T_passive1*s11+T_passive4*s12+T_passive6*s13
      a12   =T_passive1*s21+T_passive4*s22+T_passive6*s23
      a13   =T_passive1*s31+T_passive4*s32+T_passive6*s33
      a21   =T_passive4*s11+T_passive2*s12+T_passive5*s13
      a22   =T_passive4*s21+T_passive2*s22+T_passive5*s23
      a23   =T_passive4*s31+T_passive2*s32+T_passive5*s33
      a31   =T_passive6*s11+T_passive5*s12+T_passive3*s13
      a32   =T_passive6*s21+T_passive5*s22+T_passive3*s23
      a33   =T_passive6*s31+T_passive5*s32+T_passive3*s33
c
      T_passive1=    s11*a11   +s12*a21   +s13*a31
      T_passive2=    s21*a12   +s22*a22   +s23*a32
c      T_passive3=    s31*a13   +s32*a23   +s33*a33
c      T_passive4=    s11*a12   +s12*a22   +s13*a32
c      T_passive5=    s21*a13   +s22*a23   +s23*a33
c      T_passive6=    s11*a13   +s12*a23   +s13*a33
c  
c     save passive stress in material system
c      hsv(37)=T_passive1
c    
c     cauchy stress
c
      rx=1./df
c
      q1=F1*sp1+F2*sp4+F3*sp6
      q2=F2*sp2+F1*sp4+F3*sp5
      q3=F3*sp3+F1*sp6+F2*sp5
      q4=F4*sp1+F5*sp4+F6*sp6
      q5=F5*sp2+F4*sp4+F6*sp5
      q6=F6*sp3+F4*sp6+F5*sp5
      q7=F7*sp1+F8*sp4+F9*sp6
      q8=F8*sp2+F7*sp4+F9*sp5
      q9=F9*sp3+F7*sp6+F8*sp5
      ssig1=rx*(F1*q1+F2*q2+F3*q3)
      ssig2=rx*(F4*q4+F5*q5+F6*q6)
      ssig3=rx*(F7*q7+F8*q8+F9*q9)
      ssig4=rx*(F4*q1+F5*q2+F6*q3)
      ssig5=rx*(F7*q4+F8*q5+F9*q6)
      ssig6=rx*(F7*q1+F8*q2+F9*q3)
c
c     updating the stress
      sig(1)=ssig1
      sig(2)=ssig2
      sig(3)=ssig3
      sig(4)=ssig4
      sig(5)=ssig5
      sig(6)=ssig6
c
c*********************************************************************
c     stress in fiber direction
      a11   =ssig1*s11+ssig4*s12+ssig6*s13
      a12   =ssig1*s21+ssig4*s22+ssig6*s23
c      a13   =ssig1*s31+ssig4*s32+ssig6*s33
      a21   =ssig4*s11+ssig2*s12+ssig5*s13
      a22   =ssig4*s21+ssig2*s22+ssig5*s23
c      a23   =ssig4*s31+ssig2*s32+ssig5*s33
      a31   =ssig6*s11+ssig5*s12+ssig3*s13
      a32   =ssig6*s21+ssig5*s22+ssig3*s23
c      a33   =ssig6*s31+ssig5*s32+ssig3*s33
      sigf1=    s11*a11   +s12*a21   +s13*a31
      sigf2=    s21*a12   +s22*a22   +s23*a32
c      sigf3=    s31*a13   +s32*a23   +s33*a33
c      sigf4=    s11*a12   +s12*a22   +s13*a32
c      sigf5=    s21*a13   +s22*a23   +s23*a33
c      sigf6=    s11*a13   +s12*a23   +s13*a33
c
c     saving fiber and cross fiber stress and strain in hsv
c
      hsv(7)=sigf1
c      hsv(8)=sigf2
      hsv(8)=N_overlap
      hsv(9)=e1
c      hsv(10)=e2
c     save T_a1 in global system
c      hsv(10)=T_a1
c

c     Save bin populations into history variables
c     These need to be x_cb_interp? And after interpolation, update D2?

       hsv(14) = x_cb_interp(1)
       hsv(15) = x_cb_interp(2)
       hsv(16) = x_cb_interp(3)
       hsv(17) = x_cb_interp(4)
       hsv(18) = x_cb_interp(5)
       hsv(19) = x_cb_interp(6)
       hsv(20) = x_cb_interp(7)
       hsv(21) = x_cb_interp(8)
       hsv(22) = x_cb_interp(9)
       hsv(23) = x_cb_interp(10)
       hsv(24) = x_cb_interp(11)
       hsv(25) = x_cb_interp(12)
       hsv(26) = x_cb_interp(13)
       hsv(27) = x_cb_interp(14)
       hsv(28) = x_cb_interp(15)
       hsv(29) = x_cb_interp(16)
       hsv(30) = x_cb_interp(17)
       hsv(31) = x_cb_interp(18)
       hsv(32) = x_cb_interp(19)
       hsv(33) = x_cb_interp(20)
       hsv(34) = x_cb_interp(21)
       hsv(35) = x_cb_interp(22)
       hsv(36) = x_cb_interp(23)
       hsv(37) = x_cb_interp(24)
       hsv(38) = x_cb_interp(25)
       hsv(39) = x_cb_interp(26)
       hsv(40) = x_cb_interp(27)
       hsv(41) = x_cb_interp(28)
       hsv(42) = x_cb_interp(29)
       hsv(43) = x_cb_interp(30)
       hsv(44) = x_cb_interp(31)
       hsv(45) = x_cb_interp(32)
       hsv(46) = x_cb_interp(33)
       hsv(47) = x_cb_interp(34)
       hsv(48) = x_cb_interp(35)
       hsv(49) = x_cb_interp(36)
       hsv(50) = x_cb_interp(37)
       hsv(51) = x_cb_interp(38)
       hsv(52) = x_cb_interp(39)
       hsv(53) = x_cb_interp(40)
       hsv(54) = x_cb_interp(41)
       hsv(55) = x_cb(42)
       hsv(56) = x_cb(43)

       hsv(57) = N_bound
       hsv(58) = N_a_active

c     Including hsv(59) = Ca to keep record of it through time?
      hsv(59) = Ca

      return
      end



c*****************************************************************************
c*****************************************************************************

c     -------------------------------------      
c     subroutine for the Runge Kutta method
c     -------------------------------------
      subroutine RK(M, t, h, init_conds, yout)
         real*8 M(43,43), t, h, init_conds(43,1), yout(43,1)
         real*8 y_values(43,1), dyt(43,1), yt(43,1), dym(43,1),hh,h6,th
         real*8 dydt(43,1)
         y_values = init_conds
         hh = h*0.5
         h6 = h/6
         th = t+hh
         
c        M*y_values (Step 1)        
         do 15 i = 1,42
         dydt(i,1) = M(i,i)*y_values(i,1) + M(i,43)*y_values(43,1)
15       continue
         do 25 i = 1,43  
         dydt(43,1) = dydt(43,1) + M(43,i)*y_values(i,1)
25       continue
       
         do 35 i = 1,43
         yt(i,1) = y_values(i,1) + hh*dydt(i,1)
35       continue
        
c        M*yt (Step 2)
         do 45 i = 1,42
         dyt(i,1) = M(i,i)*yt(i,1) + M(i,43)*yt(43,1)
45       continue
         do 55 i = 1,43
         dyt(43,1) = dyt(43,1) + M(43,i)*yt(i,1)
55       continue

         do 65 i = 1,43
         yt(i,1) = y_values(i,1) + hh*dyt(i,1)
65       continue
         
c        (Step 3)
         do 31 i = 1,42
         dym(i,1) = M(i,i)*yt(i,1) + M(i,43)*yt(43,1)
31       continue
         do 32 i = 1,43
         dym(43,1) = dym(43,1) + M(43,i)*yt(i,1)
32       continue
         
         do 75 i = 1,43
         yt(i,1) = y_values(i,1) + h*dym(i,1)
75       continue

         do 33 i = 1,43
         dym(i,1) = dym(i,1)+dyt(i,1)
33       continue

c        (Step 4)
c        Reset dyt to zero
         do 39 i = 1,43
		    dyt(i,1) = 0.0
39       continue
         do 36 i = 1,42
         dyt(i,1) = M(i,i)*yt(i,1) + M(i,43)*yt(43,1)
36       continue
         do 37 i = 1,43
         dyt(43,1) = dyt(43,1) + M(43,i)*yt(i,1)
37       continue

         do 85 i = 1,43
         yout(i,1) = y_values(i,1) + h6*(dydt(i,1)+dyt(i,1)+2*dym(i,1))
85       continue

         return 
         
      end
c     ----------------------------------------------------------
c     Matrix multiplication subroutine
c     ----------------------------------------------------------
      subroutine matrix_mult(A,m,n,B,p,q,C)
      integer i,j,k,m,n,p,q,counter1,counter2
      real*8 A(m,n), B(p,q), C(m,q)

c     Initialize resultant matrix to zeros, else it populates strangely
      do 400 counter1 = 1,m
        do 401 counter2 = 1,q
          C(counter1,counter2) = 0.0
401     continue
400   continue
                
      do 300 i = 1,m
        do 301 j = 1,q
          do 302 k = 1,n
            C(i,j) = C(i,j) + A(i,q)*B(q,j)
302       continue
301     continue
300   continue
      
      return
      end
      
c     ----------------------------------------------------------
c     subroutine for pwl_interp_1d: piecewise linear interpolant
c     from umat45
c     ----------------------------------------------------------

      subroutine pwl_interp_1d(nd,xd,yd,ni,xi,yi)
      integer nd,ni,i,k
      real*8 xd(nd),yd(nd),xi(ni),yi(ni)
      real*8 ff
      do 26 i=1,ni
        yi(i)=0.0
26    continue 
      if (nd.eq.1) then
        do 27 i=1,ni
          yi(i)=yd(1)
27      continue      
        return
      endif
      do 29 i=1,ni
        if (xi(i).lt.xd(1).or.xi(i).gt.xd(nd)) then
           yi(i)=0.0
        else
          do 28 k=2,nd
            if (xd(k-1).le.xi(i).and.xi(i).le.xd(k)) then
              ff=(xi(i)-xd(k-1))/(xd(k)-xd(k-1))
              yi(i)=(1.0-ff)*yd(k-1)+ff*yd(k)
            endif
28        continue
        endif
29    continue 
      return
      end 

c******************************************************************
c******************************************************************
c******************************************************************

      
c      subroutine umat46 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
c     1 temper,failel,crv,cma,thhsvi,nthhsv,qmat,elsiz,idele)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c      include 'nlqparm'
c      include 'bk06.inc'
c      include 'iounits.inc'
c      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),thhsvi(*),
c     . cma(*),qmat(3,3)
c      character*5 etype
c      logical failel
c
c     Isotropic elastic material with access to history variables
c     of a corresponding thermal user material. See also
c     subroutine umat41 for information about non-thermal
c     specific information. Access to thermal history variables
c     requieres IHVE.eq.1 on *MAT_THERMAL_USER_DEFINED and
c     is only possible for a selection of shell and solid elements.
c     More information including a list of supported element
c     types is given in the comments to subroutine thusrmat.
c
c     Read only variables
c
c     cm(1)=first material constant, here young's modulus
c     cm(2)=second material constant, here poisson's ratio
c     cm(3)=.eq.0 do nothing
c           .eq.1 print mechanical and thermal history variables.
c
c     nthhsv=Number of thermal history variables.
c     thhsvi=If IHVE.eq.1 on *MAT_THERMAL_USER_DEFINED then
c            thhsvi(i) i=1...nthhsv, contains history
c            variable data from thermal user material routine
c            interpolated to current integration point using
c            shape functions for the corresponding thermal element.
c
c     compute shear modulus, g
c
c      if (ncycle.eq.1) then
c        call usermsg('mat46')
c      endif
c
c      g2 =cm(1)/(1.+cm(2))
c      g  =.5*g2
c
c      if (etype.eq.'solid'.or.etype.eq.'shl_t'.or.
c     1     etype.eq.'sld2d'.or.etype.eq.'tshel') then
c        davg=(-eps(1)-eps(2)-eps(3))/3.
c        p=-davg*cm(1)/(1.-2.*cm(2))
c        sig(1)=sig(1)+p+g2*(eps(1)+davg)
c        sig(2)=sig(2)+p+g2*(eps(2)+davg)
c        sig(3)=sig(3)+p+g2*(eps(3)+davg)
c        sig(4)=sig(4)+g*eps(4)
c        sig(5)=sig(5)+g*eps(5)
c        sig(6)=sig(6)+g*eps(6)
c
c      else if (etype.eq.'shell') then
c        gc    =capa*g
c        q1    =cm(1)*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
c        q3    =1./(q1+g2)
c        eps(3)=-q1*(eps(1)+eps(2))*q3
c        davg  =(-eps(1)-eps(2)-eps(3))/3.
c        p     =-davg*cm(1)/(1.-2.*cm(2))
c        sig(1)=sig(1)+p+g2*(eps(1)+davg)
c        sig(2)=sig(2)+p+g2*(eps(2)+davg)
c        sig(3)=0.0
c        sig(4)=sig(4)+g *eps(4)
c        sig(5)=sig(5)+gc*eps(5)
c        sig(6)=sig(6)+gc*eps(6)
cc
c      else
c       write(iotty,10) etype
c       write(iohsp,10) etype
c       write(iomsg,10) etype
c       call adios(2)
c        cerdat(1)=etype
c        call lsmsg(3,MSG_SOL+1152,ioall,ierdat,rerdat,cerdat,0)
c      endif
c
c     print history variables?
c
c      if ( nint(cm(3)).eq.1) then
c        write(iotty,*) " epsp ",epsp," thhsvi(1) ",thhsvi(1),
c     1   " thhsvi(5) ",thhsvi(5),nthhsv
c      endif
c
c10   format(/
c    1 ' *** Error element type ',a,' can not be',
c    2 '           run with the current user material model #47.')
c      return
c      end
      subroutine umat47 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,cma,qmat,elsiz,idele)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      logical failel
      character*5 etype
c
      if (ncycle.eq.1) then
        call usermsg('mat47')
      endif
c
      return
      end
      subroutine umat48 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,cma,qmat,elsiz,idele)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      logical failel
      character*5 etype
c
      if (ncycle.eq.1) then
        call usermsg('mat48')
      endif
c
      return
      end
      subroutine umat49 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,cma,qmat,elsiz,idele)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      logical failel
      character*5 etype
c
      if (ncycle.eq.1) then
        call usermsg('mat49')
      endif
c
      return
      end
      subroutine umat50 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,cma,qmat,elsiz,idele)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      logical failel
      character*5 etype
c
      if (ncycle.eq.1) then
        call usermsg('mat50')
      endif
c
      return
      end
      subroutine usermsg(msg)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     print warning message
c
      include 'iounits.inc'
      character*(*) msg
c
c     write(iotty,2100) msg
c     write(iohsp,2100) msg
c     write(iomsg,2100) msg
      cerdat(1)=msg
      call lsmsg(1,MSG_SOL+1147,ioall,ierdat,rerdat,cerdat,0)
c2100 format(/
c    1 ' *** Warning'
c    2 5x,'Default ',a,' is used',/,
c    3 5x,'Please check user material library')
c
      return
      end
      subroutine umat41v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,cma,qmat,elsizv,idelev)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*)
      dimension temps(*),crv(lq1,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      dimension idelev(*)
      logical failels(*)
      character*5 etype
c
      dimension sig(6),eps(6),hsv(NHISVAR)
c
      do i=lft,llt
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
c
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
c
        call umat41(cm,eps,sig,epsps(i),hsv,dt1siz(i),capa,etype,tt,
     1   temps(i),failels(i),crv,cma,0,elsizv(i),idelev(i))
c
        sig1(i)=sig(1)
        sig2(i)=sig(2)
        sig3(i)=sig(3)
        sig4(i)=sig(4)
        sig5(i)=sig(5)
        sig6(i)=sig(6)
        d3(i)=eps(3)
      enddo
c
      return
      end
      subroutine umat45v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,cma,qmat,elsizv,idelev)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
      common/bk36loc/index
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*)
      dimension temps(*),crv(lq1,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      dimension idelev(*)
      logical failels(*)
      character*5 etype
c
      dimension sig(6),eps(6),hsv(NHISVAR)
c
c$omp threadprivate (/bk36loc/)
      do i=lft,llt
        index=i
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
c
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
c
        hsv(1)=hsvs(i,1)
        hsv(2)=hsvs(i,2)
        hsv(3)=hsvs(i,3)
        hsv(4)=hsvs(i,4)
        hsv(5)=hsvs(i,5)
        hsv(6)=hsvs(i,6)
        hsv(7)=hsvs(i,7)
        hsv(8)=hsvs(i,8)
        hsv(9)=hsvs(i,9)
c
        call umat45(cm,eps,sig,epsps(i),hsv,dt1siz(i),capa,etype,tt,
     1   temps(i),failels(i),crv,cma,0,elsizv(i),idelev(i))
c
        sig1(i)=sig(1)
        sig2(i)=sig(2)
        sig3(i)=sig(3)
        sig4(i)=sig(4)
        sig5(i)=sig(5)
        sig6(i)=sig(6)
c
        hsvs(i,1)=hsv(1)
        hsvs(i,2)=hsv(2)
        hsvs(i,3)=hsv(3)
        hsvs(i,4)=hsv(4)
        hsvs(i,5)=hsv(5)
        hsvs(i,6)=hsv(6)
        hsvs(i,7)=hsv(7)
        hsvs(i,8)=hsv(8)
        hsvs(i,9)=hsv(9)
        d3(i)=eps(3)
      enddo
c
      return
      end
      subroutine umat44v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,eps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,cma,qmat,elsizv,idelev)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c***  isotropic plasticity with linear hardening
c
c***  updates only the deviatoric stress so that it can be used with
c     an equation of state
c
      parameter (third=1.0/3.0)
      include 'nlqparm'
c
      common/eosdloc/pc(nlq)
c
      dimension cm(*),d1(*),d2(*),d3(*),d4(*),d5(*),d6(*),
     & sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*),
     & eps(*),hsvs(nlqa,*),dt1siz(*),temps(*),crv(lq1,2,*),cma(*),
     & failels(*),qmat(nlq,3,3),elsizv(*),idelev(*)
c
      character*5 etype
      logical failels
c
c***  shear modulus, initial yield stress, hardening, and pressure cut-off
c$omp threadprivate (/eosdloc/)
      g   =cm(1)
      sy0 =cm(2)
      h   =cm(3)
      pcut=cm(4)
c
c***  plastic strain for failure
      epsfail=cm(5)
c
      ofac=1.0/(3.0*g+h)
      twog=2.0*g
c
      do i=lft,llt
c
c***    elastic deviatoric stress
        davg=third*(d1(i)+d2(i)+d3(i))
        savg=third*(sig1(i)+sig2(i)+sig3(i))
        sig1(i)=sig1(i)-savg+twog*(d1(i)-davg)
        sig2(i)=sig2(i)-savg+twog*(d2(i)-davg)
        sig3(i)=sig3(i)-savg+twog*(d3(i)-davg)
        sig4(i)=sig4(i)+g*d4(i)
        sig5(i)=sig5(i)+g*d5(i)
        sig6(i)=sig6(i)+g*d6(i)
c
c***    radial return
        aj2=sqrt(1.5*(sig1(i)**2+sig2(i)**2+sig3(i)**2)+
     &           3.0*(sig4(i)**2+sig5(i)**2+sig6(i)**2))
        sy=sy0+h*eps(i)
        eps(i)=eps(i)+ofac*max(0.0,aj2-sy)
        synew=sy0+h*eps(i)
        scale=synew/max(synew,aj2)
c
c***    scaling for radial return
        sig1(i)=scale*sig1(i)
        sig2(i)=scale*sig2(i)
        sig3(i)=scale*sig3(i)
        sig4(i)=scale*sig4(i)
        sig5(i)=scale*sig5(i)
        sig6(i)=scale*sig6(i)
c
c***    set pressure cut-off
        pc(i)=pcut
c
c***    failure due to plastic strain
        failels(i)=eps(i).gt.epsfail
c
      enddo
c
      return
      end
      subroutine umat46v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,cma,thhsv,nthhsv,qmat,elsizv,
     . idelev)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*),thhsv(nlq,*)
      dimension temps(*),crv(lq1,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      dimension idelev(*)
      logical failels(*)
      character*5 etype
c
c     For more information see comments in subroutine umat46.
c
c     nthhsv=Number of thermal history variables.
c     thhsv=If IHVE.eq.1 on *MAT_THERMAL_USER_DEFINED then
c           thhsvi(i,j) j=1...nthhsv, i=1...nlq, contains history
c           variable data from thermal user material routine
c           interpolated to current integration point using
c           shape functions for the corresponding thermal element.
c
      dimension sig(6),eps(6),hsv(NHISVAR),thhsvi(100)
c
      do i=lft,llt
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
c
        do j=1,nthhsv
          thhsvi(j)=thhsv(i,j)
        enddo
c
        call umat46(cm,eps,sig,epsps(i),hsv,dt1siz(i),capa,etype,tt,
     1   temps(i),failels(i),crv,cma,thhsvi,nthhsv,0,elsizv(i),
     2   idelev(i))
c
        do j=1,nthhsv
          thhsv(i,j)=thhsvi(j)
        enddo
c
        sig1(i)=sig(1)
        sig2(i)=sig(2)
        sig3(i)=sig(3)
        sig4(i)=sig(4)
        sig5(i)=sig(5)
        sig6(i)=sig(6)
        d3(i)=eps(3)
      enddo
c
      return
      end
      subroutine umat47v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,cma,qmat,elsizv,idelev)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*)
      dimension temps(*),crv(lq1,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      dimension idelev(*)
      logical failels(*)
      character*5 etype
c
      return
      end
c$$$      subroutine umat48v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
c$$$     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
c$$$     . etype,tt,temps,failels,nlqa,crv,cma)
c$$$      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
c$$$      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
c$$$      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*)
c$$$      dimension temps(*),crv(lq1,2,*),cma(*)
c$$$      logical failels(*)
c$$$      character*5 etype
c$$$c
c$$$      return
c$$$      end
c$$$      subroutine (cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
c$$$     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
c$$$     . etype,tt,temps,failels,nlqa,crv,cma)
c$$$      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
c$$$      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
c$$$      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*)
c$$$      dimension temps(*),crv(lq1,2,*),cma(*)
c$$$      logical failels(*)
c$$$      character*5 etype
c$$$c
c$$$      return
c$$$      end
      subroutine umat50v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,cma,qmat,elsizv,idelev)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*)
      dimension temps(*),crv(lq1,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      dimension idelev(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine loadsetud(time,lft,llt,crv,iduls,parm)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     Input (not modifiable)
c       time : analysis time
c       x    : coordinate of node or element center
c       d    : displacement of node or lement center
c       v    : velocity of node or lement center
c       temp : temperature of node or element center
c       crv  : value of LCURV at time=time
c       iduls: id of user_loading_set
c       parm : parameters from user_loading
c     Output (defined by user)
c       udl : user-defined load curve value
      include 'nlqparm'
      common/aux8loc/
     & x1(nlq),x2(nlq),x3(nlq),v1(nlq),
     & v2(nlq),v3(nlq),d1(nlq),d2(nlq),
     & d3(nlq),temp(nlq),udl(nlq),tmp(nlq,12)
c
      dimension parm(*)
c
c     sample
c     if (iduls.eq.100) then
c       do i=lft,llt
cc        to define a pressure loading proportional to sin of the angle
cc        between unit vector v,=(vx,vy,vz), and element normal
cc        take 42 x 13 as normal
c         x42=x2(i)-x4(i)
c         y42=y2(i)-y4(i) 
c         z42=z2(i)-z4(i)
c         x13=x3(i)-x1(i)
c         y13=y3(i)-y1(i) 
c         z13=z3(i)-z1(i)
c         xn = y42*z13-y13*z42
c         yn = z42*x13-z13*x42
c         zn = x42*y13-x13*y42
cc        v x normal
c         vnx= vy*zn-yn*vz
c         vny= vz*xn-zn*vx
c         vnz= vx*yn-xn*vy
c         vn = sqrt(vnx*vnx+vny*vny+vnz*vnz)
cc        abs(sin)=abs(v x n )/|n|
c         rn = 1./sqrt(xn*xn+yn*yn+zn*zn)
c         sint = vn*rn
c         fact = sint
c         udl(i)=crv*fact
c       enddo
c     else if (iduls.eq.200) then
c       do i=lft,llt
c         udl(i)=....
c       enddo
c       ....
c     endif
c$omp threadprivate (/aux8loc/)
      do i=lft,llt
       udl(i)= crv
      enddo
c
      return
      end
      subroutine loadud(fnod,dt1,time,ires,x,d,v,a,ixs,
     . numels,ixb,numelb,idrflg,tfail,isf,p,npc,fval,iob,iadd64,numelh,
     . ixh,nhex_del,nbeam_del,nshell_del,hexarray,hextim,bemarray,
     . bemtim,shlarray,shltim,parm)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     input arrays
c
c     fnod - global nodal forces
c
c     dt1 - current time step size
c     time - current problem time
c     ires - restart flag, ( 0=solution phase     )
c                          (-n=input n parameters )
c                          ( 2=restart            )
c                          ( 3=write data into dump file  )
c                          ( 4=read data from restart file)
c
c     When data is read, DUMMY arrays are passed in the call.
c     Data should be read into a local common block which is
c     written into the restart database.
c
c
c     d - displacements
c     v - velocities
c     a - acceleations
c
c     ixs - shell element connectivities (ixs(1,*)=part ID)
c                                        (ixs(2,*)=node 1)
c                                        (ixs(3,*)=node 2)
c                                        (ixs(4,*)=node 3)
c                                        (ixs(5,*)=node 4)
c
c     ixb - beam  element connectivities (ixb(1,*)=part ID)
c                                        (ixb(2,*)=node 1)
c                                        (ixb(3,*)=node 2)
c                                        (ixb(4,*)=orientation node)
c
c     ixh - shell element connectivities (ixh(1,*)=part ID)
c                                         or
c                                        (ixh(1,*)=0 implies element deleted)
c                                        (ixh(2,*)=node 1)
c                                        (ixh(3,*)=node 2)
c                                        (ixh(4,*)=node 3)
c                                        (ixh(5,*)=node 4)
c                                        (ixh(6,*)=node 5)
c                                        (ixh(7,*)=node 6)
c                                        (ixh(8,*)=node 7)
c                                        (ixh(9,*)=node 8)
c
c     numels - number of shell elements
c     numelb - number of beam  elements
c     numelh - number of solid elements
c     isf    - shell element failure flag (1=on)
c     tfail  - shell element failure time (eq.0:okay)
c                                         (ne.0:failure time)
c
c     idrflg - nonzero if dynamic relaxation phase
c     p      - load curve data pairs (abcissa,ordinate)
c     npc    - pointer into p.  (p(npc(lc)) points to the beginning
c              of load curve ID lc.  npoints=npc(lc+1)-npc(lc)=
c              number of points in the load curve.
c     fval   - fval(lc) is the value of load curve lc at t=time
c     iob    - i/o buffer
c
c
c
c     ELEMENT DELETION FROM THE USER LOADING SUBROUTINE
c
c     to delete elements the original file needs to have a
c     *DEFINE_ELEMENT_DEATH definition for one or more elements
c     of the type to be deleted.  The deletion time can be set
c     to a value that greatly exceeds the termination time for the
c     run.  With this keyword the necessary arrays are created that
c     allows elements of the same type to deleted during the run.
c
c     nhex_del if >0 element deletion option is active for solids
c     nbem_del if >0 element deletion option is active for beams
c     nshl_del if >0 element deletion option is active for shells
c     hexarray defines time to delete solid element,
c              the value should be >time
c     hextim   we check for solid element deletion when the current
c              time is greater to or equal to hextim
c     bemarray defines time to delete beam element,
c              the value should be >time
c     bemtim   we check for beam element deletion when the current
c              time is greater to or equal to bemim
c     shlarray defines time to delete shell element,
c              the value should be >time
c     shltim   we check for shell element deletion when the current
c              time is greater to or equal to shltim
c     parm     user loading parameter if ires<0
c
      include 'iounits.inc'
      include 'bigprb.inc'
      include 'txtline.inc'
c
      parameter (NPARM=1000)
c     common/usrldv/parm(NPARM)
c
      integer*8 iadd64
      real*8 x
      real*8 d
      dimension a(3,*),v(3,*),d(3,*),fnod(3,*),ixs(5,*),ixb(4,*),
     . x(3,*),tfail(*),p(*),npc(*),fval(*),iob(*),ixh(9,*),
     . hexarray(*),bemarray(*),shlarray(*),parm(*)
c     character*80 txts,mssg
c
c
      if (ires.lt.0) then
        n=abs(ires)
        write(iohsp,1030)
        mssg='reading user loading subroutine'
        if (longs) then
          do 11 i=1,n,8
          call gttxsg (txts,lcount)
          read (txts,'(8e20.0)',err=400) (parm(j),j=i,min(i+3,n))
          write(iohsp,1040) (j,parm(j),j=i,min(i+3,n))
   11     continue
        else
          do 10 i=1,n,8
          call gttxsg (txts,lcount)
          read (txts,1020,err=400) (parm(j),j=i,min(i+7,n))
          write(iohsp,1040) (j,parm(j),j=i,min(i+7,n))
   10     continue
        endif
        write(iohsp,1050)
        return
      endif
c
c     if (ires.eq.3) then
c       nloddv=160
c       call wrabsf64(iob,parm,nloddv,iadd64)
c       iadd64=iadd64+nloddv
c       return
c     endif
c
c     if (ires.eq.4) then
c       nloddv=160
c       call rdabsf64(iob,parm,nloddv,iadd64,ioerr)
c       iadd64=iadd64+nloddv
c       return
c     endif
c
c     for shells only
c
      do 20 i=1,numels
      ixs2i=ixs(2,i)
      ixs3i=ixs(3,i)
      ixs4i=ixs(4,i)
      ixs5i=ixs(5,i)
c
      xx11 =x(1,ixs2i)
      xx21 =x(2,ixs2i)
      xx31 =x(3,ixs2i)
c
      xx12 =x(1,ixs3i)
      xx22 =x(2,ixs3i)
      xx32 =x(3,ixs3i)
c
      xx13 =x(1,ixs4i)
      xx23 =x(2,ixs4i)
      xx33 =x(3,ixs4i)
c
      xx14 =x(1,ixs5i)
      xx24 =x(2,ixs5i)
      xx34 =x(3,ixs5i)
c
      fs1 =-xx11+xx12+xx13-xx14
      fs2 =-xx21+xx22+xx23-xx24
      fs3 =-xx31+xx32+xx33-xx34
c
      ft1 =-xx11-xx12+xx13+xx14
      ft2 =-xx21-xx22+xx23+xx24
      ft3 =-xx31-xx32+xx33+xx34
c
      e=fs1*fs1+fs2*fs2+fs3*fs3
      f=fs1*ft1+fs2*ft2+fs3*ft3
      g=ft1*ft1+ft2*ft2+ft3*ft3
      area=sqrt((e*g-f*f)/16.)
      tr1 =fs2*ft3-fs3*ft2
      tr2 =fs3*ft1-fs1*ft3
      tr3 =fs1*ft2-fs2*ft1
      xmg =.25/sqrt(tr1**2+tr2**2+tr3**2)
c
c     unit normal surface vector
c
c     pressure set by user here
c
      pressure=-5.e-04*fval(1)
c
      tr1 =xmg*tr1*area*pressure
      tr2 =xmg*tr2*area*pressure
      tr3 =xmg*tr3*area*pressure
c
      fnod(1,ixs2i)=fnod(1,ixs2i)-tr1
      fnod(2,ixs2i)=fnod(2,ixs2i)-tr2
      fnod(3,ixs2i)=fnod(3,ixs2i)-tr3
      fnod(1,ixs3i)=fnod(1,ixs3i)-tr1
      fnod(2,ixs3i)=fnod(2,ixs3i)-tr2
      fnod(3,ixs3i)=fnod(3,ixs3i)-tr3
      fnod(1,ixs4i)=fnod(1,ixs4i)-tr1
      fnod(2,ixs4i)=fnod(2,ixs4i)-tr2
      fnod(3,ixs4i)=fnod(3,ixs4i)-tr3
      fnod(1,ixs5i)=fnod(1,ixs5i)-tr1
      fnod(2,ixs5i)=fnod(2,ixs5i)-tr2
      fnod(3,ixs5i)=fnod(3,ixs5i)-tr3
   20 continue
c
c
      return
  400 call termin (txts,mssg,lcount,1)
 1020 format(8e10.0)
 1030 format(//' u s e r   d e f i n e d   l o a d i n g',
     .           '   p a r a m e t e r s',/)
 1040 format(
     1 5x,'   parameter number  ',i4,'=', e15.8)
 1050 format(//)
      end
      subroutine ushout(iop,sig,nzgpt,rslt,ix,x,nmtcon,nnel,mx,mte)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     rslt(1:8) array of force/moment resultants to be written to
c               the LS-TAURUS database unless changed in this user
c               routine. Components 26-33 in LS-TAURUS.
c
c     rslt(9:9) current shell thickness
c
c     iop       shell formulation flag
c               1=Hughes-Liu
c               2=Belytschko-Tsay
c               3=DKT
c               4=C0
c               5=Membrane
c               6=Hughes-Liu with full integration
c               7=Corotational Hughes-Liu with full integration
c               8=Belytschko-Leviathan
c               9=Membrane with 2 x 2 integratio
c              10=Belytschko-Wong
c              11=Corotational Hughes-Liu
c
c     sig(1:6) Stress array xx, yy, zz, xy, yz, zx
c              if iop=1 or 6 the stresses are in the global
c                            system
c              else          the stresses reference the local element
c                            system
c
c     sig(7:-) Constitutive history variables
c
c     nzgpt    Number of thru thickness integration points
c
c     ix(1:4)  element connectivities
c
c     nmtcon   length of stress+history variables at an integration poitn
c
c     nnel     internal element number
c
c     mx       part ID
c
c     mte      material type
c
      real*8 x
      dimension sig(nmtcon,*),rslt(9),ix(*),x(3,*)
      return
      end
      subroutine airusr (rbu,rbv,rba,time,dt1,dt2,param,hist,itrnon,
     . rbug,rbvg,rbag,icnv)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     user subroutine to initiate the inflation of the airbag
c
c     variables
c
c     displacements are defined at time n+1 in local system
c     velocites are defined at time n+1/2 in local system
c     accelerations are defined at time n in local system
c
c         rbu(1-3) total displacements in the local xyz directions
c         rbu(3-6) total rotations about the local xyz axes
c         rbv(1-3) velocities in the local xyz directions
c         rbv(3-6) rotational velocities about the local xyz axes
c         rba(1-3) accelerations in the local xyz directions
c         rba(3-6) rotational accelerations about the local xyz axes
c         time is the current time
c         dt1 is time step size at n-1/2
c         dt2 is time step size at n+1/2
c         param is user defined input parameters
c         hist is user defined history variables
c         itrnon is a flag to turn on the airbag inflation
c         rbug,rbvg,rbag, are similar to rbu,rbv,rba but are defined
c                         globally.
c         icnv is the airbag ID
c
c     the user subroutine sets the variable itrnon to:
c
c           itrnon=0 bag is not inflated
c           itrnon=1 bag inflation begins and this subroutine is
c                    not called again
c
      include 'iounits.inc'
      dimension rbu(6),rbv(6),rba(6),param(25),hist(25),
     . rbug(6),rbvg(6),rbag(6)
c
      itrnon=0
      ra=sqrt(rba(1)**2+rba(2)**2+rba(3)**2)
      if (ra.gt.param(1)) then
        itrnon=1
        write(iotty,100) time
        write(iohsp,100) time
        write(iomsg,100) time
      endif
 100  format (' Airbag activated at time ',1pe10.3)
c
      return
      end
      subroutine usrdmp (ioflg,iobuf,jadd64,maxsiz64,lendhu)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     dump out or read files
c
c     ioflag -
c         =0 return the length of magusr to lendhu
c         =1 read user defined data
c         =2 dump user defined data
c         useful for plot files etc.
c
      include 'bk06.inc'
      parameter (magusr=1)
      common/usrcom/usr(magusr)
      integer   iobuf
      integer*8 jadd64,maxsiz64
c
      magusrr=magusr
      if (magusrr.eq.1) return
c     if (magusr.eq.1) return
c
c     return the lenght of magusr
c
      if (ioflg.eq.0) then
        lendhu=magusr
c
c     read user's data
c
      elseif (ioflg.eq.1) then
        call rdabsf64 (iobuf,usr,magusr,jadd64,ioerr)
        jadd64=jadd64+magusr
c
c     dump user's data
c
      elseif (ioflg.eq.2) then
        call wrabsf64 (iobuf,usr,magusr,jadd64)
        jadd64=jadd64+magusr
      endif
c
      return
      end
      subroutine uctrl1 (numnp,ndof,time,dt1,dt2,prtc,pltc,frci,prto,
     . plto,frco,vt,vr,at,ar,ut,ur,xmst,xmsr,irbody,rbdyn,usrhv,
     . messag,totalm,cycle,idrint,mtype,mxrb,nrba,rbcor,x,rbv,nrbn,
     . nrb,xrb,yrb,zrb,axrb,ayrb,azrb,dtx,nmmat,rba,fvalnew,fvalold,
     . fvalmid,fvalnxt)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     user subroutine for solution control and is called at the
c     beginning of time step n+1.   The time at n+
c
c
c     note:  ls-dyna3d uses an internal numbering system to
c            accomodate arbitrary node numbering.   to access
c            information for user node n, address array location m,
c            m=lqf8(n,1). to obtain user node number, n,
c            corresponding to array address m, set n=lqfinv(m,1)
c
c     arguments:
c          numnp=number of nodal points
c          ndof =number of degrees of freedom per node
c          time =current solution time at n+1
c          dt1  =time step size between time n-1 and n
c          dt2  =time step size between time n and n+1
c          prtc =output interval for taurus time history data
c          pltc =output interval for taurus state data
c          frci =output interval for taurus interface force data
c          prto =output time for time history file
c          plto =output time for state data
c          frco =output time for force data
c          vt(3,numnp) =nodal translational velocity vector
c          vr(3,numnp) =nodal rotational velocity vector.  this
c                       array is defined if and only if ndof=6
c          at(3,numnp) =nodal translational acceleration vector
c          ar(3,numnp) =nodal rotational acceleration vector.  this
c                       array is defined if and only if ndof=6
c          ut(3,numnp) =nodal translational displacement vector
c          ur(3,numnp) =nodal rotational displacement vector.  this
c                       array is defined if and only if ndof=6
c          xmst(numnp) =reciprocal of nodal translational masses
c          xmsr(numnp) =reciprocal of nodal rotational masses.  this
c                       array is defined if and only if ndof=6
c          fvalold     =array for storing load curve values at time n
c          fvalnew     =array for setting load curve values at time n+1
c                       only load curves with 0 input points may be user
c                       defined.  When the load curve is user set, the
c                       value at time n must be stored in array fvalold.
c          fvalmid     =array for predicting load curve values at time n+3/2
c          fvalnxt     =array for predicting load curve values at time n+2
c                       for some applications it is necessary to predict the
c                       load curve values at time n+2, a time that is not known,
c                       this for instance for boundary prescribed motion. In
c                       this case the load curve values at time n+3/2 need to
c                       be predicted in fvalmid, and fvalold should be set to
c                       fvalnew and fvalnew should be set to fvalnxt. See
c                       coding below.
c
c          irbody      =0 if no rigid bodies
c          rbdyn(numnp)=flag for rigid body nodal points
c                       if deformable node then set to 1.0
c                       if rigid body node then set to 0.0
c                       defined if and only if rigid bodies are present
c                       i.e., irbody.ne.0 if no rigid bodies are
c                       present
c          usrhv(lenhv)=user defined history variables that are stored
c                       in the restart file.  lenhv=100+7*nummat where
c                       nummat is the # of materials in the problem.
c                       array usrhv is updated only in this subroutine.
c          messag      =flag for dyna3d which may be set to:
c                      'sw1.' ls-dyna3d terminates with restart file
c                      'sw3.' ls-dyna3d writes a restart file
c                     'sw4.' ls-dyna3d writes a plot state
c           totalm     =total mass in problem
c           cycle     =cycle number
c           idrint    =flag for dynamic relaxation phase
c                      .ne.0: dynamic relaxation in progress
c                      .eq.0: solution phase
c
      include 'ptimes.inc'
c
c     prtims(1-37)=output intervals for ascii files
c
c      ascii files:
c           ( 1)-cross section forces
c           ( 2)-rigid wall forces
c           ( 3)-nodal data
c           ( 4)-element data
c           ( 5)-global data
c           ( 6)-discrete elements
c           ( 7)-material energies
c           ( 8)-nodal interface forces
c           ( 9)-resultant interface forces
c           (10)-smug animator
c           (11)-spc reaction forces
c           (12)-nodal constraint resultant forces
c           (13)-airbag statistics
c           (14)-avs database
c           (15)-nodal force groups
c           (16)-output intervals for nodal boundary conditions
c           (17)-(32) unused at this time
c           (37)-auto tiebreak damage output
c
c     prtlst(32)=output times for ascii files above.  when solution time
c                exceeds the output time a print state is dumped.
c
      common/rbkeng/enrbdy,rbdyx,rbdyy,rbdyz
c
c      total rigid body energies and momentums:
c           enrbdy=rigid body kinetic energy
c           rbdyx =rigid body x-momentum
c           rbdyy =rigid body y-momentum
c           rbdyz =rigid body z-momentum
c
      common/swmke/swxmom,swymom,swzmom,swkeng
c
c      total stonewall energies and momentums:
c           swxmom=stonewall x-momentum
c           swymom=stonewall y-momentum
c           swzmom=stonewall z-momentum
c           swkeng=stonewall kinetic energy
c
      common/deengs/deeng
c
c      deeng=total discrete element energy
c
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie,selie,selke,
     . erodehg
c
c      xpe  =total internal energy in the finite elements
c
      common/sprengs/spreng
c
c      spreng=total spr energy
c
      character*(*) messag
      integer cycle
      real*8 x
      dimension vt(3,*),vr(3,*),at(3,*),ar(3,*),
     . xmst(*),xmsr(*),rbdyn(*),usrhv(*),mtype(*),mxrb(*),nrba(*),
     . rbcor(3,1),x(*),rbv(6,*),nrbn(*),nrb(*),xrb(*),yrb(*),
     . zrb(*),axrb(*),ayrb(*),azrb(*),rba(6,*),fvalnew(*),
     . fvalold(*),fvalmid(*),fvalnxt(*)
      real*8 ut(3,*),ur(3,*)
c
c     MPP special variables
c
c     remove all "cmpp" for output to a single file in MPP for UMAT routine
c
c
c     sample momentum and kinetic energy calculations
c
c     remove all comments in column 1 below to activate
      i=1
      if (i.eq.1) return
c     return
cc
cc
cc     initialize kinetic energy, xke, and x,y,z momentums.
cc
c      xke=2.*swkeng+2.*enrbdy
c      xm=swxmom+rbdyx
c      ym=swymom+rbdyy
c      zm=swzmom+rbdyz
cc
c      numnp2=numnp
c      if (ndof.eq.6) then
c        numnp2=numnp+numnp
c      endif
c      write(iotty,*)ndof
cc
cc
cc     no rigid bodies present
cc
c      if (irbody.eq.0) then
cc       note in blank comment vr follows vt.  this fact is used below.
c        do 10 n=1,numnp2
c        xmsn=1./xmst(n)
c        vn1=vt(1,n)
c        vn2=vt(2,n)
c        vn3=vt(3,n)
c        xm=xm+xmsn*vn1
c        ym=ym+xmsn*vn2
c        zm=zm+xmsn*vn3
c        xke=xke+xmsn*(vn1*vn1+vn2*vn2+vn3*vn3)
c   10   continue
cc
cc
cc     rigid bodies present
cc
c      else
cc       nodal accerations for rigid bodies
cc
c        do 12 n=1,nmmat
c        if (mtype(n).ne.20.or.mxrb(n).ne.n) go to 12
c        lrbn=nrba(n)
c        call stvlut(rbcor(1,n),x,vt,at,ar,vr,rbv(1,n),dt2,
c    .   nrbn(n),nrb(lrbn),xrb,yrb,zrb,axrb,ayrb,azrb,dtx)
c
c        rigid body nodal accelerations
c
c        if (ndof.eq.6) then
c          call rbnacc(nrbn(n),nrb(lrbn),rba(4,n),ar)
c        endif
c
c  12    continue
cc
c        do 20 n=1,numnp
c        xmsn=1./xmst(n)
c        vn1=rbdyn(n)*vt(1,n)
c        vn2=rbdyn(n)*vt(2,n)
c        vn3=rbdyn(n)*vt(3,n)
c        xm=xm+xmsn*vn1
c        ym=ym+xmsn*vn2
c        zm=zm+xmsn*vn3
c        xke=xke+xmsn*(vn1*vn1+vn2*vn2+vn3*vn3)
c   20   continue
c        if (ndof.eq.6) then
c          do 30 n=1,numnp
c          xmsn=1./xmsr(n)
c          vn1=rbdyn(n)*vr(1,n)
c          vn2=rbdyn(n)*vr(2,n)
c          vn3=rbdyn(n)*vr(3,n)
c          xm=xm+xmsn*vn1
c          ym=ym+xmsn*vn2
c          zm=zm+xmsn*vn3
c          xke=xke+xmsn*(vn1*vn1+vn2*vn2+vn3*vn3)
c   30     continue
c        endif
c      endif
cc
cc     total kinetic energy
c      xke=.5*xke
cc     total internal energy
c      xie=xpe+deeng+spreng
cc     total energy
c      xte=xke+xpe+deeng+spreng
cc     total x-rigid body velocity
c      xrbv=xm/totalm
cc     total y-rigid body velocity
c     yrbv=ym/totalm
cc     total z-rigid body velocity
c      zrbv=zm/totalm
c
c     New source to control load curve 2 in model during solution
c
c     External DYNA model load curve id
      lidext=2
c     Obtain the internal load curve id
      lidint=lqfl(lidext)
c     Put old 'new' value into fvalold
      fvalold(lidint)=fvalnew(lidint)
      fvalnew(lidint)=fvalnxt(lidint)
c     Calculate value for fvalnxt
      fvalnxt(lidint)=100.0*sin(time*3.1415926)
      fvalmid(lidint)=.5*(fvalnxt(lidint)+fvalnew(lidint))
c
c     New source to check that we can correctly obtain the nodal
c     displacements and velocities for an arbitrary node.
c     Using node 3 in this case on part A of model.
c     Obtain the displacements and velocities at t=0.3s.
c
      if (abs(time-0.3).lt.1e-6) then
c        Time is 0.3s (time step is 2.5e-6s)
         call lsopen(85,'node3.out','formatted','new',ierr,0,0)
c        External DYNA model node number
         nodext=3
c        Obtain the internal node id
         nodint=lqf8(nodext,1)
         write (85,'('' Time = '',E10.4E2,/)')time
c        Obtain and write out the nodal displacements
         write (85,'('' Node 3 x disp. = '',E10.4E2)')ut(1,nodint)
         write (85,'('' Node 3 y disp. = '',E10.4E2)')ut(2,nodint)
         write (85,'('' Node 3 z disp. = '',E10.4E2)')ut(3,nodint)
c        Obtain and write out the nodal velocities
         write (85,'('' Node 3 x velo. = '',E10.4E2)')vt(1,nodint)
         write (85,'('' Node 3 y velo. = '',E10.4E2)')vt(2,nodint)
         write (85,'('' Node 3 z velo. = '',E10.4E2)')vt(3,nodint)
         close (unit=85)
      endif
c
      return
      end
      subroutine rbnacc(nrbn,nrb,rba,ar)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      dimension nrb(*),rba(6),ar(3,*)
c
      do 20 i=1,nrbn
      ar(1,nrb(i))=rba(1)
      ar(2,nrb(i))=rba(2)
      ar(3,nrb(i))=rba(3)
   20 continue
c
      return
      end
      subroutine uctrl2(nsi,nty,time,cycle,msr,nmn,nsv,nsn,
     1 thmr,thsv,vt,xi,ut,iskip,idrint,numnp,dt2,ninput,ua,
     2 irectm,nrtm,irects,nrts)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     user subroutine for interface control
c
c     note:  ls-dyna3d uses an internal numbering system to
c            accomodate arbitrary node numbering.   to access
c            information for user node n, address array location m,
c            m=lqf8(n,1). to obtain user node number, n,
c            corresponding to array address m, set n=lqfinv(m,1)
c
c     arguments:
c          nsi        =number of sliding interface
c          nty        =interface type.
c                      .eq.4:single surface
c                      .ne.4:surface to surface
c          time       =current solution time
c          cycle      =cycle number
c          msr(nmn)   =list of master nodes numbers in internal
c                      numbering scheme
c          nmn        =number of master nodes
c          nsv(nsn)   =list of slave nodes numbers in internal
c                      numbering scheme
c          nsn        =number of slave nodes
c          thmr(nmn)  =master node thickness
c          thsv(nsn)  =slave node thickness
c          vt(3,numnp)=nodal translational velocity vector
c          xi(3,numnp)=initial coordinates at time=0
c          ut(3,numnp)=nodal translational displacement vector
c          idrint     =flag for dynamic relaxation phase
c                      .ne.0: dynamic relaxation in progress
c                      .eq.0: solution phase
c          numnp      =number of nodal points
c          dt2        =time step size at n+1/2
c          ninput     =number of variables input into ua
c          ua(*)      =users' array, first ninput locations
c                      defined by user.  the length of this
c                      array is defined on control card 10.
c                      this array is unique to interface nsi.
c          irectm(4,*)=list of master segments in internal
c                      numbering scheme
c          nrtm       =number of master segments
c          irects(4,*)=list of slave segments in internal
c                      numbering scheme
c          nrts       =number of master segments
c
c     set flag for active contact
c       iskip=0 active
c       iskip=1 inactive
c
c*******************************************************************
c
      integer cycle
      real*8 ut
      real*8 xi
      dimension msr(*),nsv(*),thmr(*),thsv(*),vt(3,*),xi(3,*),
     .          ut(3,*),ua(*),irectm(4,*),irects(4,*)
c
c     the following sample of codeing is provided to illustrate how
c     this subroutine might be used.   here we check to see if the
c     surfaces in the surface to surface contact are separated.   if
c     so the iskip=1 and the contact treatment is skipped.
c
c     if (nty.eq.4) return
c     dt2hlf=dt2/2.
c     xmins= 1.e+16
c     xmaxs=-xmins
c     ymins= 1.e+16
c     ymaxs=-ymins
c     zmins= 1.e+16
c     zmaxs=-zmins
c     xminm= 1.e+16
c     xmaxm=-xminm
c     yminm= 1.e+16
c     ymaxm=-yminm
c     zminm= 1.e+16
c     zmaxm=-zminm
c     thks=0.0
c     thkm=0.0
c     do 10 i=1,nsn
c     dsp1=ut(1,nsv(i))+dt2hlf*vt(1,nsv(i))
c     dsp2=ut(2,nsv(i))+dt2hlf*vt(2,nsv(i))
c     dsp3=ut(3,nsv(i))+dt2hlf*vt(3,nsv(i))
c     x1=xi(1,nsv(i))+dsp1
c     x2=xi(2,nsv(i))+dsp2
c     x3=xi(3,nsv(i))+dsp3
c     thks = max(thsv(i),thks)
c     xmins=min(xmins,x1)
c     xmaxs=max(xmaxs,x1)
c     ymins=min(ymins,x2)
c     ymaxs=max(ymaxs,x2)
c     zmins=min(zmins,x3)
c     zmaxs=max(zmaxs,x3)
c  10 continue
c     do 20 i=1,nmn
c     dsp1=ut(1,msr(i))+dt2hlf*vt(1,msr(i))
c     dsp2=ut(2,msr(i))+dt2hlf*vt(2,msr(i))
c     dsp3=ut(3,msr(i))+dt2hlf*vt(3,msr(i))
c     x1=xi(1,msr(i))+dsp1
c     x2=xi(2,msr(i))+dsp2
c     x3=xi(3,msr(i))+dsp3
c     thkm = max(thmr(i),thks)
c     xminm=min(xminm,x1)
c     xmaxm=max(xmaxm,x1)
c     yminm=min(yminm,x2)
c     ymaxm=max(ymaxm,x2)
c     zminm=min(zminm,x3)
c     zmaxm=max(zmaxm,x3)
c  20 continue
c
c     if thks or thkm equal zero set them to some reasonable value
c
c     if (thks.eq.0.0) then
c       e1=(xi(1,irects(1,1))-xi(1,irects(3,1)))**2
c    .    +(xi(2,irects(1,1))-xi(2,irects(3,1)))**2
c    .    +(xi(3,irects(1,1))-xi(3,irects(3,1)))**2
c       e2=(xi(1,irects(2,1))-xi(1,irects(4,1)))**2
c    .    +(xi(2,irects(2,1))-xi(2,irects(4,1)))**2
c    .    +(xi(3,irects(2,1))-xi(3,irects(4,1)))**2
c       thks=.3*sqrt(max(e1,e2))
c     endif
c     if (thkm.eq.0.0) then
c       e1=(xi(1,irectm(1,1))-xi(1,irectm(3,1)))**2
c    .    +(xi(2,irectm(1,1))-xi(2,irectm(3,1)))**2
c    .    +(xi(3,irectm(1,1))-xi(3,irectm(3,1)))**2
c       e2=(xi(1,irectm(2,1))-xi(1,irectm(4,1)))**2
c    .    +(xi(2,irectm(2,1))-xi(2,irectm(4,1)))**2
c    .    +(xi(3,irectm(2,1))-xi(3,irectm(4,1)))**2
c       thkm=.3*sqrt(max(e1,e2))
c     endif
c
c     if (xmaxs+thks.lt.xminm-thkm) go to 40
c     if (ymaxs+thks.lt.yminm-thkm) go to 40
c     if (zmaxs+thks.lt.zminm-thkm) go to 40
c     if (xmaxm+thkm.lt.xmins-thks) go to 40
c     if (ymaxm+thkm.lt.ymins-thks) go to 40
c     if (zmaxm+thkm.lt.zmins-thks) go to 40
c     iskip=0
c
c     return
c  40 iskip=1
c
      return
      end
      subroutine matusr_24(lft,llt,ietype,cm,crv,ipt,ipt_thk,nip,
     . nelbwp,lfail)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     user-defined damage subroutine for materials:
c     24, 36, 114, 123, 124, 133, 155, 182, 225, 238, 243, 251, 255
c
c     Description of Input Arguments:
c
c     sig1 - stress in x  direction, local x for shells
c     sig2 - stress in y  direction, local y for shells
c     sig3 - stress in z  direction, local z for shells
c     sig4 - stress in xy direction, local for shells
c     sig5 - stress in yz direction, local for shells
c     sig6 - stress in zx direction, local for shells
c
c     dt1siz-element time step size
c
c     incremental strains are used for shell elements
c     strain rates are used for solid elements
c
c     d1 - strain rate/increment in x  direction, local x for shells
c     d2 - strain rate/increment in y  direction, local y for shells
c     d3 - strain rate/increment in z  direction, local z for shells
c     d4 - strain rate/increment in xy direction, local for shells
c     d5 - strain rate/increment in yz direction, local for shells
c     d6 - strain rate/increment in zx direction, local for shells
c
c     wxxdt - time-step*spin rate in x  direction, only for solds       
c     wyydt - time-step*spin rate in y  direction, only for solds       
c     wzzdt - time-step*spin rate in z  direction, only for solds       
c
c     The following strain components are only available when
c     INTOUT in *DATABASE_EXTENT_BINARY is set to "STRAIN" or "ALL"
c
c     eps1 - strain in x  direction, local x for shells
c     eps2 - strain in y  direction, local y for shells
c     eps3 - strain in z  direction, local z for shells
c     eps4 - strain in xy direction, local for shells
c     eps5 - strain in yz direction, local for shells
c     eps6 - strain in zx direction, local for shells
c
c     eint - strain energy density increment
c     ep - plastic strain
c     fail - failure flag SET TO ZERO IF ELEMENT FAILS
c     n - number of elements in vector block
c     ietype - 1 - beams, 2 - shells, 3 - solids.
c
c     lft: 1st element to be processed
c     llt: last element to be processed
c
c     cm: material properties
c       -1: young's modulus
c       -2: poissons's ratio
c
c     crv: curve array
c
c     ipt: integration point number
c     ipt_thk: through-thickness integration point number
c     nip: total number of integration points
c
c     nelbwp: element id array
c
c     - solids: idele=lqfinv(nelbwp(i),2)
c     - beams : idele=lqfinv(nelbwp(i),3)
c     - shells: idele=lqfinv(nelbwp(i),4)
c
c     nintcy: current cycle
c     tt: current time
c
c     lfail: where failure flag is stored, defined by ls-dyna,
c           =0: failure is stored in fail0
c           >2: failure is stored in fail1
c
c     Be aware that in addition to failure criteria set forth in this routine,
c     failure also occurs if ep >= abs(FAIL) where FAIL is the 'failure flag'
c     provided in the mat24 input.
c
      include 'nlqparm'
      common/aux2loc/d1(nlq),d2(nlq),d3(nlq),d4(nlq),d5(nlq),d6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common/aux13loc/zeta(nlq),thick(nlq)
      common/aux14loc/
     1 sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),
     2 sig5(nlq),sig6(nlq),  ep(nlq),fail1(nlq,142)
      common/failuloc/tmp1(nlq),fail0(nlq)
      common/subtssloc/dt1siz(nlq)
      common/strnintloc/eps1(nlq),eps2(nlq),eps3(nlq),eps4(nlq),
     1 eps5(nlq),eps6(nlq)
      common/bk26/begtim,nintcy
      common/bk28/summss,xke,xpe,tt
      dimension cm(*),crv(lq1,2,*),nelbwp(*)
c
c     loop over elements
c     if (lfail.eq.0) then
c       do i=lft,llt
c        if (failure condition is met) fail0(i)=0.0
c       enddo
c     else
c       do i=lft,llt
c        if (failure condition is met) fail1(i,lfail)=0.0
c       enddo
c     endif
c$omp threadprivate (/aux2loc/)
c$omp threadprivate (/aux13loc/)
c$omp threadprivate (/aux14loc/)
c$omp threadprivate (/failuloc/)
c$omp threadprivate (/strnintloc/)
c$omp threadprivate (/subtssloc/)
      return
      end
      subroutine matusr_103 (sig1,sig2,sig3,sig4,sig5,sig6,dt1siz,
     . d1,d2,d3,d4,d5,d6,ep,fail,n,ietype,q1,q2,d1m,d2m,d3m,d4m,d5m,
     . d6m,failure_flag,s11,s12,s13,s21,s22,s23,s31,s32,s33,cm,crv,
     . parms,nparms,dmg,aux,naux,fracs,nfracs,clen,cang,dep)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     description of Input Arguments, global for solid element
c
c     sig1 - stress in x  direction, local x for shells
c     sig2 - stress in y  direction, local y for shells
c     sig3 - stress in z  direction, local z for shells
c     sig4 - stress in xy direction, local for shells
c     sig5 - stress in yz direction, local for shells
c     sig6 - stress in zx direction, local for shells
c
c     dt1siz-element time step size
c
c     incremental strains are used for shell elements
c     strain rates are used for solid elements
c
c     d1  - strain rate/increment in x  direction, local x for shells
c     d2  - strain rate/increment in y  direction, local y for shells
c     d3  - strain rate/increment in z  direction, local z for shells
c     d4  - strain rate/increment in xy direction, local for shells
c     d5  - strain rate/increment in yz direction, local for shells
c     d6  - strain rate/increment in zx direction, local for shells
c
c     d1-d6 are in the a,b,c material local system for solids
c
c     q1,q2  - direction cosines to obtain material direction (shells)
c
c     To obtain stresses in shell a,b coordinate system:
c
c           a11    =q1(i)*sig1(i)-q2(i)*sig4(i)
c           a12    =q2(i)*sig1(i)+q1(i)*sig4(i)
c           a21    =q1(i)*sig4(i)-q2(i)*sig2(i)
c           a22    =q2(i)*sig4(i)+q1(i)*sig2(i)
c           sig1(i)=q1(i)*a11-q2(i)*a21
c           sig2(i)=q2(i)*a12+q1(i)*a22
c           sig4(i)=q1(i)*a12-q2(i)*a22
c
c     To transform back to shell local
c
c           a11    = sig1(i)*q1(i)+sig4(i)*q2(i)
c           a12    =-sig1(i)*q2(i)+sig4(i)*q1(i)
c           a21    = sig4(i)*q1(i)+sig2(i)*q2(i)
c           a22    =-sig4(i)*q2(i)+sig2(i)*q1(i)
c           sig1(i)= q1(i)*a11+q2(i)*a21
c           sig2(i)=-q2(i)*a12+q1(i)*a22
c           sig4(i)= q1(i)*a12+q2(i)*a22
c
c
c     d1m - strain rate/increment in a  direction (shells)
c     d2m - strain rate/increment in b  direction (shells)
c     d3m - strain rate/increment in c  direction (shells)
c     d4m - strain rate/increment in ab direction (shells)
c     d5m - strain rate/increment in bc direction (shells)
c     d6m - strain rate/increment in ca direction (shells)
c
c     ep - plastic strain
c     fail - failure flag SET TO ZERO IF ELEMENT FAILS
c     n - number of elements in vector block
c     ietype  beams, 2 - shells, 3 - solids.
c
c     Be aware that in addition to failure criteria set forth in this routine,
c     failure also occurs if ep >= abs(FAIL) where FAIL is the 'failure flag'
c     provided in the mat24 input.
c
c     s11-s33 transformation matrix for SOLIDS only.  To transform the stresses
c     into the local system:
c
c     a11    =sig1(i)*s11(i)+sig4(i)*s12(i)+sig6(i)*s13(i)
c     a12    =sig1(i)*s21(i)+sig4(i)*s22(i)+sig6(i)*s23(i)
c     a13    =sig1(i)*s31(i)+sig4(i)*s32(i)+sig6(i)*s33(i)
c     a21    =sig4(i)*s11(i)+sig2(i)*s12(i)+sig5(i)*s13(i)
c     a22    =sig4(i)*s21(i)+sig2(i)*s22(i)+sig5(i)*s23(i)
c     a23    =sig4(i)*s31(i)+sig2(i)*s32(i)+sig5(i)*s33(i)
c     a31    =sig6(i)*s11(i)+sig5(i)*s12(i)+sig3(i)*s13(i)
c     a32    =sig6(i)*s21(i)+sig5(i)*s22(i)+sig3(i)*s23(i)
c     a33    =sig6(i)*s31(i)+sig5(i)*s32(i)+sig3(i)*s33(i)
c     sig1(i)=s11(i)*a11+s12(i)*a21+s13(i)*a31
c     sig2(i)=s21(i)*a12+s22(i)*a22+s23(i)*a32
c     sig3(i)=s31(i)*a13+s32(i)*a23+s33(i)*a33
c     sig4(i)=s11(i)*a12+s12(i)*a22+s13(i)*a32
c     sig5(i)=s21(i)*a13+s22(i)*a23+s23(i)*a33
c     sig6(i)=s11(i)*a13+s12(i)*a23+s13(i)*a33
c
c     warning: For shells the stresses must be in the shell local system
c              when returning from this subroutine.  For solids, the
c              stresses must be in the global system when returning.
c
c     Additional arguments for supporting user defined damage model
c
c     parms(nparms)      = damage parameters, input
c     dmg(i)             = damage, input/output
c     aux(i,naux)        = history variable, input/output
c     fracs(i,nfracs)    = fraction material of phases, input
c     clen(i)            = characteristic element length divided by shell
c                          thickness (for solids this is zero), input
c     cang(i)            = cosine of angle between loading and first material
c                          direction, input
c     dep(i)             = plastic strain rate, input
c
      include 'nlqparm'
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*),
     . dt1siz(*),d1(*),d2(*),d3(*),d4(*),d5(*),d6(*),ep(*),
     . fail(*),q1(*),q2(*),d1m(*),d2m(*),d3m(*),d4m(*),d5m(*),d6m(*),
     . s11(*),s12(*),s13(*),s21(*),s22(*),s23(*),s31(*),s32(*),s33(*),
     . cm(*),crv(lq1,2,*)
      dimension parms(*),dmg(*),aux(nlq,*),fracs(nlq,*),
     .     clen(*),cang(*),dep(*)
c
c     loop over elements, Put new failure criteria here.
c
       do i=1,n
c      if (failure condition is met) fail(i)=0.0
c      example
       if (ep(i).gt.abs(failure_flag)) fail(i)=0.
c
       enddo
      return
      end
      subroutine usrfrc(nosl,time,ncycle,dt2,insv,areas,xs,ys,zs,
     . lsv,ix1,ix2,ix3,ix4,aream,xx1,xx2,xx3,stfn,stf,fni,
     . dx,dy,dz,fdt2,ninput,ua,side,iisv5,niisv5,rn1,rn2,rn3,fric1,
     . fric2,fric3,fric4,bignum,fdat,iseg,fxis,fyis,fzis,ss,tt,
     . ilbsv,stfk,frc,numnp,npc,pld,lcfst,lcfdt,temp,temp_bot,
     . temp_top,isurface,nhv,ufhist,neps,ufeps,sfacn3)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     user subroutine for interface friction control (smp)
c
c     note:  ls-dyna3d uses an internal numbering system to
c            accomodate arbitrary node numbering.   to access
c            information for user node n, address array location m,
c            m=lqf8(n,1). to obtain user node number, n,
c            corresponding to array address m, set n=lqfinv(m,1)
c
c     arguments:
c
c          nosl       =number of sliding interface
c          time       =current solution time
c          ncycle     =ncycle number
c          dt2        =time step size at n+1/2
c          insv       =slave node array where the nodes are stored
c                      in ls-dyna3d internal numbering.  User numbers
c                      are given by function: lqfinv(insv(ii),1)
c                      for slave node ii.
c          areas(ii)  =slave node area (interface types 5&10 only) for
c                      slave node ii (forming contact, too)
c          xs(ii)     =x-coordinate slave node ii (projected)
c          ys(ii)     =y-coordinate slave node ii (projected)
c          zs(ii)     =z-coordinate slave node ii (projected)
c          lsv(ii)    =master segment number for slave node ii
c          ix1(ii), ix2(ii), ix3(ii), ix4(ii)
c                     =master segment nodes in ls-dyna3d internal
c                      numbering for slave node ii
c          aream(ii)  =master segment area for slave node ii.
c                      (not for forming contacts)
c          xx1(ii,4)  =x-coordinates master surface (projected) for
c                      slave node ii (not for forming contacts)
c          xx2(ii,4)  =y-coordinates master surface (projected) for
c                      slave node ii (not for forming contacts)
c          xx3(ii,4)  =z-coordinates master surface (projected) for
c                      slave node ii (not for forming contacts)
c          stfn       =slave node penalty stiffness
c          stf        =master segment penalty stiffness
c          fni        =normal force
c          dx,dy,dz   =relative x,y,z-displacement between slave node and
c                      master surface.  Multipling by fdt2 defines the
c                      relative velocity.
c          rn1,rn2,rn3=x,y, and z components of master segments normal
c                      vector
c          sfacn3     =scaling factor for implicit calculations
c
c***********************************************************************
c       frictional coefficients defined for the contact interface
c
c          fric1      =static friction coefficient
c          fric2      =dynamic friction coefficient
c          fric3      =decay constant
c          fric4      =viscous friction coefficient (setting fric4=0
c                      turns this option off)
c
c***********************************************************************
c
c          bignum     =0.0 for one way surface to surface and
c                      for surface to surface, and 1.e+10 for nodes
c                      to surface contact
c                     =0.0 for forming contact
c          ninput     =number of variables input into ua
c          ua(*)      =users' array, first ninput locations
c                      defined by user.  the length of this
c                      array is defined on control card 10.
c                      this array is unique to interface nosl.
c
c          side       ='master' for first pass.  the master
c                       surface is the surface designated in the
c                       input
c                     ='slave' for second pass after slave and
c                      master surfaces have be switched for
c                      the type 3 symmetric interface treatment.
c
c          iisv5      =an array giving the pointers to the active nodes
c                      in the arrays.
c
c          niisv5     =number of active nodes
c
c          fdat       =contact history data array containing:
c                       1 - ss (isoparametric coord x)
c                       2 - tt (isoparametric coord y)
c                       3 - old friction x-force
c                       4 - old friction y-force
c                       5 - old friction z-force
c          iseg       =contact master segment from previous step.
c          fxis       =slave node force component in global x dir.
c                      to be updated to include friction
c          fyis       =slave node force component in global y dir.
c                      to be updated to include friction
c          fzis       =slave node force component in global z dir.
c                      to be updated to include friction
c          ss(ii)     =s contact point (-1 to 1) in parametric coordinates
c                      for slave node ii.
c          tt(ii)     =t contact point (-1 to 1) in parametric coordinates
c                      for slave node ii.
c          ilbsv(ii)  =pointer for node ii into global arrays.
c          stfk(ii)   =penalty stiffness for slave node ii which was used
c                      to compute normal interface force.
c          frc(1,lsv(ii))
c                    =Coulomb friction scale factor for segment lsv(ii)
c          frc(2,lsv(ii))
c                    =viscous friction scale factor for segment lsv(ii)
c
c***********************************************************************
c       parameters for a coupled thermal-mechanical contact
c
c          numnp      = number of nodal points in the model
c          npc        = load curve pointer
c          pld        = load curve (x,y) data
c          lcfst(nosl)= load curve number for static coefficient of
c                       friction versus temperture for contact
c                       surface nosl
c          lcfdt(nosl)= load curve number for dynamic coefficient of
c                       friction versus temperture for contact
c                       surface nosl
c          temp(j)    = temperature for node point j
c          temp_bot(j)= temparature for thick thermal shell bottom
c                       surface
c          temp_top(j)= temparature for thick thermal shell top
c                       surface
c          numsh12    = number of thick thermal shells
c          itopaz(1)  = 999 ==> thermal-mechanical analysis
c          isurface   = thick thermal shell surface pointer
c
c***********************************************************************
c
c          nhv = size of user friction history data array
c          ufhist(nhv,ilbsv(ii)) = user friction history data array
c          neps= size of element history data array
c          ufeps(neps,ilbsv(ii)) = element history data array
c                                  1 - plastic strain (history var #1)
c
c***********************************************************************
      include 'nlqparm'
      integer
     1 nh01,nh02,nh03,nh04,nh05,nh06,nh07,nh08,nh09,nh10,
     2 nb01,nb02,nb03,nb04,nb05,nb06,nb07,nb08,nb09,nb10,
     3 ns01,ns02,ns03,ns04,ns05,ns06,ns07,ns08,ns09,ns10,
     4 nt01,nt02,nt03,nt04,nt05,nt06,nt07,nt08,nt09,nt10,
     5 nu01,nu02,nu03,nu04,nu05,nu06,nu07,nu08,nu09,nu10,
     6 nd01,nd02,nd03,nd04,nd05,nd06,nd07,nd08,nd09,nd10,
     7 ntbqs
      common/bk05/
     1 nh01,nh02,nh03,nh04,nh05,nh06,nh07,nh08,nh09,nh10,
     2 nb01,nb02,nb03,nb04,nb05,nb06,nb07,nb08,nb09,nb10,
     3 ns01,ns02,ns03,ns04,ns05,ns06,ns07,ns08,ns09,ns10,
     4 nt01,nt02,nt03,nt04,nt05,nt06,nt07,nt08,nt09,nt10,
     5 nu01,nu02,nu03,nu04,nu05,nu06,nu07,nu08,nu09,nu10,
     6 nd01,nd02,nd03,nd04,nd05,nd06,nd07,nd08,nd09,nd10,
     7 ntbqs
c
      integer*8 dm_hmtnum,dm_hnfegp
      integer*8 dm_bmtnum,dm_bnfegp
      integer*8 dm_smtnum,dm_snfegp
      integer*8 dm_tmtnum,dm_tnfegp
c
      common/dmbk05/
     1 dm_hmtnum,dm_hnfegp,
     2 dm_bmtnum,dm_bnfegp,
     3 dm_smtnum,dm_snfegp,
     4 dm_tmtnum,dm_tnfegp
      integer n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      common/bk07/n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      integer n4a,n4b,n4c,n4d,n4e,n4f,n4g,n4h,n7a,n7b,n7c,nusir,
     1 mpusr,mpubr,n4i,n4j,n4k,n4l
      common/bk08/n4a,n4b,n4c,n4d,n4e,n4f,n4g,n4h,n7a,n7b,n7c,nusir,
     1 mpusr,mpubr,n4i,n4j,n4k,n4l
      integer lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      common/bk13/lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      common/bmcntc/ibmcnt(101)
      common/thropt/itopaz(101),iprm_thrm(20),rprm_thrm(20)
      common/blk03/rdumy1(2),idumy6(25),numsh12,ndtot,nsl_th,lenhsv,
     . numel2,numel4,numel6,numel8,numel10
      dimension insv(*),areas(*),xs(*),ys(*),zs(*),lsv(*),ix1(*),
     . ix2(*),ix3(*),ix4(*),aream(*),xx1(nlq,*),xx2(nlq,*),xx3(nlq,*),
     . stfn(*),stf(*),fni(*),dx(*),dy(*),dz(*),ua(*),iisv5(*),
     . fdat(5,*),iseg(*),fxis(*),fyis(*),fzis(*),ss(*),tt(*),
     . ilbsv(*),stfk(*),frc(2,*),npc(*),pld(2,*),lcfst(*),lcfdt(*),
     . temp(*),temp_bot(*),temp_top(*),isurface(numnp,*),
     . ufhist(nhv,*),ufeps(neps,*)
c
      real rn1(*),rn2(*),rn3(*)
c
      integer ncycle
      character*(*) side
c
c     set coefficient of friction multipliers to 1.
      fstm=1.
      fdym=1.
c
c     loop on the number of active nodes niisv5 in contact
      do kk = 1, niisv5
c
c     get the location pointer ii for this active node
      ii = iisv5(kk)
c
c------------------------- begin thermal -------------------------------
c     recalculate coefficient of friction multipliers if
c     function of temperature
c
      if (itopaz(1) .eq. 999) then
c
c     get the node number jj
c     the user node number is obtained using node_user=lqfinv(jj,1)
      jj = insv(ii)
c
c     get the temperature of the active node.
c     use the array isurface to distinquish between a bottom,
c     midplane, or top surface node for a thick thermal shell
c      if (ibmcnt(19).ne.0) then
      if (numsh12.ne.0) then
      if (isurface(jj,nosl) .lt. 0) then
        tnode=temp_bot(jj)
      elseif (isurface(jj,nosl) .gt. 0) then
        tnode=temp_top(jj)
      else
        tnode=temp(jj)
      endif
      else
        tnode=temp(jj)
      endif
c
c     get the shape functions for this active node projected
c     onto the surface it contacts
      h1=.25*(1.-ss(ii))*(1.-tt(ii))
      h2=.25*(1.+ss(ii))*(1.-tt(ii))
      h3=.25*(1.+ss(ii))*(1.+tt(ii))
      h4=.25*(1.-ss(ii))*(1.+tt(ii))
c     get the segment temperature at the interpolation point
c     use the array isurface to distinquish between a bottom,
c     midplane, or top surface node for a thick thermal shell
c      if (ibmcnt(19).ne.0) then
      if (numsh12.ne.0) then
      if (isurface(ix1(ii),nosl) .lt. 0) then
        tseg=h1*temp_bot(ix1(ii))+h2*temp_bot(ix2(ii))
     1   +h3*temp_bot(ix3(ii))+h4*temp_bot(ix4(ii))
      elseif (isurface(ix1(ii),nosl) .gt. 0) then
        tseg=h1*temp_top(ix1(ii))+h2*temp_top(ix2(ii))
     1   +h3*temp_top(ix3(ii))+h4*temp_top(ix4(ii))
      else
        tseg=h1*temp(ix1(ii))+h2*temp(ix2(ii))
     1   +h3*temp(ix3(ii))+h4*temp(ix4(ii))
      endif
      else
        tseg=h1*temp(ix1(ii))+h2*temp(ix2(ii))
     1   +h3*temp(ix3(ii))+h4*temp(ix4(ii))
      endif
c
c     get the average contact temperature
      tavg=.5*(tnode+tseg)
c
c     get static coefficient of friction temperature multiplier fstm
      call value (nosl,lcfst,npc,pld,tavg,tavg,fstm,valp,1)
c     get dynamic coefficient of friction temperature multiplier fdym
      call value (nosl,lcfdt,npc,pld,tavg,tavg,fdym,valp,1)
      endif
c------------------------- end thermal ---------------------------------
c---------------------- begin mechanical -------------------------------
c
c     current & accumulated sliding distance
      ddist=sqrt(dx(ii)**2+dy(ii)**2+dz(ii)**2)
      ufhist(1,ilbsv(ii))=ufhist(1,ilbsv(ii))+ddist
      dist=ufhist(1,ilbsv(ii))
c
c     current sliding velocity
      vel=fdt2*ddist
c
c     current contact pressure = normal force / slave node area
      press=-fni(ii)/areas(ii)
c
c     master and slave node area
c     (needed for frict. force limit, see parameter VC in the manual)
      aream(ii)=aream(ii)*fric4*frc(2,lsv(ii))+
     .          1.e+16*(.5-sign(.5,fric4-1.e-10))
      areas(ii)=areas(ii)*fric4*frc(2,lsv(ii))+
     .          1.e+16*(.5-sign(.5,fric4-1.e-10))+bignum
      areas(ii)=areas(ii)+1.e+16*(.5-sign(.5,areas(ii)-1.e-10))
c
c     new contact force = old + current penalty force
      fdat(3,ii)=fdat(3,ii)+stfk(ii)*dx(ii)*sfacn3
      fdat(4,ii)=fdat(4,ii)+stfk(ii)*dy(ii)*sfacn3
      fdat(5,ii)=fdat(5,ii)+stfk(ii)*dz(ii)*sfacn3
c
c     tangential part = total - normal part
      proj=fdat(3,ii)*rn1(ii)+fdat(4,ii)*rn2(ii)+fdat(5,ii)*rn3(ii)
      fdat(3,ii)=fdat(3,ii)-proj*rn1(ii)
      fdat(4,ii)=fdat(4,ii)-proj*rn2(ii)
      fdat(5,ii)=fdat(5,ii)-proj*rn3(ii)
c
c     scalar tangential force = func(fric.coeff.) * normal force
      fmax=-(fdym*fric2+(fstm*fric1-fdym*fric2)*exp(-fric3*vel))
     .     *fni(ii)*frc(1,lsv(ii))
c
c     limit the friction force
      fmax=min(fmax,aream(ii),areas(ii))
      fmag=sqrt(fdat(3,ii)**2+fdat(4,ii)**2+fdat(5,ii)**2)+1.e-06
c
c     tangential force due to friction
      if (fmag.gt.fmax .and. fmag.ne.0.0) then
        sclf=fmax/fmag
        fdat(3,ii)=sclf*fdat(3,ii)
        fdat(4,ii)=sclf*fdat(4,ii)
        fdat(5,ii)=sclf*fdat(5,ii)
      endif
c
c     slave node force components in global directions
      fxis(ii)=fxis(ii)+fdat(3,ii)
      fyis(ii)=fyis(ii)+fdat(4,ii)
      fzis(ii)=fzis(ii)+fdat(5,ii)
c
      iseg(ii)  =lsv(ii)
      fdat(1,ii)=ss(ii)
      fdat(2,ii)=tt(ii)
c----------------------- end mechanical --------------------------------
c
      enddo
c
      return
      end
      subroutine useradap(eso,b,maxsiz64)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      integer
     1 nh01,nh02,nh03,nh04,nh05,nh06,nh07,nh08,nh09,nh10,
     2 nb01,nb02,nb03,nb04,nb05,nb06,nb07,nb08,nb09,nb10,
     3 ns01,ns02,ns03,ns04,ns05,ns06,ns07,ns08,ns09,ns10,
     4 nt01,nt02,nt03,nt04,nt05,nt06,nt07,nt08,nt09,nt10,
     5 nu01,nu02,nu03,nu04,nu05,nu06,nu07,nu08,nu09,nu10,
     6 nd01,nd02,nd03,nd04,nd05,nd06,nd07,nd08,nd09,nd10,
     7 ntbqs
      common/bk05/
     1 nh01,nh02,nh03,nh04,nh05,nh06,nh07,nh08,nh09,nh10,
     2 nb01,nb02,nb03,nb04,nb05,nb06,nb07,nb08,nb09,nb10,
     3 ns01,ns02,ns03,ns04,ns05,ns06,ns07,ns08,ns09,ns10,
     4 nt01,nt02,nt03,nt04,nt05,nt06,nt07,nt08,nt09,nt10,
     5 nu01,nu02,nu03,nu04,nu05,nu06,nu07,nu08,nu09,nu10,
     6 nd01,nd02,nd03,nd04,nd05,nd06,nd07,nd08,nd09,nd10,
     7 ntbqs
c
      integer*8 dm_hmtnum,dm_hnfegp
      integer*8 dm_bmtnum,dm_bnfegp
      integer*8 dm_smtnum,dm_snfegp
      integer*8 dm_tmtnum,dm_tnfegp
c
      common/dmbk05/
     1 dm_hmtnum,dm_hnfegp,
     2 dm_bmtnum,dm_bnfegp,
     3 dm_smtnum,dm_snfegp,
     4 dm_tmtnum,dm_tnfegp
      integer n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      common/bk07/n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      integer n4a,n4b,n4c,n4d,n4e,n4f,n4g,n4h,n7a,n7b,n7c,nusir,
     1 mpusr,mpubr,n4i,n4j,n4k,n4l
      common/bk08/n4a,n4b,n4c,n4d,n4e,n4f,n4g,n4h,n7a,n7b,n7c,nusir,
     1 mpusr,mpubr,n4i,n4j,n4k,n4l
      integer lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      common/bk13/lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
c
c     dynamic memory allocation stuff.  This has to be a C style 
c     include because of the use of integer*8
c
      integer i_mem(0:1)
      real r_mem(0:1)
      integer*8 i8_mem(0:1)
      real*8 r8_mem(0:1)
      real*8 real8_mem
      integer*8 m_mem(0:1)
      common/dynmem/real8_mem(0:1)
      equivalence (i_mem,r_mem,i8_mem,r8_mem,real8_mem,m_mem)
      integer*8 mem_alloc, mems_alloc,mem_realloc, memh_ptr
      integer mem_length,memh_alloc,memh_realloc,memh_length
      integer memh_enum,memh_enums,memh_attach
c      integer mem_newmark,mem_getmark,memh_getmark
      integer*8 memh_detach,memh_ralloc
      integer*8 m_to_m8
      integer*8 m_to_mm
      external mem_alloc, mems_alloc, mem_realloc, mem_length
      external memh_alloc, memh_realloc, memh_length, memh_ptr
      external memh_enum,memh_attach,memh_detach,memh_ralloc
      external memh_enums
      external m_to_m8
      external m_to_mm
      integer*8 dm_x,dm_xr,dm_v,dm_vr,dm_a,dm_ar,dm_rots,dm_rotr
      integer*8 dm_xms,dm_xmr,dm_xm2,dm_xmz,dm_xma
      integer*8 dm_xtpz,dm_disp
      integer*8 dm_me1,dm_me2,dm_me3
      integer*8 dm_masp,dm_masr
      integer*8 dm_xrb,dm_yrb,dm_zrb
      integer*8 dm_axrb,dm_ayrb,dm_azrb
      integer*8 dm_rbfx,dm_rbfy,dm_rbfz
      integer*8 dm_n45,dm_n46,dm_n47
      integer len_n45,len_n46,len_n47
      common/dynmem1/ 
     . dm_x,dm_xr,dm_v,dm_vr,dm_a,dm_ar,dm_rots,dm_rotr,
     . dm_xms,dm_xmr,dm_xm2,dm_xmz,dm_xma,
     . dm_xtpz,dm_disp,
     . dm_me1,dm_me2,dm_me3,
     . dm_masp,dm_masr,
     . dm_xrb,dm_yrb,dm_zrb,
     . dm_axrb,dm_ayrb,dm_azrb,
     . dm_rbfx,dm_rbfy,dm_rbfz,
     . dm_n45,dm_n46,dm_n47
      common/dynmem1l/len_n45,len_n46,len_n47
c
      integer*8 mem_allocm,mems_allocm,mem_reallocm
      integer memh_allocm,memh_reallocm
      external mem_allocm,mems_allocm,mem_reallocm
      external memh_allocm,memh_reallocm
      integer*8 mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      integer*8 mem_reallocmchk,mem_reallocchk
      integer memh_allocmchk,memh_allocchk
      external mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      external mem_reallocmchk,mem_reallocchk
      external memh_allocmchk,memh_allocchk
c
c NOTE: IF YOU CHANGE THIS FILE PLEASE ALSO CHANGE sorter_noloc.inc
c WHICH NEEDS TO BE IDENTICAL EXCEPT FOR nnc->znnc, lczc->zlczc
c
      integer nnc,lczc,
     &    ns11, ns12, ns13, ns14, ns15, ns16,
     &    nh11, nh12, nh13, nh14, nh15, nh16,
     &    nt11, nt12, nt13, nt14, nt15, nt16,
     &    nb11, nb12, nb13, nb14, nb15, nb16,
     &    nu11, nu12, nu13, nu14, nu15, nu16,
     &    nd11, nd12, nd13, nd14, nd15, nd16,
     &    nd17, nd18, nd19, nd20, nd21, nd22,
     &    ncf26, ncf27, ncf28,
     &    n_em_26, n_em_27, n_em_28,
     &    nscf08, nscf09, nscf10, nscf13, nscf14,
     &    nhcf11, nhcf12, nhcf13, nhcf14, nhcf15,
     &    nh_em_11, nh_em_12, nh_em_13, nh_em_14, nh_em_15,
     &    lccf1h, lccf1s, lc_em_1h
      common/sorter/nnc, lczc,
     &    ns11, ns12, ns13, ns14, ns15, ns16,
     &    nh11, nh12, nh13, nh14, nh15, nh16,
     &    nt11, nt12, nt13, nt14, nt15, nt16,
     &    nb11, nb12, nb13, nb14, nb15, nb16,
     &    nu11, nu12, nu13, nu14, nu15, nu16,
     &    nd11, nd12, nd13, nd14, nd15, nd16,
     &    nd17, nd18, nd19, nd20, nd21, nd22,
     &    ncf26, ncf27, ncf28,
     &    n_em_26, n_em_27, n_em_28,
     &    nscf08, nscf09, nscf10, nscf13, nscf14,
     &    nhcf11, nhcf12, nhcf13, nhcf14, nhcf15,
     &    nh_em_11, nh_em_12, nh_em_13, nh_em_14, nh_em_15,
     &    lccf1h, lccf1s, lc_em_1h
c
      integer*8 dm_hnsubgc,dm_hnfegc
      integer*8 dm_bnsubgc,dm_bnfegc
      integer*8 dm_snsubgc,dm_snfegc
      integer*8 dm_tnsubgc,dm_tnfegc
c
      common/dmsorter/
     1 dm_hnsubgc,dm_hnfegc,
     2 dm_bnsubgc,dm_bnfegc,
     3 dm_snsubgc,dm_snfegc,
     4 dm_tnsubgc,dm_tnfegc
c
c The size of the sorter common block.  Used in the restart and dump
c routines.
c
      integer ISORTERSIZE
      parameter (ISORTERSIZE = 68)
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
      integer*8 maxsiz64
c
c
      dimension eso(7,*),b(1)
c
      call userrefin (b(n1),b(n4f),r_mem(dm_x),i_mem(dm_smtnum),b(lc1s),
     +  b(ns06),b(n1+nmmat),b(ns13),b(ns14),eso,b(n4d+3*nmmat),maxsiz64)
c
      return
      end
      subroutine userrefin (mtype,csprop,x,mtnum,ixp,auxvec,
     + ishlfm,nshpnt,lochvs,eso,nconstp,maxsiz64)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'bk06.inc'
      include 'bk19.inc'
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
c     common/bk28/summss,xke,xpe,time
      common/daimlr/noddie,neldie,ncunk,npersv,ivtmem
      common/ellcm2/nelnt,nelet,nstart,nelbt
      common/xcntrl/kontrl(200)
c
      real*8 x
      dimension mtype(*),csprop(24,*),x(3,*),mtnum(*),ixp(5,*),
     +          auxvec(*),ishlfm(*),nshpnt(*),lochvs(*),eso(7,*),
     +          stress(7,20),nconstp(*)
      dimension output(7)
      integer   iob28
      integer*8 maxsiz64
      integer*8 izero64,jadd64
      data izero64/0/
      data init/0/
      data output/7*0.0/
c
c     Variables in the subroutine:
c
c        numnp=number of nodal points
c        numels-neldie=number of shell elements
c
c        x(1,i)=x-coordinate of nodal point i
c        x(2,i)=y-coordinate of nodal point i
c        x(3,i)=z-coordinate of nodal point i
c
c
c        ixp - shell element connectivities (ixp(1,*)=part ID)
c                                           (ixp(2,*)=node 1)
c                                           (ixp(3,*)=node 2)
c                                           (ixp(4,*)=node 3)
c                                           (ixp(5,*)=node 4)
c
c        ipt=integration point number
c        nip=number of integration points
c
c        stress(1,ipt)=local x  stress
c        stress(2,ipt)=local y  stress
c        stress(3,ipt)=local z  stress
c        stress(4,ipt)=local xy stress
c        stress(5,ipt)=local yz stress
c        stress(6,ipt)=local zx stress
c        stress(7,ipt)=effective plastic strain
c
c     Write file with data for useradaptivity:
c
c     structure of file:
c     line 1: Number of nodes, number of elements (shell)
c     part 2: Nodal points (numnp lines)
c     part 3: Material # + Element connectivities (numels lines)
c     part 4: Element materials and stresses (numels*nip lines)
c             Material #,
c             Integration point #
c             local sigma_xx
c             local sigma_yy
c             local sigma_zz
c             local sigma_xy
c             local sigma_yz
c             local sigma_zx
c             effective plastic strain
c
      call lsopen(21,'adapt.user','formatted','unknown',ierr,0,0)
      write(21,900) numnp,numels-neldie
      write(21,915) (x(1,n),x(2,n),x(3,n),n=1,numnp)
      write(21,905) (ixp(1,nshpnt(n)),ixp(2,nshpnt(n)),
     +               ixp(3,nshpnt(n)),ixp(4,nshpnt(n)),
     +               ixp(5,nshpnt(n)),n=1,numels-neldie)
c
      do 100 n=1,numels-neldie
      nel =nshpnt(n)
      lav =lochvs(nel)
      mx=ixp(1,nel)
      if (mx.eq.0) go to 100
      iop=ishlfm(mx)
      nip=csprop(2,mx)
      mte=mtype(mx)
      nmtcon=7+nconstp(mx)
c
      do 50 ipt=1,nip
      if (mte.eq.20) then
        stress(1,ipt)=0.
        stress(2,ipt)=0.
        stress(3,ipt)=0.
        stress(4,ipt)=0.
        stress(5,ipt)=0.
        stress(6,ipt)=0.
        stress(7,ipt)=0.
      else
        call blkcpy (auxvec(lav),stress(1,ipt),7)
      endif
c
      write(21,540) mx,ipt,(stress(j,ipt),j=1,7)
      lav=lav+nmtcon
   50 continue
  100 continue
      close(21)
c
c     External program for caclulating refinement indicators
c     can be called here
c
c     call my_system('External program')
c
c     Write adapt.err file for userspecified adaptivity
c
c     The value of enerden specify
c     which element which will be refined
c
c     if enerden < 1 no refinement of element
c     if enerden > 1 refinement of element
c
c
c
c     Or refinement indicators can be calculated
c     inside the next loop
c
c     call llsopen(28,'adapt.err','unformatted','unknown',ierr,0,0)
c     rewind 28
      call rwabsf64(iob28,'adapt.err',28,2**9,maxsiz64)
      nn=numels-neldie-nelet
      jadd64=0
      call wrabsf64(iob28,nn,1,jadd64)
      jadd64=1
c
      enerden=0.
      do 200 n=1,numels-neldie-nelet
c
c     calculate refinement for each element based
c     on the data in the "adapt.user" file
c
c     calculate the value of enerden here
c     and the value of enerden will be written to
c     the file "adapt.err"
c
      output(1)=enerden
      if (init .eq. 0) then
c       write(28) enerden,0.0,0.0,enerden,0.0,0.0,0.0
        output(4)=enerden
      else
c       write(28) enerden,0.0,0.0,eso(1,n),0.0,0.0,0.0
        output(4)=eso(1,n)
      endif
      call wrabsf64(iob28,output,7,jadd64)
      jadd=jadd+7
  200 continue
      kontrl(42)=jadd
      init = init+1
c     close(28)
      jadd64=0
      call rwabsf64(iob28,'keep',0,0,jadd64)
c
      return
  540 format(2i5,7(1pe15.7))
  900 format(10i8)
  905 format(5i10)
c 910 format(6(1pe12.5))
  915 format(3(1pe15.7))
      end
      subroutine usrmat (lft,llt,cm,bqs,capa,eltype,mt,ipt,
     . npc,plc,crv,rcoor,scoor,tcoor,nnm1,nip,ipt_thk)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      integer n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      common/bk07/n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      integer n4a,n4b,n4c,n4d,n4e,n4f,n4g,n4h,n7a,n7b,n7c,nusir,
     1 mpusr,mpubr,n4i,n4j,n4k,n4l
      common/bk08/n4a,n4b,n4c,n4d,n4e,n4f,n4g,n4h,n7a,n7b,n7c,nusir,
     1 mpusr,mpubr,n4i,n4j,n4k,n4l
      include 'memaia.inc'
      include 'umatss.inc'
c
c NOTE: IF YOU CHANGE THIS FILE PLEASE ALSO CHANGE sorter_noloc.inc
c WHICH NEEDS TO BE IDENTICAL EXCEPT FOR nnc->znnc, lczc->zlczc
c
      integer nnc,lczc,
     &    ns11, ns12, ns13, ns14, ns15, ns16,
     &    nh11, nh12, nh13, nh14, nh15, nh16,
     &    nt11, nt12, nt13, nt14, nt15, nt16,
     &    nb11, nb12, nb13, nb14, nb15, nb16,
     &    nu11, nu12, nu13, nu14, nu15, nu16,
     &    nd11, nd12, nd13, nd14, nd15, nd16,
     &    nd17, nd18, nd19, nd20, nd21, nd22,
     &    ncf26, ncf27, ncf28,
     &    n_em_26, n_em_27, n_em_28,
     &    nscf08, nscf09, nscf10, nscf13, nscf14,
     &    nhcf11, nhcf12, nhcf13, nhcf14, nhcf15,
     &    nh_em_11, nh_em_12, nh_em_13, nh_em_14, nh_em_15,
     &    lccf1h, lccf1s, lc_em_1h
      common/sorter/nnc, lczc,
     &    ns11, ns12, ns13, ns14, ns15, ns16,
     &    nh11, nh12, nh13, nh14, nh15, nh16,
     &    nt11, nt12, nt13, nt14, nt15, nt16,
     &    nb11, nb12, nb13, nb14, nb15, nb16,
     &    nu11, nu12, nu13, nu14, nu15, nu16,
     &    nd11, nd12, nd13, nd14, nd15, nd16,
     &    nd17, nd18, nd19, nd20, nd21, nd22,
     &    ncf26, ncf27, ncf28,
     &    n_em_26, n_em_27, n_em_28,
     &    nscf08, nscf09, nscf10, nscf13, nscf14,
     &    nhcf11, nhcf12, nhcf13, nhcf14, nhcf15,
     &    nh_em_11, nh_em_12, nh_em_13, nh_em_14, nh_em_15,
     &    lccf1h, lccf1s, lc_em_1h
c
      integer*8 dm_hnsubgc,dm_hnfegc
      integer*8 dm_bnsubgc,dm_bnfegc
      integer*8 dm_snsubgc,dm_snfegc
      integer*8 dm_tnsubgc,dm_tnfegc
c
      common/dmsorter/
     1 dm_hnsubgc,dm_hnfegc,
     2 dm_bnsubgc,dm_bnfegc,
     3 dm_snsubgc,dm_snfegc,
     4 dm_tnsubgc,dm_tnfegc
c
c The size of the sorter common block.  Used in the restart and dump
c routines.
c
      integer ISORTERSIZE
      parameter (ISORTERSIZE = 68)
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     1 ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     2 grvity,idirgv,nodspc,nspcor,numelh10,numels8,numelsx,numnpin,
     3 numelhx,nmuxct
      dimension cm(*),bqs(*),crv(lq1,2,*),npc(*),plc(*)
      character*(*) eltype
c
      if (eltype.eq.'solid'.or.eltype.eq.'shl_t'.or.
     1     eltype.eq.'sld2d'.or.eltype.eq.'tshel') then
        call urmathn (lft,llt,cm,bqs,mt,crv,npc,plc,nnm1,
     1   rcoor,scoor,tcoor,ia(n4d+3*nmmat),nip,ipt,eltype,
     2   a(nh13+numelh),a(ns13+numels),a(nt13+numelt))
      elseif (eltype.eq.'shell') then
        call urmats (lft,llt,cm,capa,mt,crv,ipt,rcoor,scoor,tcoor,
     1   nnm1,ia(n4d+3*nmmat),nip,ipt_thk,npc,plc,eltype,
     2   a(ns13+numels))
      elseif (eltype.eq.'hbeam ') then
        call urmatb (lft,llt,cm,capa,mt,crv,nnm1,ia(n4d+3*nmmat),
     1   npc,plc,eltype,a(nb13+numelb))
      elseif (eltype.eq.'dbeam') then
        call urmatd (lft,llt,cm,bqs,mt,crv,nnm1,ia(n4d+3*nmmat),
     1   npc,plc,eltype,a(nb13+numelb))
      elseif (eltype.eq.'tbeam') then
        call urmatt (lft,llt,cm,bqs,mt,crv,nnm1,ia(n4d+3*nmmat),
     1   npc,plc,eltype,a(nb13+numelb))
      endif
c
      return
      end
      subroutine urmathn(lft,llt,cm,bqs,mt,crv,npc,plc,nnm1,
     1 rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nhxbwp,nshbwp,ntsbwp)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c                                                                 |
c   modified by Xiaoyan Zhang based on Kay Sun                    |
c   7/2014                                                        |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
      integer n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      common/bk07/n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      integer lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      common/bk13/lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      include 'bk19.inc'
      include 'iounits.inc'
      include 'memaia.inc'
      include 'umatss.inc'
c
c     dynamic memory allocation stuff.  This has to be a C style 
c     include because of the use of integer*8
c
      integer i_mem(0:1)
      real r_mem(0:1)
      integer*8 i8_mem(0:1)
      real*8 r8_mem(0:1)
      real*8 real8_mem
      integer*8 m_mem(0:1)
      common/dynmem/real8_mem(0:1)
      equivalence (i_mem,r_mem,i8_mem,r8_mem,real8_mem,m_mem)
      integer*8 mem_alloc, mems_alloc,mem_realloc, memh_ptr
      integer mem_length,memh_alloc,memh_realloc,memh_length
      integer memh_enum,memh_enums,memh_attach
c      integer mem_newmark,mem_getmark,memh_getmark
      integer*8 memh_detach,memh_ralloc
      integer*8 m_to_m8
      integer*8 m_to_mm
      external mem_alloc, mems_alloc, mem_realloc, mem_length
      external memh_alloc, memh_realloc, memh_length, memh_ptr
      external memh_enum,memh_attach,memh_detach,memh_ralloc
      external memh_enums
      external m_to_m8
      external m_to_mm
      integer*8 dm_x,dm_xr,dm_v,dm_vr,dm_a,dm_ar,dm_rots,dm_rotr
      integer*8 dm_xms,dm_xmr,dm_xm2,dm_xmz,dm_xma
      integer*8 dm_xtpz,dm_disp
      integer*8 dm_me1,dm_me2,dm_me3
      integer*8 dm_masp,dm_masr
      integer*8 dm_xrb,dm_yrb,dm_zrb
      integer*8 dm_axrb,dm_ayrb,dm_azrb
      integer*8 dm_rbfx,dm_rbfy,dm_rbfz
      integer*8 dm_n45,dm_n46,dm_n47
      integer len_n45,len_n46,len_n47
      common/dynmem1/ 
     . dm_x,dm_xr,dm_v,dm_vr,dm_a,dm_ar,dm_rots,dm_rotr,
     . dm_xms,dm_xmr,dm_xm2,dm_xmz,dm_xma,
     . dm_xtpz,dm_disp,
     . dm_me1,dm_me2,dm_me3,
     . dm_masp,dm_masr,
     . dm_xrb,dm_yrb,dm_zrb,
     . dm_axrb,dm_ayrb,dm_azrb,
     . dm_rbfx,dm_rbfy,dm_rbfz,
     . dm_n45,dm_n46,dm_n47
      common/dynmem1l/len_n45,len_n46,len_n47
c
      integer*8 mem_allocm,mems_allocm,mem_reallocm
      integer memh_allocm,memh_reallocm
      external mem_allocm,mems_allocm,mem_reallocm
      external memh_allocm,memh_reallocm
      integer*8 mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      integer*8 mem_reallocmchk,mem_reallocchk
      integer memh_allocmchk,memh_allocchk
      external mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      external mem_reallocmchk,mem_reallocchk
      external memh_allocmchk,memh_allocchk
      common/bux13loc/ eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      common/aux14loc/sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),sig5(nlq),
     & sig6(nlq),epsps(nlq),hsvs(nlq,NHISVAR)
      common/aux18loc/dd(nlq),def(nlq),ddq(nlq)
      common/aux2loc/d1(nlq),d2(nlq),d3(nlq),d4(nlq),d5(nlq),d6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common/aux33loc/ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ix5(nlq),
     1 ix6(nlq),ix7(nlq),ix8(nlq),mxt(nlq)
      common/aux34loc/xioff(nlq),dx(nlq)
      common/aux35loc/rhoa(nlq),cb(nlq),davg(nlq),p(nlq)
      common/aux40loc/
     1 a11(nlq),a12(nlq),a13(nlq),a21(nlq),a22(nlq),a23(nlq),
     2 a31(nlq),a32(nlq),a33(nlq),z11(nlq),z12(nlq),z13(nlq),
     3 z21(nlq),z22(nlq),z23(nlq),z31(nlq),z32(nlq),z33(nlq)
      common/aux41loc/qq1(nlq),cbb(nlq),aj1(nlq)
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
      common/bk36loc/index
      common/soundloc/ss(nlq),sndsp(nlq),diagm(nlq),sarea(nlq),
     . dxl(nlq)
      common/subtssloc/dt1siz(nlq)
      common/vect32loc/g11(nlq),g21(nlq),g31(nlq),g12(nlq),g22(nlq),
     .     g32(nlq),g13(nlq),g23(nlq),g33(nlq),ifg(nlq)
c
      common/bk01/itherm,itemp,ntmp0,ntmp1,itempan,itempdr,itmpe(44)
      common/bk26/begtim,nintcy
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie
      logical failur
      common/failcm/failur
      common/failcmloc/ifail(nlq)
      common/failuloc/sieu(nlq),fail(nlq),ifaili(nlq),nfipts(nlq)
      common/numcpu/ncpu,ncpua,ncpub,lenvec(8)
      common/thropt/itopaz(101)
c
      dimension cm(*),bqs(*),crv(lq1,2,*),npc(*),plc(*),nconstp(*)
      character*(*) eltype
      dimension nhxbwp(*),nshbwp(*),ntsbwp(*)
c     un-commented the following "dimension" and hsv(NHISVAR) changed to hsv(100)
c      dimension eps(6),sig(6),hsv(100),temps(nlq)
      dimension q21(nlq),q22(nlq),q23(nlq)
      dimension thhsv(nlq,100),thhsvi(100)
      dimension elsizv(nlq),idelev(nlq)
c
      real mlt1,mlt2
      logical failel,failels(nlq)
c     data temps/nlq*0.0/
c     added the following 3 lines
      data temps/nlq*0.0/
      real vrat, Ebk
      dimension xM(3),xS(3),xN(3)
c
c
c$omp threadprivate (/aux2loc/)
c$omp threadprivate (/aux14loc/)
c$omp threadprivate (/aux18loc/)
c$omp threadprivate (/aux33loc/)
c$omp threadprivate (/aux34loc/)
c$omp threadprivate (/aux35loc/)
c$omp threadprivate (/aux40loc/)
c$omp threadprivate (/aux41loc/)
c$omp threadprivate (/bk36loc/)
c$omp threadprivate (/failcmloc/)
c$omp threadprivate (/failuloc/)
c$omp threadprivate (/soundloc/)
c$omp threadprivate (/subtssloc/)
c$omp threadprivate (/vect32loc/)
c$omp threadprivate (/bux13loc/)
      capa=0.0
c     added Ebk
      Ebk=0.0
c
c     processing elements
c
      mx=48*(mxt(lft)-1)
      gm=cm(mx+abs(ishrmp(mt)))
      bk=cm(mx+ibulkp(mt))
      ss(lft)=bk+4.*gm/3.
c
c     energy calculations
c
      do 20 i=lft,llt
      cb(i)=ss(lft)
      einc(i)=(d1(i)*sig1(i)+d2(i)*sig2(i)+d3(i)*sig3(i)+
     .        d4(i)*sig4(i)+d5(i)*sig5(i)+
     .        d6(i)*sig6(i)+dd(i)*bqs(i))*dt1siz(i)
 20   continue
c
c     Determining number of requested history variables
c
      no_hsvs=nconstp(mxt(lft))
      if (ivumat(mt).ne.0) no_hsvs=no_hsvs-9-6
      if (lenvec(5) .eq.1) no_hsvs=no_hsvs-6
      if (iorien(mt).ne.0) no_hsvs=no_hsvs-6
      if (ihyper(mt).ne.0) no_hsvs=no_hsvs-9
      if (iorien(mt).ne.0.and.ihyper(mt).lt.0) no_hsvs=no_hsvs-6
      if (ifipts(mt).lt.0) no_hsvs=no_hsvs-1
c
c     initialization of indeces
c
      if(ivumat(mt).ne.0) then
        if (iorien(mt).ne.0) then
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_ortho=ind_defgrad+9
          else
            ind_ortho=1+no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_ortho
            ind_ortho=ind_ortho+1
          endif
          ind_q11=ind_ortho
          ind_q12=ind_ortho+1
          ind_q13=ind_ortho+2
          ind_q31=ind_ortho+3
          ind_q32=ind_ortho+4
          ind_q33=ind_ortho+5
          ind_last=ind_q33
          if (ihyper(mt).lt.0) then
            ind_qstore=ind_last+1
            ind_last=ind_last+6
          endif
        else
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_last=ind_defgrad+8
          else
            ind_last=no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_last+1
            ind_last=ind_last+1
          endif
        endif
        ind_vumat=ind_last+1
        ind_last=ind_vumat+9+6-1       
c        ind_last=ind_last+1
c        write(*,'(A,i5,A,i5,A,i5,a,i5,a,i5,a,i5)')'NHSV:',no_hsvs,
c     $       'DF',ind_defgrad,'ORTH',ind_ortho,'VMAT',ind_vumat,
c     $       'LST',ind_last,'QSTO',ind_qstore
      else
        if (iorien(mt).ne.0) then
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_ortho=ind_defgrad+9
          else
            ind_ortho=1+no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_ortho
            ind_ortho=ind_ortho+1
          endif
          ind_q11=ind_ortho
          ind_q12=ind_ortho+1
          ind_q13=ind_ortho+2
          ind_q31=ind_ortho+3
          ind_q32=ind_ortho+4
          ind_q33=ind_ortho+5
          ind_last=ind_q33
          if (ihyper(mt).lt.0) then
            ind_qstore=ind_last+1
            ind_last=ind_last+6
          endif
        else
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_last=ind_defgrad+8
          else
            ind_last=no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_last+1
            ind_last=ind_last+1
          endif
        endif
      endif
c
      if (ind_last.gt.NHISVAR) then
c       write(iotty,900) lqfmiv(mxt(lft))
c       write(iohsp,900) lqfmiv(mxt(lft))
c       write(iomsg,900) lqfmiv(mxt(lft))
c       call adios(2)
        ierdat(1)=lqfmiv(mxt(lft))
        call lsmsg(3,MSG_SOL+1149,ioall,ierdat,rerdat,cerdat,0)
      endif
c
c$$$      if (tt.eq.begtim.and.nintcy.eq.0) then
c$$$c
c$$$c       move transformation matrix in history variables array
c$$$c
c$$$        if (iorien(mt).ne.0) then
c$$$          do i=lft,llt
c$$$            hsvs(i,ind_q33)=hsvs(i,6)
c$$$            hsvs(i,ind_q32)=hsvs(i,5)
c$$$            hsvs(i,ind_q31)=hsvs(i,4)
c$$$            hsvs(i,ind_q13)=hsvs(i,3)
c$$$            hsvs(i,ind_q12)=hsvs(i,2)
c$$$            hsvs(i,ind_q11)=hsvs(i,1)
c$$$          enddo
c$$$c
c$$$          if (no_hsvs.gt.0) then
c$$$            do j=1,no_hsvs
c$$$              do i=lft,llt
c$$$                hsvs(i,j)=0.0
c$$$              enddo
c$$$            enddo
c$$$          endif
c$$$        endif
c$$$      endif
c
      if (ipt.eq.1.and.ifipts(mt).ne.0) then
        do i=lft,llt
          nfipts(i)=0
        enddo
      endif
c
c     storage of deformation gradient
c
      if (ivumat(mt).ne.0) then
        ifo11=ind_vumat
        ifo21=ind_vumat+1
        ifo31=ind_vumat+2
        ifo12=ind_vumat+3
        ifo22=ind_vumat+4
        ifo32=ind_vumat+5
        ifo13=ind_vumat+6
        ifo23=ind_vumat+7
        ifo33=ind_vumat+8
      endif
      if (ihyper(mt).ne.0) then
        if11=ind_defgrad
        if21=ind_defgrad+1
        if31=ind_defgrad+2
        if12=ind_defgrad+3
        if22=ind_defgrad+4
        if32=ind_defgrad+5
        if13=ind_defgrad+6
        if23=ind_defgrad+7
        if33=ind_defgrad+8
c
c
c       deformation gradient in global system at time = 0
c
        if (tt.eq.begtim.and.nintcy.eq.0) then
          do i=lft,llt
            hsvs(i,if11)=1.0
            hsvs(i,if21)=0.0
            hsvs(i,if31)=0.0
            hsvs(i,if12)=0.0
            hsvs(i,if22)=1.0
            hsvs(i,if32)=0.0
            hsvs(i,if13)=0.0
            hsvs(i,if23)=0.0
            hsvs(i,if33)=1.0
          enddo
          if (ivumat(mt).ne.0) then
            do i=lft,llt
              hsvs(i,ifo11)=1.0
              hsvs(i,ifo21)=0.0
              hsvs(i,ifo31)=0.0
              hsvs(i,ifo12)=0.0
              hsvs(i,ifo22)=1.0
              hsvs(i,ifo32)=0.0
              hsvs(i,ifo13)=0.0
              hsvs(i,ifo23)=0.0
              hsvs(i,ifo33)=1.0
            enddo
          endif
        endif
c
c       deformation gradient calculated exactly
c
        if (ifg(lft).eq.1) then
          do i=lft,llt
            hsvs(i,if11)=g11(i)
            hsvs(i,if21)=g21(i)
            hsvs(i,if31)=g31(i)
            hsvs(i,if12)=g12(i)
            hsvs(i,if22)=g22(i)
            hsvs(i,if32)=g32(i)
            hsvs(i,if13)=g13(i)
            hsvs(i,if23)=g23(i)
            hsvs(i,if33)=g33(i)
          enddo
c
c
c       update deformation gradient incrementally in global system
c       at time = n+1
c
        else
          call integf(hsvs(1,if11),hsvs(1,if21),hsvs(1,if31),
     &     hsvs(1,if12),hsvs(1,if22),hsvs(1,if32),
     &     hsvs(1,if13),hsvs(1,if23),hsvs(1,if33),lft,llt)
        endif
      endif
c
      if (iorien(mt).ne.0) then
c
c       set up global to element transformation matrix a11,...,a33
c
        if (ihyper(mt).ge.0.and.eltype.ne.'tshel') then
c        if (ihyper(mt).eq.1 .or. ihyper(mt).eq.0) then
          if (eltype.ne.'sld2d') then
            call lcsm22(lft,llt)
          else
            call cs2d_g2l(lft,llt)
          endif
        endif
c     
        do 10 i=lft,llt
          q21(i)=hsvs(i,ind_q32)*hsvs(i,ind_q13)-
     1         hsvs(i,ind_q12)*hsvs(i,ind_q33)
          q22(i)=hsvs(i,ind_q33)*hsvs(i,ind_q11)-
     1         hsvs(i,ind_q13)*hsvs(i,ind_q31)
          q23(i)=hsvs(i,ind_q31)*hsvs(i,ind_q12)-
     1         hsvs(i,ind_q11)*hsvs(i,ind_q32)
 10     continue
c        write(*,*)'Q21,22,23',q21(1),q22(1),q23(1)
c     
c       set up global to material transformation matrix z11,...,z33
c
        if (ihyper(mt).ge.0) then
c     if (ihyper(mt).eq.1 .or. ihyper(mt).eq.0) then
          if (eltype.ne.'tshel') then
            do i=lft,llt
              z11(i)=a11(i)*hsvs(i,ind_q11)+a21(i)*hsvs(i,ind_q12)+
     1             a31(i)*hsvs(i,ind_q13)
              z12(i)=a12(i)*hsvs(i,ind_q11)+a22(i)*hsvs(i,ind_q12)+
     1             a32(i)*hsvs(i,ind_q13)
              z13(i)=a13(i)*hsvs(i,ind_q11)+a23(i)*hsvs(i,ind_q12)+
     1             a33(i)*hsvs(i,ind_q13)
              z21(i)=a11(i)*q21(i)+a21(i)*q22(i)+a31(i)*q23(i)
              z22(i)=a12(i)*q21(i)+a22(i)*q22(i)+a32(i)*q23(i)
              z23(i)=a13(i)*q21(i)+a23(i)*q22(i)+a33(i)*q23(i)
              z31(i)=a11(i)*hsvs(i,ind_q31)+a21(i)*hsvs(i,ind_q32)+
     1             a31(i)*hsvs(i,ind_q33)
              z32(i)=a12(i)*hsvs(i,ind_q31)+a22(i)*hsvs(i,ind_q32)+
     1             a32(i)*hsvs(i,ind_q33)
              z33(i)=a13(i)*hsvs(i,ind_q31)+a23(i)*hsvs(i,ind_q32)+
     1             a33(i)*hsvs(i,ind_q33)
            enddo
          else
            do i=lft,llt
              z11(i)=hsvs(i,ind_q11)
              z12(i)=hsvs(i,ind_q12)
              z13(i)=hsvs(i,ind_q13)
              z21(i)=q21(i)
              z22(i)=q22(i)
              z23(i)=q23(i)
              z31(i)=hsvs(i,ind_q31)
              z32(i)=hsvs(i,ind_q32)
              z33(i)=hsvs(i,ind_q33)
            enddo
          endif
c     
c     transform stresses and strains to local system
c     
          call ldsm22(lft,llt)
          call ldem22(lft,llt)
c     
c     else
c        endif
      elseif (ihyper(mt).lt.0) then
c     if ( (ihyper(mt).lt.0) .or. (ihyper(mt).eq.2) ) then
c     
          if (tt.eq.begtim.and.nintcy.eq.0) then
c     
c           set up global to element transformation matrix a11,...,a33
c
            if (eltype.ne.'tshel') then
              if (eltype.ne.'sld2d') then
                call lcsm22(lft,llt)
              else
                call cs2d_g2l(lft,llt)
              endif
c
c             set up global to material transformation matrix z11,...,z33
c
              do i=lft,llt
                z11(i)=a11(i)*hsvs(i,ind_q11)+a21(i)*hsvs(i,ind_q12)+
     1               a31(i)*hsvs(i,ind_q13)
                z12(i)=a12(i)*hsvs(i,ind_q11)+a22(i)*hsvs(i,ind_q12)+
     1               a32(i)*hsvs(i,ind_q13)
                z13(i)=a13(i)*hsvs(i,ind_q11)+a23(i)*hsvs(i,ind_q12)+
     1               a33(i)*hsvs(i,ind_q13)
                z21(i)=a11(i)*q21(i)+a21(i)*q22(i)+a31(i)*q23(i)
                z22(i)=a12(i)*q21(i)+a22(i)*q22(i)+a32(i)*q23(i)
                z23(i)=a13(i)*q21(i)+a23(i)*q22(i)+a33(i)*q23(i)
                z31(i)=a11(i)*hsvs(i,ind_q31)+a21(i)*hsvs(i,ind_q32)+
     1               a31(i)*hsvs(i,ind_q33)
                z32(i)=a12(i)*hsvs(i,ind_q31)+a22(i)*hsvs(i,ind_q32)+
     1               a32(i)*hsvs(i,ind_q33)
                z33(i)=a13(i)*hsvs(i,ind_q31)+a23(i)*hsvs(i,ind_q32)+
     1               a33(i)*hsvs(i,ind_q33)
              enddo
            else
              do i=lft,llt
                z11(i)=hsvs(i,ind_q11)
                z12(i)=hsvs(i,ind_q12)
                z13(i)=hsvs(i,ind_q13)
                z21(i)=q21(i)
                z22(i)=q22(i)
                z23(i)=q23(i)
                z31(i)=hsvs(i,ind_q31)
                z32(i)=hsvs(i,ind_q32)
                z33(i)=hsvs(i,ind_q33)
              enddo
            endif
c            write(*,*)'A1:',a11(1),a22(1),a33(1)
c            write(*,*)'Z1:',z11(1),z22(1),z33(1)
c     
c     Store correct transformation used for transforming
c     the deformation gradient
c     
            do i=lft,llt
              hsvs(i,ind_qstore  )=hsvs(i,ind_q11)
              hsvs(i,ind_qstore+1)=hsvs(i,ind_q12)
              hsvs(i,ind_qstore+2)=hsvs(i,ind_q13)
              hsvs(i,ind_qstore+3)=hsvs(i,ind_q31)
              hsvs(i,ind_qstore+4)=hsvs(i,ind_q32)
              hsvs(i,ind_qstore+5)=hsvs(i,ind_q33)
              hsvs(i,ind_q11)=z11(i)
              hsvs(i,ind_q12)=z12(i)
              hsvs(i,ind_q13)=z13(i)
              hsvs(i,ind_q31)=z31(i)
              hsvs(i,ind_q32)=z32(i)
              hsvs(i,ind_q33)=z33(i)
            enddo
c     
c     set up element to material transformation matrix z11,...,z33
c     
          else
            do i=lft,llt
              z11(i)=hsvs(i,ind_q11)
              z12(i)=hsvs(i,ind_q12)
              z13(i)=hsvs(i,ind_q13)
              z21(i)=q21(i)
              z22(i)=q22(i)
              z23(i)=q23(i)
              z31(i)=hsvs(i,ind_q31)
              z32(i)=hsvs(i,ind_q32)
              z33(i)=hsvs(i,ind_q33)
            enddo
          endif
        endif
c     
c        deformation gradient transformed to local system
c        as F_bar=Q^T*F
c
         if (ihyper(mt).gt.0) then
c        if (ihyper(mt).eq.1) then
          do i=lft,llt
            f1=z11(i)*hsvs(i,if11)+
     1           z12(i)*hsvs(i,if21)+z13(i)*hsvs(i,if31)
            f2=z21(i)*hsvs(i,if11)+
     1           z22(i)*hsvs(i,if21)+z23(i)*hsvs(i,if31)
            hsvs(i,if31)=z31(i)*hsvs(i,if11)+
     1           z32(i)*hsvs(i,if21)+z33(i)*hsvs(i,if31)
            hsvs(i,if11)=f1
            hsvs(i,if21)=f2
            f1=z11(i)*hsvs(i,if12)+
     1           z12(i)*hsvs(i,if22)+z13(i)*hsvs(i,if32)
            f2=z21(i)*hsvs(i,if12)+
     1           z22(i)*hsvs(i,if22)+z23(i)*hsvs(i,if32)
            hsvs(i,if32)=z31(i)*hsvs(i,if12)+
     1           z32(i)*hsvs(i,if22)+z33(i)*hsvs(i,if32)
            hsvs(i,if12)=f1
            hsvs(i,if22)=f2
            f1=z11(i)*hsvs(i,if13)+
     1           z12(i)*hsvs(i,if23)+z13(i)*hsvs(i,if33)
            f2=z21(i)*hsvs(i,if13)+
     1           z22(i)*hsvs(i,if23)+z23(i)*hsvs(i,if33)
            hsvs(i,if33)=z31(i)*hsvs(i,if13)+
     1           z32(i)*hsvs(i,if23)+z33(i)*hsvs(i,if33)
            hsvs(i,if13)=f1
            hsvs(i,if23)=f2
          enddo
c
c     deformation gradient transformed to co-rotated coordinate system
c     as F_bar=Q^T*F*Q (suggested by Richard Hamm, P&G 2011)
        else if (ihyper(mt).eq.-2) then
          do i=lft,llt
            f11=z11(i)*hsvs(i,if11)+
     1           z12(i)*hsvs(i,if21)+z13(i)*hsvs(i,if31)
            f21=z21(i)*hsvs(i,if11)+
     1           z22(i)*hsvs(i,if21)+z23(i)*hsvs(i,if31)
            f31=z31(i)*hsvs(i,if11)+
     1           z32(i)*hsvs(i,if21)+z33(i)*hsvs(i,if31)
            f12=z11(i)*hsvs(i,if12)+
     1           z12(i)*hsvs(i,if22)+z13(i)*hsvs(i,if32)
            f22=z21(i)*hsvs(i,if12)+
     1           z22(i)*hsvs(i,if22)+z23(i)*hsvs(i,if32)
            f32=z31(i)*hsvs(i,if12)+
     1           z32(i)*hsvs(i,if22)+z33(i)*hsvs(i,if32)
            f13=z11(i)*hsvs(i,if13)+
     1           z12(i)*hsvs(i,if23)+z13(i)*hsvs(i,if33)
            f23=z21(i)*hsvs(i,if13)+
     1           z22(i)*hsvs(i,if23)+z23(i)*hsvs(i,if33)
            f33=z31(i)*hsvs(i,if13)+
     1           z32(i)*hsvs(i,if23)+z33(i)*hsvs(i,if33)
            hsvs(i,if11)=f11*z11(i)+f12*z12(i)+f13*z13(i)
            hsvs(i,if12)=f11*z21(i)+f12*z22(i)+f13*z23(i)
            hsvs(i,if13)=f11*z31(i)+f12*z32(i)+f13*z33(i)
            hsvs(i,if21)=f21*z11(i)+f22*z12(i)+f23*z13(i)
            hsvs(i,if22)=f21*z21(i)+f22*z22(i)+f23*z23(i)
            hsvs(i,if23)=f21*z31(i)+f22*z32(i)+f23*z33(i)
            hsvs(i,if31)=f31*z11(i)+f32*z12(i)+f33*z13(i)
            hsvs(i,if32)=f31*z21(i)+f32*z22(i)+f33*z23(i)
            hsvs(i,if33)=f31*z31(i)+f32*z32(i)+f33*z33(i)
          enddo
c     
c
c        deformation gradient transformed to local system
c        as F_bar=F*Q
c
c        else if (ihyper(mt).lt.0) then
        else if (ihyper(mt).eq.-1) then
          do i=lft,llt
            f1=z11(i)*hsvs(i,if11)+
     1           z12(i)*hsvs(i,if12)+z13(i)*hsvs(i,if13)
            f2=z21(i)*hsvs(i,if11)+
     1           z22(i)*hsvs(i,if12)+z23(i)*hsvs(i,if13)
            hsvs(i,if13)=z31(i)*hsvs(i,if11)+
     1           z32(i)*hsvs(i,if12)+z33(i)*hsvs(i,if13)
            hsvs(i,if11)=f1
            hsvs(i,if12)=f2
            f1=z11(i)*hsvs(i,if21)+
     1           z12(i)*hsvs(i,if22)+z13(i)*hsvs(i,if23)
            f2=z21(i)*hsvs(i,if21)+
     1           z22(i)*hsvs(i,if22)+z23(i)*hsvs(i,if23)
            hsvs(i,if23)=z31(i)*hsvs(i,if21)+
     1           z32(i)*hsvs(i,if22)+z33(i)*hsvs(i,if23)
            hsvs(i,if21)=f1
            hsvs(i,if22)=f2
            f1=z11(i)*hsvs(i,if31)+
     1           z12(i)*hsvs(i,if32)+z13(i)*hsvs(i,if33)
            f2=z21(i)*hsvs(i,if31)+
     1           z22(i)*hsvs(i,if32)+z23(i)*hsvs(i,if33)
            hsvs(i,if33)=z31(i)*hsvs(i,if31)+
     1           z32(i)*hsvs(i,if32)+z33(i)*hsvs(i,if33)
            hsvs(i,if31)=f1
            hsvs(i,if32)=f2
          enddo
        endif
      endif
c     
c     if ishrmp(mt).lt.0 then the flag for temperature is active
c     temperatures are stored in the temps array.
c
c        write(*,*)'Z2:',z11(1),z22(1),z33(1)
c        write(*,*)'dyn21_2 F:',
c     $       hsvs(1,if11),hsvs(1,if22),hsvs(1,if33)
      num_nods=8
      if (eltype.eq.'shl_t') num_nods=4
      if (eltype.eq.'sld2d') num_nods=4
      if (ishrmp(mt).lt.0) then
        call usr_temps (a(ntmp0+1),a(n19),itemp,temps,lft,llt,ix1,ix2,
     .   ix3,ix4,ix5,ix6,ix7,ix8,num_nods,itopaz,0,0.,0.,0.)
      else
        do i=lft,llt
          temps(i)=0.0
        enddo
      endif
c
c     convert rates to increments
c
      do i=lft,llt
        d1(i)=d1(i)*dt1siz(i)
        d2(i)=d2(i)*dt1siz(i)
        if (eltype.eq.'sld2d') then
          d3(i)=0.0
        else
          d3(i)=d3(i)*dt1siz(i)
        endif
        d4(i)=d4(i)*dt1siz(i)
        d5(i)=d5(i)*dt1siz(i)
        d6(i)=d6(i)*dt1siz(i)
      enddo
c
c     Get thermal history variables from thermal
c     user material
c
      call get_thhsv(nthhsv,thhsv,itemp,lft,llt,
     1 num_nods,itopaz,nnm1,rcoor,scoor,tcoor)
c
c     characteristic element size
      do i=lft,llt
        elsizv(i)=dx(i)
      enddo
c
c     element id
      if (eltype.eq.'shl_t'.or.eltype.eq.'sld2d') then
        do i=lft,llt
          ielm=nshbwp(nnm1+i)
          if (ielm.gt.0) idelev(i)=lqfinv(ielm,4)
        enddo
      elseif (eltype.eq.'tshel') then
        do i=lft,llt
          ielm=ntsbwp(nnm1+i)
          if (ielm.gt.0) idelev(i)=lqfinv(ielm,5)
        enddo
      else
        do i=lft,llt
          ielm=nhxbwp(nnm1+i)
          if (ielm.gt.0) idelev(i)=lqfinv(ielm,2)
        enddo
      endif
c
c     memory pointer for extra memory data
      lcma=ncma+nmmat+ia(ncma+mxt(lft)-1)-1
c
c     vector parameter
c
      ivect=ivectr(mt)
c
      mte40=mt-40
c
c     scalar umat
c
      if (ivect.eq.0) then
        do 90 i=lft,llt
c
c       gathering of variables
c       added the following lines
c       volume ratio
        vrat=def(i)
c       add in transformation matrix
        xM(1)=z11(i)
        xM(2)=z12(i)
        xM(3)=z13(i)
        xS(1)=z21(i)
        xS(2)=z22(i)
        xS(3)=z23(i)
        xN(1)=z31(i)
        xN(2)=z32(i)
        xN(3)=z33(i)
c
        index =i
        dt1   =dt1siz(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
        epsp=epsps(i)
        temper=temps(i)
        failel=.false.
        if (ifipts(mt).lt.0) failel=hsvs(i,ind_fipts).ne.0.
        if (no_hsvs.gt.0) then
          do 30 j=1,no_hsvs
          hsv(j)=hsvs(i,j)
 30       continue
        endif
        if (ihyper(mt).ne.0) then
          do 40 j=1,9
          hsv(no_hsvs+j)=hsvs(i,ind_defgrad-1+j)
 40       continue
        endif
        if (nthhsv.gt.0) then
          do j=1,nthhsv
            thhsvi(j)=thhsv(i,j)
          enddo
        endif
        elsiz=elsizv(i)
        idele=idelev(i)
c
c       call user developed subroutines here
c       changed the call to "umat45"
c
        go to (41,42,43,44,45,46,47,48,49,50), mte40
 41     call umat41 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
 42     call umat42 (cm(mx+1),eps,sig,hsv,dt1,
     1   capa,eltype,tt,crv,a(lcma),temper,0,elsiz,idele)
        go to 60
 43     call umat43 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
 44     call umat44 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
 45     call umat45 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1       capa,'solid',tt,temper,failel,crv,xM,xS,xN,
     1       vrat,Eb1)
c 45     call umat45 (cm(mx+1),eps,sig,epsp,hsv,dt1,
c     1   capa,eltype,tt,temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
 46     call umat46 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,a(lcma),thhsvi,nthhsv,0,
     2   elsiz,idele)
        go to 60
 47     call umat47 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
 48     call umat48 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
 49     call umat49 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,a(lcma),0,elsiz,idele)
        go to 60
 50     call umat50 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,a(lcma),0,elsiz,idele)
c
c       scattering of variables
c
 60     if (no_hsvs.gt.0) then
          do j=1,no_hsvs
            hsvs(i,j)=hsv(j)
          enddo
        endif
        if (ihyper(mt).ne.0) then
          do j=1,9
            hsvs(i,ind_defgrad-1+j)=hsv(no_hsvs+j)
          enddo
        endif
c
        if (ifipts(mt).eq.0) failel=.false.
c
        if (failel) then
          sig1(i)=0.
          sig2(i)=0.
          sig3(i)=0.
          sig4(i)=0.
          sig5(i)=0.
          sig6(i)=0.
          if (ifipts(mt).lt.0) hsvs(i,ind_fipts)=1.
          nfipts(i)=nfipts(i)+1
          if (ifipts(mt).ge.0) then
          failur=.true.
          ifail(i)=1
          fail(i)=0.
          elseif (nfipts(i).ge.nint(cm(mx-ifipts(mt)))) then
          failur=.true.
          ifail(i)=1
          fail(i)=0.
          endif
        else
          sig1(i)=sig(1)
          sig2(i)=sig(2)
          sig3(i)=sig(3)
          sig4(i)=sig(4)
          sig5(i)=sig(5)
          sig6(i)=sig(6)
        endif
        epsps(i)=epsp
c
C *******************************************************************   
c lines added by Daniel Einstein for sound speed calculation
c
c        if (Eb1.gt.Ebk) then
          Ebk=Eb1
c        endif
c
 90   continue
c     
c
      ss(lft)=Ebk
      do 91 i=lft,llt
        cb(i)=ss(lft)
 91   continue    
C *******************************************************************
c
c     vector umat
c
      else
c
c       assume no failed elements
c
        if (ifipts(mt).ge.0) then
        do i=lft,llt
          failels(i)=.false.
        enddo
        else
          do i=lft,llt
            failels(i)=hsvs(i,ind_fipts).ne.0.
          enddo
        endif
c
c       call user defined subroutines here
c       first history variable is the one requested
c
        if(ivumat(mt).eq.0) then
c
        if (mt.eq.41) then
          call umat41v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,a(lcma),
     .     0,elsizv,idelev)
        elseif (mt.eq.42) then
          call umat42v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,a(lcma),
     .     0,elsizv,idelev)
        elseif (mt.eq.43) then
          call umat43v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,a(lcma),
     .     0,elsizv,idelev)
        elseif (mt.eq.44) then
          call umat44v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,a(lcma),
     .     0,elsizv,idelev)
        elseif (mt.eq.45) then
          call umat45v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,a(lcma),
     .     0,elsizv,idelev)
        elseif (mt.eq.46) then
          call umat46v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,a(lcma),
     .     thhsv,nthhsv,0,elsizv,idelev)
        elseif (mt.eq.47) then
          call umat47v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,a(lcma),
     .     0,elsizv,idelev)
        elseif (mt.eq.48) then
c         ipt=1
          call umat48v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,ipt,a(n8),a(n9),a(lcma),0,
     .     elsizv,idelev)
        elseif (mt.eq.49) then
c         ipt=1
          call umat49v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,ipt,temps,a(lcma),0,elsizv,
     .     idelev)
        elseif (mt.eq.50) then
          call umat50v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,a(lcma),
     .     0,elsizv,idelev)
        elseif (mt.eq.281) then
c         ipt=1
          call mat281h(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
        elseif (mt.eq.282) then
c         ipt=1
          call mat282h(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
c       elseif (mt.eq.283) then
c         call umat283h(cm(mx+1),d1,d2,d3,d4,d5,d6,
c    .     sig1,sig2,sig3,sig4,
c    .     sig5,sig6,epsps,hsvs,
c    .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,
c    .     a(lcma))
        elseif (mt.eq.284) then
c         ipt=1
          call mat284h(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
c       elseif (mt.eq.285) then
c         call umat285h(cm(mx+1),d1,d2,d3,d4,d5,d6,
c    .     sig1,sig2,sig3,sig4,
c    .     sig5,sig6,epsps,hsvs,
c    .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv)
        elseif (mt.eq.286) then
c         ipt=1
          call mat286h(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
        elseif (mt.eq.287) then
c         ipt=1
          call mat287h(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
        endif
c
        endif
        if (ivumat(mt).ne.0) then 
          ndi=3
          nsh=3
          nblock=llt-lft+1
          curveid=ivumat(mt)
          icid=int(curveid)
          if (icid.lt.0) then
            iicid=lcids(-icid)
            k2=npc(iicid)
            k3=npc(iicid+1)
            nprops=nint((k3-k2)/2.0)
          else
            nprops=40+icmadd(mt)
          endif
c          write(*,*)'INDICES DF:',
c     $         if11,if22,if33,if12,if21,if23,if32,if13,if31
c          write(*,*)'INDICES VMAT:',
c     $         ifo11,ifo22,ifo33,ifo12,ifo21,ifo23,ifo32,ifo13,ifo31
          call ext_vecumat(lft,llt,cm,a(lcma),bqs,mt,crv,npc,plc,nnm1,
     $         curveid,rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,
     $         nprops,ndi,nsh,
     $         if11,if22,if33,if12,if21,if23,if32,if13,if31,no_hsvs,
     $         ifo11,ifo22,ifo33,ifo12,ifo21,ifo23,ifo32,ifo13,ifo31)
        endif
c     
        if (ifipts(mt).ne.0) then
          do i=lft,llt
            if (failels(i)) then
              sig1(i)=0.
              sig2(i)=0.
              sig3(i)=0.
              sig4(i)=0.
              sig5(i)=0.
              sig6(i)=0.
              if (ifipts(mt).lt.0) hsvs(i,ind_fipts)=1.
              nfipts(i)=nfipts(i)+1
              if (ifipts(mt).ge.0) then
                failur=.true.
                ifail(i)=1
                fail(i)=0.
              elseif (nfipts(i).ge.nint(cm(mx-ifipts(mt)))) then
                failur=.true.
                ifail(i)=1
                fail(i)=0.
              endif
            endif
          enddo
        endif
      endif
c
c     energy calculations
c
      do 100 i=lft,llt
      einc(i)=(d1(i)*sig1(i)+d2(i)*sig2(i)+d3(i)*sig3(i)+
     . d4(i)*sig4(i)+d5(i)*sig5(i)+d6(i)*sig6(i))+einc(i)
 100  continue
c
      if (iorien(mt).ne.0) then
c
c       transform stresses to global system
c     
        if (ihyper(mt).ge.0) then
c     if ((ihyper(mt).eq.0) .or. (ihyper(mt).eq.1)) then
          call gblm22 (lft,llt)
        endif
c     
c     if ((ihyper(mt).eq.1) .or. (ihyper(mt).eq.0) ) then
        if (ihyper(mt).gt.0) then
c     
c     deformation gradient stored in global system
c
          if (ifg(lft).eq.1) then
            do i=lft,llt
              hsvs(i,if11)=g11(i)
              hsvs(i,if21)=g21(i)
              hsvs(i,if31)=g31(i)
              hsvs(i,if12)=g12(i)
              hsvs(i,if22)=g22(i)
              hsvs(i,if32)=g32(i)
              hsvs(i,if13)=g13(i)
              hsvs(i,if23)=g23(i)
              hsvs(i,if33)=g33(i)
            enddo
c     
c     deformation gradient transformed back to global system
c
          else
            do i=lft,llt
              f1=z11(i)*hsvs(i,if11)+
     1             z21(i)*hsvs(i,if21)+z31(i)*hsvs(i,if31)
              f2=z12(i)*hsvs(i,if11)+
     1             z22(i)*hsvs(i,if21)+z32(i)*hsvs(i,if31)
              hsvs(i,if31)=z13(i)*hsvs(i,if11)+
     1             z23(i)*hsvs(i,if21)+z33(i)*hsvs(i,if31)
              hsvs(i,if11)=f1
              hsvs(i,if21)=f2
              f1=z11(i)*hsvs(i,if12)+
     1             z21(i)*hsvs(i,if22)+z31(i)*hsvs(i,if32)
              f2=z12(i)*hsvs(i,if12)+
     1             z22(i)*hsvs(i,if22)+z32(i)*hsvs(i,if32)
              hsvs(i,if32)=z13(i)*hsvs(i,if12)+
     1             z23(i)*hsvs(i,if22)+z33(i)*hsvs(i,if32)
              hsvs(i,if12)=f1
              hsvs(i,if22)=f2
              f1=z11(i)*hsvs(i,if13)+
     1             z21(i)*hsvs(i,if23)+z31(i)*hsvs(i,if33)
              f2=z12(i)*hsvs(i,if13)+
     1             z22(i)*hsvs(i,if23)+z32(i)*hsvs(i,if33)
              hsvs(i,if33)=z13(i)*hsvs(i,if13)+
     1             z23(i)*hsvs(i,if23)+z33(i)*hsvs(i,if33)
              hsvs(i,if13)=f1
              hsvs(i,if23)=f2
            enddo
          endif
c     
        elseif (ihyper(mt).eq.-2) then
c     
c     deformation gradient stored in global system
c     
          if (ifg(lft).eq.1) then
            do i=lft,llt
              hsvs(i,if11)=g11(i)
              hsvs(i,if21)=g21(i)
              hsvs(i,if31)=g31(i)
              hsvs(i,if12)=g12(i)
              hsvs(i,if22)=g22(i)
              hsvs(i,if32)=g32(i)
              hsvs(i,if13)=g13(i)
              hsvs(i,if23)=g23(i)
              hsvs(i,if33)=g33(i)
            enddo
c     
c     deformation gradient transformed back to global system
c     
          else
            do i=lft,llt
              f11=z11(i)*hsvs(i,if11)+
     1             z21(i)*hsvs(i,if21)+z31(i)*hsvs(i,if31)
              f21=z12(i)*hsvs(i,if11)+
     1             z22(i)*hsvs(i,if21)+z32(i)*hsvs(i,if31)
              f31=z13(i)*hsvs(i,if11)+
     1             z23(i)*hsvs(i,if21)+z33(i)*hsvs(i,if31)
c     hsvs(i,if11)=f1
c     hsvs(i,if21)=f2
              f12=z11(i)*hsvs(i,if12)+
     1             z21(i)*hsvs(i,if22)+z31(i)*hsvs(i,if32)
              f22=z12(i)*hsvs(i,if12)+
     1             z22(i)*hsvs(i,if22)+z32(i)*hsvs(i,if32)
              f32=z13(i)*hsvs(i,if12)+
     1             z23(i)*hsvs(i,if22)+z33(i)*hsvs(i,if32)
c     hsvs(i,if12)=f1
c     hsvs(i,if22)=f2
              f13=z11(i)*hsvs(i,if13)+
     1             z21(i)*hsvs(i,if23)+z31(i)*hsvs(i,if33)
              f23=z12(i)*hsvs(i,if13)+
     1             z22(i)*hsvs(i,if23)+z32(i)*hsvs(i,if33)
              f33=z13(i)*hsvs(i,if13)+
     1             z23(i)*hsvs(i,if23)+z33(i)*hsvs(i,if33)
c     hsvs(i,if13)=f1
c     hsvs(i,if23)=f2
              hsvs(i,if11)=f11*z11(i)+f12*z21(i)+f13*z31(i)
              hsvs(i,if12)=f11*z12(i)+f12*z22(i)+f13*z32(i)
              hsvs(i,if13)=f11*z13(i)+f12*z23(i)+f13*z33(i)
              hsvs(i,if21)=f21*z11(i)+f22*z21(i)+f23*z31(i)
              hsvs(i,if22)=f21*z12(i)+f22*z22(i)+f23*z32(i)
              hsvs(i,if23)=f21*z13(i)+f22*z23(i)+f23*z33(i)
              hsvs(i,if31)=f31*z11(i)+f32*z21(i)+f33*z31(i)
              hsvs(i,if32)=f31*z12(i)+f32*z22(i)+f33*z32(i)
              hsvs(i,if33)=f31*z13(i)+f32*z23(i)+f33*z33(i)
            enddo
          endif
c     elseif (ihyper(mt).lt.0) then
        elseif (ihyper(mt).eq.-1) then
c     
c     deformation gradient stored in global system
c     
          if (ifg(lft).eq.1) then
            do i=lft,llt
              hsvs(i,if11)=g11(i)
              hsvs(i,if21)=g21(i)
              hsvs(i,if31)=g31(i)
              hsvs(i,if12)=g12(i)
              hsvs(i,if22)=g22(i)
              hsvs(i,if32)=g32(i)
              hsvs(i,if13)=g13(i)
              hsvs(i,if23)=g23(i)
              hsvs(i,if33)=g33(i)
            enddo
c     
c     deformation gradient transformed back to global system
c     
          else
            do i=lft,llt
              f1=z11(i)*hsvs(i,if11)+
     1             z21(i)*hsvs(i,if12)+z31(i)*hsvs(i,if13)
              f2=z12(i)*hsvs(i,if11)+
     1             z22(i)*hsvs(i,if12)+z32(i)*hsvs(i,if13)
              hsvs(i,if13)=z13(i)*hsvs(i,if11)+
     1             z23(i)*hsvs(i,if12)+z33(i)*hsvs(i,if13)
              hsvs(i,if11)=f1
              hsvs(i,if12)=f2
              f1=z11(i)*hsvs(i,if21)+
     1             z21(i)*hsvs(i,if22)+z31(i)*hsvs(i,if23)
              f2=z12(i)*hsvs(i,if21)+
     1             z22(i)*hsvs(i,if22)+z32(i)*hsvs(i,if23)
              hsvs(i,if23)=z13(i)*hsvs(i,if21)+
     1             z23(i)*hsvs(i,if22)+z33(i)*hsvs(i,if23)
              hsvs(i,if21)=f1
              hsvs(i,if22)=f2
              f1=z11(i)*hsvs(i,if31)+
     1             z21(i)*hsvs(i,if32)+z31(i)*hsvs(i,if33)
              f2=z12(i)*hsvs(i,if31)+
     1             z22(i)*hsvs(i,if32)+z32(i)*hsvs(i,if33)
              hsvs(i,if33)=z13(i)*hsvs(i,if31)+
     1             z23(i)*hsvs(i,if32)+z33(i)*hsvs(i,if33)
              hsvs(i,if31)=f1
              hsvs(i,if32)=f2
            enddo
          endif
        endif
      endif
c     
      return
c     900  format(/
c    1 ' *** Error User defined material in part',i12,
c    2 ' exceeds the current limit'/
c    3 10x,'of history variables for solids.')
      end
      subroutine usr_temps (tnew,fval,itemp,ctmp,lft,llt,ix1,ix2,
     . ix3,ix4,ix5,ix6,ix7,ix8,num_nods,itopaz,nnm1,rcoor,scoor,tcoor)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      common/blk03/sigma,tolb,nummat,numnpt,numel,numsh4,iunits,iband,
     1 nbs,nsl,nslvt,nmsrt,numelst,nprof,irtyp,itmaxb,numeltt,iphase,
     2 igenl,igenm,igene,igend,isotr,iem,isoln,igeom,mcut,numsh12,
     3 ndtot,nsl_th,lenhsv,numel2,numel4,numel6,numel8,numel10
      dimension tnew(*),fval(*),ctmp(*),ix1(*),ix2(*),ix3(*),ix4(*),
     .          ix5(*),ix6(*),ix7(*),ix8(*),itopaz(*)
c
c    temperatures at element center
c
      if (itemp.gt.0) then
        tempr=fval(itemp)
        do i=lft,llt
          ctmp(i)=tempr
        enddo
c
      elseif (itemp.lt.0 .or. itopaz(1).eq.999) then
        if (num_nods.eq.8) then
          do i=lft,llt
            ctmp(i) = .125*(tnew(ix1(i))+tnew(ix2(i))+tnew(ix3(i))
     .                      +tnew(ix4(i))+tnew(ix5(i))+tnew(ix6(i))
     .                      +tnew(ix7(i))+tnew(ix8(i)))
          enddo
c
        elseif (num_nods.eq.4) then
          if (numsh12.ne.0) then
           call get_tmp_tpz (ctmp,lft,llt,nnm1,rcoor,scoor,tcoor)
          else
            do i=lft,llt
            ctmp(i)=.25*(tnew(ix1(i))+tnew(ix2(i))
     .                  +tnew(ix3(i))+tnew(ix4(i)))
            enddo
          endif
c
        else
          do i=lft,llt
            ctmp(i)=.50*(tnew(ix1(i))+tnew(ix2(i)))
          enddo
        endif
      endif
c
      return
      end
      subroutine urtanh(cm,mt,crv,lft,llt,nconstp,nnm1,eltype)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     assume that F and all coordinate transformations have been
c     performed by this time, and everything is already loaded in
c     common blocks
c
      include 'nlqparm'
      include 'nhisparm.inc'
      integer n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      common/bk07/n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      integer lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      common/bk13/lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      include 'bk19.inc'
      include 'memaia.inc'
      include 'umatss.inc'
      common/bux13loc/ eps(6),sig(6),hsv(NHISVAR),temps(nlq),es(6,6)
      common/aux14loc/sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),sig5(nlq),
     & sig6(nlq),epsps(nlq),hsvs(nlq,NHISVAR)
      common/aux18loc/dd(nlq),def(nlq),ddq(nlq)
      common/aux33loc/
     1 ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ixs(nlq,4),mxt(nlq)
      common/aux40loc/
     1 a11(nlq),a12(nlq),a13(nlq),a21(nlq),a22(nlq),a23(nlq),
     2 a31(nlq),a32(nlq),a33(nlq),z11(nlq),z12(nlq),z13(nlq),
     3 z21(nlq),z22(nlq),z23(nlq),z31(nlq),z32(nlq),z33(nlq)
      common/aux35loc/rhoa(nlq),cb(nlq),davg(nlq),p(nlq)
      common/aux41loc/qq1(nlq),cbb(nlq),aj1(nlq)
      common/bk36loc/index
      common/soundloc/ss(nlq),sndsp(nlq),diagm(nlq),sarea(nlq),
     . dxl(nlq)
      common/subtssloc/dt1siz(nlq)
      common/vect8loc/dsave(nlq,6,6)
c
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
      common/bk01/itherm,itemp,ntmp0,ntmp1,itempan,itempdr,itmpe(44)
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie
      common/aux2loc/d1(nlq),d2(nlq),d3(nlq),d4(nlq),d5(nlq),d6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common/bk26/begtim,nintcy
      common/numcpu/ncpu,ncpua,ncpub,lenvec(8)
      common/thropt/itopaz(101)
      common/tsbsis/h(8,5,5),pr(8,5,5),ps(8,5,5),pt(8,5,5),ipt,
     1 nip,wgts(10,10),szeta(10,10)
c
      dimension cm(*),crv(lq1,2,*),nconstp(*)
      character*(*) eltype
c     dimension eps(6),sig(6),hsv(NHISVAR),temps(nlq),es(6,6)
      dimension sigmaold(nlq,6)
      dimension dfgold(nlq,9)
      logical failel,failels(nlq)
c     data temps/nlq*0.0/
c
c$omp threadprivate (/aux2loc/)
c$omp threadprivate (/aux14loc/)
c$omp threadprivate (/aux18loc/)
c$omp threadprivate (/aux33loc/)
c$omp threadprivate (/aux35loc/)
c$omp threadprivate (/aux40loc/)
c$omp threadprivate (/aux41loc/)
c$omp threadprivate (/bk36loc/)
c$omp threadprivate (/soundloc/)
c$omp threadprivate (/subtssloc/)
c$omp threadprivate (/vect8loc/)
c$omp threadprivate (/bux13loc/)
      capa=0.0
      mx=48*(mxt(lft)-1)
c
c     Determining number of requested history variables
c
      no_hsvs=nconstp(mxt(lft))
      if (lenvec(5) .eq.1) no_hsvs=no_hsvs-6
      if (iorien(mt).ne.0) no_hsvs=no_hsvs-6
      if (ihyper(mt).ne.0) no_hsvs=no_hsvs-9
      if (iorien(mt).ne.0.and.ihyper(mt).lt.0) no_hsvs=no_hsvs-6
      if (ifipts(mt).lt.0) no_hsvs=no_hsvs-1
c
c     initialization of indeces
c
      if (iorien(mt).ne.0) then
        if (ihyper(mt).ne.0) then
          ind_defgrad=1+no_hsvs
          ind_ortho=ind_defgrad+9
        else
          ind_ortho=1+no_hsvs
        endif
        if (ifipts(mt).lt.0) then
          ind_fipts=ind_ortho
          ind_ortho=ind_ortho+1
        endif
        ind_q11=ind_ortho
        ind_q12=ind_ortho+1
        ind_q13=ind_ortho+2
        ind_q31=ind_ortho+3
        ind_q32=ind_ortho+4
        ind_q33=ind_ortho+5
        ind_last=ind_q33
        if (ihyper(mt).lt.0) then
          ind_qstore=ind_last+1
          ind_last=ind_last+6
        endif
      else
        if (ihyper(mt).ne.0) then
          ind_defgrad=1+no_hsvs
          ind_last=ind_defgrad+8
        else
          ind_last=no_hsvs
        endif
        if (ifipts(mt).lt.0) then
          ind_fipts=ind_last+1
          ind_last=ind_last+1
        endif
      endif
c
      if (iorien(mt).ne.0) then
c
c       store global stresses
c
        if (ihyper(mt).ge.0) then
          do i=lft,llt
            sigmaold(i,1)=sig1(i)
            sigmaold(i,2)=sig2(i)
            sigmaold(i,3)=sig3(i)
            sigmaold(i,4)=sig4(i)
            sigmaold(i,5)=sig5(i)
            sigmaold(i,6)=sig6(i)
          enddo
c
c         transform stresses to local system
c
          call ldsm22(lft,llt)
        endif
c
c       storage of deformation gradient
c
        if (ihyper(mt).gt.0) then
          if11=ind_defgrad
          if21=ind_defgrad+1
          if31=ind_defgrad+2
          if12=ind_defgrad+3
          if22=ind_defgrad+4
          if32=ind_defgrad+5
          if13=ind_defgrad+6
          if23=ind_defgrad+7
          if33=ind_defgrad+8
c
c         store deformation gradient in global system
c
          do i=lft,llt
            dfgold(i,1)=hsvs(i,if11)
            dfgold(i,2)=hsvs(i,if21)
            dfgold(i,3)=hsvs(i,if31)
            dfgold(i,4)=hsvs(i,if12)
            dfgold(i,5)=hsvs(i,if22)
            dfgold(i,6)=hsvs(i,if32)
            dfgold(i,7)=hsvs(i,if13)
            dfgold(i,8)=hsvs(i,if23)
            dfgold(i,9)=hsvs(i,if33)
          enddo
c
c         transform deformation gradient to local system
c         as F_bar=Q^T*F
c
          do i=lft,llt
            f1=z11(i)*hsvs(i,if11)+
     1         z12(i)*hsvs(i,if21)+z13(i)*hsvs(i,if31)
            f2=z21(i)*hsvs(i,if11)+
     1         z22(i)*hsvs(i,if21)+z23(i)*hsvs(i,if31)
            hsvs(i,if31)=z31(i)*hsvs(i,if11)+
     1         z32(i)*hsvs(i,if21)+z33(i)*hsvs(i,if31)
            hsvs(i,if11)=f1
            hsvs(i,if21)=f2
            f1=z11(i)*hsvs(i,if12)+
     1         z12(i)*hsvs(i,if22)+z13(i)*hsvs(i,if32)
            f2=z21(i)*hsvs(i,if12)+
     1         z22(i)*hsvs(i,if22)+z23(i)*hsvs(i,if32)
            hsvs(i,if32)=z31(i)*hsvs(i,if12)+
     1         z32(i)*hsvs(i,if22)+z33(i)*hsvs(i,if32)
            hsvs(i,if12)=f1
            hsvs(i,if22)=f2
            f1=z11(i)*hsvs(i,if13)+
     1         z12(i)*hsvs(i,if23)+z13(i)*hsvs(i,if33)
            f2=z21(i)*hsvs(i,if13)+
     1         z22(i)*hsvs(i,if23)+z23(i)*hsvs(i,if33)
            hsvs(i,if33)=z31(i)*hsvs(i,if13)+
     1         z32(i)*hsvs(i,if23)+z33(i)*hsvs(i,if33)
            hsvs(i,if13)=f1
            hsvs(i,if23)=f2
          enddo
c
c       storage of deformation gradient
c
        else if (ihyper(mt).lt.0) then
          if11=ind_defgrad
          if21=ind_defgrad+1
          if31=ind_defgrad+2
          if12=ind_defgrad+3
          if22=ind_defgrad+4
          if32=ind_defgrad+5
          if13=ind_defgrad+6
          if23=ind_defgrad+7
          if33=ind_defgrad+8
c
c         store deformation gradient in global system
c
          do i=lft,llt
            dfgold(i,1)=hsvs(i,if11)
            dfgold(i,2)=hsvs(i,if21)
            dfgold(i,3)=hsvs(i,if31)
            dfgold(i,4)=hsvs(i,if12)
            dfgold(i,5)=hsvs(i,if22)
            dfgold(i,6)=hsvs(i,if32)
            dfgold(i,7)=hsvs(i,if13)
            dfgold(i,8)=hsvs(i,if23)
            dfgold(i,9)=hsvs(i,if33)
          enddo
c
c         transform deformation gradient to local system
c         as F_bar=F*Q
c
          do i=lft,llt
            f1=z11(i)*hsvs(i,if11)+
     1         z12(i)*hsvs(i,if12)+z13(i)*hsvs(i,if13)
            f2=z21(i)*hsvs(i,if11)+
     1         z22(i)*hsvs(i,if12)+z23(i)*hsvs(i,if13)
            hsvs(i,if13)=z31(i)*hsvs(i,if11)+
     1         z32(i)*hsvs(i,if12)+z33(i)*hsvs(i,if13)
            hsvs(i,if11)=f1
            hsvs(i,if12)=f2
            f1=z11(i)*hsvs(i,if21)+
     1         z12(i)*hsvs(i,if22)+z13(i)*hsvs(i,if23)
            f2=z21(i)*hsvs(i,if21)+
     1         z22(i)*hsvs(i,if22)+z23(i)*hsvs(i,if23)
            hsvs(i,if23)=z31(i)*hsvs(i,if21)+
     1         z32(i)*hsvs(i,if22)+z33(i)*hsvs(i,if23)
            hsvs(i,if21)=f1
            hsvs(i,if22)=f2
            f1=z11(i)*hsvs(i,if31)+
     1         z12(i)*hsvs(i,if32)+z13(i)*hsvs(i,if33)
            f2=z21(i)*hsvs(i,if31)+
     1         z22(i)*hsvs(i,if32)+z23(i)*hsvs(i,if33)
            hsvs(i,if33)=z31(i)*hsvs(i,if31)+
     1         z32(i)*hsvs(i,if32)+z33(i)*hsvs(i,if33)
            hsvs(i,if31)=f1
            hsvs(i,if32)=f2
          enddo
        endif
      endif
c
c     if ishrmp(mt).lt.0 then the flag for temperature is active
c     temperatures are stored in the temps array.
c
      num_nods=8
      if (eltype.eq.'shl_t') num_nods=4
      if (eltype.eq.'sld2d') num_nods=4
      if (ishrmp(mt).lt.0) then
        call usr_temps (a(ntmp0+1),a(n19),itemp,temps,lft,llt,ix1,ix2,
     .   ix3,ix4,ixs(1,1),ixs(1,2),ixs(1,3),ixs(1,4),num_nods,
     .   itopaz,0,0.,0.,0.)
      else
        do i=lft,llt
          temps(i)=0.0
        enddo
      endif
c
c     memory pointer for extra memory data
      lcma=ncma+nmmat+ia(ncma+mxt(lft)-1)-1
c
c     vector flag
c
      ivect=ivectr(mt)
c
      mte40=mt-40
c
c     scalar utan
c
      if (ivect.eq.0) then
        do 90 i=lft,llt
c
c       gathering of variables
c
        index =i
        dt1   =dt1siz(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
        epsp=epsps(i)
        temper=temps(i)
        failel=.false.
        if (ifipts(mt).lt.0) failel=hsvs(i,ind_fipts).ne.0.
        if (no_hsvs.gt.0) then
          do 30 j=1,no_hsvs
          hsv(j)=hsvs(i,j)
 30       continue
        endif
        if (ihyper(mt).ne.0) then
          do 40 j=1,9
          hsv(no_hsvs+j)=hsvs(i,ind_defgrad-1+j)
 40       continue
        endif
c
c       setting material tangent modulus to = 0
c
        do k=1,6
          do j=1,6
            es(j,k)=0.
          enddo
        enddo
c
c       call user developed subroutines here
c
        go to (41,42,43,44,45,46,47,48,49,50), mte40
   41   call utan41 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,failel,a(lcma),0)
        go to 60
   42   call utan42 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,failel,a(lcma),0)
        go to 60
   43   call utan43 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,failel,a(lcma),0)
        go to 60
   44   call utan44 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,failel,a(lcma),0)
        go to 60
   45   call utan45 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,failel,a(lcma),0)
        go to 60
   46   call utan46 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,failel,a(lcma),0)
        go to 60
   47   call utan47 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,failel,a(lcma),0)
        go to 60
   48   call utan48 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,failel,a(lcma),0)
        go to 60
   49   call utan49 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,failel,a(lcma),0)
        go to 60
 50     call utan50 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,failel,a(lcma),0)
c
c       scattering tangent modulus
c
 60     do k=1,6
          do j=1,6
           dsave(i,j,k)=es(j,k)
          enddo
        enddo
 90    continue
c
c
c     vector umat
c
      else
c
c       setting material tangent modulus to = 0
c
        do k=1,6
          do j=1,6
            do i=lft,llt
              dsave(i,j,k)=0.
            enddo
          enddo
        enddo
c
        if (ifipts(mt).ge.0) then
          do i=lft,llt
            failels(i)=.false.          
          enddo
        else
          do i=lft,llt
            failels(i)=hsvs(i,ind_fipts).ne.0.
          enddo
        endif
c
c       call user developed routines here
c       first history variable is the one requested
c
        if (mt.eq.41) then
          call utan41v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,failels,
     .     a(lcma),0)
        elseif (mt.eq.42) then
          call utan42v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,failels,
     .     a(lcma),0)
        elseif (mt.eq.43) then
          call utan43v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,failels,
     .     a(lcma),0)
        elseif (mt.eq.44) then
          call utan44v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,failels,
     .     a(lcma),0)
        elseif (mt.eq.45) then
          call utan45v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,failels,
     .     a(lcma),0)
        elseif (mt.eq.46) then
          call utan46v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,failels,
     .     a(lcma),0)
        elseif (mt.eq.47) then
          call utan47v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,failels,
     .     a(lcma),0)
        elseif (mt.eq.48) then
          call utan48v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,failels,
     .     a(lcma),0)
        elseif (mt.eq.49) then
          call utan49v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,failels,
     .     a(lcma),0)
        elseif (mt.eq.50) then
          call utan50v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,failels,
     .     a(lcma),0)
        elseif (mt.eq.286) then
          call utan286v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,failels,
     .     a(lcma))
        elseif (mt.eq.287) then
          call utan287v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,failels,
     .     a(lcma))
        endif
      endif
c
      if (iorien(mt).ne.0) then
c
c       transform tangent stiffness to global system
c
        if (ihyper(mt).ge.0) then
          call transform44(z11,z12,z13,z21,z22,z23,z31,z32,z33,
     1     lft,llt)
c
c         restore stresses in global system
c
          do i=lft,llt
            sig1(i)=sigmaold(i,1)
            sig2(i)=sigmaold(i,2)
            sig3(i)=sigmaold(i,3)
            sig4(i)=sigmaold(i,4)
            sig5(i)=sigmaold(i,5)
            sig6(i)=sigmaold(i,6)
          enddo
        endif
c
c       restore deformation gradient in global system
c
        if (ihyper(mt).ne.0) then
          do i=lft,llt
            hsvs(i,if11)=dfgold(i,1)
            hsvs(i,if21)=dfgold(i,2)
            hsvs(i,if31)=dfgold(i,3)
            hsvs(i,if12)=dfgold(i,4)
            hsvs(i,if22)=dfgold(i,5)
            hsvs(i,if32)=dfgold(i,6)
            hsvs(i,if13)=dfgold(i,7)
            hsvs(i,if23)=dfgold(i,8)
            hsvs(i,if33)=dfgold(i,9)
          enddo
        endif
      endif
c
      return
      end
      subroutine urtans(cm,mt,capa,lft,llt,nconstp,nnm1)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
      integer n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      common/bk07/n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      integer lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      common/bk13/lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1   lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2   lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
      include 'bk19.inc'
      include 'memaia.inc'
      include 'umatss.inc'
c
c     dynamic memory allocation stuff.  This has to be a C style 
c     include because of the use of integer*8
c
      integer i_mem(0:1)
      real r_mem(0:1)
      integer*8 i8_mem(0:1)
      real*8 r8_mem(0:1)
      real*8 real8_mem
      integer*8 m_mem(0:1)
      common/dynmem/real8_mem(0:1)
      equivalence (i_mem,r_mem,i8_mem,r8_mem,real8_mem,m_mem)
      integer*8 mem_alloc, mems_alloc,mem_realloc, memh_ptr
      integer mem_length,memh_alloc,memh_realloc,memh_length
      integer memh_enum,memh_enums,memh_attach
c      integer mem_newmark,mem_getmark,memh_getmark
      integer*8 memh_detach,memh_ralloc
      integer*8 m_to_m8
      integer*8 m_to_mm
      external mem_alloc, mems_alloc, mem_realloc, mem_length
      external memh_alloc, memh_realloc, memh_length, memh_ptr
      external memh_enum,memh_attach,memh_detach,memh_ralloc
      external memh_enums
      external m_to_m8
      external m_to_mm
      integer*8 dm_x,dm_xr,dm_v,dm_vr,dm_a,dm_ar,dm_rots,dm_rotr
      integer*8 dm_xms,dm_xmr,dm_xm2,dm_xmz,dm_xma
      integer*8 dm_xtpz,dm_disp
      integer*8 dm_me1,dm_me2,dm_me3
      integer*8 dm_masp,dm_masr
      integer*8 dm_xrb,dm_yrb,dm_zrb
      integer*8 dm_axrb,dm_ayrb,dm_azrb
      integer*8 dm_rbfx,dm_rbfy,dm_rbfz
      integer*8 dm_n45,dm_n46,dm_n47
      integer len_n45,len_n46,len_n47
      common/dynmem1/ 
     . dm_x,dm_xr,dm_v,dm_vr,dm_a,dm_ar,dm_rots,dm_rotr,
     . dm_xms,dm_xmr,dm_xm2,dm_xmz,dm_xma,
     . dm_xtpz,dm_disp,
     . dm_me1,dm_me2,dm_me3,
     . dm_masp,dm_masr,
     . dm_xrb,dm_yrb,dm_zrb,
     . dm_axrb,dm_ayrb,dm_azrb,
     . dm_rbfx,dm_rbfy,dm_rbfz,
     . dm_n45,dm_n46,dm_n47
      common/dynmem1l/len_n45,len_n46,len_n47
c
      integer*8 mem_allocm,mems_allocm,mem_reallocm
      integer memh_allocm,memh_reallocm
      external mem_allocm,mems_allocm,mem_reallocm
      external memh_allocm,memh_reallocm
      integer*8 mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      integer*8 mem_reallocmchk,mem_reallocchk
      integer memh_allocmchk,memh_allocchk
      external mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      external mem_reallocmchk,mem_reallocchk
      external memh_allocmchk,memh_allocchk
c
      common/frozloc/yhat(nlq,3,4),qloc(nlq,3,3),qmat(nlq,3,3)
      common/bux13loc/ eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      common/aux2loc/e1(nlq),e2(nlq),d3(nlq),e4(nlq),e5(nlq),e6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common /aux8cloc/
     & d1(nlq),d2(nlq),d4(nlq),d5(nlq),d6(nlq),
     & dfg11(nlq),dfg21(nlq),dfg31(nlq),dfg12(nlq),dfg22(nlq),
     & dfg32(nlq),dfg13(nlq),dfg23(nlq),dfg33(nlq),stg5(nlq),
     & stg6(nlq),a11(nlq),a12(nlq),a21(nlq),a22(nlq)
      common/aux14loc/
     & sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),
     & sig5(nlq),sig6(nlq),epsps(nlq),hsvs(nlq,NHISVAR)
      common/aux19loc/
     1 sign0(nlq),sign1(nlq),sign2(nlq),sign3(nlq),sign4(nlq),
     2 sign5(nlq),sign6(nlq)
      common/aux33loc/
     1 ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ixs(nlq,4),mxt(nlq)
      common/bk36loc/index
      common/failuloc/sieu(nlq),fail(nlq)
      common/hourgloc/ym(nlq),gm(nlq),ifsv(nlq)
      common/shlopt/istrn,istupd,ibelyts,miter
      common/shloptloc/ibelyt
      common/soundloc/ss(nlq),sndsp(nlq),diagm(nlq),sarea(nlq),
     . dxl(nlq)
      common/vect8loc/dsave(nlq,6,6)
c
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
      common/bk01/itherm,itemp,ntmp0,ntmp1,itempan,itempdr,itmpe(44)
      common/bk26/begtim,nintcy
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie
      common/numcpu/ncpu,ncpua,ncpub,lenvec(8)
      common/slcntr/islcnt(21)
      common/subtssloc/dt1siz(nlq)
      common/thropt/itopaz(101)
c
      dimension cm(*),nconstp(*)
c     dimension eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      dimension es(6,6)
      dimension sigmaold(nlq,6)
      logical failel,failels(nlq)
      dimension qmats(3,3)
c     data temps/nlq*0.0/
c
c     Note that no temperatures are calculated here
c
c$omp threadprivate (/aux2loc/)
c$omp threadprivate (/aux8cloc/)
c$omp threadprivate (/aux14loc/)
c$omp threadprivate (/aux19loc/)
c$omp threadprivate (/aux33loc/)
c$omp threadprivate (/bk36loc/)
c$omp threadprivate (/failuloc/)
c$omp threadprivate (/hourgloc/)
c$omp threadprivate (/shloptloc/)
c$omp threadprivate (/soundloc/)
c$omp threadprivate (/subtssloc/)
c$omp threadprivate (/vect8loc/)
c$omp threadprivate (/bux13loc/)
c$omp threadprivate (/frozloc/)
      mx=48*(mxt(lft)-1)
c
c     Determining number of requested history variables
c
      no_hsvs=nconstp(mxt(lft))
      if (lenvec(5) .eq.1) no_hsvs=no_hsvs-6
      if (iorien(mt).ne.0) no_hsvs=no_hsvs-2
      if (ihyper(mt).ne.0) no_hsvs=no_hsvs-9
      if (ifipts(mt).lt.0) no_hsvs=no_hsvs-1
c
c     initializing indeces
c
      if (iorien(mt).ne.0) then
        if (ihyper(mt).ne.0) then
          ind_defgrad=1+no_hsvs
          ind_ortho=ind_defgrad+9
        else
          ind_ortho=1+no_hsvs
        endif
        if (ifipts(mt).lt.0) then
          ind_fipts=ind_ortho
          ind_ortho=ind_ortho+1
        endif
        ind_q1=ind_ortho
        ind_q2=ind_ortho+1
        ind_last=ind_q2
      else
        if (ihyper(mt).ne.0) then
          ind_defgrad=1+no_hsvs
          ind_last=ind_defgrad+8
        else
          ind_last=no_hsvs
        endif
        if (ifipts(mt).lt.0) then
          ind_fipts=ind_last+1
          ind_last=ind_last+1
        endif
      endif
c
      if (iorien(mt).ne.0) then
c
c       store element local stresses
c
        if (ihyper(mt).ge.0) then
          do i=lft,llt
            sigmaold(i,1)=sig1(i)
            sigmaold(i,2)=sig2(i)
            sigmaold(i,3)=sig3(i)
            sigmaold(i,4)=sig4(i)
            sigmaold(i,5)=sig5(i)
            sigmaold(i,6)=sig6(i)
          enddo
c
c         transform stresses to material system
c
          do 20 i=lft,llt
          stg5(i)=sig5(i)
          stg6(i)=sig6(i)
          a11(i) =hsvs(i,ind_q1)*sig1(i)-hsvs(i,ind_q2)*sig4(i)
          a12(i) =hsvs(i,ind_q2)*sig1(i)+hsvs(i,ind_q1)*sig4(i)
          a21(i) =hsvs(i,ind_q1)*sig4(i)-hsvs(i,ind_q2)*sig2(i)
          a22(i) =hsvs(i,ind_q2)*sig4(i)+hsvs(i,ind_q1)*sig2(i)
          sig1(i)=hsvs(i,ind_q1)*a11(i)-hsvs(i,ind_q2)*a21(i)
          sig2(i)=hsvs(i,ind_q2)*a12(i)+hsvs(i,ind_q1)*a22(i)
          sig4(i)=hsvs(i,ind_q1)*a12(i)-hsvs(i,ind_q2)*a22(i)
          sig5(i)=hsvs(i,ind_q2)*stg6(i)+hsvs(i,ind_q1)*stg5(i)
          sig6(i)=hsvs(i,ind_q1)*stg6(i)-hsvs(i,ind_q2)*stg5(i)
 20       continue
        endif
c
c       storage of deformation gradient
c
        if (ihyper(mt).gt.0) then
          if11=ind_defgrad
          if21=ind_defgrad+1
          if31=ind_defgrad+2
          if12=ind_defgrad+3
          if22=ind_defgrad+4
          if32=ind_defgrad+5
          if13=ind_defgrad+6
          if23=ind_defgrad+7
          if33=ind_defgrad+8
c
c         transform deformation gradient to material system
c         as F_bar=Q^T*F
c
          do i=lft,llt
c
c           store deformation gradient in element system
c
            dfg11(i)=hsvs(i,if11)
            dfg21(i)=hsvs(i,if21)
            dfg12(i)=hsvs(i,if12)
            dfg22(i)=hsvs(i,if22)
            dfg13(i)=hsvs(i,if13)
            dfg23(i)=hsvs(i,if23)
c
            ftemp=hsvs(i,ind_q1)*hsvs(i,if11)-
     1            hsvs(i,ind_q2)*hsvs(i,if21)
            hsvs(i,if21)=hsvs(i,ind_q2)*hsvs(i,if11)+
     1            hsvs(i,ind_q1)*hsvs(i,if21)
            hsvs(i,if11)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if12)-
     1            hsvs(i,ind_q2)*hsvs(i,if22)
            hsvs(i,if22)=hsvs(i,ind_q2)*hsvs(i,if12)+
     1            hsvs(i,ind_q1)*hsvs(i,if22)
            hsvs(i,if12)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if13)-
     1            hsvs(i,ind_q2)*hsvs(i,if23)
            hsvs(i,if23)=hsvs(i,ind_q2)*hsvs(i,if13)+
     1            hsvs(i,ind_q1)*hsvs(i,if23)
            hsvs(i,if13)=ftemp
          enddo
c
c       storage of deformation gradient
c
        else if (ihyper(mt).lt.0) then
          if11=ind_defgrad
          if21=ind_defgrad+1
          if31=ind_defgrad+2
          if12=ind_defgrad+3
          if22=ind_defgrad+4
          if32=ind_defgrad+5
          if13=ind_defgrad+6
          if23=ind_defgrad+7
          if33=ind_defgrad+8
c
c         transform deformation gradient to material system
c         as F_bar=F*Q
c
          do i=lft,llt
c
c           store deformation gradient in element system
c
            dfg11(i)=hsvs(i,if11)
            dfg21(i)=hsvs(i,if21)
            dfg31(i)=hsvs(i,if31)
            dfg12(i)=hsvs(i,if12)
            dfg22(i)=hsvs(i,if22)
            dfg32(i)=hsvs(i,if32)
c
            ftemp=hsvs(i,ind_q1)*hsvs(i,if11)-
     1            hsvs(i,ind_q2)*hsvs(i,if12)
            hsvs(i,if12)=hsvs(i,ind_q2)*hsvs(i,if11)+
     1            hsvs(i,ind_q1)*hsvs(i,if12)
            hsvs(i,if11)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if21)-
     1            hsvs(i,ind_q2)*hsvs(i,if22)
            hsvs(i,if22)=hsvs(i,ind_q2)*hsvs(i,if21)+
     1            hsvs(i,ind_q1)*hsvs(i,if22)
            hsvs(i,if21)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if31)-
     1            hsvs(i,ind_q2)*hsvs(i,if32)
            hsvs(i,if32)=hsvs(i,ind_q2)*hsvs(i,if31)+
     1            hsvs(i,ind_q1)*hsvs(i,if32)
            hsvs(i,if31)=ftemp
          enddo
        endif
      endif
c
c     memory pointer for extra memory data
      lcma=ncma+nmmat+ia(ncma+mxt(lft)-1)-1
c
c     vector parameter
c
      ivect=ivectr(mt)
c
      mte40=mt-40
c
      if (ivect.eq.0) then
c
c       scalar utan
c
        do 120 i=lft,llt
c
c       gathering of variables
c
        index=i
        dt1   =dt1siz(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
        epsp=epsps(i)
        temper=temps(i)
        failel=.false.
        if (ifipts(mt).lt.0) failel=hsvs(i,ind_fipts).ne.0.
        if (no_hsvs.gt.0) then
          do j=1,no_hsvs
            hsv(j)=hsvs(i,j)
          enddo
        endif
        if (ihyper(mt).ne.0) then
          do j=1,9
            hsv(no_hsvs+j)=hsvs(i,ind_defgrad-1+j)
          enddo
          if (iabs(ihyper(mt)).eq.3) then
             do k=1,3
                do j=1,3
                   qmats(j,k)=qmat(i,j,k)
                enddo
             enddo
          endif
        endif
c
c       setting material tangent modulus to = 0
c
        do k=1,6
          do j=1,6
            es(j,k)=0.
          enddo
        enddo
c
c       call user developed subroutines here
c
        go to (41,42,43,44,45,46,47,48,49,50), mte40
 41     call utan41 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),failel,a(lcma),qmats)
        go to 60
 42     call utan42 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),failel,a(lcma),qmats)
        go to 60
 43     call utan43 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),failel,a(lcma),qmats)
        go to 60
 44     call utan44 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),failel,a(lcma),qmats)
        go to 60
 45     call utan45 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),failel,a(lcma),qmats)
        go to 60
 46     call utan46 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),failel,a(lcma),qmats)
        go to 60
 47     call utan47 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),failel,a(lcma),qmats)
        go to 60
 48     call utan48 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),failel,a(lcma),qmats)
        go to 60
 49     call utan49 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),failel,a(lcma),qmats)
        go to 60
 50     call utan50 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),failel,a(lcma),qmats)
c
c       scattering material tangent modulus
c
 60     do k=1,6
          do j=1,6
           dsave(i,j,k)=es(j,k)
         enddo
        enddo
 120    continue
c
c     vector utan
c
      else
c
c       setting material tangent modulus to = 0
c
        do k=1,6
          do j=1,6
            do i=lft,llt
              dsave(i,j,k)=0.
            enddo
          enddo
        enddo
c
        if (ifipts(mt).ge.0) then
          do i=lft,llt
            failels(i)=.false.          
          enddo
        else
          do i=lft,llt
            failels(i)=hsvs(i,ind_fipts).ne.0.
          enddo
        endif
c
c       call user developed routines here
c       first history variable is the one requested
c
        if (mt.eq.41) then
          call utan41v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.42) then
          call utan42v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.43) then
          call utan43v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.44) then
          call utan44v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.45) then
          call utan45v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.46) then
          call utan46v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.47) then
          call utan47v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.48) then
          call utan48v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.49) then
          call utan49v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.50) then
          call utan50v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.286) then
          call utan286v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),failels,
     .     a(lcma))
        elseif (mt.eq.287) then
          call utan287v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),failels,
     .     a(lcma))
        endif
      endif
c
      if (iorien(mt).ne.0) then
c
c       transform tangent stiffness to element system
c
        if (ihyper(mt).ge.0) then
          call transform44s(hsvs(1,ind_q1),hsvs(1,ind_q2),lft,llt)
c
c         restore element local stresses
c
          do i=lft,llt
            sig1(i)=sigmaold(i,1)
            sig2(i)=sigmaold(i,2)
            sig3(i)=sigmaold(i,3)
            sig4(i)=sigmaold(i,4)
            sig5(i)=sigmaold(i,5)
            sig6(i)=sigmaold(i,6)
          enddo
        endif
c
c       restore deformation gradient in element system
c
        if (ihyper(mt).ne.0) then
          do i=lft,llt
            hsvs(i,if11)=dfg11(i)
            hsvs(i,if21)=dfg21(i)
            hsvs(i,if12)=dfg12(i)
            hsvs(i,if22)=dfg22(i)
            hsvs(i,if13)=dfg13(i)
            hsvs(i,if23)=dfg23(i)
          enddo
        endif
      endif
c
      return
      end
      subroutine utan50v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan49v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan48v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan47v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan46v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan42v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan44v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan43v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan45v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      logical failels(*)
      character*5 etype
c
      dimension sig(6),eps(6),hsv(NHISVAR),es(6,6)
      logical failel
c
      do i=lft,llt
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
c
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
c
        hsv(1)=hsvs(i,1)
        hsv(2)=hsvs(i,2)
        hsv(3)=hsvs(i,3)
        hsv(4)=hsvs(i,4)
        hsv(5)=hsvs(i,5)
        hsv(6)=hsvs(i,6)
        hsv(7)=hsvs(i,7)
        hsv(8)=hsvs(i,8)
        hsv(9)=hsvs(i,9)
c
        do k=1,6
          do j=1,6
            es(j,k)=0.
          enddo
        enddo
c
        failel=failels(i)
c
        call utan45(cm,eps,sig,epsps(i),hsv,dt1siz(i),capa,etype,tt,
     1   temps(i),es,crv,failel,cma,0)
c
        do k=1,6
          do j=1,6
            dsave(i,j,k)=es(j,k)
          enddo
        enddo
      enddo
c
      return
      end
      subroutine utan41v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      logical failels(*)
      character*5 etype
c
      dimension sig(6),eps(6),hsv(NHISVAR),es(6,6)
      logical failel
c
      do i=lft,llt
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
c
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
c
        do k=1,6
          do j=1,6
            es(j,k)=0.
          enddo
        enddo
c
        failel=failels(i)
c
        call utan41(cm,eps,sig,epsps(i),hsv,dt1siz(i),capa,etype,tt,
     1   temps(i),es,crv,failel,cma,0)
c
        do k=1,6
          do j=1,6
            dsave(i,j,k)=es(j,k)
          enddo
        enddo
      enddo
c
      return
      end
      subroutine utan41(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      factor=1.
      if (failel) factor=1.e-8
c
      g=factor*.5*abs(cm(1))/(1.+cm(2))
      b=factor*abs(cm(1))/3./(1.-2.*cm(2))
      bg23=b-2.*g/3.
      bg43=b+4.*g/3.
c
      es(1,1)=bg43
      es(2,2)=bg43
      es(3,3)=bg43
      es(2,1)=bg23
      es(3,1)=bg23
      es(3,2)=bg23
      es(1,2)=es(2,1)
      es(1,3)=es(3,1)
      es(2,3)=es(3,2)
      es(4,4)=g
      es(5,5)=g
      es(6,6)=g
c
      return
      end
      subroutine utan45(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*)
      dimension es(6,*),crv(lq1,2,*),cma(*),qmat(3,3)
      logical failel
      character*5 etype
c
c     no history variables, NHV=0
c     deformation gradient stored in hsv(1),...,hsv(9)
c
c     compute jacobian
c
      detf=hsv(1)*(hsv(5)*hsv(9)-hsv(6)*hsv(8))
     1     -hsv(2)*(hsv(4)*hsv(9)-hsv(6)*hsv(7))
     2     +hsv(3)*(hsv(4)*hsv(8)-hsv(5)*hsv(7))
c
c     compute lame parameters
c
      xlambda=cm(1)*cm(2)/((1.+cm(2))*(1.-2.*cm(2)))
      xmu=.5*cm(1)/(1.+cm(2))
c
c     compute tangent stiffness
c     same for both shells and solids
c
      detf=max(detf,1.e-8)
      detfinv=1./detf
      dmu=xmu-xlambda*log(detf)
      es(1,1)=detfinv*(xlambda+2.*dmu)
      es(2,2)=detfinv*(xlambda+2.*dmu)
      es(3,3)=detfinv*(xlambda+2.*dmu)
      es(4,4)=detfinv*dmu
      es(5,5)=detfinv*dmu
      es(6,6)=detfinv*dmu
      es(2,1)=detfinv*xlambda
      es(3,2)=detfinv*xlambda
      es(3,1)=detfinv*xlambda
      es(1,2)=es(2,1)
      es(2,3)=es(3,2)
      es(1,3)=es(3,1)
c
      return
      end
      subroutine utan43(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      return
      end
      subroutine utan44(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      return
      end
      subroutine utan42(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      return
      end
      subroutine utan46(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      return
      end
      subroutine utan47(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      return
      end
      subroutine utan48(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      return
      end
      subroutine utan49(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      return
      end
      subroutine utan50(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
      return
      end
      subroutine uweldfail(idweld,strr,fail,fibl,cm,tt,lft,llt)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     spotweld failure routine
c
c***  local coordinates: x is tangent to beam, y & z are normal
c
c***  variables
c          idweld ---- weld id number
c          strr ------ stress resultants
c                      (1) axial (x direction) force
c                      (2) y shear force
c                      (3) z shear force
c                      (4) moment about z
c                      (5) moment about y
c                      (6) torsional resultant
c          fail ------ failure flag
c                      =0 not failed
c                      =1 fail element
c          fibl ------ fibl(1,i)=diameter of weld
c          cm -------- 6 constants supplied by user
c          tt -------- current simulation time
c          lft,llt --- do-loop range for strr
c
      dimension idweld(*),strr(8,*),fail(*),cm(*),fibl(5,*)
c
c***  demonstration failure model
c
c***  geometry
c
c***  failure constants
      if (cm(1).gt.0.) then
        oamax=1./cm(1)
      else
        oamax=0.
      endif
      if (cm(2).gt.0.) then
        osmax=1./cm(2)
      else
        osmax=0.
      endif
c
c***  evaluate failure criterion:
c       axial stress = axial force/area + moment/sectionmod
c       shear stress = shear force/area + torsion/(2.*sectionmod)
c       failure surface = (axial stress/cm(1))**2 +
c                         (shear stress/cm(2))**2
      do i=lft,llt
      Z       =3.14159265359*fibl(1,i)**3/32
      osctmod =1.0/Z
      hosctmod=.5*osctmod
      oarea   =4./(3.14159265359*fibl(1,i)**2)
        axial=max(0.,strr(1,i))*oarea+
     &          sqrt(strr(4,i)**2+strr(5,i)**2)*osctmod
        shear=sqrt(strr(2,i)**2+strr(3,i)**2)*oarea+
     &        abs(strr(6,i))*hosctmod
        if ((axial*oamax)**2+(shear*osmax)**2.gt.1.0) fail(i)=1.
      enddo
c
      return
      end
      subroutine uweldfail22(idweld,strr,rscale,fail,fibl,cm,tt,lft,llt)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     sample user defined spotweld failure routine with 22 parameters
c
c
c***  local coordinates: x is tangent to beam, y & z are normal
c
c***  variables
c          idweld ---- weld id number
c          strr stress resultants at beam ends 1 and 2
c                      (1) axial (x direction) force
c                      (2) y shear force
c                      (3) z shear force
c                      (4) moment about z
c                      (5) moment about y
c                      (6) torsional resultant
c          fail ------ failure flag
c                      =0 not failed
c                      =1 fail element
c          fibl ------ fibl(1,i)=diameter of weld
c                      fibl(2,i)=pointer to rscale data
c          rscale ---- resultant force scale factor for
c                      *CONTROL_SPOTWELD_BEAM effects
c                      (1) node 1 tensile factor
c                      (2) node 2 tensile factor
c                      (3) node 1 shear factor
c                      (4) node 2 shear factor
c          cm -------- up to 22 constants supplied by user
c          tt -------- current simulation time
c          lft,llt --- do-loop range for strr
c
c     *** NOTE *** for solid element welds, the rscale factors are not
c                  useful and are set to 1.0
c
      dimension idweld(*),strr(8,*),rscale(4,*),fail(*),cm(*),fibl(5,*)
c
c***  demonstration failure model
c
c***  failure constants
      if (cm(1).gt.0.) then
        oamax=1./cm(1)
      else
        oamax=0.
      endif
      if (cm(2).gt.0.) then
        osmax=1./cm(2)
      else
        osmax=0.
      endif
c
c***  evaluate failure criterion:
c       axial stress = axial force/area + moment/sectionmod
c       shear stress = shear force/area + torsion/(2.*sectionmod)
c       failure surface = (axial stress/cm(1))**2 +
c                         (shear stress/cm(2))**2
      do i=lft,llt
      Z       =3.14159265359*fibl(1,i)**3/32
      osctmod =1.0/Z
      hosctmod=.5*osctmod
      oarea   =4./(3.14159265359*fibl(1,i)**2)
      axial=max(0.,strr(1,i))*oarea+
     &        sqrt(strr(4,i)**2+strr(5,i)**2)*osctmod
      shear=sqrt(strr(2,i)**2+strr(3,i)**2)*oarea+
     &      abs(strr(6,i))*hosctmod
      j=nint(fibl(2,i))
      axial1=axial*rscale(1,j)
      axial2=axial*rscale(2,j)
      shear1=shear*rscale(3,j)
      shear2=shear*rscale(4,j)
      if ((axial1*oamax)**2+(shear1*osmax)**2.gt.1.0) fail(i)=1.
      if ((axial2*oamax)**2+(shear2*osmax)**2.gt.1.0) fail(i)=1.
      enddo
c
      return
      end
      subroutine uweldfail12(eps_to_fail,beta,dmg,uhisv,num_uhisv,
     . idweld,iwtype,strr,eps,esr,npt,diameter,cm,idmgopt,tt,ncycle,
     . lft,llt)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     Sample user defined spotweld damage and failure routine with 22
c     parameters; this routine can be used to damage or fail beam welds
c     and hex weld assemblies.
c
c     There are two ways to use it.  The first way is to do a damage
c     initiation calculation and set eps_to_fail and beta when damage
c     should initiate.  Based on the values set (see below) either
c     immediate failure, or linear or exponential damage and failure
c     based on effective plastic strain will occur.  The second way is
c     for the user to do their own dmaage calculation and update the
c     dmg variable thoughout the calculation.
c           
c    
c     subroutine output variables (may be updated by user)
c
c        eps_to_fail  Plastic strain from damage initiation until
c                     failure.  A zero value means that damage has not
c                     yet begun.  Any negative value will cause
c                     immediate failure.  In the cycle when eps_to_fail
c                     is set, the current effective plastic strain will
c                     be stored and used as the strain at which damage
c                     inititiates (eps_di).  Note that eps_to_fail and
c                     beta can only be set only one time for each
c                     integration point and subsequent changes will be
c                     ignored.
c
c        beta ....... Set to positive number for exponential
c                     damage.  Set to zero for linear damage.
c
c                     For beta=0, the decay function is
c
c                        dmg=(eps-eps_di)/eps_to_fail
c
c                      For beta>0, the decay function is
c
c                        dmg=1.-exp(-beta*(eps-eps_di)/eps_to_fail)
c
c                     For exponential damage, the damage will grow
c                     according to the exponential equation above until
c                     the plastic strain reaches eps_to_fail at which
c                     time the weld will be failed.
c
c        uhisv ...... Weld history values.  These values will be
c                     retained, but will not be used for the
c                     calculation.
c
c     subroutine input/output variable (may be either calculated by
c     ls-dyna or updated by the user)
c     
c
c        dmg ........ Current damage at the integration point.
c                     Initially, dmg=0.  When dmg=1.0, the integration
c                     point will fail.  When all integration points
c                     fail, the element or assembly will be deleted.
c                     If the user wants ls-dyna to calculate damage
c                     using eps_to_fail and beta, then damage should
c                     not be changed.  Optionally, the user may update
c                     the dmg variable throughout the calculation and
c                     the user's value will be applied.
c                    
c
c     subroutine input variables
c
c        num_uhisv .. number of weld history variables (will be
c                     equal to FVAL on card 3 of *MAT_SPOTWELD)
c        idweld ..... weld id number
c
c        iwtype ..... weld type
c                      1=beam element
c                      2=solid element assembly
c
c        strr ....... stress resultants
c                       (1) axial (x direction) force
c                       (2) y shear force
c                       (3) z shear force
c                       (4) moment about z
c                       (5) moment about y
c                       (6) torsional mement
c         
c        eps ........ effective plastic strain at the integration point
c                     (or hex weld of an assembly)
c
c        esr ........ effective strain rate at the integration point
c                     (or hex weld of an assembly)
c
c        npt ........ number of integration points in the weld; for hex
c                     weld assemblies, npt is the number of elements
c                     in the assembly since each element has 1 integration
c                     point 
c
c        diameter ... diameter of weld
c        cm ......... up to 22 constants supplied by user
c        idmgopt .... damage option parameter DMGOPT
c        tt ......... current simulation time
c        ncycle ..... current cycle number
c        lft,llt .... do-loop range
c
c
c     For hex assmeblies, some variables may be dumped to d3plot using
c     NEIPH on *DATABASE_EXTENT_BINARY.  The table below shows the
c     variable number for each.
c
c        variable name    extra history variable #
c
c             dmg                   8
c            eps_di                10
c          eps_to_fail             11
c             beta                 12
c             esr                  13
c
      dimension eps_to_fail(*),beta(*),dmg(npt,*),uhisv(num_uhisv,*),
     . idweld(*),strr(8,*),eps(npt,*),esr(npt,*),diameter(*),cm(*)
c
c***  demonstration failure model
c
c***  failure constants
      if (cm(1).gt.0.) then
        oamax=1./cm(1)
      else
        oamax=0.
      endif
      if (cm(2).gt.0.) then
        osmax=1./cm(2)
      else
        osmax=0.
      endif
c
c***  evaluate damage initiation criterion:
c       axial stress = axial force/area + moment/sectionmod
c       shear stress = shear force/area + torsion/(2.*sectionmod)
c       failure surface = (axial stress/cm(1))**2 +
c                         (shear stress/cm(2))**2
      do i=lft,llt
      if (eps_to_fail(i).eq.0.) then
        Z       =3.14159265359*diameter(i)**3/32
        osctmod =1.0/Z
        hosctmod=.5*osctmod
        oarea   =4./(3.14159265359*diameter(i)**2)
c
        axial=max(0.,strr(1,i))*oarea+
     &         sqrt(strr(4,i)**2+strr(5,i)**2)*osctmod
        shear=sqrt(strr(2,i)**2+strr(3,i)**2)*oarea+
     &         abs(strr(6,i))*hosctmod
c
        if ((axial*oamax)**2+(shear*osmax)**2.gt.1.0) then
          eps_to_fail(i)=cm(3)
          beta(i)=cm(4)
        endif
      endif
      enddo
      return
      end
      subroutine hdot(lft,llt,hmin,hmax,iactive,e,h,dt,const,s)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     user integration in time of the smoothing length
c
c     lft      : left indice of the loop
c     llt      : right indice of the loop
c     hmin     : Minimum smoothing length defined in the corresponding
c                *SECTION_SPH
c     hmax     : Maximum smoothing length defined in the corresponding
c                *SECTION_SPH
c     iactive  : variable that indicates if the particle is active or not
c                = 1, particle is active
c                = 0, particle is not active
c     e        : tensor of gradient of deformation
c     h        : smoothing length
c     dt       : current time step
c     const    : constant = dt/dimension of the problem
c     s        : stress tensor
c
      dimension h(1),iactive(1),e(9,1),ro(1),smthi(1),matsph(1)
      dimension s(6,1)
c
c$omp paralleldo default(none)
c$omp& shared(lft,llt,iactive,e,h,const,dt1,hmin,hmax)
c$omp& private(i,divv,hnew)
c
      do 200 i=lft,llt
      divv = e(1,i) + e(2,i) + e(3,i)
c
      hnew = (1.+const*divv)*h(i)
      h(i) = max(hmin,hnew)
      h(i) = min(h(i),hmax)
 200  continue
c
      return
      end
      subroutine usrhcon(h,uc,nc,d,p,xs,xm,vs,vm,ts,tm,
     . istria,imtria,crv)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     User subroutine for interface conductivity control (mpp/smp)
c
c     To define thermal contact conductivity
c
c     Variables:
c
c     h      = thermal contact conductivity (output)
c     uc(nc) = array of user defined parameters (input)
c     nc     = number of user defined conductivity parameters (input)
c     d      = distance (gap) between slave and master segment (input)
c     p      = interface pressure (input)
c     xs     = slave segment coordinates (input)
c     xm     = master segment coordinates (input)
c     vs     = slave segment velocities (input)
c     vm     = master segment velocities (input)
c     ts     = slave segment temperatures (input)
c     tm     = master segment temperatures (input)
c     istria = 1 iff slave segment is a triangle (input)
c     imtria = 1 iff master segment is a triangle (input)
c     crv    = curve array (input)
c
c     Note that usrvrel can be called to compute relative sliding velocity
c
      include 'nlqparm'
      dimension uc(nc)
      dimension xs(3,4),xm(3,4)
      dimension vs(3,4),vm(3,4)
      dimension crv(lq1,2,*)
c
c     Define conductivity as function of pressure via load curve
c
      call crvval(crv,uc(1),p,h,dh)
c
      return
      end
      subroutine usrvrel(vrel,xs,vs,vm,istria,imtria)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     compute sliding velocity vrel and that can be called from usrhcon.
c
c     Variables:
c     vrel    = relative sliding velocity (output)
c     xs      = slave segment coordinates (input)
c     vs      = slave segment velocities (input)
c     vm      = master segment velocities (input)
c     istria  = 1 iff slave segment is a triangle (input)
c     imtria  = 1 iff master segment is a triangle (input)
c
      dimension xs(3,4),vs(3,4),vm(3,4)
      dimension xnrm(3),vr(3)
c
      jmax=4
      if (istria.eq.1) jmax=3
      do i=1,3
         vr(i)=0.
         do j=1,jmax
            vr(i)=vr(i)+vs(i,j)
         enddo
         vr(i)=vr(i)/jmax
      enddo
      jmax=4
      if (imtria.eq.1) jmax=3
      do i=1,3
         vrl=0.
         do j=1,jmax
            vrl=vrl+vm(i,j)
         enddo
         vr(i)=vr(i)-vrl/jmax
      enddo
      x31=xs(1,3)-xs(1,1)
      y31=xs(2,3)-xs(2,1)
      z31=xs(3,3)-xs(3,1)
      x42=xs(1,4)-xs(1,2)
      y42=xs(2,4)-xs(2,2)
      z42=xs(3,4)-xs(3,2)
      xnrm(1)=y31*z42-z31*y42
      xnrm(2)=z31*x42-x31*z42
      xnrm(3)=x31*y42-y31*x42
      xlen=sqrt(xnrm(1)**2+xnrm(2)**2+xnrm(3)**2)
      xnrm(1)=xnrm(1)/xlen
      xnrm(2)=xnrm(2)/xlen
      xnrm(3)=xnrm(3)/xlen
      vnrm=vr(1)*xnrm(1)+vr(2)*xnrm(2)+vr(3)*xnrm(3)
      vr(1)=vr(1)-vnrm*xnrm(1)
      vr(2)=vr(2)-vnrm*xnrm(2)
      vr(3)=vr(3)-vnrm*xnrm(3)
      vrel=sqrt(vr(1)**2+vr(2)**2+vr(3)**2)
      return
      end
      subroutine usrflux(fl,flp,x,tnpl,tnl,nodes,
     . alpha,atime,atemp,dt,time,fhsv,nfhsv,crv,v)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     User subroutine for boundary thermal flux
c
c     To define thermal flux parameter (heat per surface area and time)
c
c     Variables:
c
c     fl              = flux intensity (output)
c     flp             = flux intensity derivative wrt atemp (output)
c     x(3,nodes)      = global segment coordinates (input)
c     v(3,nodes)      = global segment velocities (input)
c     tnpl(nodes)     = temperatures at time time (input)
c     tnl(nodes)      = temperatures at time time-dt (input)
c     nodes           = number of nodes in segment (3,4 or 6) (input)
c     alpha           = time integration parameter (input)
c     atime           = time+(alpha-1)*dt
c     atemp           = average segment temperature at time atime
c     dt              = time step size (input)
c     time            = time at which the new temperature is sought (input)
c     fhsv(nfhsv)     = flux history variables (input/output)
c     nfhsv           = number of flux history variables for this segment (input)
c     crv             = curve array (input)
c
      include 'nlqparm'
      dimension x(3,*),tnpl(*),tnl(*),v(3,*)
      dimension fhsv(*),crv(lq1,2,*),func(9)
c---------------------------------------------------------------------
c Example 1
c     Define flux by linear convection
c     that optionally decays (in an ad-hoc way) as power
c     dissipates from surface
c
c     fhsv(1) = convection coefficient
c     fhsv(2) = ambient temperature
c     fhsv(3) = total amount of energy per surface area available
c     fhsv(4) = dissipated energy per surface area at current
c
c---- uncomment following lines for example 1
c     hcon=fhsv(1)
c     tinf=fhsv(2)
c     flin=hcon*(tinf-atemp)
c     if (nfhsv.gt.2) then
c        q=(1.-fhsv(4)/fhsv(3))/
c    .        (1.+.5*dt*flin/fhsv(3))
c        flp=-q*hcon
c        if (q.gt.1.) then
c           q=1.
c           flp=-hcon
c        elseif (q.lt.0.) then
c           q=0.
c           flp=0.
c        endif
c        fl=q*flin
c        fhsv(4)=fhsv(4)+dt*.5*fl
c        fhsv(4)=min(fhsv(3),fhsv(4))
c     else
c        fl=flin
c        flp=-hcon
c     endif
c     fl=-fl
c     flp=-flp
c---------------------------------------------------------------------
c Example 2
c     Define flux using a function ID
c
c     fhsv(1) =  function ID
c
      fl=0.
      flp=0.
c-- get segment centroid coordinates
      func(1)=(x(1,1)+x(1,2)+x(1,3)+x(1,4))/4.
      func(2)=(x(2,1)+x(2,2)+x(2,3)+x(2,4))/4.
      func(3)=(x(3,1)+x(3,2)+x(3,3)+x(3,4))/4.
c-- get segment average velocity
      func(4)=(v(1,1)+v(1,2)+v(1,3)+v(1,4))/4.
      func(5)=(v(2,1)+v(2,2)+v(2,3)+v(2,4))/4.
      func(6)=(v(3,1)+v(3,2)+v(3,3)+v(3,4))/4.
c-- get segment average temperature
      func(7)=(tnl(1)+tnl(2)+tnl(3)+tnl(4))/4.
c-- T_infinity is not a typical input for flux, set it to zero
      func(8)=0.
c-- get time
      func(9)=atime
c-- the first history variable for this segment is a function id
c-- make it an integer
      idff=int(fhsv(1))
c-- convert user function ID to the internal numbering system
      internal_f=ind_gf(idff)
c-- check to see if internal_f is valid. If so, then evaluate function
      if (internal_f .gt. 0) then
        fl=eval_function(internal_f,func,ierr)
        flp=0
      endif
c
      return
      end
      subroutine ujntfrc(idfrc,idjnt,const,t,dt,statevar,disp,vel,frc)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c***  user subroutine for applying forces to joint degrees of freedom
c
c***  variables
c          idfrc ----- ID for the joint force
c          idjnt ----- ID for the joint
c          const ----- array of user supplied constants
c          t --------- current time
c          dt -------- time step size
c          statevar -- array of user integrated state variables
c          disp ------ array of relative joint displacements
c          vel ------- array of relative joint velocities
c          frc ------- array of forces and torques calculated by user
c
      dimension const(*),statevar(*),disp(*),vel(*),frc(*)
c
c***  example for a spring and damper for a revolute joint
c     torque =  c*angular_velocity + k*(angle-offset)
      frc(1)=const(1)*vel(1)+const(2)*(disp(1)-const(3))
c
c***  statevariable integration: statevar(1)=0.5*t**2
      statevar(1)=statevar(1)+dt*t
c
      return
      end
      subroutine ujntfdrv(njntf,idat,kdat,idatjf,kdatjf,const,
     &                    statevar,q,rbv,rbf,x)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'iounits.inc'
c
c***  driver for the user joint force subroutines
c
c***  variables
c          njntf ----- number of joint forces
c          idatjf ---- integer user joint force data
c                      1) ID
c                      2) joint ID
c                      3) number of state variables
c                      4) joint index
c                      5) rigid body 1
c                      6) rigid body 2
c                      7) pointer into state variables
c          idat ------ 1) joint type
c                      2:7) nodes
c          const ----- user constants
c          statevar -- user state variables
c          q --------- integrated displacements
c          rbv ------- rigid body velocities
c          rbf ------- rigid body forces
c
c***  overall structure
c
c       qdot = B*{rigid body velocities}
c       q    = q + dt*qdot
c       f    = user_function(q,qdot,t,dt)
c       f_rb = f_rb - B^T*f
c       ener = ener + dt*f*qdot
c
      common/bk02/iburn,dt1,dt2,isdo,tet10jtol
      common/bk28/summss,xke,xpe,tt,xte0
c
      real*8 x
      character*19 atype(14)
      dimension idat(kdat,*),idatjf(kdatjf,*),const(48,*),
     &          statevar(*),q(6,*),rbv(6,*),rbf(6,*),x(3,*)
c
      data atype/'spherical          ',
     2           'revolute           ',
     3           'cylindrical        ',
     4           'planar             ',
     5           'universal          ',
     6           'translational      ',
     7           'locking            ',
     8           'translational motor',
     9           'rotational motor   ',
     .           'gears              ',
     1           'rack and pinion    ',
     2           'constant velocity  ',
     3           'pulley             ',
     4           'screw              '/
c
c***  over joint forces
      do 1000 ijntf=1,njntf
c
c***    the basics
        idjntf=idatjf(1,ijntf)
        idjnt =idatjf(2,ijntf)
        indjnt=idatjf(4,ijntf)
        irb1  =idatjf(5,ijntf)
        irb2  =idatjf(6,ijntf)
        ipsv  =idatjf(7,ijntf)
        jnttyp=idat(indjnt,1)
c
c***    branch on joint type
        go to (  10,  20,  30,  40,  50,  60,  70,
     &           80,  90, 100, 110, 120, 130, 140), jnttyp
c
c***    spherical
   10   go to 911
c
c***    revolute
   20   call rjufrc(idjntf,idjnt,const(1,ijntf),tt,dt1,
     &              statevar(ipsv),q(1,ijntf),
     &              rbv(1,irb1),rbf(1,irb1),rbv(1,irb2),rbf(1,irb2),
     &              x,idat,kdat,indjnt)
        go to 1000
c
c***    cylindrical
   30   go to 911
c
c***    planar
   40   go to 911
c
c***    universal
   50   go to 911
c
c***    translational
   60   go to 911
c
c***    locking
   70   go to 911
c
c***    tranlational motor
   80   go to 911
c
c***    rotational motor
   90   go to 911
c
c***    gears
  100   go to 911
c
c***    rack and pinion
  110   go to 911
c
c***    constant velocity
  120   go to 911
c
c***    pulley
  130   go to 911
c
c***    screw
  140   go to 911
c
c***    error message -- should never get here
c 911   write(iohsp,1010) idjntf,atype(jnttyp),idjnt
c       write(iomsg,1010) idjntf,atype(jnttyp),idjnt
c       write(iotty,1010) idjntf,atype(jnttyp),idjnt
  911   continue
        ierdat(1)=idjntf
        cerdat(1)=atype(jnttyp)
        ierdat(2)=idjnt
        call lsmsg(3,MSG_SOL+1153,ioall,ierdat,rerdat,cerdat,0)
c
 1000 continue
c
      return
c1010 format(1x/
c    &  '*** Error user joint force ID',i12,' for ',a19,i12/
c    &  '          this type of joint is currently not supported',
c    &  '          for user joint forces.')
      end
      subroutine rjufrc(idfrc,idjnt,const,t,dt,statevar,q,
     &                  rbv1,rbf1,rbv2,rbf2,x,idat,kdat,indjnt)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c***  revolute joint user force driver
c
c***  variables
c          idfrc ----- id of user force
c          idjnt ----- id of joint
c          const ----- user constants
c          t --------- current time
c          dt -------- time step size
c          statevar -- user state variables
c          q --------- displacements
c          rbv1,2 ---- rigid body velocities
c          rbf1,2 ---- rigid body forces
c          x --------- coordinates
c          idat ------ list of joint nodes
c          indjnt ---- joint indext
c
      real*8 x
      dimension const(*),statevar(*),q(*),
     &          rbv1(*),rbf1(*),rbv2(*),rbf2(*),x(3,*),idat(kdat,*)
c
c***  nodes
      n1=idat(indjnt,2)
      n2=idat(indjnt,3)
      n3=idat(indjnt,4)
      n4=idat(indjnt,5)
c
c***  angular velocity
      dx11=x(1,n3)-x(1,n1)
      dx21=x(2,n3)-x(2,n1)
      dx31=x(3,n3)-x(3,n1)
      omag=1.0/sqrt(dx11**2+dx21**2+dx31**2)
      dx11=omag*dx11
      dx21=omag*dx21
      dx31=omag*dx31
c
      dx12=x(1,n4)-x(1,n2)
      dx22=x(2,n4)-x(2,n2)
      dx32=x(3,n4)-x(3,n2)
      omag=1.0/sqrt(dx12**2+dx22**2+dx32**2)
      dx12=omag*dx12
      dx22=omag*dx22
      dx32=omag*dx32
c
      angv=rbv2(4)*dx12+rbv2(5)*dx22+rbv2(6)*dx32
     &    -rbv1(4)*dx11-rbv1(5)*dx21-rbv1(6)*dx31
c
c***  angular displacement integrated in time
c     copy so user can't overwrite our history variable value
      q(1)=q(1)+dt*angv
      ang=q(1)
c
c***  call user routine
      call ujntfrc(idfrc,idjnt,const,t,dt,statevar,ang,angv,torq)
c
c***  add torq to rhs
      rbf2(4)=rbf2(4)-dx12*torq
      rbf2(5)=rbf2(5)-dx22*torq
      rbf2(6)=rbf2(6)-dx32*torq
      rbf1(4)=rbf1(4)+dx11*torq
      rbf1(5)=rbf1(5)+dx21*torq
      rbf1(6)=rbf1(6)+dx31*torq
c
      return
      end
      subroutine matfailusercontrol(user)
      dimension user(*)
c
c     user interface to for user failure subroutine
c     when iuserfail equals to 1, user has the option to access stress, 
c     and history information for each element
c     Then user can make any change with 'subroutine matuserfail'
c     
c
      iuserfail=0
      if(iuserfail.eq.1) then
      call matfailusercontrol1(user)
      end if
      return
      end
      subroutine matuserfail(stress,shisv,nip,nplane,nhis)
      dimension stres(8,4,9),shisv(90,4,9)
c
c     user can access stress and history information now 
c     so, user's failure model can be defined here.
c
c     no more than 9 integration point in the thickness
c
c     stress(1,*,*): zeta (local coordinate) in the thickness direction
c     stress(8,*,*): effective strain
c
c     shisv(n,*,*) the nth history variable in one integration
c
c     Stress and history can be obtained in the following format
c     do i = 1, nplane
c     do j = 1, nip
c       write(*,'(t1,8e10.3)') (stres(k,i,j),k=1,8)
c       do m = 1, nhis,8
c         n = min(nhis,m+7)
c         write(*,'(t1,8e10.3)') (shisv(k,i,j),k=m,n)
c       end do
c     end do
c     end do
      return
      end
      subroutine rwumat(iphase,txt,lmc,ipnt,lmca,ipnta,idfoff,
     . fctmas,fcttim,fctlen,fcttem,incout)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'memaia.inc'
      include 'iodeg.inc'
c
c     In this routine, material parameter cards are read from
c     *MAT_USER_DEFINED_MATERIAL_MODELS in phase 1. 
c     In phase 2, the parameters are written as real numbers
c     to the structured input (dyna.str).
c
c     Main variables:
c       iphase : phase flag
c       txt    : incoming text
c       lmc    : number of material constants
c       ipnt   : pointer to material constants
c       lmca   : number of additional material constants
c       ipnta  : pointer to additional material constants
c
c     Offset and transformation factors from *INCLUDE_TRANSFORM 
c     can be used in phase 1:
c       idfoff : func id offset number, i.e., load curves and tables
c       fctmas : mass transformation factor
c       fcttim : time transformation factor
c       fctlen : length transformation factor
c       fcttem : temperature transformation factor
c       incout : the transformed included file will be written to dyna.inc
c                =0, no output file
c                =1, dyna.inc includes the complete transformed input
c
      character*200 txt
      character*80 mssg
      character*10 out10(8)
      dimension prop(48)
c
      isexample=0
c
c---  PHASE 1: read keyword input
c
c     In phase 1, the read format can be defined arbitrarily.
c     Parameters can be stored as real a(...) or integer ia(...)
c     as shown in the example. 
c
      if (iphase.eq.1) then
        mssg='readin *MAT_USER_DEFINED_MATERIAL_MODELS'
c
c       default read material constants: LMC
        if (isexample.eq.0) then
          do i=1,lmc,8
            iend=min(i+7,lmc)
            call gttxln (txt,lcount,2,8)
            read (txt,'(8e10.0)',err=900) (a(ipnt+k),k=i,iend)
          enddo
c
c       example: read 5 parameters
c       p1  - real with no unit 
c       p2  - real with time unit
c       p3  - real with stress unit
c       ip4 - integer with no unit
c       ip5 - load curve id (stress vs. strain)
        elseif (isexample.eq.1) then
          call gttxln (txt,lcount,2,8)
          read (txt,'(3e10.0,2i10)',err=900) p1,p2,p3,ip4,ip5
          p2=p2*fcttim
          fctfor=fctmas*fctlen/(fcttim*fcttim)
          fctpre=fctfor/(fctlen*fctlen)
          p3=p3*fctpre
          if (ip5.ne.0) then
            ip5=ip5+idfoff
c           adtbmdu : Add TaBle or Curve to be MoDified User-defined
c           sfa = scale factor for abscissa value
c           sf0 = scale factor for ordinate value
            sfa=1.0
            sfo=fctpre
            call adtbmdu(ip5,sfa,sfo,ipnt,ipnta,lmca)
          endif
          a(ipnt+1)=p1
          a(ipnt+2)=p2
          a(ipnt+3)=p3
          ia(ipnt+4)=ip4
          ia(ipnt+5)=ip5
        endif
c
c       read extra material constants: LMCA
        if (lmca.gt.0) then
          do i=1,lmca,8
            iend=min(i+7,lmca)
            call gttxln (txt,lcount,2,8)
            read (txt,'(8e10.0)',err=900) (a(ipnta+k),k=i,iend)
          enddo
        endif
c
c---  PHASE 2: write structured input
c
c     In phase 2, parameters must be transformed to real numbers,
c     before writing them to the structured input. Therefore,
c     external load curve ids should be changed to internal ids
c     as shown in the example. 
c
      elseif (iphase.eq.2) then
c
c       default write 48 parameters (LMC max)
        if (isexample.eq.0) then
          do i=1,48,8
            ii=min(48-i+1,8)
            do j=1,ii
              call form10(out10(j),a(ipnt+i+j-1))
            enddo
            write(io,'(8a10)') (out10(j),j=1,ii)
            call wttxsg(1)
          enddo
c
c       example: write 5 parameters and 43 zeros (sum=48!)
c       p1  - real with no unit 
c       p2  - real with time unit
c       p3  - real with stress unit
c       ip4 - integer with no unit
c       ip5 - load curve id (stress vs. strain)
        elseif (isexample.eq.1) then
          do i=1,48
            prop(i)=a(ipnt+i)
          enddo
          prop(4)=dble(ia(ipnt+4))
c         change load curve id to internal number with lcidi
          lc=ia(ipnt+5)
          prop(5)=dble(lcidi(lc))
          do i=1,48,8
            ii=min(48-i+1,8)
            do j=1,ii
c             change to special output format with form10
              call form10(out10(j),prop(i+j-1))
            enddo
            write(io,'(8a10)') (out10(j),j=1,ii)
            call wttxsg(1)
          enddo
        endif
c
c       write extra material constants: LMCA
        if (lmca.gt.0) then
          do i=1,lmca,8
            ii=min(lmca-i+1,8)
            write(io,'(1p8e10.3)') (a(ipnta+i+j-1),j=1,ii)
            call wttxsg(1)
          enddo
        endif
c
      endif
c
      return
  900 call terminn(txt,mssg,lcount,1)
      end
      subroutine usermatfpert ( x, y, z, fail )
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     user specification of material failure
c
c     Inputs:
c     x, y, z: coordinates
c
c     Returned:
c     fail:  failure criterion (stress, yield)
c
      fail = x
c
      return
      end
      subroutine chknhue(nhsv)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'iounits.inc'
      include 'nhisparm.inc'
c
      if (nhsv.gt.NHISVUE) then
        write(iotty,100) NHISVUE
        write(iohsp,100) NHISVUE
        write(iomsg,100) NHISVUE
      endif
      nhsv=min(NHISVUE,nhsv)
  100 format(/' *** Warning number of user element history variables',
     .       /'             is limited by NHISVUE in nhisparm.inc:',i10)
c
      return
      end
      subroutine chkusercomm
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     check the length of common block size
c
      include 'nlqparm'
      include 'nhisparm.inc'
      include 'iounits.inc'
c
      parameter (NBEL8A=max(100,2*(NHISVUE+45)))
      parameter (NBEL7 =max(2601,
     . 1180+576*NXDOFUE+64*NXDOFUE*NXDOFUE+12*5))
      common/aux14loc/ax14zz1(nlq,7+NHISVAR)
      common/bel8loc/bel8lz1(nlq,2*(NHISVUE+45))
      common/bel8aloc/bel8lz2(nlq,NBEL8A)
      common/bel7loc/bel7lz(nlq,NBEL7)
      common/usercom/iusercomm(4)
c
c$omp threadprivate (/aux14loc/)
c$omp threadprivate (/bel8loc/)
c$omp threadprivate (/bel8aloc/)
c$omp threadprivate (/bel7loc/)
      l14=nlq*(7+NHISVAR)
      l8 =nlq*2*(NHISVUE+45)
      l8a=nlq*NBEL8A
      l7 =nlq*NBEL7
c
      if (l14.gt.iusercomm(1)) then
        ierdat(1)=l14
        ierdat(2)=iusercomm(1)
        cerdat(1)='/aux14loc/'
        call lsmsg(2,MSG_OTH+81,iotty,ierdat,rerdat,cerdat,0)
        call cstop('E r r o r   t e r m i n a t i o n')
      endif
      if (l8.gt.iusercomm(2)) then
        ierdat(1)=l8
        ierdat(2)=iusercomm(2)
        cerdat(1)='/bel8loc/'
        call lsmsg(2,MSG_OTH+81,iotty,ierdat,rerdat,cerdat,0)
	call cstop('E r r o r   t e r m i n a t i o n')
      endif
      if (l8a.gt.iusercomm(3)) then
        ierdat(1)=l8a
        ierdat(2)=iusercomm(3)
        cerdat(1)='/bel8aloc/'
        call lsmsg(2,MSG_OTH+81,iotty,ierdat,rerdat,cerdat,0)
        call cstop('E r r o r   t e r m i n a t i o n')
      endif
      if (l7.gt.iusercomm(4)) then
        ierdat(1)=l7
        ierdat(2)=iusercomm(4)
        cerdat(1)='/bel7loc/'
        call lsmsg(2,MSG_OTH+81,iotty,ierdat,rerdat,cerdat,0)
        call cstop('E r r o r   t e r m i n a t i o n')
      endif
c
      return
      end

