      program prob_bd
c Likelihood code to be run on a mock galaxy. Input:catalogue of PNe and input file parameters and fi file. Check if table.dat is clipping as the confidence level you desire (default 0.04). 


      parameter(m=100000,e=2.718281828)
      parameter(pi=3.141592654)

      double precision id(m),ra(m),dec(m),vel(m),gi(m),cat(m),fb(m)   
      double precision l(m),v(m),sv(m),sr(m),ss(m),rmax(m),rs(m)
      double precision xp,yp,xss,yss,ps,Rc,Rmc,Rsc,Dc,Dmc,Dsc,pa,incl
      double precision vned, i_rad,lcut(m),ppm(m),xm(m),ym(m),sz_c(m)  
      double precision cos_phi(m),sin_phi(m),xs(m),ys(m),ysi(m)
      double precision vvs(m),sigmas(m),v_halos(m),bd(m),bb(m),xsi(m)
      double precision mm(m),mxl(m),mxv(m),mxs(m),mxsv(m),mxh(m)
      double precision mxf(m),rgal(m),vhel_av(m),rgalgc(m),b(m),mxvh(m)  
      double precision likes,loglike,likes_av,loglikes_av,sg_c(m)
      double precision xs_c(m),ys_c(m),ysi_c(m),vel_c(m),bin_c(m)
      double precision fb_c(m),gi_c(m),cat_c(m),likes_c(m),comp(m)
      double precision mmm(m),vhel_all(m),vel_all_av(m),bi(m),ri(m)
      double precision xs_ca(m),ys_ca(m),ysi_ca(m),vel_ca(m),sg(m)
      double precision fb_ca(m),bba(m),bda(m),gi_ca(m),cat_ca(m),ftot(m)
      double precision bin_ca(m),likes_ca(m),i_ca(m),i_c(m),sz(m)
      double precision logbda(m),logbba(m),cos_phica(m),comp_c(m)
      double precision ram(m),ras(m),decm(m),decs(m),rah(m),dech(m)
      double precision errb(m),errr(m),sm(m),sm_c(m),ra_c(m),dec_c(m)
      double precision tot_red,tot_red_disk,tot_red_bulge,tot_blue
      double precision tot_blue_disk,tot_blue_bulge,ll(m)
      

c INPUT & OUTPUT FILES
      
c     catalogue (without :) created with set_catalogue.sm
       open(11,file="2768GC_complete_clean_dec.dat",status="old")
c     f
       open(18,file="N2768_f.GC.dat",status="old")
c     input created with create_input.f
       open(21,file="input.N2768.gc.dat",status="old")
c     rejection 0.04 nearly 2 sigma
       open(56,file="table.dat",status="unknown")
c    likelihood from pne
c       open(10,file="likelihood.bulge.dat",status="unknown")
c bin
c       open(15,file="bin.dat",status="unknown")


c     coordinate (xs,ys,x,y,ysi)
c       open(44,file="coord.dat",status="unknown")
c     rgal(h), h, raggio di ogni bin 
       open(15,file="bin_2768.dat",status="unknown")
c     coordinate and vel v_av, xs, ysi
c       open(12,file="vprofile.dat",status="unknown")
c     likelihood
       open(10,file="likelihood_2768.dat",status="unknown")
       
c     data used in the program before and after sigma-clipping
c       open(66,file="cleancatalogue.dat",status="unknown")
c     errors
c       open(77,file='error.dat',status='unknown')
c       open(57,file="fdisk.dat",status="unknown")

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     reading parameters

        xp=0
        yp=0
        Rc=0
        Rmc=0
        Rsc=0
        Dc=0
        Dmc=0
        Dsc=0
        pa=0
        incl=0
        vned=0    
        Racc=0
        Decc=0
        nbin=0

      read(21,*)xp,yp,xss,yss,ps,Rc,Rmc,Rsc,Dc,Dmc,Dsc,pa,incl,vned,nbin

      ai=incl
 
      pa_rad=(pa)*pi/180      

       i_rad=ai*pi/180

       write(*,*)' i_rad=',i_rad,' i_deg=',ai, 'pa= ', pa
       write(*,*) xp,yp,Rc,Rmc,Rsc,Dc,Dmc,Dsc

       cos_i=COS(i_rad)
       sin_i=SIN(i_rad)

       write(*,*)'cos_i=',cos_i,' sin_i=',sin_i

       write(*,*)'log(e)', log(e)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     probability

       do j=1,10
          read(56,*)lcut(j)
          write(*,*)lcut(j)
       enddo

 56    close(56)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       write(*,*)'starting main galaxy'


       n=0

       do i=1,m

c          Read the file containin the fi values


          read(18,*,end=118)fb(i)

c   Read coordinate ra and dec in arcsec

         read(11,*,end=1111)comp(i),ra(i),dec(i),vel(i),b(i),ri(i)

          n=n+1
          ppm(i)=n

        gi(i)=b(i)-ri(i)
        write(*,*)'gi',gi(i)

     

          Racc=Rc*3600 + Rmc*60 + Rsc
          Decc=Dc*3600 + Dmc*60 + Dsc

c          write(*,*)Racc, Decc

         xm(i)=0
c         ra(i)=(rah(i)+ram(i)/60+ras(i)/3600)*15
          xm(i)=( ra(i)/15*3600 )-Racc
c           xm(i)=( ra(i))-Racc
    
         ym(i)=0
c         dec(i)=dech(i)+decm(i)/60+decs(i)/3600
          ym(i)=dec(i)*3600 - Decc
c         ym(i)=dec(i) - Decc
      
          
          y_rad=Decc*pi/(3600*180)
         

          xm(i)=xm(i)*15*cos(y_rad)
c            xm(i)=xm(i)*15
c            open(56,file="1023_PNe_f.cat",status="unknown")
c            write(56,556)ra(i),dec(i),vel(i),fb(i)
            enddo
c 556        format(4(1x,f8.3))

 118   close(18)
 1111  close(11)
       close(56)



c     rotating the galaxi in order to have the major axis on the y axis and check the rotation...

      do i=1,n

c     minor axis
          xs(i)=xm(i)*COS(pa_rad) - ym(i)*SIN(pa_rad)
c     major axis
          ys(i)=xm(i)*SIN(pa_rad) + ym(i)*COS(pa_rad)


c     average velocity & sigma
        
          av=av+vel(i)
          av2=av2+vel(i)**2
        
       enddo

      av=av/n
      
     
       av2=av2/n
       sigma_los=sqrt(av2-av**2)

       write(*,*) ' av=',av,' av2=',av2,' sig=',sigma_los,' n=',n
       write(*,*)'v_ned=',vned

c     obtaining the inclination angle and correct coordinates

       do i=1,n
          
          
          xsi(i)=xs(i)*sqrt(cos_i)
          ysi(i)=ys(i)/sqrt(cos_i)
          rs(i)=sqrt(xsi(i)**2+ysi(i)**2)

c     obtain azimuthal angle in the galaxy plain

          cos_phi(i)=xsi(i)/rs(i)
          sin_phi(i)=ysi(i)/rs(i)
          vhel_av(i)=vel(i)-vned
          b(i)=rs(i)
       enddo

       do i=1,n
          do j=i+1,n
             if(b(i).gt.b(j))then
                a=b(i)
                b(i)=b(j)
                b(j)=a
             endif
          enddo                           
       enddo
                                                  
c       do i=1,n
c          write(*,*)b(i),i,xs(i),ys(i)
c       enddo

   
c GC obtaining prob for each GC
 
           open(17,file="unliky.dat",status="unknown")
           open(19,file="all.dat",status="unknown")
           open(20,file="for_lodo.dat",status="unknown")
           open(22,file="bin_vb.dat",status="unknown")

           do h=1,nbin
c          do i=1,n
             rgalgc(h)=b(int((h*n/nbin)))
c            enddo
              write(*,*)'rgalgc', rgalgc(h)
              write(22,*) rgalgc(h)
             enddo


               likes=0
               tot_gc=0
              
               do i=1,n
                  vvs(i)=0
               enddo

               do i=1,n
                  sigmas(i)=0
                  v_halos(i)=0
                  bd(i)=0
                  bb(i)=0   
                  xs_c(i)=0
                  ys_c(i)=0
                  ysi_c(i)=0
                  vel_c(i)=0
                  fb_c(i)=0
                  gi_c(i)=0
c                  cat_c(i)=0
                  likes_c(i)=0
                  bin_c(i)=0
                  comp_c(i)=0
                  sg_c(i)=0
                  sz_c(i)=0
                  sm_c(i)=0
                  ra_c(i)=0
                  dec_c(i)=0
               enddo




                 
               cont=0
               cont1=0

               do h=1,nbin

                  read(10,*)mxl(h),mxv(h),mxs(h),mxsv(h),mxh(h),mxf(h),
     &            s,s,mm(h),s,mxvh(h)

                  read(15,*)rgal(h)
                  write(*,*)'rgal',rgal(h)

                  if(rgal(nbin).le.rgalgc(nbin))rgal(nbin)=rgalgc(nbin)

                  loglike=0
                  likes=0
                  likes_av=0
                  mm(h)=0

                  do i=1,n

                       loglike=0
c                       write(*,*)'rs',rs(i)

                     if((int(rgal(h-1)).lt.int(rs(i))).and.
     &                    (int(rs(i)).le.int(rgal(h))))then

                        cont=cont+1

                        mm(h)=mm(h)+1

                        sigmas(i)=sqrt((((mxs(h)*sin_i)**2)*
     &                       (sin_phi(i)**2))+
     &                       ((((mxsv(h)*sin_i))**2)*
     &                       (cos_phi(i)**2))) 

                      


                        vvs(i)=(((vhel_av(i)-(mxv(h)*sin_i
     &                       *(cos_phi(i))))**2))/(2*(sigmas(i)**2))

                        
                        v_halos(i)=((vhel_av(i)-(mxvh(h)*
     &                               sin_i*cos_phi(i)))**2)/
     &                               (2*(mxh(h)**2))

c                        v_halos(i)=((vhel_av(i)**2)/
c     &                              (2*(mxh(h)**2)))


                           likes=(((1-fb(i))*exp(-(vvs(i)))
     &                          /sigmas(i))+(fb(i)*
     &                          exp(-(v_halos(i)))/mxh(h)))

                           loglike=log(((1-fb(i))*exp(-(vvs(i)))
     &                          /sigmas(i))+(fb(i)*
     &                          exp(-(v_halos(i)))/mxh(h)))

                           bb(cont)=log(fb(i)*
     &                          exp(-(v_halos(i)))/mxh(h))

                           bd(cont)=log((1-fb(i))*exp(-(vvs(i)))
     &                          /sigmas(i))

                           xs_c(cont)=-xs(i)
                           ys_c(cont)=ys(i)
                           ysi_c(cont)=ysi(i)
                           vel_c(cont)=vel(i)
                           fb_c(cont)=fb(i)
                           gi_c(cont)=gi(i)
c                           cat_c(cont)=cat(i)
                           likes_c(cont)=loglike
                           bin_c(cont)=h
                           i_c(cont)=i
                           comp_c(cont)=comp(i)
                           sg_c(cont)=sg(i)
                           sz_c(cont)=sz(i)
                           sm_c(cont)=sm(i)
                           ra_c(cont)=ra(i)
                           dec_c(cont)=dec(i)

c                           do j=400,1200

c                              write(*,*)'i',i
c
                                          
c                  bda(j)=0
c                  bba(j)=0   
c                  vel_ca(j)=0
c                  fb_ca(j)=0
c                  likes_ca(i)=0
c                  bin_ca(i)=0

              


c                              vhel_all(j)=j
c                              vel_all_av(j)=vhel_all(j)-vned
c                              loglike=0
c                              cont1=cont1+1
c                              mmm(h)=mmm(h)+1

c                        sigmas(j)=sqrt((((mxs(h)*sin_i)**2)*
c     &                       (sin_phi(i)**2))+
c     &                       ((((mxsv(h)*sin_i))**2)*
c     &                       (cos_phi(i)**2))) 

                      


c                        vvs(j)=(((vel_all_av(j)-(mxv(h)*sin_i
c     &                       *(cos_phi(i))))**2))/(2*(sigmas(i)**2))

                       

c                        v_halos(j)=((vel_all_av(j)**2)/
c     &                              (2*(mxh(h)**2)))


c                           likes=(((1-fb(i))*exp(-(vvs(j)))
c     &                          /sigmas(j))+(fb(i)*
c     &                          exp(-(v_halos(j)))/mxh(h)))

c                           loglike=log(((1-fb(i))*exp(-(vvs(j)))
c     &                          /sigmas(j))+(fb(i)*
c     &                          exp(-(v_halos(j)))/mxh(h)))

c                           bba(cont1)=(fb(i)*
c     &                          exp(-(v_halos(j)))/mxh(h))/likes
c                           logbba(cont1)=log(fb(i)*
c     &                          exp(-(v_halos(j)))/mxh(h))
                         


c                           bda(cont1)=((1-fb(i))*exp(-(vvs(j)))
c     &                          /sigmas(j))/likes
c                           logbda(cont1)=log((1-fb(i))*exp(-(vvs(j)))
c     &                          /sigmas(j))

c
c                           vel_ca(cont1)=vhel_all(j)
c                           cos_phica(cont1)=cos_phi(i)
c                           fb_ca(cont1)=fb(i)
c                           likes_ca(cont1)=loglike
c                           bin_ca(cont1)=h
c                          i_ca(cont1)=i
c                           enddo

 
                              endif


                             
               
               enddo

             
 
                 

              
               
 
                write(*,*)mm(h),mmm(h)
                tot_gc=mm(h)+tot_gc

               enddo

c               write(*,*)'here'

               write(*,*)n,tot_gc,cont,cont1

               open(66,file="bin_rej.dat",status="unknown")

               do h=1,nbin
                  write(*,*)rgal(h),rgalgc(h),mm(h)
                  read(66,*)ll(h)
              
                  enddo

                     tot_red_disk=0
                     tot_red=0
                     tot_red_bulge=0
                     tot_blue_disk=0
                     tot_blue=0
                     tot_blue_bulge=0

                     ntot_gc=0
                     n_rej=0

                   
                 do i=1,n
c                     write(*,*)gi_c(i)
                    write(17,771)ra_c(i),dec_c(i),xs_c(i),ys_c(i),
     &               ysi_c(i),vel_c(i),
     &              fb_c(i),bb(i),bd(i),gi_c(i),likes_c(i),sg_c(i),
     &              sz_c(i),comp_c(i),bin_c(i),i_c(i)


                   bb(i)=e**(bb(i))
                   bd(i)=e**(bd(i))
                   ftot(i)=(bb(i)/(bb(i)+bd(i)))
                   
                   
                   if(comp_c(i).ge.0) ntot_gc=ntot_gc+1

              do h=1,nbin
                 
                


                     if((int(rgal(h-1)).lt.int(rs(i))).and.
     &                    (int(rs(i)).le.int(rgal(h))))then

              if(likes_c(i).le.ll(h))comp_c(i)=5
c.and.(comp_c(i).eq.1))comp_c(i)=5

                  if(comp_c(i).eq.5)n_rej=n_rej+1

                   if(comp_c(i).ne.5)then
                      write(*,*)'gi_c',gi_c(i)
                      if(gi_c(i)>=0.55)then
                      tot_red_disk=tot_red_disk+(1-ftot(i))
                      tot_red=tot_red+1
                      tot_red_bulge=tot_red_bulge+ftot(i)
                      endif
                         if(gi_c(i)<=0.55)then
                          tot_blue_disk=tot_blue_disk+(1-ftot(i))
                          tot_blue=tot_blue+1
                          tot_blue_bulge=tot_blue_bulge+ftot(i)
                          endif
                          endif
                          endif

c                   write(*,*)likes_c(i),comp_c(i)

c                   if(likes_c(i).le.-6.85)comp_c(i)=5


                     write(20,772)ra_c(i),dec_c(i),xs_c(i),ys_c(i),
     &              vel_c(i),
     &              fb_c(i),ftot(i),gi_c(i),comp_c(i),bin_c(i),i_c(i),
     &               likes_c(i),sm_c(i)

                    enddo
                    enddo

                   write(*,*)'red=',tot_red,'red_d=',tot_red_disk,
     &                       'red_bulge=',tot_red_bulge,
     &                        'blue=',tot_blue,'blue_d=',tot_blue_disk,
     &                       'blue_bulge=',tot_blue_bulge,
     &                   'tot=',ntot_gc,'rej=',n_rej,
     &                    ntot_gc-n_rej,tot_red+tot_blue
c                     do i=1,n
c                   write(17,771)xs_c(i),ys_c(i),ysi_c(i),vel_c(i),
c     &              fb_c(i),bb(i),bd(i),likes_c(i),bin_c(i),i_c(i)
c                    enddo


c                  do i=1,cont1
c                    write(19,778)vel_ca(i),fb_ca(i),likes_ca(i),
c     &              bin_ca(i),i_ca(i),logbba(i),logbda(i),cos_phica(i)
c                    enddo

                    






 771           format(2(1x,f8.4),14(1x,f8.2))
 772          format(2(1x,f8.4),11(1x,f8.2))
c 778           format(8(1x,f8.2))
               
               close(10)
               close(15)
               close(19)
               close(17)

              stop
              end
 
