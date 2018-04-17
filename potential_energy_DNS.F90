Program potential_energy
         use variables
         implicit none
        
    integer     NY_GRID
    character*35   file_th_xyz
    character*38   file_th_1d
!        real*8        TH_3D(1:NX,1:NZ,1:NY_TOT) ,KAPPA_T_3D_TOT(1:NX,1:NZ,1:NY_TOT),&
!                      TH_1D(1:NY_3D_TOT_TEMP),TH_TOT_3D(1:NX,1:NZ,1:NY_TEMP), &
!                      tmp1,tmp2,DRHO_DZSTAR,tmp3,DELRHO2_TEMP,tmp4,KAPPA_T_TEMP
                 
REAL*8         TEMP1(1:((NY_3D_TOT_TEMP+1)/2)),TEMP2(1:((NY_3D_TOT_TEMP+1)/2)), &
               TEMP3(1:((NY_3D_TOT_TEMP+1)/2)),TEMP4(1:((NY_3D_TOT_TEMP+1)/2))

! REAL*8    GYF(0:NY_TOT+1), GY(0:NY_TOT+1),  &
!              DYF(0:NY_TOT+1), DY(0:NY_TOT+1), &
!              DV(1:NX,1:NZ,1:NY_TEMP),    &
!              DV_1D(1:NY_3D_TOT_TEMP),        &
!              Z_STAR_1D(1:NY_3D_TOT_TEMP),    &
!              Z_1D(1:NY_3D_TOT_TEMP),   &
!              DISP(1:NY_3D_TOT_TEMP)

REAL*8  DELRHOX2(1:NX,1:NZ,1:NY_TEMP),DELRHOZ2(1:NX,1:NZ,1:NY_TEMP), &
        DELRHOY2(1:NX,1:NZ,1:NY_TEMP), DELRHO2(1:NX,1:NZ,1:NY_TEMP), &
        DELRHO2_1D(1:NY_3D_TOT_TEMP),PHI_D,PHI_D_TEMP,PHI_D_NEW,PHI_D_TEMP_NEW,ZSTAR_TEMP_2,TH_TEMP, &
!        KAPPA_T_1D(1:NY_3D_TOT_TEMP)
        FLUXTERM_X(1:NX,1:NZ,1:NY_TEMP),FLUXTERM_Z(1:NX,1:NZ,1:NY_TEMP), &
        FLUXTERM_Y(1:NX,1:NZ,1:NY_TEMP),FLUXTERM(1:NX,1:NZ,1:NY_TEMP), &
        FLUXTERM_1D(1:NY_3D_TOT_TEMP)

        CHARACTER*47 filn_vel_yz,filn_vel_xy
   
         pi = 4.0d0 * ATAN(1.0d0)
         write(6,*) 'PI = ', pi,'NX =', NX, ' NY = ', NY_TOT, 'NZ= ',NZ

        angle_beta  = beta*pi/180.0d0
        angle_alpha = alpha*pi/180.0d0
        N_infty= omega0/sin(angle_alpha)
        drho_dz= -N_infty*N_infty*rho_0/gravity_g
!        drho_dz = -1.0d0


        LAMBDA_X = 1.63287d0

         K_X_R = (2.0d0*PI/ LAMBDA_X)*(cos(angle_beta)+ &
             sin(angle_beta)/tan(angle_alpha))

 !      LX =  2.0d0*(2.0d0*PI/K_X_R)
        LX = 5.D0

        NU = 1.0E-6
        PR = 7.0D0

        LZ = 1.0d0
        DX = LX/NX
        DZ = LZ/NZ

         write(6,*) 'LX =', LX , 'Ri_code', N_infty**2.0d0, (NX+1)*DX


        open(21,file='ygrid.txt',form='formatted',status='old')
        READ(21,*) NY_GRID
        IF(NY_GRID .EQ. NY_TOT) THEN
         write(6,*) 'Grids matched'
        END IF 

         DO J=1,NY_TOT+1
           READ(21,*) GY(j)
         END DO
         DO J=1,NY_TOT
           READ(21,*) GYF(j)
         END DO
         CLOSE(21)


! Define ghost cells
         GYF(0)=2.d0*GYF(1)-GYF(2)
         GYF(NY_TOT+1)=2.d0*GYF(NY_TOT)-GYF(NY_TOT-1)
         GY(0)=2.d0*GY(1)-GY(2)

! Define grid spacing
         DO J=1,NY_TOT+1
           DY(J)=(GYF(J)-GYF(J-1))
         END DO
         DO J=0,NY_TOT
           DYF(J)=(GY(J+1)-GY(J))
         END DO

         DO J=1,NY_TOT
          DO I=1,NX
           xpoint(I,J) = DX*(i-1)
           ypoint(I,J) = gyf(j)
          ENDDO
          DO k=1,NZ
           z_point(k,J) = DZ*(k-1)
           y_point(k,J) = gyf(j)
          ENDDO
         ENDDO
 
!        DO J=1,NY_TOT
!          DO K=1,NZ
!           DO I=1,NX
!            ypoint3D(I,K,J) = gyf(j)
!            xpoint3D(I,K,J) = DX*(i-1)
!            zpoint3D(I,K,J) = DZ*(k-1)
!          ENDDO
!         ENDDO
!        ENDDO



        DO J=1,NY_TOT
         DO I=1,NX
           xpoint_m(I,J) = xpoint(I,J)*cos(angle_beta) - ypoint(I,J)*sin(angle_beta)
           ypoint_m(I,J) = xpoint(I,J)*sin(angle_beta) + ypoint(I,J)*cos(angle_beta)
         ENDDO
        ENDDO


!        DO K=1,NZ
!         DO J=1,NY_TOT
!          DO I=1,NX
!          xpoint3D_m(I,K,J) = xpoint_m(I,J)
!          ypoint3D_m(I,K,J) = ypoint_m(I,J)
!          ENDDO
!         ENDDO
!        ENDDO

        COUNTER = 0

        DO T_S=T_S_bin,TIME_STEP,100


         COUNTER = COUNTER + 100
         write(*,*) 'Now beginning timestep', T_S

            file_th_xyz = 'density_data/data_th_'  &
              //CHAR(MOD(COUNTER,10000000)/1000000+48)  &
              //CHAR(MOD(COUNTER,1000000)/100000+48)    &
              //CHAR(MOD(COUNTER,100000)/10000+48)      &
              //CHAR(MOD(COUNTER,10000)/1000+48)        &
              //CHAR(MOD(COUNTER,1000)/100+48)          &
              //CHAR(MOD(COUNTER,100)/10+48)            &
              //CHAR(MOD(COUNTER,10)+48) //             &
              '.pln'

              file_th_1d = 'density_data_1d/data_th_'  &
              //CHAR(MOD(COUNTER,10000000)/1000000+48)  &
              //CHAR(MOD(COUNTER,1000000)/100000+48)    &
              //CHAR(MOD(COUNTER,100000)/10000+48)      &
              //CHAR(MOD(COUNTER,10000)/1000+48)        &
              //CHAR(MOD(COUNTER,1000)/100+48)          &
              //CHAR(MOD(COUNTER,100)/10+48)            &
              //CHAR(MOD(COUNTER,10)+48) //             &
              '.pln'

         open(1011,file= file_th_xyz,status='old',form='unformatted')


         read(1011) TIME
         read(1011)(((TH_3D(I,K,J),i=1,NX),k=1,NZ),j=1,NY_TOT)


         close(1011) 
             

!!!!!!!!!!!!!!!  READING KAPPA_T !!!!!!!!!!!!!!!!!

!         filn_vel_xy = 'add_pln_xy/data_vel_xy_'  &
!              //CHAR(MOD(COUNTER/100,10000000)/1000000+48)  &
!              //CHAR(MOD(COUNTER/100,1000000)/100000+48)    &
!              //CHAR(MOD(COUNTER/100,100000)/10000+48)      &
!              //CHAR(MOD(COUNTER/100,10000)/1000+48)        &
!              //CHAR(MOD(COUNTER/100,1000)/100+48)          &
!              //CHAR(MOD(COUNTER/100,100)/10+48)            &
!              //CHAR(MOD(COUNTER/100,10)+48) //             &
!              '.pln'

!         open(401,file= filn_vel_xy,status='old',form='unformatted')


!         read(401) TIME_TEMP
!         read(401) (TEMP_2(j),j=1,NY_TOT),             &
!                   ((TEMP_1(i,j),i=1,NX),j=1,NY_TOT),   &
!                   ((TEMP_1(i,j),i=1,NX),j=1,NY_TOT),   &
!                   ((TEMP_1(i,j),i=1,NX),j=1,NY_TOT),   &
!                   ((TEMP_1(i,j),i=1,NX),j=1,NY_TOT),   &
!                   ((TEMP_1(i,j),i=1,NX),j=1,NY_TOT),   &
!                   (TEMP_2(j),j=1,NY_TOT),   &
!                   (TEMP_2(j),j=1,NY_TOT),   &
!                   (TEMP_2(j),j=1,NY_TOT),   &
!                   (TEMP_2(j),j=1,NY_TOT),   &
!                   (TEMP_2(j),j=1,NY_TOT),   &
!                   (KAPPA_T_MEAN(j),j=1,NY_TOT),   &
!                   (TEMP_2(j),j=1,NY_TOT),   &
!                   (TEMP_2(j),j=1,NY_TOT),   &
!                   (TEMP_2(j),j=1,NY_TOT),   &
!                   (TEMP_2(j),j=1,NY_TOT),   &
!                   (TEMP_2(j),j=1,NY_TOT),   &
!                   (TEMP_2(j),j=1,NY_TOT),   &
!                   (TEMP_2(j),j=1,NY_TOT),   &
!                   (TEMP_2(j),j=1,NY_TOT),   &
!                   (TEMP_2(j),j=1,NY_TOT), &
!                   (TEMP_2(j),j=1,NY_TOT), &
!                   (TEMP_2(j),j=1,NY_TOT), &
!                   (TEMP_2(j),j=1,NY_TOT), &
!                   (TEMP_2(j),j=1,NY_TOT), &
!                   (TEMP_2(j),j=1,NY_TOT), &
!                   (TEMP_2(j),j=1,NY_TOT), &
!                   (TEMP_2(j),j=1,NY_TOT), &
!                   (TEMP_2(j),j=1,NY_TOT),  &
                  ! (TEMP_2(j),j=1,NY_TOT) 

!              IF (TIME .NE. TIME_TEMP) THEN
!                  write(*,*) 'Error detected'
!                     stop 
!              ENDIF



         DO J=1,NY_TOT
          DO I=1,NX
          TH_B(i,j) = 1000.0 + (ypoint_m(i,j)-ypoint_m(NX,NY_TOT))*drho_dz
!          TH_B(i,j) = (ypoint_m(i,j)-ypoint_m(NX,NY_TOT))*drho_dz
         ENDDO
        ENDDO
         
         DO K=1,NZ
          DO J=1,NY_TEMP
           DO I=1,NX
            TH_TOT_3D(I,K,J)  = TH_B(i,j) + TH_3D(I,K,J)
!            TH_TOT_3D(I,K,J)  = TH_3D(I,K,J) 
         ENDDO
        ENDDO
       ENDDO          
        

 
 
          DO k=1,NZ
           DO J=1,NY_TEMP
            DO i=1,NX
              DV(I,K,J) = DYF(J)*DX*DZ

               IF(I .EQ. 1) THEN
                DELRHOX2(I,K,J) = ((TH_TOT_3D(I+1,K,J)-TH_TOT_3D(I,K,J))/DX)**2.0D0 
               ELSEIF(I .EQ. NX) THEN
                DELRHOX2(I,K,J) = ((TH_TOT_3D(I,K,J)-TH_TOT_3D(I-1,K,J))/DX)**2.0D0
               ELSE
                DELRHOX2(I,K,J) = ((TH_TOT_3D(I+1,K,J)-TH_TOT_3D(I-1,K,J))/(2.0D0*DX))**2.0D0
               ENDIF
  
          
               IF(K .EQ. 1) THEN
                DELRHOZ2(I,K,J) = ((TH_TOT_3D(I,K+1,J)-TH_TOT_3D(I,K,J))/DZ)**2.0D0  
               ELSEIF(K .EQ. NZ) THEN
                DELRHOZ2(I,K,J) = ((TH_TOT_3D(I,K,J)-TH_TOT_3D(I,K-1,J))/DZ)**2.0D0
               ELSE
                DELRHOZ2(I,K,J) = ((TH_TOT_3D(I,K+1,J)-TH_TOT_3D(I,K-1,J))/(2.0D0*DZ))**2.0D0
               ENDIF

               IF(J .EQ. 1) THEN
                DELRHOY2(I,K,J) = ((TH_TOT_3D(I,K,J+1)-TH_TOT_3D(I,K,J))/(GYF(J+1)-GYF(J)))**2.0D0
               ELSEIF(J .EQ. NY_TEMP) THEN
                DELRHOY2(I,K,J) = ((TH_TOT_3D(I,K,J)-TH_TOT_3D(I,K,J-1))/(GYF(J)-GYF(J-1)))**2.0D0
               ELSE
                DELRHOY2(I,K,J) = ((TH_TOT_3D(I,K,J+1)-TH_TOT_3D(I,K,J-1))/(GYF(J+1)-GYF(J-1)))**2.0D0
               ENDIF
          
               DELRHO2(I,K,J) = DELRHOX2(I,K,J) + DELRHOZ2(I,K,J) + DELRHOY2(I,K,J) 
              ENDDO
             ENDDO
            ENDDO



          DO k=1,NZ
           DO J=1,NY_TEMP
            DO i=1,NX

!!!!!!!!!!!!!!!!!!!!!!  CALCULATING THE TERM {DEL. (KAPPA_T DEL RHO)} !!!!!!!!!!!!!

               IF(I .EQ. 1) THEN
                FLUXTERM_X(I,K,J) = ((NU/PR)*(TH_TOT_3D(I+2,K,J)-2.0*TH_TOT_3D(I+1,K,J)+TH_TOT_3D(I,K,J))/(DX**2.0)) !+ &
!                                    (((TH_TOT_3D(I+1,K,J)-TH_TOT_3D(I,K,J))/(DX))*((KAPPA_T_3D_TOT(I+1,K,J)-KAPPA_T_3D_TOT(I,K,J))/(DX)))
               ELSE IF(I .EQ. 2) THEN
                FLUXTERM_X(I,K,J) = ((NU/PR)*(TH_TOT_3D(I+2,K,J)-2.0*TH_TOT_3D(I+1,K,J)+TH_TOT_3D(I,K,J))/(DX**2.0)) !+ &
!                                    (((TH_TOT_3D(I+1,K,J)-TH_TOT_3D(I-1,K,J))/(2.0D0*DX))*((KAPPA_T_3D_TOT(I+1,K,J)-KAPPA_T_3D_TOT(I-1,K,J))/(2.0D0*DX)))
               ELSEIF(I .EQ. NX) THEN
                FLUXTERM_X(I,K,J) = ((NU/PR)*(TH_TOT_3D(I,K,J)-2.0*TH_TOT_3D(I-1,K,J)+TH_TOT_3D(I-2,K,J))/(DX**2.0)) !+ &
!                                    (((TH_TOT_3D(I,K,J)-TH_TOT_3D(I-1,K,J))/(1.0D0*DX))*((KAPPA_T_3D_TOT(I,K,J)-KAPPA_T_3D_TOT(I-1,K,J))/(1.0D0*DX)))
               ELSEIF(I .EQ. NX-1) THEN
                FLUXTERM_X(I,K,J) = ((NU/PR)*(TH_TOT_3D(I,K,J)-2.0*TH_TOT_3D(I-1,K,J)+TH_TOT_3D(I-2,K,J))/(DX**2.0)) !+ &
!                                    (((TH_TOT_3D(I+1,K,J)-TH_TOT_3D(I-1,K,J))/(2.0D0*DX))*((KAPPA_T_3D_TOT(I+1,K,J)-KAPPA_T_3D_TOT(I-1,K,J))/(2.0D0*DX)))
               ELSE
                FLUXTERM_X(I,K,J) = ((NU/PR)*(TH_TOT_3D(I+1,K,J)-2.0*TH_TOT_3D(I,K,J)+TH_TOT_3D(I-1,K,J))/(1.0D0*(DX**2.0))) !+ &
!                                    (((TH_TOT_3D(I+1,K,J)-TH_TOT_3D(I-1,K,J))/(2.0D0*DX))*((KAPPA_T_3D_TOT(I+1,K,J)-KAPPA_T_3D_TOT(I-1,K,J))/(2.0D0*DX)))
               ENDIF
      
               IF(K .EQ. 1) THEN
                FLUXTERM_Z(I,K,J) = ((NU/PR )*(TH_TOT_3D(I,K+2,J)-2.0*TH_TOT_3D(I,K+1,J)+TH_TOT_3D(I,K,J))/(DZ**2.0)) !+ &
!                                    (((TH_TOT_3D(I,K+1,J)-TH_TOT_3D(I,K,J))/(DZ))*((KAPPA_T_3D_TOT(I,K+1,J)-KAPPA_T_3D_TOT(I,K,J))/(DZ)))
               ELSE IF(K .EQ. 2) THEN
                FLUXTERM_Z(I,K,J) = ((NU/PR )*(TH_TOT_3D(I,K+2,J)-2.0*TH_TOT_3D(I,K+1,J)+TH_TOT_3D(I,K,J))/(DZ**2.0)) !+ &
!                                    (((TH_TOT_3D(I,K+1,J)-TH_TOT_3D(I,K-1,J))/(2.0D0*DZ))*((KAPPA_T_3D_TOT(I,K+1,J)-KAPPA_T_3D_TOT(I,K-1,J))/(2.0D0*DZ)))
               ELSEIF(K .EQ. NZ) THEN
                FLUXTERM_Z(I,K,J) = ((NU/PR )*(TH_TOT_3D(I,K,J)-2.0*TH_TOT_3D(I,K-1,J)+TH_TOT_3D(I,K-2,J))/(DZ**2.0)) !+ &
!                                    (((TH_TOT_3D(I,K,J)-TH_TOT_3D(I,K-1,J))/(1.0D0*DZ))*((KAPPA_T_3D_TOT(I,K,J)-KAPPA_T_3D_TOT(I,K-1,J))/(1.0D0*DZ)))
               ELSEIF(K .EQ. NZ-1) THEN
                FLUXTERM_Z(I,K,J) = ((NU/PR )*(TH_TOT_3D(I,K,J)-2.0*TH_TOT_3D(I,K-1,J)+TH_TOT_3D(I,K-2,J))/(DZ**2.0)) !+ &
!                                    (((TH_TOT_3D(I,K+1,J)-TH_TOT_3D(I,K-1,J))/(2.0D0*DZ))*((KAPPA_T_3D_TOT(I,K+1,J)-KAPPA_T_3D_TOT(I,K-1,J))/(2.0D0*DZ)))
               ELSE
                FLUXTERM_Z(I,K,J) = ((NU/PR )*(TH_TOT_3D(I,K+1,J)-2.0*TH_TOT_3D(I,K,J)+TH_TOT_3D(I,K-1,J))/(1.0D0*(DZ**2.0))) !+ &
!                                    (((TH_TOT_3D(I,K+1,J)-TH_TOT_3D(I,K-1,J))/(2.0D0*DZ))*((KAPPA_T_3D_TOT(I,K+1,J)-KAPPA_T_3D_TOT(I,K-1,J))/(2.0D0*DZ)))
               ENDIF


             IF(J .EQ. 1) THEN
                FLUXTERM_Y(I,K,J) = (  (NU/PR )*  &
                                    (((TH_TOT_3D(I,K,J+2)-TH_TOT_3D(I,K,J+1))/(DYF(J+1)*DY(J+2))) - ((TH_TOT_3D(I,K,J+1)-TH_TOT_3D(I,K,J))/(DYF(J+1)*DY(J+1))))  ) !+ &
!                                    (((TH_TOT_3D(I,K,J+1)-TH_TOT_3D(I,K,J))/(DY(J+1)))*((KAPPA_T_3D_TOT(I,K,J+1)-KAPPA_T_3D_TOT(I,K,J))/(DY(J+1))))
             ELSE IF(J .EQ. 2) THEN
                FLUXTERM_Y(I,K,J) = (  (NU/PR)*  &
                                    (((TH_TOT_3D(I,K,J+2)-TH_TOT_3D(I,K,J+1))/(DYF(J+1)*DY(J+2))) - ((TH_TOT_3D(I,K,J+1)-TH_TOT_3D(I,K,J))/(DYF(J+1)*DY(J+1))))  ) !+ &
!                                    (((TH_TOT_3D(I,K,J+1)-TH_TOT_3D(I,K,J-1))/(GYF(J+1)-GYF(J-1)))*((KAPPA_T_3D_TOT(I,K,J+1)-KAPPA_T_3D_TOT(I,K,J-1))/(GYF(J+1)-GYF(J-1))))
             ELSE IF(J .EQ. NY_TEMP-1) THEN
                FLUXTERM_Y(I,K,J) = (  (NU/PR)*  &
                                    (((TH_TOT_3D(I,K,J)-TH_TOT_3D(I,K,J-1))/(DY(J)*DYF(J-1))) - ((TH_TOT_3D(I,K,J-1)-TH_TOT_3D(I,K,J-2))/(DY(J-1)*DYF(J-1))))  ) !+ &
!                                    (((TH_TOT_3D(I,K,J+1)-TH_TOT_3D(I,K,J-1))/(GYF(J+1)-GYF(J-1)))*((KAPPA_T_3D_TOT(I,K,J+1)-KAPPA_T_3D_TOT(I,K,J-1))/(GYF(J+1)-GYF(J-1))))
             ELSE IF(J .EQ. NY_TEMP) THEN
              FLUXTERM_Y(I,K,J) = (  (NU/PR)*    &
                                    (((TH_TOT_3D(I,K,J)-TH_TOT_3D(I,K,J-1))/(DY(J)*DYF(J-1))) - ((TH_TOT_3D(I,K,J-1)-TH_TOT_3D(I,K,J-2))/(DY(J-1)*DYF(J-1))))  ) !+ &
!                                    (((TH_TOT_3D(I,K,J)-TH_TOT_3D(I,K,J-1))/(DY(J)))*((KAPPA_T_3D_TOT(I,K,J)-KAPPA_T_3D_TOT(I,K,J-1))/(DY(J))))
            ELSE
              FLUXTERM_Y(I,K,J) = (  (NU/PR )*  &
                                    (((TH_TOT_3D(I,K,J+1)-TH_TOT_3D(I,K,J))/(DY(J+1)*DYF(J))) - ((TH_TOT_3D(I,K,J)-TH_TOT_3D(I,K,J-1))/(DY(J)*DYF(J))))  )
            ENDIF


              FLUXTERM(I,K,J) = FLUXTERM_X(I,K,J) + FLUXTERM_Z(I,K,J) + FLUXTERM_Y(I,K,J)

           END DO
          END DO
         END DO 
      
         TPE = 0  
         m = 0

          DO k=1,NZ
           DO J=1,NY_TEMP
            DO i=1,NX
!             TPE = TPE + TH_TOT_3D(I,K,J)*gravity_g*ypoint_m(i,j)*DYF(J)*DX*DZ !DV(I,K,J)
              TPE = TPE + TH_TOT_3D(I,K,J)*gravity_g*(ypoint_m(i,j)+((DX/2)*sin(angle_beta)+(GYF(1)-GY(1))*cos(angle_beta)))*DYF(J)*DX*DZ !DV(I,K,J)

            ENDDO
           ENDDO
          ENDDO


           DO j=1,NY_TEMP
             DO i=1,NX

             IF(MOD(I,2) .eq. 0) THEN
              DO k=NZ,1,-1
               m = m+1
               DV_1D(m) = DYF(J)*DX*DZ !DV(I,K,J)
               TH_1D(m) = TH_TOT_3D(I,K,J)
               DELRHO2_1D(m)=DELRHO2(I,K,J)
               FLUXTERM_1D(m) = FLUXTERM(I,K,J)
             END DO
 
             ELSE
              DO k=1,NZ
               m = m+1
               DV_1D(m) = DYF(J)*DX*DZ !DV(I,K,J)
               TH_1D(m) = TH_TOT_3D(I,K,J)
               DELRHO2_1D(m)=DELRHO2(I,K,J)
               FLUXTERM_1D(m) = FLUXTERM(I,K,J)
             END DO
            ENDIF

           END DO
          END DO

       DO i = 1,NY_3D_TOT_TEMP
          TH_1D_INIT(i) = TH_1D(i)
       ENDDO

       call  CALCULATE_Z_AND_ZSTAR_INITIAL


!          write(6,*) 'TOTAL VOLUME COMPARISON',(DX*NX*DZ*NZ*SUM(DYF(1:NY_TEMP))), SUM(DV_1D(1:NY_3D_TOT_TEMP)),m, ypoint_m(1,1),GYF(1)-GY(1)

!!!!!!!!!!!!!!  IMPORTANT TO NOTE THAT WE NEED TO KEEP TRACK OF THE DV FOR EACH ELEMENT !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!  WHEN SORTING THE DENSITY FIELD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         write(*,*) 'Starting merge sort '
!!!!!!!!!!!!!!!!!!  MERGE SORTING   !!!!!!!!!!!!!!!!!!!!!!!!

    call MergeSort(TH_1D,NY_3D_TOT_TEMP,TEMP1,DV_1D,TEMP2,DELRHO2_1D,TEMP3,FLUXTERM_1D,TEMP4)
!   write(*,'(A,/,10I3)')'Sorted array :',A

!         write(*,*) 'Sorting done' 
!!!!!!!!!!!!!!!!!!! END MERGE SORTING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO i = 1,NY_3D_TOT_TEMP/2
       tmp1 = TH_1D(I)
       TH_1D(I)=TH_1D(NY_3D_TOT_TEMP-(I-1))
       TH_1D(NY_3D_TOT_TEMP-(I-1)) = tmp1

       tmp2 = DV_1D(I)
       DV_1D(I)=DV_1D(NY_3D_TOT_TEMP-(I-1))
       DV_1D(NY_3D_TOT_TEMP-(I-1)) = tmp2

       tmp3 = DELRHO2_1D(I)
       DELRHO2_1D(I)=DELRHO2_1D(NY_3D_TOT_TEMP-(I-1))
       DELRHO2_1D(NY_3D_TOT_TEMP-(I-1)) = tmp3

       tmp4 = FLUXTERM_1D(I)
       FLUXTERM_1D(I)=FLUXTERM_1D(NY_3D_TOT_TEMP-(I-1))
       FLUXTERM_1D(NY_3D_TOT_TEMP-(I-1)) = tmp4
      ENDDO 


       BPE_TOT=0.D0
       DBPE_DT = 0.d0
       PHI_D = 0.d0
       PHI_D_NEW = 0.D0

   call  CALCULATE_Z_AND_ZSTAR

!    write(6,*) 'Z_1D MIN MAX', MINVAL(Z_STAR_1D_INITIAL), MAXVAL(Z_STAR_1D_INITIAL),SUM(DV_1D(1:NY_3D_TOT_TEMP))


      i=2

       DO WHILE ( I .LE. (NY_3D_TOT_TEMP)) 
        matchcount=1

!        IF(I .GT. 1) THEN
!          DZSTAR_DRHO(I) = (Z_STAR_1D(I)-Z_STAR_1D(I-1))/(TH_1D(I)-TH_1D(I-1))
!        ELSE
!          DZSTAR_DRHO(I) = (Z_STAR_1D(I+1)-Z_STAR_1D(I))/(TH_1D(I+1)-TH_1D(I))
!        ENDIF

        ZSTAR_TEMP_2 = Z_STAR_1D(I-1)
        ZSTAR_TEMP = Z_STAR_1D(I)*DV_1D(I)
        TH_TEMP = TH_1D(I-1)
        DV_TEMP = DV_1D(I)
        DELRHO2_TEMP = DELRHO2_1D(I)*DV_1D(I)
        FLUXTERM_TEMP = FLUXTERM_1D(I)*DV_1D(I)

!!!!!!!!!!! CHECKING IF BOTH VOLUME ELEMENTS HAVE SAME DENSITY !!!!!!!!!

       IF((I .NE. NY_3D_TOT_TEMP) .AND. (TH_1D(I) .EQ. TH_1D(I+1))) THEN

         DO WHILE (TH_1D(I) .EQ. TH_1D(I+1))
           ZSTAR_TEMP =  ZSTAR_TEMP + Z_STAR_1D(I+1)*DV_1D(I+1)
           DV_TEMP = DV_TEMP + DV_1D(I+1)
           DELRHO2_TEMP = DELRHO2_TEMP + DELRHO2_1D(I+1)*DV_1D(I+1)
           FLUXTERM_TEMP = FLUXTERM_TEMP + FLUXTERM_1D(I+1)*DV_1D(I+1)
           I = I+1
            matchcount = matchcount + 1

         END DO        

          BPE_TEMP = TH_1D(I)*gravity_g*(ZSTAR_TEMP/DV_TEMP)*DV_TEMP
!          PHI_D_TEMP =  -(NU/PR)*gravity_g*((ZSTAR_TEMP/DV_TEMP-ZSTAR_TEMP_2)/(TH_1D(I)-TH_TEMP))* &
!          (DELRHO2_TEMP/DV_TEMP)*DV_TEMP

           PHI_D_TEMP =  -(NU/PR)*gravity_g*((Z_STAR_1D(I+1)-ZSTAR_TEMP_2)/(TH_1D(I+1)-TH_TEMP))* &
          (DELRHO2_TEMP/DV_TEMP)*DV_TEMP

!          - (1.d0*(KAPPA_T_TEMP/matchcount)*gravity_g*((ZSTAR_TEMP/matchcount-ZSTAR_TEMP_2)/(TH_1D(I)-TH_TEMP))* &
!         DV_TEMP*(DELRHO2_TEMP/matchcount))
           PHI_D_TEMP_NEW = gravity_g*(ZSTAR_TEMP/DV_TEMP)*(FLUXTERM_TEMP/DV_TEMP)*DV_TEMP

          I = I+1
       ELSE
!         TH_MATCHED=.FALSE.
         IF(I .EQ. NY_3D_TOT_TEMP) THEN
         write(*,*) 'MAX I REACHED',NY_3D_TOT_TEMP
         ENDIF

          BPE_TEMP = TH_1D(I)*gravity_g*(ZSTAR_TEMP/DV_TEMP)*DV_TEMP
          IF( I .NE. NY_3D_TOT_TEMP) THEN
          PHI_D_TEMP =  -(NU/PR)*gravity_g*((Z_STAR_1D(I+1)-Z_PREVIOUS_TEMP)/(TH_1D(I+1)-TH_1D(I-1)))*DV_TEMP*(DELRHO2_TEMP/DV_TEMP)
          ELSE
          PHI_D_TEMP =  -(NU/PR)*gravity_g*((Z_STAR_1D(I)-Z_STAR_1D(I-1))/(TH_1D(I)-TH_1D(I-1)))*DV_TEMP*(DELRHO2_TEMP/DV_TEMP)
          ENDIF
!-  (1.d0*KAPPA_T_TEMP*gravity_g*((Z_STAR_1D(I+1)-Z_STAR_1D(I-1))/(TH_1D(I+1)-TH_1D(I-1)))*DV_TEMP*DELRHO2_TEMP)
          PHI_D_TEMP_NEW = gravity_g*(ZSTAR_TEMP/DV_TEMP)*(FLUXTERM_TEMP/DV_TEMP)*DV_TEMP

          I = I+1
       ENDIF 

       Z_PREVIOUS_TEMP = ZSTAR_TEMP/DV_TEMP
       BPE_TOT = BPE_TOT + BPE_TEMP
       PHI_D = PHI_D + PHI_D_TEMP
       PHI_D_NEW = PHI_D_NEW + PHI_D_TEMP_NEW

       IF( (BPE_TEMP-5) .EQ. (BPE_TEMP+5)) THEN
        write(*,*) 'Nan detected', i, Z_1D(I)
        stop
       END IF  

      END DO
 
      DIFF_FLUX = 0
 
      DO I =1,NX  
       DO K=1,NZ
        DIFF_FLUX = DIFF_FLUX + (NU/PR)*gravity_g*ypoint_m(I,NY_TEMP)*(TH_TOT_3D(I,K,NY_TEMP)-TH_TOT_3D(I,K,NY_TEMP-1))/(ypoint(I,NY_TEMP)-ypoint(I,NY_TEMP-1))*DX*DZ
      ENDDO
     ENDDO

      IF(COUNTER .GT. 1) THEN
       DBPE_DT = (BPE_TOT-BPE_TOT_OLD)/(TIME-TIME_OLD)
      END IF 


    DO JJ=1,2      
 
      L_T(JJ) = 0.D0 
      N_ZSTAR_AVG(JJ) = 0.D0
      DISP_SUM = 0.D0
      COUNT1=0
      COUNT2=0
      INDEX_1=2
      INDEX_2=1000
      TP = 2*pi/omega0
      TIME_RATIO = TIME/TP

     IF(JJ .EQ. 1) THEN     !!!!!!!!!!!!!!! FIRST INVERSION

      DO
        IF(Z_STAR_1D(INDEX_1) .LT. 2.0) THEN
          INDEX_1=INDEX_1+1
        ELSE
          exit
         ENDIF
      END DO


      DO
        
        IF(Z_STAR_1D(INDEX_2) .LT. 8.0) THEN
          INDEX_2=INDEX_2+1
        ELSE
          exit
        ENDIF
 
      ENDDO
    
      write(6,*) 'JJ indexes', JJ,INDEX_1, INDEX_2
 
     ELSE       !!!!! SECOND INVERSION
 
       INDEX_1=2

     DO        
        IF(Z_STAR_1D(INDEX_2) .LT. 2.0) THEN
          INDEX_2=INDEX_2+1
        ELSE
          exit
        ENDIF
      
      ENDDO

       write(6,*) 'JJ indexes', JJ,INDEX_1, INDEX_2

     ENDIF



      DO I=INDEX_1,INDEX_2
!  
         DISP(I) = Z_STAR_1D(I) - Z_STAR_1D_INITIAL(I)
         DRHO_DZSTAR(1)=DRHO_DZSTAR(2)
         N_ZSTAR(I) = sqrt((-10.0/1000.0)*DRHO_DZSTAR(i))
        
          IF( (N_ZSTAR(i)-N_ZSTAR(i)).NE. 0) THEN
!           write(6,*) 'error index', i,DRHO_DZSTAR(i)
           stop
          ENDIF

!         IF(DISP(I) .LT. 0) THEN
!          COUNT1=COUNT1+1
!         ELSEIF(DISP(I) .GT. 0) THEN
!          COUNT2=COUNT2+1
!         ENDIF

         L_T(JJ) = L_T(JJ) + (DISP(I)**2.0)*DV_1D(I)
         N_ZSTAR_AVG(JJ) = N_ZSTAR_AVG(JJ) + N_ZSTAR(I)*DV_1D(I)
         DISP_SUM = DISP_SUM + DISP(I)
      ENDDO 

      N_ZSTAR_AVG(JJ) = N_ZSTAR_AVG(JJ)/(SUM(DV_1D(INDEX_1:INDEX_2)))
      L_T(JJ) = sqrt(L_T(JJ)/(SUM(DV_1D(INDEX_1:INDEX_2))) )
      L_O = 0.8*L_T(JJ)
      EPS_INFERRED(JJ) = (L_O**2)*(N_ZSTAR_AVG(JJ)**3.0)
   ENDDO     !!!!!!!!!END OF JJ (INVERSION) LOOP 


!      open(1012,file= file_th_1d,status='unknown',form='unformatted')


!         write(1012) TIME
!         write(1012)(TH_1D_INIT(i),i=1,NY_3D_TOT_TEMP),  &
!                   (TH_1D(i),i=1,NY_3D_TOT_TEMP),        &
!                   (Z_STAR_1D_INITIAL(i),i=1,NY_3D_TOT_TEMP),  &
!                   (Z_STAR_1D(i),i=1,NY_3D_TOT_TEMP)           
                   

!         close(1012)


!      write(*,*) 'Writing statistics to a file', DBPE_DT,TPE-BPE_TOT,PHI_D
      write(*,*) 'Thorpe scale', TIME_RATIO,L_T(1),L_T(2), N_ZSTAR_AVG(1),N_ZSTAR_AVG(2), EPS_INFERRED(1), EPS_INFERRED(2), PHI_D, PHI_D_NEW,DIFF_FLUX
 
      open(40,file='BPE.dat',form='formatted',status='unknown',position='append') 
       
      write(40,500) TIME, BPE_TOT, DBPE_DT, TPE, TPE-BPE_TOT,  PHI_D, PHI_D_NEW, L_T(1),L_T(2), N_ZSTAR_AVG(1),N_ZSTAR_AVG(2),EPS_INFERRED(1),EPS_INFERRED(2)
500   format(F20.8,12F35.16)
      

       TIME_OLD = TIME
       BPE_TOT_OLD = BPE_TOT
 
     
      END DO
!      ENDING TIME LOOP

       stop
      end




subroutine Merge(A,NA,B,NB,C,NC,D,E,F,X,Y,Z,P,Q,R)

   integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
   real*8, intent(in out) :: A(NA),D(NA),X(NA),P(NA)       ! B overlays C(NA+1:NC)
   real*8, intent(in)     :: B(NB),E(NB),Y(NB),Q(NB)
   real*8, intent(in out) :: C(NC),F(NC),Z(NC),R(NC)

   integer :: I,J,K

   I = 1; J = 1; K = 1;
   do while(I <= NA .and. J <= NB)
      if (A(I) <= B(J)) then
         C(K) = A(I)
         F(K) = D(I)
         Z(K) = X(I)
         R(K) = P(I)
         I = I+1
      else
         C(K) = B(J)
         F(K) = E(J)
         Z(K) = Y(J)
         R(K) = Q(J)
         J = J+1
      endif
      K = K + 1
   enddo
   do while (I <= NA)
      C(K) = A(I)
      F(K) = D(I)
      Z(K) = X(I)
      R(K) = P(I)
      I = I + 1
      K = K + 1
   enddo
   return
end subroutine merge

recursive subroutine MergeSort(A,N,TEMP1,AT,TEMP2,ATT,TEMP3,ATTT,TEMP4)

   integer, intent(in) :: N
   real*8, dimension(N), intent(in out) :: A
   real*8, dimension((N+1)/2), intent (out) :: TEMP1, TEMP2, TEMP3,TEMP4
   REAL*8, dimension(N), intent(in out) :: AT, ATT,ATTT

   integer :: NA,NB
   real*8  :: V

   if (N < 2) return
   if (N == 2) then
      if (A(1) > A(2)) then
         V = A(1)
         A(1) = A(2)
         A(2) = V

         V = AT(1)
         AT(1)=AT(2)
         AT(2) = V

         V = ATT(1)
         ATT(1)=ATT(2)
         ATT(2) = V

         V = ATTT(1)
         ATTT(1)=ATTT(2)
         ATTT(2) = V

!         V = ATTT(1)
!         ATTT(1)=ATTT(2)
!         ATTT(2) = V
      endif
      return
   endif
   NA=(N+1)/2
   NB=N-NA

   call MergeSort(A,NA,TEMP1,AT,TEMP2,ATT,TEMP3,ATTT,TEMP4)
   call MergeSort(A(NA+1),NB,TEMP1,AT(NA+1),TEMP2,ATT(NA+1),TEMP3,ATTT(NA+1),TEMP4)

   if (A(NA) > A(NA+1)) then
      TEMP1(1:NA)=A(1:NA)
      TEMP2(1:NA)=AT(1:NA)
      TEMP3(1:NA)=ATT(1:NA)  
      TEMP4(1:NA)=ATTT(1:NA)
      call Merge(TEMP1,NA,A(NA+1),NB,A,N,TEMP2,AT(NA+1),AT,TEMP3,ATT(NA+1),ATT,TEMP4,ATTT(NA+1),ATTT)
   endif
   return

end subroutine MergeSort




subroutine CALCULATE_Z_AND_ZSTAR

       use variables
!       implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  CALCULATING THE BPE USING THE SORTED FIELD !!!!!!!!!!!!!!!!

       METHOD1 = .TRUE.
       METHOD2 = .FALSE. 
       METHOD3 = .FALSE. 

       DV_1D_TOT =0

       DO i = 1,NY_3D_TOT_TEMP
         DV_1D_TOT = DV_1D_TOT + DV_1D(I)

         IF(METHOD1) THEN
          Z_1D(I)= ABS(sqrt(DV_1D_TOT*sin(2*angle_beta)/LZ))
            
           IF((Z_1D(I)/(sin(angle_beta))) .LE. LX) THEN
 !           DO NOTHING
           ELSE
             I1 = I
!             write(*,*)'exceeded the dimension in x- direction'
             METHOD1 = .FALSE.
             METHOD2 = .TRUE.
             GOTO 232
           END IF     
         
         ELSEIF(METHOD2) THEN
232         Z_1D(I)= cos(angle_beta)*(DV_1D_TOT + (0.5d0*(LX**2.0d0)*LZ*tan(angle_beta)))/(LX*LZ)

           IF((Z_1D(I)/(cos(angle_beta))) .LE. (GY(NY_TEMP+1)-GY(1))) THEN
 !           DO NOTHING
           ELSE
             I2 = I
!             write(*,*)'exceeded the dimension in z- direction'
             METHOD2 = .FALSE.
             METHOD3 = .TRUE.
             GOTO 246
           END IF
 
           ELSEIF(METHOD3) THEN

           IF(I .NE. NY_3D_TOT_TEMP) THEN
246         Z_1D(I)= ABS(SQRT((LX*(GY(NY_TEMP+1)-GY(1))*sin(2*angle_beta))-DV_1D_TOT*sin(2*angle_beta)/LZ)) 
            Z_1D(I) = LX*sin(angle_beta) + (GY(NY_TEMP+1)-GY(1))*cos(angle_beta) - Z_1D(I)

           ELSE
            Z_1D(I) = LX*sin(angle_beta) + (GY(NY_TEMP+1)-GY(1))*cos(angle_beta)
           END IF

           END IF
 
       
           IF(Z_1D(I) .LE. 0) THEN
!            write(*,*)'Error detected',i
            stop
           ENDIF
           
       END DO


!!!!!!!!!!!!!!! CALCULATING Z STAR !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!  Z STAR IS THE CENTER OF COLUME/MASS FOR EACH VOLUME ELEMENT, THIS HIEGHT IS USED TO CALCULATE THE BACKGROUND POTENTIAL ENERGY

       i=1
       DO WHILE ( I .LE. NY_3D_TOT_TEMP) 

        IF(I .EQ. 1) THEN
         Z_STAR_1D(I)= (2.0d0/3.0d0)*(Z_1D(I))
!         Z_STAR_1D(I)= (1.0d0/2.0d0)*(Z_1D(I))
        ELSEIF(I .EQ. NY_3D_TOT_TEMP) THEN
         Z_STAR_1D(I)= (Z_1D(I)+2.0D0*Z_1D(I-1))/3
        ELSE
!         Z_STAR_1D(I)= 0.5d0*(Z_1D(I-1) + Z_1D(I))
         !!!!!!!!!!  Centroid of a trapezoid with a shorter side and b the longer side !!!!!!!!!!!!!!!!!!! 
         IF( I .LT. I1) THEN

          a = Z_1D(I-1)*tan(angle_beta)+Z_1D(I-1)*(1.0d0/tan(angle_beta))
          b = Z_1D(I)*tan(angle_beta)+Z_1D(I)*(1.0d0/tan(angle_beta))
          Z_STAR_1D(I)=Z_1D(I)-(1.0D0/3.0D0)*(Z_1D(I)-Z_1D(I-1))*(2*a+b)/(a+b)
 
         ELSE IF( I .LT. I2) THEN
            Z_STAR_1D(I)= 0.5d0*(Z_1D(I-1) + Z_1D(I))
        ELSE
          Z0 =  LX*sin(angle_beta) + (GY(NY_TEMP+1)-GY(1))*cos(angle_beta)
          b = (Z0-Z_1D(I-1))*tan(angle_beta)+(Z0-Z_1D(I-1))*(1.0d0/tan(angle_beta))
          a = (Z0-Z_1D(I))*tan(angle_beta)+(Z0-Z_1D(I))*(1.0d0/tan(angle_beta))
          Z_STAR_1D(I)=Z_1D(I-1)+(1.0D0/3.0D0)*(Z_1D(I)-Z_1D(I-1))*(2*a+b)/(a+b)
         ENDIF  

        END IF

!!!!!!!!!!! CHECKING IF THERE IS ANY ANOMALY IN DRHO_DZSTAR !!!!!!!!!!
       IF(I .GT. 1) THEN
        DRHO_DZSTAR(i) = (TH_1D(I)-TH_1D(I-1))/(Z_STAR_1D(I)-Z_STAR_1D(I-1))
!         DRHO_DZSTAR = (TH_1D(I)-TH_1D(I-1))
        IF((DRHO_DZSTAR(i) .GT. 0)) THEN
          write(*,*) 'Error detected drhodz* is positive' , i
          stop
        END IF
      ENDIF

        

        I = I+1
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine CALCULATE_Z_AND_ZSTAR


subroutine CALCULATE_Z_AND_ZSTAR_INITIAL

       use variables
!       implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  CALCULATING THE BPE USING THE SORTED FIELD !!!!!!!!!!!!!!!!

       METHOD1 = .TRUE.
       METHOD2 = .FALSE. 
       METHOD3 = .FALSE. 

       DV_1D_TOT =0

       DO i = 1,NY_3D_TOT_TEMP
         DV_1D_TOT = DV_1D_TOT + DV_1D(I)

         IF(METHOD1) THEN
          Z_1D_INITIAL(I)= ABS(sqrt(DV_1D_TOT*sin(2*angle_beta)/LZ))
            
           IF((Z_1D_INITIAL(I)/(sin(angle_beta))) .LE. LX) THEN
 !           DO NOTHING
           ELSE
             I1 = I
!             write(*,*)'exceeded the dimension in x- direction'
             METHOD1 = .FALSE.
             METHOD2 = .TRUE.
             GOTO 232
           END IF     
         
         ELSEIF(METHOD2) THEN
232         Z_1D_INITIAL(I)= cos(angle_beta)*(DV_1D_TOT + (0.5d0*(LX**2.0d0)*LZ*tan(angle_beta)))/(LX*LZ)

           IF((Z_1D_INITIAL(I)/(cos(angle_beta))) .LE. (GY(NY_TEMP+1)-GY(1))) THEN
 !           DO NOTHING
           ELSE
             I2 = I
!             write(*,*)'exceeded the dimension in z- direction'
             METHOD2 = .FALSE.
             METHOD3 = .TRUE.
             GOTO 246
           END IF
 
           ELSEIF(METHOD3) THEN

           IF(I .NE. NY_3D_TOT_TEMP) THEN
246         Z_1D_INITIAL(I)= ABS(SQRT((LX*(GY(NY_TEMP+1)-GY(1))*sin(2*angle_beta))-DV_1D_TOT*sin(2*angle_beta)/LZ)) 
            Z_1D_INITIAL(I) = LX*sin(angle_beta) + (GY(NY_TEMP+1)-GY(1))*cos(angle_beta) - Z_1D_INITIAL(I)

           ELSE
            Z_1D_INITIAL(I) = LX*sin(angle_beta) + (GY(NY_TEMP+1)-GY(1))*cos(angle_beta)
           END IF

           END IF
 
       
           IF(Z_1D_INITIAL(I) .LE. 0) THEN
            write(*,*)'Error detected',i
            stop
           ENDIF
           
       END DO


!!!!!!!!!!!!!!! CALCULATING Z STAR !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!  Z STAR IS THE CENTER OF COLUME/MASS FOR EACH VOLUME ELEMENT, THIS HIEGHT IS USED TO CALCULATE THE BACKGROUND POTENTIAL ENERGY

       i=1
       DO WHILE ( I .LE. NY_3D_TOT_TEMP) 

        IF(I .EQ. 1) THEN
         Z_STAR_1D_INITIAL(I)= (2.0d0/3.0d0)*(Z_1D_INITIAL(I))
!         Z_STAR_1D(I)= (1.0d0/2.0d0)*(Z_1D(I))
        ELSEIF(I .EQ. NY_3D_TOT_TEMP) THEN
         Z_STAR_1D_INITIAL(I)= (Z_1D_INITIAL(I)+2.0D0*Z_1D_INITIAL(I-1))/3
        ELSE
!         Z_STAR_1D(I)= 0.5d0*(Z_1D(I-1) + Z_1D(I))
         !!!!!!!!!!  Centroid of a trapezoid with a shorter side and b the longer side !!!!!!!!!!!!!!!!!!! 
         IF( I .LT. I1) THEN

          a = Z_1D_INITIAL(I-1)*tan(angle_beta)+Z_1D_INITIAL(I-1)*(1.0d0/tan(angle_beta))
          b = Z_1D_INITIAL(I)*tan(angle_beta)+Z_1D_INITIAL(I)*(1.0d0/tan(angle_beta))
          Z_STAR_1D_INITIAL(I)=Z_1D_INITIAL(I)-(1.0D0/3.0D0)*(Z_1D_INITIAL(I)-Z_1D_INITIAL(I-1))*(2*a+b)/(a+b)
 
         ELSE IF( I .LT. I2) THEN
            Z_STAR_1D_INITIAL(I)= 0.5d0*(Z_1D_INITIAL(I-1) + Z_1D_INITIAL(I))
        ELSE
          Z0 =  LX*sin(angle_beta) + (GY(NY_TEMP+1)-GY(1))*cos(angle_beta)
          b = (Z0-Z_1D_INITIAL(I-1))*tan(angle_beta)+(Z0-Z_1D_INITIAL(I-1))*(1.0d0/tan(angle_beta))
          a = (Z0-Z_1D_INITIAL(I))*tan(angle_beta)+(Z0-Z_1D_INITIAL(I))*(1.0d0/tan(angle_beta))
          Z_STAR_1D_INITIAL(I)=Z_1D_INITIAL(I-1)+(1.0D0/3.0D0)*(Z_1D_INITIAL(I)-Z_1D_INITIAL(I-1))*(2*a+b)/(a+b)
         ENDIF  

        END IF

!!!!!!!!!!! CHECKING IF THERE IS ANY ANOMALY IN DRHO_DZSTAR !!!!!!!!!!
!        IF(I .GT. 1) THEN
!         DRHO_DZSTAR = (TH_1D(I)-TH_1D(I-1))/(Z_STAR_1D_INITIAL(I)-Z_STAR_1D_INITIAL(I-1))
! !         DRHO_DZSTAR = (TH_1D(I)-TH_1D(I-1))
!         IF((DRHO_DZSTAR .GT. 0)) THEN
!           write(*,*) 'Error detected drhodz* is positive' , i
!           stop
!         END IF
!        ENDIF

        I = I+1
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine CALCULATE_Z_AND_ZSTAR_INITIAL

!SUBROUTINE WRITE_TECPLOT_3D_DENSITY_DATA
!       use variables
!       implicit none
!  CHARACTER*32 filn_tec_xyz,flow_data

!!!!!!!!!!!!!!!!!!!! Tecplot variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       integer  debug,ier,itot,ymax
!       integer  tecini,tecdat,teczne,tecnod,tecfil,tecend
!       integer  visdouble,disdouble, INDEXES
!       character*1 nulchar


!!!!!!!!! WRITING 3D DENSITY DATA TO TECPLOT !!!!!!!!!!!!!!!!!!!!!!!!!!
!         nulchar = char(0)
!        debug   = 0
!        visdouble = 0
!        disdouble = 1
!        imax = NX
!        kmax = NZ
!        jmax = NY_TOT

!        filn_tec_xyz = 'density_data/data_th_' &
!              //CHAR(MOD(COUNTER,10000000)/1000000+48)  &
!              //CHAR(MOD(COUNTER,1000000)/100000+48)    &
!              //CHAR(MOD(COUNTER,100000)/10000+48)      &
!              //CHAR(MOD(COUNTER,10000)/1000+48)        &
!              //CHAR(MOD(COUNTER,1000)/100+48)          &
!              //CHAR(MOD(COUNTER,100)/10+48)            &
!              //CHAR(MOD(COUNTER,10)+48) //             &
!              '.plt'



!       write(6,*)COUNTER,filn_tec_xyz,TIME
!       Write( flow_data,111) 'Time=',TIME
!111    format(A,f10.5)

!       ier = tecini(flow_data//nulchar,&
!                  &'x,y,z,x_T,y_T,rho_3D'&
!                  &//nulchar,&
!                  &filn_tec_xyz//nulchar,&
!                  &'.'//nulchar,&
!                  &debug,visdouble)
!!       endif
!      write(6,*) 'strat writing in plt format'
      !
      ! Write the zone header information.
      !
!      ier = teczne('flow_data'//nulchar, &
!                  imax,kmax,jmax,        &
       !
      ! Write out the field data.
!       itot = imax*kmax*jmax
!      ier = tecdat(itot,xpoint3D,disdouble)
!      ier = tecdat(itot,ypoint3D,disdouble)
!      ier = tecdat(itot,zpoint3D,disdouble)
!      ier = tecdat(itot,xpoint3D_m,disdouble)
!      ier = tecdat(itot,ypoint3D_m,disdouble)
!      ier = tecdat(itot,TH_3D,disdouble)
!      ier = tecend()

!     return 
!    end subroutine WRITE_TECPLOT_3D_DENSITY_DATA
