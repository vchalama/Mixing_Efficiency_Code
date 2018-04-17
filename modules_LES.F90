module variables

integer        i,j,m,NY,NY_min,NY_max,COUNTER,STEP,I1,I2,matchcount, COUNT1,COUNT2,JJ,   &
           TIME_STEP,T_S,T_S_bin,NX,n_pc,NY_TOT,k,NZ,imax,jmax,kmax, NY_3D_TOT,NY_3D_TOT_TEMP,NY_TEMP, INDEX_1, INDEX_2, tmp6, tmp_index


        REAL*8        TIME,drho_dz, angle_beta,omega0,N_infty, a,b,Z0, &
                      gravity_g, rho_0,dx,LX,LZ,dz,lw_rho,pi,beta, alpha,angle_alpha, &
                      Max_vel, Vel_max,K_X_R,LAMBDA_X,BPE_TOT,TIME_OLD,BPE_TOT_OLD,DBPE_DT,TPE, &
                      BPE_TEMP,ZSTAR_TEMP,DV_TEMP,DIFF_FLUX,NU,PR,TIME_TEMP,L_T(1:2),DISP_SUM(1:2),DISP_SUM_TEMP,MIN_DISP,MAX_DISP,DV_1D_TOT, TP, TIME_RATIO, TEMP_SUM


        LOGICAL       LES,METHOD1,METHOD2,METHOD3,TH_MATCHED
        parameter     (n_pc = 32,  NY =20, NX=128, NZ=64, &
                      NY_TOT = NY*n_pc-(n_pc-1),          &
                      TIME_STEP=500000,                   &
                      T_S_bin=30000, NY_3D_TOT=NX*NZ*NY_TOT, &
                      NY_TEMP_MIN=300, NY_TEMP_MAX=490, NY_TEMP=NY_TEMP_MAX-NY_TEMP_MIN+1, &
                      NY_3D_TOT_TEMP=NX*NZ*(NY_TEMP_MAX-NY_TEMP_MIN+1))    

        parameter     (beta = 5.0d0, alpha=5.0d0,omega0=0.00014075, &
                      gravity_g = 10.0d0, rho_0=1000.0d0  )


REAL*8               xpoint(1:NX,1:NY_TOT), &
                      ypoint(1:NX,1:NY_TOT),xpoint_m(1:NX,1:NY_TOT),    &
                      ypoint_m(1:NX,1:NY_TOT),UME(1:NY_TOT),            &
                      z_point(1:NZ,1:NY_TOT),y_point(1:NZ,1:NY_TOT),    &
                      TH_B(1:NX,1:NY_TOT),TEMP_1(1:NX,1:NY_TOT),TEMP_2(1:NY_TOT)


REAL*8    GYF(0:NY_TOT+1), GY(0:NY_TOT+1),  &
             DYF(0:NY_TOT+1), DY(0:NY_TOT+1), &
             DV(1:NX,1:NZ,NY_TEMP_MIN:NY_TEMP_MAX),    &
             DV_1D(1:NY_3D_TOT_TEMP),        &
             Z_STAR_1D(1:NY_3D_TOT_TEMP),    &
             Z_1D(1:NY_3D_TOT_TEMP),   &
             Z_STAR_1D_INITIAL(1:NY_3D_TOT_TEMP),    &
             Z_1D_INITIAL(1:NY_3D_TOT_TEMP),   &
             DISP(1:NY_3D_TOT_TEMP),DISP_TEMP(1:NY_3D_TOT_TEMP)

real*8        TH_3D(1:NX,1:NZ,1:NY_TOT) ,KAPPA_T_3D_TOT(1:NX,1:NZ,1:NY_TOT),&
                      TH_1D(1:NY_3D_TOT_TEMP),TH_TOT_3D(1:NX,1:NZ,NY_TEMP_MIN:NY_TEMP_MAX), TH_1D_INIT(1:NY_3D_TOT_TEMP), DV_INIT(1:NY_3D_TOT_TEMP),&
                      tmp1,tmp2,DRHO_DZSTAR(1:NY_3D_TOT_TEMP),tmp3,DELRHO2_TEMP,tmp4,KAPPA_T_TEMP,tmp5,FLUXTERM_TEMP, &
                      N_ZSTAR(1:NY_3D_TOT_TEMP),EPS_INFERRED(1:2),N_ZSTAR_AVG(1:2),L_O,ZSTAR_TEMP_2,KAPPA_T_TEMP_2,Z_PREVIOUS_TEMP
end module variables
