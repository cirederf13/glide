MODULE definitions 
      
      type parm
         integer n,m_max,napt,nzirc,narar,dummy,contr
         integer iflag_error
         double precision edot_mean,deltat,xL,sigma2,t_total,lat1,lat2,lon1
         double precision bot,lon2,dlon,kappa,hp,Ts,Tb,zl,angle,aspect
         character :: run*5,topofile*100
      end type

      type matr
         integer,dimension(:),allocatable :: isys,B,nsystems
         double precision,dimension(:),allocatable :: edot,edot_pr,zc,zcp,eps,elev,eps_dat,eps_res,work,syn_age,elev_true
         double precision,dimension(:),allocatable :: x_true,y_true,x_dum_true,y_dum_true 
         double precision,dimension(:),allocatable :: BB,edot_dat,edot_pr_dum,Cmm,eps_dum,x_dum,y_dum,misfits
         double precision,dimension(:),allocatable :: R,kernel,sf
         double precision,dimension(:),allocatable :: ta,a_error,x,y,cov_eps,depth,depth1,zz,tsteps,tsteps_sum
         double precision,dimension(:,:),allocatable :: A,G,cov,cpost,Y1,Y2,Y3,H,II,distm,ages,Y9,Cee_full,Cee_pr
      end type

end MODULE definitions

