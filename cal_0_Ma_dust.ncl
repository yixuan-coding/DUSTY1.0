;;；To calculate the 0 Ma dust emissions and scale the area-weighted global mean emission to observation
;;; The factorisation is carried on simultaneously
;;;do not include the specific exp names and directories
;;;
;;;all variables set to "double"


begin

     ;print("start")  

;;;--- switches
write_dust=1 ;;write nc file at 1
write_factor=1
do_scaling=1 ;; scale to obs at 1


expname=" " ;;

lsm_root_dir=" "
exp_dir=" "
n_month=12
monthname=new(n_month,string)
monthname=(/"ja","fb","mr","ar","my","jn","jl","ag","sp","ot","nv","dc"/)
lsname=" " ;; 

nlon=96
nlat=73
n_hour_total=12*30*720 ;;how many hours are calculated [12*30*720]

hour_l=new((/nlat,nlon/),double)
hour_lb=new((/nlat,nlon/),double)
hour_lbm=new((/nlat,nlon/),double)
hour_lbs=new((/nlat,nlon/),double)
hour_lbU=new((/nlat,nlon/),double)
hour_lbms=new((/nlat,nlon/),double)
hour_lbmU=new((/nlat,nlon/),double)
hour_lbsU=new((/nlat,nlon/),double)
hour_lbmsU=new((/nlat,nlon/),double)

scaled_hour_l=new((/nlat,nlon/),double)
scaled_hour_lb=new((/nlat,nlon/),double)
scaled_hour_lbm=new((/nlat,nlon/),double)
scaled_hour_lbs=new((/nlat,nlon/),double)
scaled_hour_lbU=new((/nlat,nlon/),double)
scaled_hour_lbms=new((/nlat,nlon/),double)
scaled_hour_lbmU=new((/nlat,nlon/),double)
scaled_hour_lbsU=new((/nlat,nlon/),double)
scaled_hour_lbmsU=new((/nlat,nlon/),double)

dhour_b=new((/nlat,nlon/),double)
dhour_m=new((/nlat,nlon/),double)
dhour_s=new((/nlat,nlon/),double)
dhour_U=new((/nlat,nlon/),double)
dhour_all_l=new((/nlat,nlon/),double)
dhour_all_lb=new((/nlat,nlon/),double)

dD_m=new((/nlat,nlon/),double)
dD_m!0="latitude"
dD_m!1="longitude"
latitude=fspan(90.0,-90.0,73)
longitude=fspan(0.0,356.25,96)
latitude@units="degrees_north"
longitude@units="degrees_east"
dD_m&longitude=longitude
dD_m&latitude=latitude

emi_rate_lbmsU=new((/nlat,nlon/),double)
emi_rate_lbmsU!0="latitude"
emi_rate_lbmsU!1="longitude"
latitude=fspan(90.0,-90.0,73)
longitude=fspan(0.0,356.25,96)
latitude@units="degrees_north"
longitude@units="degrees_east"
emi_rate_lbmsU&longitude=longitude
emi_rate_lbmsU&latitude=latitude

C1_total=new(1,double)
C1_avg=new(1,double)
C2_total=new(1,double)
C2_avg=new(1,double)
C3_total=new(1,double)
C3_avg=new(1,double)
C4_total=new(1,double)
C4_avg=new(1,double)
C5_total=new(1,double)
C5_avg=new(1,double)

;;;---dust model parameters
t_baresoil_array=(/0.263293, 0.276331, 0.273549, 0.251742, 0.262899/)
t_smoisture_array=(/22.657526, 21.608974, 23.838291, 23.155424, 22.88455/)
t_snow_array=(/31.750608, 47.453629, 42.489281, 42.004654, 49.7556/)
t_wind_min_aray=(/1.609415, 1.974174, 1.559718, 1.500514, 0.983392/)
t_wind_max_array=(/2.831963, 2.723186, 3.437541, 2.668725, 2.520874/)


model_number = 1

t_bsoil = t_baresoil_array(model_number-1)
t_smoisture = t_smoisture_array(model_number-1)
t_snow = t_snow_array(model_number-1)
t_wind_min = t_wind_min_aray(model_number-1)
t_wind_max = t_wind_max_array(model_number-1)

print("model_number="+model_number+", t_bsoil="+t_bsoil+", t_smoisture="+t_smoisture+" t_snow="+t_snow+" t_wind_min="+t_wind_min+"t_wind_max="+t_wind_max)

;;;---the main process

f_lsm=addfile(lsm_root_dir+lsname+"/inidata/"+lsname+".qrparm.mask.nc","r")


land=(/f_lsm->lsm/) 
var_land=flt2dble(land(0,0,:,:))
lat=f_lsm->latitude

yearhead="0000000"

;;;--scalling stuff---
D_obs = 1.96474*10^(-10) ;;kg m-2 s-1
rad = 4.0*atan(1.0)/180.0
clat = cos(lat*rad)

do y=0,29 ;;30 in total
print("y="+(y+1))
    do m=0,11 ;; 12 in total
    ;print("y="+(y+1)+",m="+(m+1))
   
        month=monthname(m)
        
        if (y .lt. 9) then yearname="0"+(y+1)
        else 
        yearname=(y+1)
        end if

        filename_pi=expname+"a#pi"+yearhead+yearname+month+"+.nc"
        file_pi=exp_dir+filename_pi
        f_pi=addfile(file_pi,"r")

        bsoil=(/f_pi->fracPFTs_mm_srf/)
        var_bsoil=flt2dble(bsoil(0,7,:,:))
     
        filename_pa=expname+"a#pa"+yearhead+yearname+month+"+.nc"    
        file_pa=exp_dir+filename_pa
        f_pa=addfile(file_pa,"r")
        

        smoisture=(/f_pa->sm_soil/)
        snow=(/f_pa->snowdepth_srf/)
 

        file_wind=exp_dir+"regrided_wind/regrided_wind_"+filename_pa
        f_wind=addfile(file_wind,"r")
        wind_u=(/f_wind->u_10m/)
        wind_v=(/f_wind->v_10m/)

do h=0,719 ;;720 in total
    ;print("h="+(h+1))
    var_smoisture=flt2dble(smoisture(h,0,:,:))
    var_snow=flt2dble(snow(h,0,:,:))
    var_wind_u=flt2dble(wind_u(h,0,:,:))
    var_wind_v=flt2dble(wind_v(h,0,:,:))   

    D_l = where(var_land .eq. 1, var_land, 0)  

    D_b = where(var_bsoil .lt. t_bsoil, 0, var_bsoil)
    D_b = where(ismissing(D_b),0,D_b)

    D_m = where(var_smoisture .gt. t_smoisture, 0, (t_smoisture-var_smoisture)/t_smoisture)
    D_m = where(ismissing(D_m),0,D_m)

    D_s = where(var_snow .gt. t_snow, 0, (t_snow-var_snow)/t_snow)
    D_s = where(ismissing(D_s),0,D_s)
    
    U = sqrt(var_wind_u*var_wind_u + var_wind_v*var_wind_v)

    F_1 = where(U .ge. t_wind_min, 1, 0)
    F_2 = where(U .gt. t_wind_max, 0, 1)
    F_3 = where(U .gt. t_wind_max, 1, 0)
    F_u = (U - t_wind_min) * U^2
    F_u2 = (t_wind_max - t_wind_min) * t_wind_max^2
    D_U = F_1*F_2*F_u + F_3*F_u2


        n_hour = y*12*720+m*720+h

        hour_l = D_l
        hour_lb = D_l * D_b
        hour_lbm = D_l * D_b * D_m
        hour_lbs = D_l * D_b * D_s
        hour_lbU = D_l * D_b * D_U
        hour_lbms = D_l * D_b * D_m * D_s    
        hour_lbsU = D_l * D_b * D_s * D_U
        hour_lbmU = D_l * D_b * D_m * D_U  
        hour_lbmsU = D_l * D_b * D_m * D_s * D_U
        
        ;;---scaling to observation
        if (do_scaling .eq. 1) then
        
        C1 = D_obs/wgt_areaave_Wrap(hour_lbms, clat, 1.0, 1)
       
        C2 = C1*(wgt_areaave_Wrap(hour_lbms, clat, 1.0, 1)/wgt_areaave_Wrap(hour_lbmsU, clat, 1.0, 1))
        C3 = C1*(wgt_areaave_Wrap(hour_lb, clat, 1.0, 1)/wgt_areaave_Wrap(hour_lbU, clat, 1.0, 1))
        C4 = C1*(wgt_areaave_Wrap(hour_lbm, clat, 1.0, 1)/wgt_areaave_Wrap(hour_lbmU, clat, 1.0, 1))
        C5 = C1*(wgt_areaave_Wrap(hour_lbs, clat, 1.0, 1)/wgt_areaave_Wrap(hour_lbsU, clat, 1.0, 1))
        else 
        C1 = 1 
        C2 = 1 
        C3 = 1 
        C4 = 1 
        C5 = 1 
    
        end if 

     
        scaled_hour_l = C1 * D_l
        scaled_hour_lb = C1 * D_l * D_b
        scaled_hour_lbm = C1 * D_l * D_b * D_m
        scaled_hour_lbs = C1 * D_l * D_b * D_s
        scaled_hour_lbms = C1 * D_l * D_b * D_m * D_s  
        scaled_hour_lbU = C3 * D_l * D_b * D_U          
        scaled_hour_lbsU = C5 * D_l * D_b * D_s * D_U
        scaled_hour_lbmU = C4 * D_l * D_b * D_m * D_U       
        scaled_hour_lbmsU = C2 * D_l * D_b * D_m * D_s * D_U

        ;;---linear-sum factorisation
        dhour_all_l = scaled_hour_lbmsU - scaled_hour_l
        dhour_all_lb = scaled_hour_lbmsU - scaled_hour_lb
        dhour_b = scaled_hour_lb - scaled_hour_l
        dhour_m = 1.0/6.0*(2*(scaled_hour_lbm-scaled_hour_lb)+(scaled_hour_lbmU-scaled_hour_lbU)+(scaled_hour_lbms-scaled_hour_lbs)+2*(scaled_hour_lbmsU-scaled_hour_lbsU))
        dhour_s = 1.0/6.0*(2*(scaled_hour_lbs-scaled_hour_lb)+(scaled_hour_lbms-scaled_hour_lbm)+(scaled_hour_lbsU-scaled_hour_lbU)+2*(scaled_hour_lbmsU-scaled_hour_lbmU))
        dhour_U = 1.0/6.0*(2*(scaled_hour_lbU-scaled_hour_lb)+(scaled_hour_lbmU-scaled_hour_lbm)+(scaled_hour_lbsU-scaled_hour_lbs)+2*(scaled_hour_lbmsU-scaled_hour_lbms))


        ;;---sum up in loops
        if (y .eq. 0 .and. m .eq. 0 .and. h .eq. 0) then 
            sum_scaled_hour_l = scaled_hour_l
            sum_scaled_hour_lb = scaled_hour_lb 
            sum_scaled_hour_lbm = scaled_hour_lbm 
            sum_scaled_hour_lbs = scaled_hour_lbs 
            sum_scaled_hour_lbU = scaled_hour_lbU 
            sum_scaled_hour_lbms = scaled_hour_lbms 
            sum_scaled_hour_lbsU = scaled_hour_lbsU 
            sum_scaled_hour_lbmU = scaled_hour_lbmU 
            sum_scaled_hour_lbmsU = scaled_hour_lbmsU

            sum_dhour_all_l = dhour_all_l 
            sum_dhour_all_lb = dhour_all_lb
            sum_dhour_b = dhour_b 
            sum_dhour_m = dhour_m 
            sum_dhour_s = dhour_s 
            sum_dhour_U = dhour_U
        else 
            sum_scaled_hour_l = sum_scaled_hour_l + scaled_hour_l
            sum_scaled_hour_lb = sum_scaled_hour_lb + scaled_hour_lb 
            sum_scaled_hour_lbm = sum_scaled_hour_lbm + scaled_hour_lbm 
            sum_scaled_hour_lbs = sum_scaled_hour_lbs + scaled_hour_lbs 
            sum_scaled_hour_lbU = sum_scaled_hour_lbU + scaled_hour_lbU 
            sum_scaled_hour_lbms = sum_scaled_hour_lbms + scaled_hour_lbms 
            sum_scaled_hour_lbsU = sum_scaled_hour_lbsU + scaled_hour_lbsU 
            sum_scaled_hour_lbmU = sum_scaled_hour_lbmU + scaled_hour_lbmU 
            sum_scaled_hour_lbmsU = sum_scaled_hour_lbmsU + scaled_hour_lbmsU
    
            sum_dhour_all_l = sum_dhour_all_l + dhour_all_l 
            sum_dhour_all_lb = sum_dhour_all_lb + dhour_all_lb
            sum_dhour_b = sum_dhour_b + dhour_b 
            sum_dhour_m = sum_dhour_m + dhour_m 
            sum_dhour_s = sum_dhour_s + dhour_s 
            sum_dhour_U = sum_dhour_U + dhour_U

        end if 

        ;;---average the constants---
        if (h .eq. 0 .and. m .eq. 0 .and. y .eq. 0)
        n = 0
        C1_total = C1
        C1_avg = C1
        C2_total = C2
        C2_avg = C2
        C3_total = C3
        C3_avg = C3
        C4_total = C4
        C4_avg = C4
        C5_total = C5
        C5_avg = C5
        else
        n=n+1
        C1_total = C1_total + C1
        C1_avg = C1_total / (n+1)
        C2_total = C2_total + C2
        C2_avg = C2_total / (n+1)
        C3_total = C3_total + C3
        C3_avg = C3_total / (n+1)
        C4_total = C4_total + C4
        C4_avg = C4_total / (n+1)
        C5_total = C5_total + C5
        C5_avg = C5_total / (n+1)        
        end if
 
        
    

end do  ;; end h
 
end do ;; end m         
    
end do  ;; end y
print("n = "+n+", C1_avg = "+C1_avg+", C2_avg = "+C2_avg+", C3_avg = "+C3_avg+", C4_avg = "+C4_avg+", C5_avg = "+C5_avg)
;con = new(5,double)
con = (/C1_avg, C2_avg, C3_avg, C4_avg, C5_avg/)
;asciiwrite("./constants_from_model_"+model_number+".txt", con)

;;;---do time average and output netcdf files

emi_rate_l = sum_scaled_hour_l / n_hour_total
emi_rate_lb = sum_scaled_hour_lb / n_hour_total
emi_rate_lbm = sum_scaled_hour_lbm / n_hour_total
emi_rate_lbs = sum_scaled_hour_lbs / n_hour_total
emi_rate_lbU = sum_scaled_hour_lbU / n_hour_total
emi_rate_lbms = sum_scaled_hour_lbms / n_hour_total
emi_rate_lbsU = sum_scaled_hour_lbsU / n_hour_total
emi_rate_lbmU = sum_scaled_hour_lbmU / n_hour_total
emi_rate_lbmsU = sum_scaled_hour_lbmsU / n_hour_total

copy_VarMeta(emi_rate_lbmsU,emi_rate_lb)
copy_VarMeta(emi_rate_lbmsU,emi_rate_lbm)
copy_VarMeta(emi_rate_lbmsU,emi_rate_lbs)
copy_VarMeta(emi_rate_lbmsU,emi_rate_lbU)
copy_VarMeta(emi_rate_lbmsU,emi_rate_lbms)
copy_VarMeta(emi_rate_lbmsU,emi_rate_lbsU)
copy_VarMeta(emi_rate_lbmsU,emi_rate_lbmU)
copy_VarMeta(emi_rate_lbmsU,emi_rate_l)

if (write_dust .eq. 1) then

dir_o = "./emission/model_"+model_number+"/"                     ; Output directory

fil_o = "dust_of_"+expname+".nc"         ; Output file
system("/bin/rm -f " + dir_o + fil_o)  ;delete exist file
f_out=addfile(dir_o + fil_o, "c")

f_out->emidust_l = emi_rate_l
f_out->emidust_lb = emi_rate_lb
f_out->emidust_lbm = emi_rate_lbm 
f_out->emidust_lbs = emi_rate_lbs
f_out->emidust_lbU = emi_rate_lbU
f_out->emidust_lbms = emi_rate_lbms
f_out->emidust_lbsU = emi_rate_lbsU
f_out->emidust_lbmU = emi_rate_lbmU
f_out->emidust_lbmsU = emi_rate_lbmsU

end if 

if (write_factor .eq. 1) then
dD_m = sum_dhour_m / n_hour_total
dD_s = sum_dhour_s / n_hour_total
dD_U = sum_dhour_U / n_hour_total
dD_b = sum_dhour_b / n_hour_total
dD_all_l = sum_dhour_all_l / n_hour_total
dD_all_lb = sum_dhour_all_lb / n_hour_total

copy_VarMeta(dD_m, dD_s)
copy_VarMeta(dD_m, dD_U)
copy_VarMeta(dD_m, dD_b)
copy_VarMeta(dD_m, dD_all_l)
copy_VarMeta(dD_m, dD_all_lb)

dir_o = "./factor/model_"+model_number+"/"                     ; Output directory
fil_o = "factors_of_"+expname+".nc" 
system("/bin/rm -f " + dir_o + fil_o)  ;delete exist file
f_out=addfile(dir_o + fil_o, "c")
f_out->dust_factor_all_l=dD_all_l
f_out->dust_factor_all_lb=dD_all_lb
f_out->dust_factor_l=emi_rate_l
f_out->dust_factor_b=dD_b
f_out->dust_factor_m=dD_m
f_out->dust_factor_s=dD_s
f_out->dust_factor_U=dD_U

end if 



end
