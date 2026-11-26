
#include <iostream>
#include <string>

#include "file_IO.h"
#include "parameters.h"
#include "dispersion_properties.h"
// #include "transport_properties_SPHS.h"
// #include "transport_properties_SID.h"
// #include "osmotic_pressure_CS.h"
#include "ported_lib.h"
#include "boost_matrix.h"


#include "matrix.h"

#define True 1
#define False 0


using namespace std;



int
help()
{
  cout << "Usage: mbla case.inp [next_case.inp ...]\nSee docs/SCIENTIFIC_GUIDE.org for SI units, output columns and validation limits.\n";
  return 0;

}

bool
print_intermediate_report(long n)
{
  if (n<=50 && n%10==0)
    return TRUE;
  else if (n<=100 && n%20==0)
    return TRUE;
  else if (n<=500 && n%50==0)
    return TRUE;
  // else if (n<=1000 && n%100==0)
  else if (n%100==0)
    return TRUE;
  // else if (n<=10000 && n%1000=0)
  //   return TRUE;
  return FALSE;
}


int
report_details(const string& base_fout, AUX_MATRIX& AUX_h1, AUX_MATRIX& AUX_h2, DISP_TYPE& DISP, PARAMETER_SET& PARS, MATRIX& z_div_L_arr, MATRIX& yt_arr)
{
  /*
    When DETAILED_REPORT=TRUE in the .inp file
    This will report all the dimensionless quantity.
    The default fout name is provided by base_fout, which indicates the xxx.dat
    Now, the report of F values and outlet values will write xxx.dat.Rout and xxx.dat.Fval    
   */
  for (long nn=0; nn<PARS.Ny; nn++)
    {
      AUX_h1.report_outlet_profile(nn, 0) = yt_arr(nn)/PARS.R_channel;
      // divided by R_channel is introduced BY MISTAKE. However, it does not change anything except the base unit. In order to avoid confusion, for the time beging, this will survive here.
      AUX_h1.report_outlet_profile(nn, 1) = AUX_h2.phi_yt_div_phib_arr(PARS.Nz - 1, nn)*PARS.phi_bulk;
      AUX_h1.report_outlet_profile(nn, 2) = AUX_h1.shear_rate_nodim(PARS.Nz - 1, nn);
      AUX_h1.report_outlet_profile(nn, 3) = AUX_h1.shear_stress_nodim(PARS.Nz - 1, nn);
      AUX_h1.report_outlet_profile(nn, 4) = DISP.fcn_Dc_div_D0(AUX_h2.phi_yt_div_phib_arr(PARS.Nz - 1, nn)*PARS.phi_bulk, PARS, AUX_h1.shear_rate_nodim(PARS.Nz - 1, nn));
      AUX_h1.report_outlet_profile(nn, 5) = DISP.fcn_eta_div_eta0(AUX_h2.phi_yt_div_phib_arr(PARS.Nz - 1, nn)*PARS.phi_bulk, PARS, AUX_h1.shear_stress_nodim(PARS.Nz - 1, nn));
					   
    }
  for (long nn=0; nn<PARS.Nz; nn++)
    {
      AUX_h1.F_info_arr(nn, 0) = z_div_L_arr(nn);
    }
  ofstream f_outlet_report((base_fout + ".Rout"), ofstream::out);
  if (!f_outlet_report) throw runtime_error("cannot create outlet report");
  AUX_h1.report_outlet_profile.fprint_basic(f_outlet_report);
  finish_output(f_outlet_report, base_fout + ".Rout");

  ofstream f_F_report((base_fout +  ".Fval"), ofstream::out);
  if (!f_F_report) throw runtime_error("cannot create integral report");
  AUX_h1.F_info_arr.fprint_basic(f_F_report);
  finish_output(f_F_report, base_fout + ".Fval");

  return 0;
}


int
run_cases
(int argc, char* argv[])
{
  if (argc==1 || (argc==2 && string(argv[1])=="--help"))
    {
      help();
      return 0;
    }
  else
    {
      for(long inp_cnt = 1; inp_cnt < argc; inp_cnt++)
	{
          /*
            This for-loop is deisgn to perform sequence of jobs using the given input file.
            The sequence is exactly the same with the given index for arguments.
            Note that the all the controls should be properly set with the given input files, which means the inheritance from the previous computation in the given input file is exactly matched with the output file name from the previous input file.
          */
	  COND given_condition(argv[inp_cnt]);
	  given_condition.cond_print();

	  PARAMETER_SET PARS;
	  PARS.init(given_condition);
	  PARS.print();	  

	  DISP_TYPE DISP(PARS);

	  MATRIX z_div_L_arr(PARS.Nz, 1, 0.);
	  gen_z_div_L_arr(z_div_L_arr, PARS);

	  MATRIX phiw_set_1(PARS.Nz, 1, PARS.phiw_init/PARS.phi_bulk); // set initial value 
	  phiw_set_1(0) = PARS.phi_bulk/PARS.phi_bulk; // force to set inlet boundary condition

	  if (given_condition("CONT")=="TRUE")
	    {
	      // when the initial wall concentration is provided
	      // this must be come with "FN_CONT", which states the data file name used as the initial condition
	      MATRIX dat_CONT(given_condition("FN_CONT").c_str());
              if (dat_CONT.rows != PARS.Nz || dat_CONT.cols < 2)
                throw runtime_error("continuation table requires Nz rows and at least two columns");
	      for (int k=1; k<PARS.Nz; k++)
		{
		  phiw_set_1(k) = dat_CONT(k, 1)/PARS.phi_bulk;
		  cout << k << '\t' << dat_CONT(k, 1) << endl;
		}
	    }
	  // phiw_set_1.copy_obj(phiw_CONT);
	  
	  MATRIX phiw_set_2(phiw_set_1);

	  MATRIX Rcake_div_Rmem_arr(PARS.Nz, 1, 0.);
	  MATRIX Rcake_div_Rmem_new_arr(PARS.Nz, 1, 0.);
	  MATRIX Pi_div_DLP_arr(PARS.Nz, 1, 0.);

	  MATRIX yt_arr;
	  gen_yt_arr(yt_arr, PARS); // generate yt_arr with adaptive step size

	  bool IDENT_brake = False;
	  double err_delta = 0.;

	  AUX_MATRIX AUX_h1(PARS), AUX_h2(PARS);//, AUX_h3(PARS);
	  // auxiliary matrix objects
	  // AUX_h1 will be used for general aspect
	  // AUX_h2 will be used for intermediate purpose, especially with the CP layer
	  // AUX_h3 will be used for intermediate purpose, especially with the cake layer
	  shear_rate_initialization(AUX_h1.shear_rate_nodim, yt_arr, PARS); // initialization with PARS.shear_rate_ast for the wall shear rate. The default setting is 2u0/R_channel 
	  shear_stress_initialization(AUX_h1.shear_stress_nodim, AUX_h1.shear_rate_nodim, yt_arr, PARS, DISP);
	  
	  long ind_onset_cake_pre = PARS.Nz + 1; // define there is no cake layer (i.e., index is outside Nz)
	  long ind_onset_cake = ind_onset_cake_pre;
	  long total_ITER = 0;
	  for (long n=0; n<PARS.N_iter; n++)
	    {
	      ind_onset_cake_pre = ind_onset_cake;
	      phiw_set_1.copy_obj(phiw_set_2); // deep copy from phiw_arr_1

	      Rcake_div_Rmem_arr.copy_obj(Rcake_div_Rmem_new_arr); // deep copy from previous Rcake_div_Rmem calculation
	      AUX_h1.shear_rate_nodim_pre.copy_obj(AUX_h1.shear_rate_nodim); // deep copy for backup the previous shear rate calculations
	      AUX_h1.shear_stress_nodim_pre.copy_obj(AUX_h1.shear_stress_nodim); // deep copy from the previously calculated shear stress
	      update_Pi_div_DLP_arr(Pi_div_DLP_arr, phiw_set_1, AUX_h1.shear_rate_nodim_pre, DISP, PARS); // update osmotic pressure
	      ind_onset_cake = phiw_update(phiw_set_2, Rcake_div_Rmem_new_arr, z_div_L_arr, phiw_set_1, ind_onset_cake_pre, Rcake_div_Rmem_arr, Pi_div_DLP_arr, yt_arr, n, DISP, PARS, AUX_h1, AUX_h2); // upd	      
	      
	      for (long i=0; i<PARS.Nz; ++i)
        if (!std::isfinite(phiw_set_2(i)) || !std::isfinite(Rcake_div_Rmem_new_arr(i)))
          throw runtime_error("nonfinite fixed-point iterate; refine/review the numerical case");
      long peak_index_phiw = phiw_set_2.argmax_1d(); // check the maximum peak index of phiw
	      double peak_value_phiw = phiw_set_2(peak_index_phiw); // check the maximum peak value of phiw
	      double err = ABS_ERR_data(phiw_set_2, phiw_set_1); // relative error between phiw_set_2 and phiw_set_1

	      double delta_f_last_1 = Rcake_div_Rmem_new_arr(PARS.Nz-1);
	      // identification delta_f_last_1 > 0 is better than ind_onset_cake
	      // because ind_onset_cake does not depends on the past history with weight


	      if (delta_f_last_1 > 0.)
	      	{
	      	  err_delta = fabs(1. - Rcake_div_Rmem_new_arr(PARS.Nz - 1)/Rcake_div_Rmem_arr(PARS.Nz - 1));
		  
	      	  if (err < PARS.TOL_phiw && err_delta < PARS.TOL_delta_f)
	      	    {
	      	      IDENT_brake = True;
	      	      printf("the iteration is converaged (with cake) at iter=%04ld\n", n);
	      	    }
	      	}
	      else
	      	{

	      	  if (err < PARS.TOL_phiw) // this is the case when phiw is converged
	      	    {
	      	      IDENT_brake = True;
	      	      printf("the iteration is converaged (w/o cake) at iter=%04ld\n", n);
	      	    }
	      	}

	      if (IDENT_brake && n>30)
	      	{
		  if (PARS.DETAILED_REPORT)
		    {
		      report_details(given_condition("fout"), AUX_h1, AUX_h2, DISP, PARS, z_div_L_arr, yt_arr);
		    }
	      	  break;
	      	}

	      if(print_intermediate_report(n))
		{
		  cout << "\nn_ITER: " << n << endl;
		  // cout << "uz(0)/u0 = " << uz_div_u0_z0 << '\t' << "vw(0)/vw0 = " << vw_div_vw0_z0 << endl;
		  cout << "cond_GT: Pin=" << PARS.Pin << ", Pout=" << PARS.Pout << endl;
		  printf("n=%ld\n[CP] phiw(0)=%4.4f, phiw(L)=%4.4f, max:(phiw(%4.3f))=%4.4f\n", n, phiw_set_2(0)*PARS.phi_bulk, phiw_set_2(PARS.Nz - 1)*PARS.phi_bulk, z_div_L_arr(peak_index_phiw), phiw_set_2(peak_index_phiw)*PARS.phi_bulk);
		  printf("norm(1. - p2/p1) : %4.3e, weight : %4.3e, weight_f : %4.3e\n",err, PARS.weight, PARS.weight_f);
		  if(PARS.str_Dc_div_D0=="SID")
		    {
		      double val_1 = SID::Dc_short_div_D0(phiw_set_2(PARS.Nz-1)*PARS.phi_bulk, PARS, AUX_h1.shear_rate_nodim(PARS.Nz-1, 0));
		      double val_2 = SPHS::Dc_short_div_D0(phiw_set_2(PARS.Nz-1)*PARS.phi_bulk, PARS);		      
		      printf("[SID] at outlet wall: (actual) Pe_a_wL = %4.3e (Pa_a_ast = %4.3e, Pe_a_wL/Pe_a_ast = %4.3e), D_SID = %4.3e, D_HS = %4.3e, D_SID/D_HS = %4.3e\n", AUX_h1.shear_rate_nodim(PARS.Nz-1, 0), PARS.Pe_a, AUX_h1.shear_rate_nodim(PARS.Nz-1, 0)/PARS.Pe_a, val_1, val_2, val_1/val_2);
		      if(PARS.str_eta_div_eta0=="TEST")
			{
			  double val_1_etar = SID::eta_div_eta0_test(phiw_set_2(PARS.Nz-1)*PARS.phi_bulk, PARS, AUX_h1.shear_stress_nodim(PARS.Nz-1, 0));
			  double val_2_etar = SID::eta_div_eta0_test(phiw_set_2(PARS.Nz-1)*PARS.phi_bulk, PARS, 0.);		      
			  printf("                      tau_wL/tau_ast  = %4.3e (tau_ast = %4.3e Pa), etar_SID = %4.3e, etar_SID/etar_EQ = %4.3e\n", AUX_h1.shear_stress_nodim(PARS.Nz-1, 0), PARS.tau_ast, val_1_etar, val_1_etar/val_2_etar);
			}
		    }
		  if (ind_onset_cake <= PARS.Nz)
		    {
		      printf("[CAKE] (z/L)[%ld] = %0.3f\n", ind_onset_cake, z_div_L_arr(ind_onset_cake));
		      printf("\tind=%ld, Rcake(L)/Rmem=%4.3e, abs(1-dR_2/dR_1)=%4.3e\n", ind_onset_cake, Rcake_div_Rmem_arr(PARS.Nz-1), err_delta);
		      printf("\td_cake/R=%4.3e, uL_div_u0=%4.3e\n", PARS.CAKE_convert_unit*Rcake_div_Rmem_arr(PARS.Nz-1), AUX_h1.uz_div_u0_arr(PARS.Nz-1));
		    }
		  
		}
	      total_ITER = n; 
	    }

	  // test
	  
	  printf("\nFINAL CHECK POINT (total iteration number = %ld):\n", total_ITER);
	  printf("COUNT force replace phiw : %ld\n", PARS.count_negative_phiw);
	  fprint_analysis_with_fouling(z_div_L_arr, phiw_set_2, Rcake_div_Rmem_new_arr, yt_arr, DISP, PARS, given_condition("fout").c_str(), AUX_h1);
          cout << "MBLA_STATUS=" << (IDENT_brake ? "converged" : "iteration_limit") << endl;
          if (!IDENT_brake) return 3;

	  
	}
      
    }
  return 0;
}

// Report invalid input/numerical failures to batch callers instead of GSL aborts.
int main(int argc, char* argv[]) {
  gsl_set_error_handler_off();
  try { return run_cases(argc, argv); }
  catch (const std::exception& error) {
    std::cerr << "mbla: " << error.what() << '\n';
    return 2;
  }
}
