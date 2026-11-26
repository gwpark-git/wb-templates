
#include "parameters.h"

// SI inputs become the inherited HF scales: R and L in metres, pressures in
// Pa, eta0 in Pa s, Lp in m/(Pa s). Preserve get_Pper's linear HP convention.

double get_Pin(double DLP, double Pout)
{
  return DLP + Pout;
}

double get_Pper(double DLP, double DTP, double k, double Pout)
{
  /*
    For a small K-value which satisfies K^2\sim \epsilon_\delta, it is safe to ignore the pure-solvent correction for DTP with "2 tanh(k/2.)/k" which, in series expansion, shows 1 - K^2/12 + O(K^4).  In such case, the permeate pressure prediction now is based on the linear approximation.
   */
  // return (1./2.)*DLP - DTP/(2.*tanh(k/2.)/k) + Pout;
  return (1./2.)*DLP - DTP + Pout; 
}

int
PARAMETER_SET::
print()
{
  cout << setw(50) << "N_iter, weight, weight_f, R_channel, L_channel: " << N_iter << '\t' << weight << '\t' << weight_f << '\t' << R_channel << '\t' << L_channel << endl;
  cout << setw(50) << "Lp, eta0, T: " << Lp << '\t' << eta0 << '\t' << T << endl;

  cout << setw(50) << "k, preU: " << k << '\t' << preU << '\t' << endl;
  cout << setw(50)<< "DLP, DTP_HP, DTP: " << DLP << '\t' << DTP_HP << '\t' << DTP << endl;
  cout << setw(50)<< "Pin, Pout, Pper: " << Pin << '\t' << Pout << '\t' << Pper << endl;
  cout << setw(50)<< "a, a_H: " << a << '\t' << a_H << endl;
  cout << setw(50)<< "Va, kT, D0: " << Va << '\t' << kT << '\t' << D0 << endl;

  cout << setw(50)<< "u0, vw0: " << u0 << '\t' << vw0 << endl;
  cout << setw(50)<< "alpha_ast, beta_ast: " << alpha_ast << '\t' << beta_ast << endl;

  cout << setw(50)<< "Nz, Nr_bulk, Ny: " << Nz << '\t' << Nr_bulk << '\t' << Ny << endl;
  cout << setw(50)<< "dz, dr: " << dz << '\t' << dr << endl;
  cout << setw(50)<< "gamma: "<< gamma << endl;

  cout << setw(50)<< "phi_freeze, phi_bulk, phiw_init: " << phi_freeze << '\t' << phi_bulk << '\t' << phiw_init << endl;
  cout << setw(50)<< "epsilon_d: " << epsilon_d << endl;
  // cout << setw(50)<< "delta_f_ast, porosity, zeta_f_KC: " << delta_f_ast << '\t' << porosity << '\t' << zeta_f_KC << endl;
  cout << setw(50)<< "Pi_ast: " <<  Pi_ast << endl;
  cout << setw(50)<< "TOL_phiw, TOL_delta_f: " << TOL_phiw << '\t' << TOL_delta_f << endl;
  cout << setw(50)<< "Type of dispersion: " << endl;
  cout << setw(50)<< "                    Dc_div_D0 = " << str_Dc_div_D0 << endl;
  cout << setw(50)<< "                    eta_div_eta0 = " << str_eta_div_eta0 << endl;
  cout << setw(50)<< "                    Pi = " << str_Pi << endl;
  cout << setw(50)<< "                    Pe_a = " << Pe_a << endl;
  return 0;
}

int
PARAMETER_SET::
init(COND& given_cond)
  {
    // linear_P_profile=FALSE;
    // if (given_cond("P_profile")=="LINEAR")
    //   {
    // 	cout << "\tWARNING: P_profile is forced to be LINEAR\n";
    // 	linear_P_profile=TRUE;
    //   }
      
    // COND has already checked the complete finite integral value and its
    // <=100000 bound. Consume that value, including 2.4e3 or 0x960, not atoi's
    // integer prefix (2 or 0 for those valid spellings).
    N_iter = static_cast<long>(std::stod(given_cond("N_iter")));
    weight = atof(given_cond("weight").c_str());
    // weight_f = atof(given_cond("weight_f").c_str());
    weight_f = weight; // currently, the weight_f is set with weight.
    // weight_f = weight/10.;

    R_channel = atof(given_cond("R_channel").c_str());
    L_channel = atof(given_cond("L_channel").c_str());
    
    Nz = atof(given_cond("Nz").c_str());
    Nr_bulk = atof(given_cond("Nr").c_str());
    Ny = 0; // this will be over-written when we generate yt_arr

    dr = (1./double(Nr_bulk))*R_channel; // it is not a strictly necessary at this moment
    dz = (1./double(Nz-1))*L_channel; // it makes keeps the starting value as zero and the last value as L_channel.


    Lp = atof(given_cond("Lp").c_str());
    eta0 = atof(given_cond("eta0").c_str());
    T = atof(given_cond("T").c_str());
    k =  4.*sqrt(pow(L_channel, 2.0) * Lp * eta0 /pow(R_channel, 3.0));


    // check point inlet-pressure or inlet-velocity
    Pout = atof(given_cond("ref_Pout").c_str());
    // inlet_BC = given_cond("inlet_BC").c_str();
    inlet_BC = given_cond("inlet_BC");
    legacy_inlet_velocity = given_cond("inlet_velocity_order") == "1";
    
    // cout << "INLET BC = " << inlet_BC << endl;
    if (inlet_BC == "velocity")
      {
	u0 = atof(given_cond("u_inlet").c_str());
	double DLP_HP =  u0*4.*eta0*L_channel/(R_channel*R_channel); // HP prediction
	DLP = DLP_HP; 
	
	Pin = DLP + Pout; // HP prediction. Pin is the subject to be changed.
	Pper = get_Pper(DLP, atof(given_cond("ref_DTP").c_str()), k, Pout);
	DTP_HP = (1/2.)*(Pin + Pout) - Pper;
	DTP = DTP_HP; // this does not change anything since the inlet-velocity condition is used

      }
    else
      {
	DLP =  atof(given_cond("DLP").c_str());
	
	Pin = get_Pin(DLP, Pout);
	Pper = get_Pper(DLP, atof(given_cond("ref_DTP").c_str()), k, Pout);
	u0 = pow(R_channel, 2.0)*DLP/(4.*eta0*L_channel);
	DTP_HP = (1/2.)*(Pin + Pout) - Pper;
	DTP = DTP_HP*2.*tanh(k/2.)/k;
      }
    
    // Preserve the constant used by the original edition. GSL 2.8 uses
    // 1.3806504e-23 J/K, which differs from current exact SI k_B. Updating it
    // would be a separately labelled numerical input change.
    kT = GSL_CONST_MKSA_BOLTZMANN*T;
    preU = sqrt(Lp*R_channel/eta0);
    
    vw0 =  Lp*DTP_HP;

    dg_ast = 2.* u0 / R_channel ;
    shear_rate_ast = dg_ast;
    
    a = atof(given_cond("a_particle").c_str());
    gamma = atof(given_cond("gamma").c_str());
    a_H = a*gamma;
    D0 = kT/(6.* M_PI * eta0 * a_H);
    Va = (4./3.)* M_PI *pow(a, 3.0);
    epsilon_d = D0/(R_channel * vw0);

    Pi_ast = kT/Va; // cofficient for osmotic pressure
    alpha_ast = DTP_HP/DLP;
    beta_ast = 4.* (L_channel / R_channel)* vw0 / u0; // this uses the cylindrical geometry

    tau_ast = kT/pow(a, 3.0);

    
    phi_freeze =  atof(given_cond("phi_freeze").c_str());
    phi_bulk =  atof(given_cond("phi_bulk").c_str());
    phiw_init = atof(given_cond("phiw_init").c_str());


    TOL_phiw = atof(given_cond("TOL_phiw").c_str());
    TOL_delta_f = atof(given_cond("TOL_delta_f").c_str());

    // note that the following parameters will stay with string
    // then .c_str() will be used later to convert into char*
    str_Dc_div_D0 = given_cond("fcn_Dc_div_D0");
    str_eta_div_eta0 = given_cond("fcn_eta_div_eta0");
    str_Pi = given_cond("fcn_Pi");

    // related with cake layer
    R_mem = 1./(eta0 * Lp); // clean membrane resistance


    Pe_a = dg_ast * a_H * a_H / D0;

    if (given_cond("DETAILED_REPORT") == "TRUE")
      {
	DETAILED_REPORT = TRUE;
      }

    CAKE_convert_unit = 0.;
    h_mem = 0.;
    if (given_cond("ACCOUNT_CAKE_THICKNESS") == "TRUE")
      {
	CAKE_THICKNESS = TRUE;
	h_mem = atof(given_cond("membrane_thickness").c_str());
	double kappa_cake = (4./180.)*a*a*pow(1. - phi_cake, 3)/(phi_cake*phi_cake); // Carman-Kozeny relation based on Mulder et al.
	double kappa_mem = Lp*eta0*h_mem; // assuming h/R is small enough. 
	CAKE_convert_unit = (h_mem/R_channel)*(kappa_cake/kappa_mem);
      }
    
    
    return 0;
  }
