/*
SNr Network Model from Whalen et al. 2021, in press
Original code by Ryan Phillips, edited by Tim Whalen, edited by John Parker, last edited Jan 2023
Contact: timcwhalen@gmail.com or owner of repo, or authors of corresponding paper

Simulates a conductance-based (or quadratic integrate & fire) SNr network model with inputs from synthetic GPe spike trains with jitter to study propagation of delta rhythms in Parkinsonism

----------

To generate results from the paper, run python script run_snr_jitter_sims.py

Running this file will produce several output text files in ./data
- .prp: properties file, which lists connection matrices, synaptic strengths and tonic excitation
- .hst: two column /t-separated of spike times and index of spiked neuron
- .gpe: like .hst, but for the synthetic GPe spike trains
- .stn: like .gpe for synthetic STN neurons (but none are used in the paper, so this file will be blank)
- .sp: voltages of each SNr neuron over time; each row is a timestep, first column is t, each subsequent column is a /t-separated voltage for every cell

See companion MATLAB code and Python code to use these outputs to generate the figures from Whalen et al. 2021
*/

#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <random>
#include <iterator>
#include <functional>
#include <cmath>

double dt;

// Equation for gate variables:
inline double x_inf(double v,double V12,double k) { return 1./(1.+exp(-(v-V12)/k)); }
inline double x_inf2(double v,double V12,double k,double min) { return (1.-min)/(1.+exp(-(v-V12)/k)); }
inline double tau_inf(double tau,double v,double tauV12,double tauk) { return tau/(cosh((v-tauV12)/tauk)); }
inline double tau_inf2(double tau0,double v,double tau1,double phi,double sigma0,double sigma1) { return tau0 + (tau1-tau0)/(exp((phi-v)/sigma0)+exp((phi-v)/sigma1)); }
inline double tau_inf3(double tau0,double v,double tau1,double tauV12,double tauk) { return tau0 + tau1/( 1 + exp((v-tauV12)/tauk) ); }
inline double Gr(double v) { return 4.34e-5*exp(-0.0211539274*v); }
inline double Gd1(double v) { return 0.075+0.043*tanh((v+20)-20); }
inline double Gv(double v) { return (10.6408- 14.6408*exp(-v/42.7671))/v; }

const double c_som=100.0, c_dend=0.4*c_som, c_total = c_som+c_dend;
double cl_coupling = 0; // on =1 off = 0

//NaF
const double mV12 = -30.2, mk = 6.2;  // from Ding 2011
const double mtau0 = 0.05, mtau1 = 0.05,mphi = 1, msigma0 = 1, msigma1 = 1;
const double hV12 = -63.3, hk = -8.1;
const double htau0 = 0.59, htau1 = 35.1, hphi = -43, hsigma0 = 10, hsigma1 = -5;
const double sV12 = -30, sk = -0.4;
const double stau0 = 10, stau1 = 50, sphi = -40, ssigma0 = 18.3, ssigma1 = -10, smin = 0.15;
double gNaF =35*c_som;

//NaP
double NAP_BLK = 1; // 0 = blocked, 1 = standard
const double napmV12 = -50, napmk = 3;
const double napmtau0 = 0.03, napmtau1 = 0.146, napmphi = -42.6, napmsigma0 = 14.4, napmsigma1 = -14.4, napmmin = 0.0;
const double naphV12 = -57, naphk = -4.0;
const double naphtau0 = 10, naphtau1 = 17, naphphi = -34, naphsigma0 = 26, naphsigma1 = -31.9, naphmin = 0.154;

//CaH (high voltage activated Ca (parameters are from Elsen 1998)
const double mcV12 = -27.5, mck = 3.0;
const double mctau = 0.5;
const double hcV12 = -52.4, hck = -5.2;
const double hctau = 18;

//CaL (Low voltage activated Ca (parameters are from Elsen 1998)
const double mcl_V12 = -59.15, mcl_k = 2.36;
const double mcl_tau = 3.2;
const double hcl_V12 = -82.4, hcl_k = -5.31;
const double hcl_tau = 11.6;

//Ca Pump
double Caout=4, alphaCa=1e-6/c_som, tauCa=250, Cain0=0.00000005; // units of Ca are in mM
double Kpump=1e-3,Vpump=Kpump/tauCa;
double ry_f = 5./1000.; // frequency of spontaneious ca release
double I_ry = 0,tau_ry = 500.0;

//Kdr
const double kmV12 = -26, kmk = 7.8;
const double kmtau0 = 0.1, kmtau1 = 14, kmphi = -26, kmsigma0 = 13, kmsigma1 = -12;
const double khV12 = -20, khk = -10;
const double khtau0 = 5.0, khtau1 = 20, khphi = 0, khsigma0 = 10, khsigma1 = -10, khmin = 0.6;
double gKdr = 50*c_som;

//CAN
double K_CAN = 1.15*1e-6, nc = 3;
double ECAN=-18.0;
double canh_V12 = -95.0, canh_k = -14.0;
double can_tau0 = 300, can_tau1 = 350, can_tau_V12 = -60, can_tau_k = 3.0;

const double hcanC12 = .00015, hcank = -.0002;
const double hcantau = 100;

//SK
double K_SK = 0.0004, nc_SK = 4, tauSK_inf =.1; // Parameters from Xia 1998 (units of ca in simulations ire in mM)

//TRPC3
const double ETRPC3 = -37;
double TRPC3_BLOCK=1;

//HCN
const double mhcn_V12 = -76.4, mhcn_k = -3.3;
const double mhcn_tau0 = 0.0, mhcn_tau1 = 3625 ,mhcn_phi = -76.4, mhcn_sigma0 = 6.56, mhcn_sigma1 = -7.48;

//NMDA
const double eta = 0.33, Mg = 1.3, gammanmda = 0.05;
const double hnmda12=0.0004, knmda = -0.00005, hnmda_tau = 500;
double Pnmda = .00;
double ENMDA = -20, gNMDA=0.000;

//ChR2
double Echr2 = 0;
double gma_chr2=0.1, eps1=0.8535, eps2=0.14;
double sig_ret=10e-20, Gd2=0.05, tau_chr2=1.3;

//Arch
double Earch = -145;

//Cl- dynamics for GABA_A
double tauCL=180, CLin0=0.0037, CLout = 120e-3;
double rho_som = 0.5, rho_dend = 1, r_dend = 2e-6, Vol_som = pow((c_som),3.0/2.0);
double alphaCL= 8.85e-5/(rho_som*pow((c_som),3.0/2.0));
double alphaCL_dend = 1.77e-11/(rho_dend*c_dend*r_dend);
double tauCL_dend=18, CLin0_dend = 0.0017;
double CsCd=c_som/(c_som+c_dend);
double gaba_fixed = 0,EGABA_fixed=-64,EGABA_dend_fixed=-80;
double hco3in =11.8e-3, hco3out=25.0e-3, phco3=1, pcl=4; // Doyon 2016 frontiers in cellular neurosci
double EHCO3=26.54*log(hco3out/hco3in);
double gcltonic = 0.0,gcltonic_dend = 0.0 ;

double vclamp=0;
double GPe_f_0 = 0, GPe_f = GPe_f_0, F_OFF = 0, D_OFF = 0, GPe_strength = 0.2, Str_strength=0.55, str_freq = 1.1;;

//Somatic and Dendritic facilitation and depression parameters

//facilitation
double f_F = 0.125, F_0 = 0.145, tau_F = 1000; // parameters for facilitation from connelly 2010
//decay
double f_D = 0.565, D_0 = 1, tau_D = 1000;    // parameters for decay: tau is taken from connelly 2010
double pD = 0.67, pDr = 1-pD;// proportion able to depress

//PV, LhX6, STN Frequency parameters
double f_pv=100.0/1000.0, f_lhx6=100.0/1000.0, f_stn = 100.0/1000.0;
double dpv = 0.0, dlhx6=0.0;
double w_pv = 0.3, w_lhx6 = 0.1, w_stn = 0.03;

//g_tonic_STN
double gSTN_ton0 = 0, gSTN_ton1 = 0;
double ginh_ton0 = 0, ginh_ton1 = 0;

double ENa=50 ,EK=-90,EHCN=-30;

const double Esyn=-64; //-Higgs 2016
const double Eexc=0.0;
const double EGABA_ton=-70; // for tonic inhibition, only EGABA not changed by fixedGABA flag
const double Edrive=-10;

double rnd() {return double(rand())/RAND_MAX; }

double tauexc=5;
double tausyn= 3; // ~2.5 from connelly 2010 for an animal at P20
double tausyn_dend=7.2;
double syncstart = 1;

bool QIF = 0; // switch from conductance-based to simple QIF model of SNr, use input flag qif to modify. TODO: make this OOP
double qif_noise_sigma = 0;
std::normal_distribution<double> qif_wnoise_dist(0,1); // standard normal, multiply by sqrt(dt) for each timestep noise
std::default_random_engine qif_generator(149084257);

class Neuron
{
	public:

// parameters
	double gCa,gCaL,gNaP,gCAN,gChR2,gArch,gSK, gHCN, gTRPC3, gSD, gKCC2, gKCC2_dend, ID;
	double E_leak,g_L;

// dynamical variables
	double v,dv,v_dend, dv_dend;
	double m,h,s; // INaF
	double napm,naph; // INaP
	double mc,hc; // ICaHigh
	double mcl,hcl; // ICaLow
	double n; // IKdr
	double km,kh; // INaP
	double mSK;
	double m_hcn;
	double h_can;
	double h_nmda, m_nmda;
	double Cain;
	double t_lastsp;
	double gsyn,gsyn_dend,gexc,ton_dend_dr,gton_stn,gton_inh;
	double tlastsp;
	double CLin, CLin_dend;
    double ECL, ECL_dend,ECL_SOMdend, ECL_dendSOM, ECa, EGABA, EGABA_dend, ISOMdend, IDENDsom, I_CL_SOMdend, I_CL_dendSOM;
	double OP1,OP2,CL1, CL2, Pchr2; //IChR2
	double Oarch;//arch channel states
	double Inoise, ICa, ISK;
	double pv_sp, lhx6_sp, stn_sp,T;
	double PVFREQ, LHX6FREQ, STNFREQ;
	double F_pre, D_pre;
	double gdendsyn[100],fdendsyn[100], gDENDSYN, fDENDSYN, gsomsyn[100],dsomsyn[100], gSOMSYN, dSOMSYN, gexc_dend[100], gEXCDEND;

	// for QIF
	double tau_qif = 1;
	double a0_qif = 1;
	double vrest_qif = -65;
	double vreset_qif = -70;

	Neuron();
	void init(double,double,double,double,double,double,double,double,double,double); // seems like this is never used
	void step(double dt,double drive,double fcan, double fsynCa, double irr, double iapplied, double ext_inh, double ext_exc);
	int spike(double vth,double dt,double drive,double fcan, double fsynCa, double irr, double iapplied, double ext_inh, double ext_exc);
};

struct connection
{
	int source,target;
	double weight;
};

struct stepReturn // hack to let Population.step to return a few arrays
{
	int sp;
	std::vector<double> ext_inh, int_inh;
};

class Population
{
	public:

	Neuron* net;
	int size;

	connection *w;
	connection *wg;
	connection *wstn;
	int n_snr_syn;
	int n_gpe_syn;
	int n_stn_syn;

	double vth;

	std::vector<std::vector<bool>> connect;
	std::vector<std::vector<bool>> connect_gpe;
	std::vector<std::vector<bool>> connect_stn;

	Population(int n,double gca0,double gca1,double gcal0,double gcal1,double gl0,double gl1,double gnap0,double gnap1,
		double gcan0,double gcan1,double gchr20,double gchr21,double garch0,double garch1,
		double snr_wght0, double snr_wght1, double gpe_wght0, double gpe_wght1, double stn_wght0, double stn_wght1,
		double stimnum,double gsk0,double gsk1,double ghcn0,double ghcn1,double gtrpc30,double gtrpc31, double gsd, double gkcc2, double gkcc2_dend,
		int n_ext_gpe, int n_ext_stn, std::vector<std::vector<bool>> conn, std::vector<std::vector<bool>> conn_gpe, std::vector<std::vector<bool>> conn_stn, std::ofstream* poutpoint,
		bool gauss_rand, int seed)
	{
		connect = conn;
		connect_gpe = conn_gpe;
		connect_stn = conn_stn;

		if (QIF) {vth = -20;} else {vth = -35;}
		size=n;
		net=new Neuron [size];

		*poutpoint<<"SNr Neuron Tonic Exc. Drives:"<<std::endl;
		std::default_random_engine stnton_generator(seed+644);
		std::normal_distribution<double> stnton_dist((gSTN_ton0+gSTN_ton1)/2,(gSTN_ton1-gSTN_ton0)/4); // only used if gauss_norm==1
		for(int i=0;i<size;i++)
		{
			net[i].gCa=c_som*gca0+c_som*(gca1-gca0)*rnd();
			net[i].gCaL=c_som*gcal0+c_som*(gcal1-gcal0)*rnd();
			net[i].gCAN=c_som*gcan0+c_som*(gcan1-gcan0)*rnd();
			net[i].gSK=c_som*gsk0+c_som*(gsk1-gsk0)*rnd();
			net[i].gNaP=c_som*gnap0+c_som*(gnap1-gnap0)*rnd();
			net[i].gHCN=c_som*ghcn0+c_som*(ghcn1-ghcn0)*rnd();
			net[i].gTRPC3=c_dend*gtrpc30+c_som*(gtrpc31-gtrpc30)*rnd();
			net[i].gSD = c_som*gsd;
			net[i].gKCC2 = c_som*gkcc2;
			net[i].gKCC2_dend = c_dend*gkcc2_dend;
			net[i].gton_inh=c_som*ginh_ton0+c_som*(ginh_ton1-ginh_ton0)*rnd();

			double wght;
			if (gauss_rand) {
				double temp_wght = -999;
				while (temp_wght < gSTN_ton0 || temp_wght > gSTN_ton1) { // outside of bounds
					temp_wght = stnton_dist(stnton_generator);
				}
				wght = temp_wght;
			} else {
				wght = gSTN_ton0+(gSTN_ton1-gSTN_ton0)*rnd();
			}
			net[i].gton_stn=c_dend*wght;

			*poutpoint<<net[i].gton_stn<<" ";

			net[i].ton_dend_dr = 0;
			net[i].g_L=c_som*gl0+(gl1-gl0)*rnd();
			net[i].gChR2=c_som*0.0;
			net[i].gArch=c_som*0.0;
			net[i].t_lastsp = 0;
			net[i].ID = i;
			if (i <= stimnum){

			net[i].gArch=c_som*garch0+(garch1-garch0)*rnd();
			net[i].gChR2=c_som*gchr20+(gchr21-gchr20)*rnd();
			}
		}
		*poutpoint<<std::endl;

		// SNr sparse connection matrix, here we change binary to weights
		std::default_random_engine generator(seed+177);
		std::normal_distribution<double> wght_dist_snr((snr_wght0+snr_wght1)/2,(snr_wght1-snr_wght0)/4); // only used if gauss_norm==1
		w=new connection [size*size];
		*poutpoint<<std::endl<<"SNr-SNr Synapse Strengths:"<<std::endl;//std::endl;

		n_snr_syn=0;
		for(int j=0;j<size;j++) { for(int i=0;i<size;i++)
		{
			if (!connect.at(i).at(j)) continue;
			w[n_snr_syn].source=i;
			w[n_snr_syn].target=j;
			double wght;
			if (gauss_rand) {
				double temp_wght = -999;
				while (temp_wght < snr_wght0 || temp_wght > snr_wght1) { // outside of bounds
					temp_wght = wght_dist_snr(generator);
				}
				wght = temp_wght;
			} else {
				wght = snr_wght0+rnd()*snr_wght1;
				//wght = 0;
			}
			w[n_snr_syn].weight=c_som*wght*connect.at(i).at(j); // TODO: delete *connect.at(i).at(j) as if already determiens these are 1's, but test to make sure it's actually binary
			*poutpoint<<c_som*wght<<" ";
			n_snr_syn++;
		} *poutpoint<<std::endl; }

		// GPe-SNr sparse connection matrix
		*poutpoint<<std::endl<<"SNr-GPe Synapse Strengths:"<<std::endl;
		std::normal_distribution<double> wght_dist_gpe((gpe_wght0+gpe_wght1)/2,(gpe_wght1-gpe_wght0)/4); // only used if gauss_norm==1
		wg=new connection [size*n_ext_gpe];
		n_gpe_syn=0;
		for(int j=0;j<size;j++) { for(int i=0;i<n_ext_gpe;i++)
		{
			if (!connect_gpe.at(i).at(j)) continue;
			wg[n_gpe_syn].source=i;
			wg[n_gpe_syn].target=j;
			double wght;
			if (gauss_rand) {
				double temp_wght = -999;
				while (temp_wght < gpe_wght0 || temp_wght > gpe_wght1) { // outside of bounds
					temp_wght = wght_dist_gpe(generator);
				}
				wght = temp_wght;
			} else {
				wght = gpe_wght0+rnd()*gpe_wght1;
			}
			wg[n_gpe_syn].weight=c_som*wght; // TODO: change all GPe_wght's to not *c_som until here
			*poutpoint<<c_som*wght<<" ";
			n_gpe_syn++;
		} *poutpoint<<std::endl; }

		// STN-SNr sparse connection matrix
		*poutpoint<<std::endl<<"SNr-STN Synapse Strengths:"<<std::endl;
		std::normal_distribution<double> wght_dist_stn((stn_wght0+stn_wght1)/2,(stn_wght1-stn_wght0)/4); // only used if gauss_norm==1
		wstn=new connection [size*n_ext_stn];
		n_stn_syn=0;
		for(int j=0;j<size;j++) { for(int i=0;i<n_ext_stn;i++)
		{
			if (!connect_stn.at(i).at(j)) continue;
			wstn[n_stn_syn].source=i;
			wstn[n_stn_syn].target=j;
			double wght;
			if (gauss_rand) {
				double temp_wght = -999;
				while (temp_wght < stn_wght0 || temp_wght > stn_wght1) { // outside of bounds
					temp_wght = wght_dist_stn(generator);
				}
				wght = temp_wght;
			} else {
				wght = stn_wght0+rnd()*stn_wght1;
			}
			wstn[n_stn_syn].weight=c_dend*wght; // TODO: change all STN_wght's to not *c_som until here
			*poutpoint<<c_dend*wght<<" ";
			n_stn_syn++;
		} *poutpoint<<std::endl; }
	}

	~Population() { delete w; delete net; }

	stepReturn step(double dt,double drive,int* spk,double fcan,double fw, double fsynCa, double irr, double iapplied, bool ext_gpe_spike[], int ext_gpe_size, bool ext_stn_spike[], int ext_stn_size)
	{
		stepReturn ret;
		int sp=0;
		double ext_inh[size] = {0}; // stim for each post neuron, all initialized to 0
		double ext_exc[size] = {0};
		double int_inh[size] = {0}; // only for returning and writing purposes

		for(int i=0;i<n_gpe_syn;i++) if(ext_gpe_spike[wg[i].source]) ext_inh[wg[i].target]+=wg[i].weight;
		for(int i=0;i<n_stn_syn;i++) if(ext_stn_spike[wstn[i].source]) ext_exc[wstn[i].target]+=wstn[i].weight;

		for(int i=0;i<size;i++) {
			spk[i]=net[i].spike(vth,dt,drive,fcan,fsynCa,irr,iapplied, ext_inh[i], ext_exc[i]); sp+=spk[i];
		}
		for(int i=0;i<n_snr_syn;i++) if(spk[w[i].source]) {
			net[w[i].target].gsyn+=w[i].weight*fw;
			int_inh[w[i].target]+=w[i].weight*fw;
		}
		for(int i=0;i<size;i++) if(spk[i]) net[i].t_lastsp = 0;

		ret.sp = sp;
		ret.ext_inh = std::vector<double>(ext_inh, ext_inh + sizeof ext_inh / sizeof ext_inh[0]);
		ret.int_inh = std::vector<double>(int_inh, int_inh + sizeof int_inh / sizeof int_inh[0]);
		return ret;
	}
};

Neuron::Neuron()
{
	v=-60+rnd()*15*syncstart;
	if(QIF) {
		tau_qif = 1;
		a0_qif = 1;
		vrest_qif = -60;
		vreset_qif = -70;
	} else {
		dv=0.0;v_dend=-60+rnd()*2*syncstart; dv_dend=0.0; m=.1; h=.9; s=.9; km=.01; kh=.9; napm=.01; naph=.04; n=.9; Cain=.00025; mc=.001; hc=.001; mcl=.1; hcl=.1; m_hcn=0.01; h_can=0.05; h_nmda=0.99; m_nmda=0; CLin = 0.006;CLin_dend = 0.006;
		gCa=.01;
		gCaL=.01;
		gNaP= 0.1;
		gCAN=0.01;
		gSK=.4;
		gKCC2=c_som*1.0;
		gKCC2_dend=c_dend*1.0;
		gHCN=0.1;
		gTRPC3=0.1;
		gChR2=0.04;
		gArch=0.04;
		E_leak=-60; g_L=0.068;
		OP1 =0.00; OP2=0.00; CL1=0.99; CL2=0.01; Pchr2=0.1;
		Oarch = 0.0;
		D_pre=D_0;F_pre=F_0;
		gDENDSYN=0; for(int i=0;i<100;i++) fdendsyn[i]=F_0; for(int i=0;i<100;i++) dsomsyn[i]=D_0;
	}
}

int Neuron::spike(double vth,double dt,double drive,double fcan, double fsynCa, double irr, double iapplied, double ext_inh, double ext_exc) // ext_inh is external inh. stim by TW
{
	double vpre=v;
	step(dt,drive,fcan,fsynCa, irr,iapplied, ext_inh, ext_exc);
	bool spiked = v>=vth;
	if (QIF && spiked) {
		v = vreset_qif;
	}
	return (vpre<vth && spiked);
}

void Neuron::init(double gca,double gcal,double glk,double gnap,double gcan,double gchr2,double garch,double gsk,double ghcn ,double gtrpc3)
{
	gCa=gca;
	gCaL=gcal;
	g_L=glk;
	gNaP=gnap;
	gCAN=gcan;
	gSK=gsk;
	gHCN=ghcn;
	gTRPC3=gtrpc3;
	gChR2=gchr2;
	gArch=garch;
	m=rnd();
	h=rnd();
	km=rnd();
	kh=rnd();
	s=rnd();
	napm=rnd();
	naph=rnd();
	mc=rnd();
	mSK=rnd();
	h_can=0.05;
	hc=rnd();
	n=rnd();
	Cain=5e-5*rnd();
	CLin = 12e-3;
	OP1=0.00; OP2=0.00; CL1=0.99; CL2=0.01; Pchr2=0.1;
	Oarch = 0.0;
	D_pre=D_0;F_pre=F_0;
	for(int i=0;i<100;i++) fdendsyn[i]=F_0;; for(int i=0;i<100;i++) dsomsyn[i]=D_0;
}

void Neuron::step(double dt,double drive,double fcan,double fsynCa, double irr, double iapplied, double ext_inh, double ext_exc) // iapplied is built-in applied current, ext_inh is external inh. stim by TW
{
	T+=dt;
	if(!QIF) {
		// Normalization: PV,LhX6, STN spikes
		pv_sp=0;
		lhx6_sp=0;
		stn_sp=0;
		// Poisson PV,LhX6, STN spikes spikes
		PVFREQ=f_pv+dpv;
		LHX6FREQ=f_lhx6-dpv+dlhx6;
		STNFREQ=f_stn-1*f_stn*sin(0.005*T)*((f_pv-dpv)/f_pv)*((f_lhx6-0.25*dlhx6)/f_lhx6);
		if(dlhx6<=0){STNFREQ=f_stn-1*f_stn*sin(0.005*T)*((f_pv-dpv)/f_pv);}

		if(PVFREQ*dt>rnd()){pv_sp=1;}
		if(LHX6FREQ*dt>rnd()){lhx6_sp=1;}
		if(STNFREQ*dt>rnd()){stn_sp=1;}

		ECa=26.54*log(Caout/Cain)/2;
		ECL=-26.54*log(CLout/CLin);
		ECL_SOMdend=-26.54*log(CLin_dend/CLin);
		ECL_dendSOM=-26.54*log(CLin/CLin_dend);
		ECL_dend=-26.54*log(CLout/CLin_dend);
		EGABA=-26.54*log((pcl*CLout +phco3*hco3out)/(pcl*CLin +phco3*hco3in));
		EGABA_dend=-26.54*log((pcl*CLout +phco3*hco3out)/(pcl*CLin_dend +phco3*hco3in));
		if(gaba_fixed==1){ECL=EGABA_fixed;EGABA=EGABA_fixed;EGABA_dend=EGABA_dend_fixed;}
		double INaF = gNaF*m*m*m*h*s;
		double INaP = gNaP*napm*napm*napm*naph*NAP_BLK;
		ICa = gCa*mc*hc*(v-ECa);
		double ICaL = gCaL*mcl*hcl;
		double IHCN = gHCN*m_hcn*(v-EHCN);
		double ITRPC3 = 0.0;
		double ITRPC3_dend = gTRPC3*(v_dend-ETRPC3)*TRPC3_BLOCK;
		double GNMDA = gNMDA*1*h_nmda/(1+eta*Mg*exp(-gammanmda*v));
		m_nmda = h_nmda/(1+eta*Mg*exp(-gammanmda*v));
		double IKdr = gKdr*km*km*km*km*kh;
		ISK = gSK/(1.+pow(K_SK/Cain,nc_SK));
		double ICAN = gCAN*fcan/(1.+pow(K_CAN/Cain,nc))*h_can;
		double IChR2 = gChR2*Gv(v)*(OP1+gma_chr2*OP2)*(v-Echr2);
		double IArch = gArch*Oarch*(v-Earch);
		ISOMdend = gSD/(CsCd)*(v_dend-v); // current from som to dend
		IDENDsom = gSD/(1+CsCd)*(v-v_dend); // current from dend to som

		I_CL_SOMdend = gSD/(CsCd)*((v_dend-v)+ECL_SOMdend);
		I_CL_dendSOM = gSD/(1+CsCd)*((v-v_dend)+ECL_dendSOM);
		double Iexc = gexc*(v-Eexc);

		Inoise = 0.0;
		dv=(-INaF*(v-ENa)-INaP*(v-ENa)-IKdr*(v-EK)-ICAN*(v-ECAN)-ISK*(v-EK)-IHCN-ITRPC3-ICa-ICaL*(v-ECa)-GNMDA*(v-ENMDA)-IChR2-IArch-g_L*(v-E_leak)-(gsyn)*(v-EGABA)-gton_inh*(v-EGABA_ton)-(drive)*(v-Edrive)+iapplied +Inoise-IDENDsom-Iexc-gSOMSYN*(v-EGABA))/c_som;
		v+=dv*dt;
		dv_dend=(-(gsyn_dend)*(v_dend-EGABA_dend) -gDENDSYN*(v_dend-EGABA_dend) -ISOMdend-ITRPC3_dend -gEXCDEND*(v_dend-Eexc) - gton_stn*(v_dend-Eexc) +ton_dend_dr)/c_dend;
		v_dend+=dv_dend*dt;
		if(vclamp == 1){v=-60;}

		m+=(x_inf(v,mV12,mk)-m)*(1.-exp(-dt/tau_inf2(mtau0,v,mtau1,mphi,msigma0,msigma1)));
		h+=(x_inf(v,hV12,hk)-h)*(1.-exp(-dt/tau_inf2(htau0,v,htau1,hphi,hsigma0,hsigma1)));
		s+=(x_inf2(v,sV12,sk,smin)-s)*(1.-exp(-dt/tau_inf2(stau0,v,stau1,sphi,ssigma0,ssigma1)));

		//fast potassium activation/inactivation
		km+=(x_inf(v,kmV12,kmk)-km)*(1.-exp(-dt/tau_inf2(kmtau0,v,kmtau1,kmphi,kmsigma0,kmsigma1)));
		kh+=(x_inf2(v,khV12,khk,khmin)-kh)*(1.-exp(-dt/tau_inf2(khtau0,v,khtau1,khphi,khsigma0,khsigma1)));

		//INaP activation/inactivation
		napm+=(x_inf(v,napmV12,napmk)-napm)*(1.-exp(-dt/tau_inf2(napmtau0,v,napmtau1,napmphi,napmsigma0,napmsigma1)));
		naph+=(x_inf2(v,naphV12,naphk,naphmin)-naph)*(1.-exp(-dt/tau_inf2(naphtau0,v,naphtau1,naphphi,naphsigma0,naphsigma1)));

		//High voltage Ca activation & inactivation
		mc+=(x_inf(v,mcV12,mck)-mc)*(1.-exp(-dt/mctau));
		hc+=(x_inf(v,hcV12,hck)-hc)*(1.-exp(-dt/hctau));

		//Low voltage Ca activation & inactivaiton
		mcl+=(x_inf(v,mcl_V12,mcl_k)-mcl)*(1.-exp(-dt/mcl_tau));
		hcl+=(x_inf(v,hcl_V12,hcl_k)-hcl)*(1.-exp(-dt/hcl_tau));

		//NMDA Inactivation
		h_nmda+=(x_inf(Cain,hnmda12,knmda)-h_nmda)*(1.-exp(-dt/hnmda_tau));

		//I_HCN
		m_hcn+=(x_inf(v,mhcn_V12,mhcn_k)-m_hcn)*(1.-exp(-dt/tau_inf2(mhcn_tau0,v,mhcn_tau1,mhcn_phi,mhcn_sigma0,mhcn_sigma1)));

		//ISK
		mSK+=(1/(1.+pow(K_SK/Cain,nc_SK)) -mSK)*(1.-exp(-dt/tauSK_inf));

		//CAN Channel Activation
		h_can+=(x_inf(v,canh_V12,canh_k)-h_can)*(1.-exp(-dt/tau_inf3(can_tau0,v,can_tau1,can_tau_V12,can_tau_k)));

		//Ca2+ and Cl- dynamics
		double ry_ca = 0;
		double ry_p = 0;

		if(t_lastsp>=0){ry_p = (t_lastsp/1750)*ry_f*dt;}
		if(t_lastsp>=1750){ry_p = ry_f*dt;}
		if(ry_p > rnd()){ry_ca =I_ry;}

		Cain+=(-alphaCa*ICa -alphaCa*ICaL*(v-ECa) + ry_ca - 5*alphaCa*Pnmda*GNMDA*(v-ECa) + (Cain0-Cain)/tauCa)*dt; // units are in mM or 10-3M

		//KCC2 pump
		double x_som = (EHCO3-EGABA)/(EHCO3-ECL);
		double x_dend = (EHCO3-EGABA_dend)/(EHCO3-ECL_dend);
		CLin+= (-alphaCL*(gKCC2*(ECL-EK)-x_som*gsyn*(v-ECL)-x_som*gcltonic*(v-ECL)) - cl_coupling*(CLin-CLin_dend)/(100*2))*dt;
		CLin_dend+= (-alphaCL_dend*(gKCC2_dend*(ECL_dend-EK)-x_dend*gsyn_dend*(v_dend-ECL_dend)-x_dend*gDENDSYN*(v_dend-ECL_dend)-x_dend*gcltonic_dend*(v-ECL_dend))+ cl_coupling*(CLin-CLin_dend)/(40*2))*dt;

		t_lastsp +=dt;

		gEXCDEND =0;
		gDENDSYN =0;
		fDENDSYN =0;
		gSOMSYN =0;
		dSOMSYN =0;

		for(int i=0;i<100;i++) {
			//Poison spike train
			if(D_OFF==1){dsomsyn[i]=1;}
			if(GPe_f*dt/1000>=rnd()){gsomsyn[i] += (pD*dsomsyn[i]+pDr)*c_som*GPe_strength/100;}
			gSOMSYN += gsomsyn[i];
			dSOMSYN += dsomsyn[i]/100;
			if(F_OFF==1){fdendsyn[i]=1;}
			if(str_freq*dt/1000>rnd()){gdendsyn[i] += c_dend*fdendsyn[i]*Str_strength/100; fdendsyn[i]+= f_F*(1-fdendsyn[i]); }
			gDENDSYN += gdendsyn[i];
			fDENDSYN += fdendsyn[i]/100;
			}

		for(int i=0;i<100;i++) {
			fdendsyn[i]+=((F_0-fdendsyn[i])/tau_F)*dt;
			dsomsyn[i]+=((D_0-dsomsyn[i])/tau_D)*dt;
			gdendsyn[i]*=exp(-dt/tausyn_dend);
			gsomsyn[i]*=exp(-dt/tausyn);
			gexc_dend[i]*=exp(-dt/tauexc);
			}

		//ChR2 gating variables
		if (gChR2 != 0) {
			double chr2irr = irr;
			if (irr <=0){chr2irr = 0;}
			double F = sig_ret*chr2irr*3.072772229e18; // last term = lamba/(w_loss*hc)
			double k1_chr2 = eps1*F*Pchr2;
			double k2_chr2 = eps2*F*Pchr2;
			double e12 = 0.011+ 0.005*log(1+chr2irr/0.024);
			double e21 = 0.008 + 0.004*log(1+chr2irr/0.024);
			double S0 = 0.5*(1+tanh(120*(100*chr2irr-0.1)));
			double dOP1 = (k1_chr2*CL1 + e21*OP2 - Gd1(v)*OP1 - e12*OP1 );
			double dOP2 = (k2_chr2*CL2 + e12*OP1 - Gd2*OP2 - e21*OP2);
			double dCL1 = (Gr(v)*CL2 + Gd1(v)*OP1 - k1_chr2*CL1);
			double dCL2 = (Gd2*OP2 - k2_chr2*CL2 - Gr(v)*CL2);
			OP1 += dt*dOP1;
			OP2 += dt*dOP2;
			CL1 += dt*dCL1;
			CL2 += dt*dCL2;
			Pchr2 += dt*((S0-Pchr2)/tau_chr2);
			}

		//Arch gating variables
		if (gArch != 0) {
			double archirr = irr;
			if (irr >=0){archirr = 0;}
			double dOarch = tanh(0.095*-1*archirr);
			Oarch += (dOarch-Oarch)*(1.-exp(-dt/25)); // open state
		}

	} else { // if QIF
		double noise = 0;
		if (qif_noise_sigma) {
			noise = sqrt(dt) * qif_noise_sigma * qif_wnoise_dist(qif_generator);
		}
		dv= (a0_qif*(v-vrest_qif) * (v-vrest_qif))/tau_qif + (-gton_stn*(v-Eexc) -gsyn*(v-EGABA_fixed))/c_som;
		v+=dv*dt + noise;
	}

	D_pre+=((D_0-D_pre)/tau_D)*dt;
	F_pre+=((F_0-F_pre)/tau_F)*dt;

	gsyn+=ext_inh*(pD*D_pre+pDr); // external inhibition
	D_pre += -f_D*D_pre;
	gsyn*=exp(-dt/tausyn);
	gsyn_dend*=exp(-dt/tausyn_dend);
	gexc+=ext_exc;
	gexc*=exp(-dt/tauexc);
}

using namespace std;

ostream& operator <<(ostream& os,Neuron& N)
{
	return (os<<N.v<<'\t'<<N.m<<'\t'<<N.h<<'\t'<<N.napm<<'\t'<<N.naph<<'\t'<<N.n<<'\t'<<N.Cain);
}

istream& operator >>(istream& is,Neuron& N)
{
    return (is>>N.v>>N.m>>N.h>>N.napm>>N.naph>>N.n>>N.Cain);
}

ostream& operator <<(ostream& os,Population& p)
{
	for(int i=0;i<p.size;i++) os<<p.net[i]<<'\t';
	return os;
}

istream& operator >>(istream& is,Population& p)
{
	for(int i=0;i<p.size;i++) is>>p.net[i];
	return is;
}

class StimulusNeuron
{
public:
	double start_stim;
	double stop_stim;
	std::function<double(const double&)> rate_func; // MADE PUBLIC TO INTRODUCE JITTER


	// returns 1 if spiked, 0 otherwise
	virtual bool step(double t){}
};

class RegularNeuron: public StimulusNeuron // regular as in pacemaking
{
	double period_stim; // period in ms
	double last_stim;
	int cellnum; // index in total GPe population
	std::ofstream* goutpoint; // output spikes, g originally meant GPe but can be any output file (e.g. for STN)
	public:
		RegularNeuron(){};
		RegularNeuron(double start, double stop, double p, int gcellnum, std::ofstream* outpoint){
			start_stim = start, stop_stim = stop, period_stim = p, last_stim=0;
			cellnum = gcellnum; goutpoint = outpoint;
		}

	bool step(double t) {
		if (t>=start_stim && t<stop_stim) {
			last_stim+=dt;
			if(last_stim>period_stim) {
				last_stim = 4*rnd()-2; // +-2ms in ISI
				*goutpoint<<(t/1000)<<'\t'<<cellnum<<'\t'<<endl;
				return 1;
			}
		}
		return 0;
	}
};

class PoissonNeuron: public StimulusNeuron
{
	double lambda;
	double refrac;
	double last_stim;
	double next_stim;
	int cellnum; // index in total GPe population
	std::ofstream* goutpoint;

	double poissonWaitTime(double lambda, double refrac){
		return refrac + log(1-rnd())/(-lambda/1000);
	}

	public:
		PoissonNeuron(){};
		PoissonNeuron(double start, double stop, double lamb, double ref, int gcellnum, std::ofstream* outpoint){
			start_stim = start, stop_stim = stop, lambda = lamb, refrac = ref;
			cellnum = gcellnum; goutpoint = outpoint;
			last_stim=0;
			next_stim = poissonWaitTime(lambda,refrac);
		}

	bool step(double t) {
		if (t>=start_stim && t<stop_stim) {
			last_stim+=dt;
			if(last_stim>next_stim) {
				last_stim = 0;
				next_stim = poissonWaitTime(lambda, refrac);
				*goutpoint<<(t/1000)<<'\t'<<cellnum<<'\t'<<endl;
				return 1;
			}
		}
		return 0;
	}
};

class PoissonOscillatingNeuron: public StimulusNeuron
{
	double lambda; // center
	double osc_lambda; // actual lambda(t) i.e. output of sine
	double refrac;
	double osc_freq;
	double osc_mod;
	double osc_delay; // not phase! measured in miliseconds, and integrated independent of freq
	double last_stim;
	int cellnum; // index in total GPe population
	std::ofstream* goutpoint;

	public:
		PoissonOscillatingNeuron(){};
		PoissonOscillatingNeuron(double start, double stop, double lamb, double ref, double freq, double mod, double delay, int gcellnum, std::ofstream* outpoint){
			start_stim = start, stop_stim = stop, lambda = lamb, refrac = ref, osc_freq = freq, osc_mod = mod, osc_delay = delay;
			cellnum = gcellnum; goutpoint = outpoint;
			last_stim=0;
		}

	bool step(double t) {
		if (t>=start_stim && t<stop_stim) {
			last_stim+=dt;
			if (last_stim>refrac){
				double osc_lambda = (lambda + osc_mod*sin(2*M_PI*osc_freq/1000*(t+osc_delay)))/1000; // rate param as f(t)
				if (rnd()<osc_lambda*dt) {
					last_stim = 0;
					*goutpoint<<(t/1000)<<'\t'<<cellnum<<'\t'<<endl;
					return 1;
				}
			}
		}
		return 0;
	}
};

// oscillation functions for PoissonFunctionNeuron
double sineWave (double t, double center, double osc_mod, double osc_freq, double osc_delay)
{
	double lambda = (center + osc_mod*sin(2*M_PI*osc_freq/1000*(t+osc_delay)))/1000;
	return lambda;
};

double squareWave (double t, double center, double osc_mod, double osc_freq, double osc_delay, double frac_up)
{
	// square wave with frac_up portion in up state, rest in down state, starting with up if osc_delay==0

	//double rnd_shift = ((double) rand()) / (double) RAND_MAX;
	//osc_freq = osc_freq * (0.5*rnd_shift + 0.75);

	//double period = 1000/(osc_freq*(0.5*rnd_shift+0.75)); // total os cperiod in ms
	double period = 1000/(osc_freq); // total os cperiod in ms
	double lambda = (center-osc_mod + 2*osc_mod*(fmod(t+osc_delay, period)/period < frac_up))/1000;
	return lambda;
};

class PoissonFunctionNeuron: public StimulusNeuron // Poisson neuron with lambda determined by arbitrary function of time
{
	// note that lambda here always refers to rate of process, not lmabda function "func"
	double osc_lambda; // actual lambda(t)
	double refrac;
	double last_stim;
	int cellnum; // index in total GPe population
	std::ofstream* goutpoint;


	public:
		std::function<double(const double&)> rate_func; // MADE PUBLIC TO INTRODUCE JITTER
		PoissonFunctionNeuron(){};
		PoissonFunctionNeuron(double start, double stop, double ref, std::function<double(const double&)> func, int gcellnum, std::ofstream* outpoint){
			start_stim = start, stop_stim = stop, refrac = ref;
			cellnum = gcellnum; goutpoint = outpoint; rate_func = func;
			last_stim=0;
		}

	bool step(double t) {
		if (t>=start_stim && t<stop_stim) {
			last_stim+=dt;
			if (last_stim>refrac){
				double osc_lambda = rate_func(t); // rate param as f(t)
				if (rnd()<osc_lambda*dt) {
					last_stim = 0;
					*goutpoint<<(t/1000)<<'\t'<<cellnum<<'\t'<<endl;
					return 1;
				}
			}
		}
		return 0;
	}
};


class ImportedNeuron: public StimulusNeuron
{
	std::vector<double> spikes;
	double tminus; // subtract from all spike times, pretend it's t=0, then add start time to get on same time as simulation
	int index = 0; // spike index for walking through vector

	public:
		ImportedNeuron(){};
		ImportedNeuron(double start, double stop, float t0, std::string spfile){
			start_stim = start, stop_stim = stop, tminus = t0;
			// not including gout for imported neurons

			// populate spikes vector from file
			std::ifstream in(spfile.c_str());
			if(!in){
				std::cerr << "Cannot open file: "<<spfile<<std::endl;
			}
			std::string str; // each line of file should be one spike time
			double strf;
			while (std::getline(in, str)){
				if(str.size() > 0) {
					strf = stof(str)-tminus;
					if (strf>0){ // greater than chosen starting point from train
						if (strf>start+stop) { // greater than needed
							break;
						}
						spikes.push_back(strf+start);
					}
				}
			}
		spikes.push_back(9999999.0); // hack to quickly avoid overflow when checking next
		in.close();
	}

	bool step(double t) {
		if (t>=spikes.at(index)) {
			index++;
			return 1;
		}
		return 0;
	}
};

bool doubleEqual(double a, double b)
{
    return fabs(a - b) < .0000001;
}

int main(int argc,char** argv)
{
	dt=.05;
	double T=1000;
	int size=1;
	double gca[2]={0.7,.7},gcal[2]={0.0, 0.0},gtrpc3[2]={0.015, 0.015},gleak[2]={0.0295, 0.0295},gnap[2]={0.175,0.175},gcan[2]={0.0,0.0},gsk[2]={0.039, 0.039}, snr_wght[2]={0.1,0.1}, gpe_wght[2]={0.2, 0.2}, stn_wght[2]={0.025, 0.025}, ghcn[2]={0.0,0.0},gsd=0.65;
	double dr[2]={0.0,0.0},fc[2]={1,1},fcw[2]={1,1}, fsca[2]={0.00,0.00}, tauDrug = -5000, gchr2[2] = {0.0,0.0}, garch[2] = {0.0,0.0}, gaba[2]={0,0}, gkcc2 = 0.005,gkcc2_dend = 0.1;
	double stimnum= 100, seed=0, iapp[2]={0, 0}, iapplied=0.0, mean_v=0, mean_gaba_som=0, mean_gaba_dend=0;
	double tpause=0, tpause1=0, tpause2=0, tpause3=0, t_ramp=0, last_net_sp = 0, beta_F = 0, beta_A=0, beta_p=0, kcc2BLK = 1;
	string fname;
	string propfile;
	int output=0;

	int arch = 0;
	float osc_freq = 1; // Hz, frequency of forced Poisson oscillation from GPe
	float prob_gpe_cross = 0.5; // probability of GPe crossing to opposite population (0 and 1 are total segregation, 0.5 is no segregation)
	float prob_snr_cross = 0.5;
	float prob_stn_cross = 0.5;
	double fracup_gpe = 0.5; // fraction up-state for GPe, only with rectangular oscillations
	double fracup_gpe2 = 0.5; // for 2nd GPe oscillating pop, if pop2_type==1 (does nothing if 0)
	double fracup_stn = 0.5;
	std::string osc_shape_gpe = "sine"; // string defining osc shape for GPe, either "sine" or "rect" for now
	std::string osc_shape_stn = "sine";
	float osc_cent_gpe = 25; // % center rate for oscillating GPe neurons
	float osc_cent_gpe2 = 25; // % for 2nd GPe oscillating pop, if pop2_type==1 (does nothing if 0)
	float osc_cent_stn = 15; // % center rate for oscillating STN neurons, 15+-10 gives mean ~12 Hz for 35% fracup
	float osc_mod_range_gpe[] = {20, 20}; // range to randomly choose osc amplitude
	float osc_mod_range_gpe2[] = {20, 20}; // for 2nd GPe oscillating pop, if pop2_type==1 (does nothing if 0)
	float osc_mod_range_stn[] = {10, 10}; // range to randomly choose osc amplitude
	float reg_rate_range_gpe[] = {25, 25}; // range of rates to randomly choose regular GPe neurons
	int n_delays_gpe = 0; // number of GPe neurons to have same delay, 0 indicates all same delay
	int n_delays_stn = 0;
	int n_gpe = 0;

	// which distirbution to pick delays from
	// if gaussian, need delay_std
	// if uniform, need delay_uniform_max_gpe (min is 0)
	std::string phase_dist_type_gpe = "gaussian";
	float delay_std_gpe = 34.6164; // derived from real GPe data in sim_sdf_xcorr.m, can be overwritten in flags
	float delay_uniform_max_gpe = 1000/osc_freq; // max of uniform distribution, default gives whole phase circle but must be changed if osc_freq is given in flags

	std::string phase_dist_type_stn = "gaussian";
	float delay_std_stn = 34.6164; // NOT DERIVED FROM REAL DATA YET this is just GPe's value
	float delay_uniform_max_stn = 1000/osc_freq; // max of uniform distribution, default gives whole phase circle but must be changed if osc_freq is given in flags
	float bip_snr = 0.5; // for arch 2, probability that synapse is SNr synapse (p for binomial distirbution)

	bool write_syn = 0; // will write exti, inti and gsyn files at each timestep
	bool write_ca = 0; // will write ca and ksk files at each timestep

	bool control_model = 0; // run with regular neurons instead of oscillating, only for arch 21+ TODO
	bool gauss_rand = 0; // to choose random parameters with truncated Gaussian (1st 2 std's) instead of uniform

	int implement_jitter = 0;
	int jitt1 = 5000;
	int jitt2 = 10000;
	double shift = 0.01;
	// TODO: gauss_rand currently only for synaptic weights and stnton
	for(int i=1;i<argc;i++)
	{
		if(strcmp(argv[i],"-o")==0) output=1;
		else if(strcmp(argv[i],"-dt")==0) dt=atof(argv[++i]);
		else if(strcmp(argv[i],"-n_pop")==0) size=atoi(argv[++i]);
		else if(strcmp(argv[i],"-T")==0) T=atof(argv[++i])*1000;
		else if(strcmp(argv[i],"-d")==0) { dr[0]=atof(argv[++i]); dr[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-w")==0) { snr_wght[0]=atof(argv[++i]); snr_wght[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-wg")==0) { gpe_wght[0]=atof(argv[++i]); gpe_wght[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-wstn")==0) { stn_wght[0]=atof(argv[++i]); stn_wght[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-tgs")==0) tausyn=atof(argv[++i]);
		else if(strcmp(argv[i],"-nap")==0) { gnap[0]=atof(argv[++i]); gnap[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-na")==0) gNaF=atof(argv[++i]);
		else if(strcmp(argv[i],"-kdr")==0) gKdr=atof(argv[++i]);
		else if(strcmp(argv[i],"-cahigh")==0) { gca[0]=atof(argv[++i]); gca[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-calow")==0) { gcal[0]=atof(argv[++i]); gcal[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-glk")==0) { gleak[0]=atof(argv[++i]); gleak[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-tca")==0) tauCa=atof(argv[++i]);
		else if(strcmp(argv[i],"-kp")==0) Kpump=atof(argv[++i]);
		else if(strcmp(argv[i],"-vp")==0) Vpump=atof(argv[++i]);
		else if(strcmp(argv[i],"-can")==0) { gcan[0]=atof(argv[++i]); gcan[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-sk")==0) { gsk[0]=atof(argv[++i]); gsk[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-hcn")==0) { ghcn[0]=atof(argv[++i]); ghcn[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-fc")==0) { fc[0]=atof(argv[++i]); fc[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-fw")==0) { fcw[0]=atof(argv[++i]); fcw[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-fsca")==0) { fsca[0]=atof(argv[++i]); fsca[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-ca0")==0) Cain0=atof(argv[++i]);
		else if(strcmp(argv[i],"-gchr2")==0) { gchr2[0]=atof(argv[++i]); gchr2[1]=atof(argv[++i]); }//chr2 channel conductance
		else if(strcmp(argv[i],"-garch")==0) { garch[0]=atof(argv[++i]); garch[1]=atof(argv[++i]); }//Arch channel conductance
		else if(strcmp(argv[i],"-stimnum")==0) stimnum=atof(argv[++i]); //number of neurons that are stimulated by opto stimulation (arch or chr2)
		else if(strcmp(argv[i],"-sd")==0) seed=atof(argv[++i]);
		else if(strcmp(argv[i],"-gnmda")==0) gNMDA=atof(argv[++i]);
		else if(strcmp(argv[i],"-CLin0")==0) CLin0=atof(argv[++i]);
		else if(strcmp(argv[i],"-CLin0_dend")==0) CLin0_dend=atof(argv[++i]);
		else if(strcmp(argv[i],"-iapp")==0) { iapp[0]=c_som*atof(argv[++i]); iapp[1]=c_som*atof(argv[++i]); }
		else if(strcmp(argv[i],"-c3")==0) { gtrpc3[0]=atof(argv[++i]); gtrpc3[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-pnmda")==0) Pnmda=atof(argv[++i]);
		else if(strcmp(argv[i],"-gaba")==0) { gaba[0]=atof(argv[++i]); gaba[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-iry")==0) I_ry=atof(argv[++i]);
		else if(strcmp(argv[i],"-DD")==0)  gtrpc3[0]=atof(argv[++i]);
		else if(strcmp(argv[i],"-fixedGABA")==0) { gaba_fixed=1; EGABA_fixed=atof(argv[++i]); }
		else if(strcmp(argv[i],"-fixedGABAdend")==0) { gaba_fixed=1; EGABA_dend_fixed=atof(argv[++i]); }
		else if(strcmp(argv[i],"-tauCL")==0)  tauCL=atof(argv[++i]);
		else if(strcmp(argv[i],"-tauCL_dend")==0)  tauCL_dend=atof(argv[++i]);
		else if(strcmp(argv[i],"-vclamp")==0)  vclamp=atof(argv[++i]);
		else if(strcmp(argv[i],"-gkcc2")==0)  gkcc2=atof(argv[++i]);
		else if(strcmp(argv[i],"-gkcc2_d")==0)  gkcc2_dend=atof(argv[++i]);
		else if(strcmp(argv[i],"-alphacl")==0)  alphaCL=atof(argv[++i]);
		else if(strcmp(argv[i],"-fD")==0)  f_D=atof(argv[++i]);
		else if(strcmp(argv[i],"-gclton")==0)  gcltonic =c_som*atof(argv[++i]);
		else if(strcmp(argv[i],"-gclton_d")==0)  gcltonic_dend =c_dend*atof(argv[++i]);
		else if(strcmp(argv[i],"-gsd")==0)  gsd =atof(argv[++i]);
		else if(strcmp(argv[i],"-doff")==0)  D_OFF =atof(argv[++i]);
		else if(strcmp(argv[i],"-foff")==0)  F_OFF =atof(argv[++i]);
		else if(strcmp(argv[i],"-gpe_w")==0)  GPe_strength  =atof(argv[++i]);
		else if(strcmp(argv[i],"-str_w")==0)  Str_strength =atof(argv[++i]);
		else if(strcmp(argv[i],"-bA")==0)  beta_A =atof(argv[++i]);
		else if(strcmp(argv[i],"-bF")==0)  beta_F =atof(argv[++i]);
		else if(strcmp(argv[i],"-gpeF")==0)  GPe_f_0 =atof(argv[++i]);
		else if(strcmp(argv[i],"-cl_coupling")==0)  cl_coupling =atof(argv[++i]);
		else if(strcmp(argv[i],"-kcc2blk")==0)  kcc2BLK =atof(argv[++i]);
		else if(strcmp(argv[i],"-syncstart")==0)  syncstart =atof(argv[++i]);
		else if(strcmp(argv[i],"-dep_off")==0) f_D=0;
		else if(strcmp(argv[i],"-fac_off")==0) f_F=0;
		else if(strcmp(argv[i],"-gstnton")==0) { gSTN_ton0=atof(argv[++i]); gSTN_ton1=atof(argv[++i]); }
		else if(strcmp(argv[i],"-arch")==0) arch=atoi(argv[++i]); // architecture - 1 for pathways, 2 for competitive
		else if(strcmp(argv[i],"-pg_cross")==0) prob_gpe_cross=atof(argv[++i]); // for partial segregation architectures
		else if(strcmp(argv[i],"-ps_cross")==0) prob_snr_cross=atof(argv[++i]);
		else if(strcmp(argv[i],"-pstn_cross")==0) prob_stn_cross=atof(argv[++i]);
		else if(strcmp(argv[i],"-fracup_gpe")==0) fracup_gpe=atof(argv[++i]); // fraction up-state for rectangular oscillations
		else if(strcmp(argv[i],"-fracup_gpe2")==0) fracup_gpe2=atof(argv[++i]);
		else if(strcmp(argv[i],"-fracup_stn")==0) fracup_stn=atof(argv[++i]);
		else if(strcmp(argv[i],"-osc_freq")==0) osc_freq=atof(argv[++i]);
		else if(strcmp(argv[i],"-osc_shape_gpe")==0) osc_shape_gpe=argv[++i];
		else if(strcmp(argv[i],"-osc_shape_stn")==0) osc_shape_stn=argv[++i];
		else if(strcmp(argv[i],"-osc_mod_gpe")==0) { osc_mod_range_gpe[0]=atof(argv[++i]); osc_mod_range_gpe[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-osc_mod_gpe2")==0) { osc_mod_range_gpe2[0]=atof(argv[++i]); osc_mod_range_gpe2[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-osc_mod_stn")==0) { osc_mod_range_stn[0]=atof(argv[++i]); osc_mod_range_stn[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-osc_cent_gpe")==0) osc_cent_gpe=atof(argv[++i]);
		else if(strcmp(argv[i],"-osc_cent_gpe2")==0) osc_cent_gpe2=atof(argv[++i]);
		else if(strcmp(argv[i],"-osc_cent_stn")==0) osc_cent_stn=atof(argv[++i]);
		else if(strcmp(argv[i],"-reg_rate_gpe")==0) { reg_rate_range_gpe[0]=atof(argv[++i]); reg_rate_range_gpe[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-n_delays_gpe")==0) n_delays_gpe=atoi(argv[++i]);
		else if(strcmp(argv[i],"-n_delays_stn")==0) n_delays_stn=atoi(argv[++i]);
		else if(strcmp(argv[i],"-n_gpe")==0) n_gpe=atoi(argv[++i]); // # of simulated GPe neurons, currently only used in arch 18
		else if(strcmp(argv[i],"-delay_std_gpe")==0) delay_std_gpe=atof(argv[++i]);
		else if(strcmp(argv[i],"-delay_std_stn")==0) delay_std_stn=atof(argv[++i]);
		else if(strcmp(argv[i],"-delay_uniform_max_gpe")==0) delay_uniform_max_gpe=atof(argv[++i]);
		else if(strcmp(argv[i],"-delay_uniform_max_stn")==0) delay_uniform_max_stn=atof(argv[++i]);
		else if(strcmp(argv[i],"-phase_dist_type_gpe")==0) phase_dist_type_gpe=argv[++i];
		else if(strcmp(argv[i],"-phase_dist_type_stn")==0) phase_dist_type_stn=argv[++i];
		else if(strcmp(argv[i],"-bip_snr")==0) bip_snr=atof(argv[++i]);
		else if(strcmp(argv[i],"-control")==0) control_model=atoi(argv[++i]);
		else if(strcmp(argv[i],"-gauss_rand")==0) gauss_rand=atoi(argv[++i]);
		else if(strcmp(argv[i],"-write_syn")==0) write_syn=atoi(argv[++i]);
		else if(strcmp(argv[i],"-write_ca")==0) write_ca=atoi(argv[++i]);
		else if(strcmp(argv[i],"-qif")==0) QIF=atoi(argv[++i]);
		else if(strcmp(argv[i],"-qif_sigma")==0) qif_noise_sigma=atof(argv[++i]);
		else if(strcmp(argv[i],"-jitt1")==0) jitt1=atof(argv[++i]);
		else if(strcmp(argv[i],"-jitt2")==0) jitt2=atof(argv[++i]);
		else if(strcmp(argv[i],"-shift")==0) shift=atof(argv[++i]);
		else if(strcmp(argv[i],"-implement_jitter")==0) implement_jitter=atof(argv[++i]);
		else fname=argv[i];
	}
	std::ofstream pout(fname+".prp"); // outfile for various network properties
	std::ofstream gout(fname+".gpe"); // outfile for gpe spike times (like SNr's .hst)
	std::ofstream stout(fname+".stn"); // outfile for stn spike times (like SNr's .hst)
	std::ofstream extiout; // outfile for each cell's external (gpe) inhibition at each step
	std::ofstream intiout; // outfile for each cell's internal (snr) inhibition at each stepc
	std::ofstream gsynout; // outfile for each cell's gsyn
	std::ofstream caout; // outfile for each cell's ca current
	std::ofstream kskout; // outfile for each cell's KSK current
	if (write_syn) {
		extiout.open(fname+".exti");
		intiout.open(fname+".inti");
		gsynout.open(fname+".gsyn");
	}
	if (write_ca) {
		caout.open(fname+".ca");
		kskout.open(fname+".ksk");
	}

	int n_import_gpe = 0;
	int size_pop1 = 0; // # of SNr and GPe neurons in first population, typically oscillating pop
	int size_pop2 = 0; // for 2nd population, typically Poisson / Poisson-receiving pop
	int size_pop_reg = 0; // for a regularly firing population of GPe neurons (irrelevant to SNr) TODO: implement past arch 21
	bool pop2_type = 0; // if 0, pop2 GPe neurons are Poisson. If 1, another oscillating pop
	int n_import_stn = 0;
	int n_pois_stn = 0;
	int n_osc_stn = 0;
	std::vector<std::vector<bool>> conn_gpe; // GPe-SNr stim
	std::vector<std::vector<bool>> conn_stn; // STN-SNr stim
	std::vector<std::vector<bool>> conn; // intra-SNr
	bool rand_phases_gpe = 1; // pick random phases for each set of same_delay_modulo_gpe GPe neurons
	bool rand_stn_phases = 1; // pick random phases for each set of same_delay_modulo_gpe GPe neurons
	std::vector<float> phases_gpe; // only need if rand_phases_gpe==0
	std::vector<float> phases_stn; // only need if rand_stn_phases==0
	int same_delay_modulo_gpe = 1; // set to >1 e.g. for arch 5, if you want the same delay for a set of consecutive neurons
	int same_delay_modulo_stn = 1;
	switch(arch) {
		case 1: { // two pathways, with full or partial segregation
			// prob_gpe_cross and prob_snr_cross set in input flags
			rand_phases_gpe = 1; // same_delay_modulo_gpe makes them all the same
    		size_pop1 = size/2; // if segregated pathways, first n_osc GPe go to first n_osc SNr
    		size_pop2 = size-size_pop1; // these define both GPe and SNr pop sizes
    		int n_conn_g2s = 4; // if larger than size/2, infinite loops will happen while trying to avoid double synapses
			int n_conn_s2s = 4;

			if (n_delays_gpe) { same_delay_modulo_gpe = n_delays_gpe;} // input flag
    		else { same_delay_modulo_gpe = size; } // all same delay

			std::vector<bool> cs1(size, 0); // rows = presynaptic GPe, columns postsynaptic SNr
			for (int i = 0; i < (size_pop1+size_pop2); ++i) {
				conn_gpe.push_back(cs1); // initialize conn mat of zeros
				conn.push_back(cs1);
			}

			// GPe->SNr connections
			pout<<"SNr neurons receiving GPe input from:"<<std::endl;
			for (int s = 0; s < size; ++s) { // s = SNr neuron index
				for (int i = 0; i < n_conn_g2s; ++i) {
					int g; // GPe neuron index
					if (rnd() > (prob_gpe_cross*(s<size_pop1) + (1-prob_gpe_cross)*(s>=size_pop1))) { // don't cross if pop 1, do cross if pop2
						while(true) { // avoid double synapses
							g = floor(rnd()*size_pop1);
							if (!(conn_gpe.at(g).at(s))) { break; } // new synapse
						}
					}
					else {
						while(true) {
							g = floor(rnd()*size_pop2) + size_pop1;
							if (!(conn_gpe.at(g).at(s))) { break; } // new synapse
						}
					}
					conn_gpe.at(g).at(s) = 1;
					pout<<g<<" ";
				}
				pout<<std::endl;
			}

			// SNr->SNr connections   separated only to print to .prp in separate sections...
			pout<<"SNr neurons receiving SNr input from:"<<std::endl;
			for (int s = 0; s < size; ++s) { // s = post neuron index
				for (int i = 0; i < n_conn_s2s; ++i) {
					int g; // pre neuron index
					if (rnd() > (prob_snr_cross*(s<size_pop1) + (1-prob_snr_cross)*(s>=size_pop1))) {
						while(true) { // avoid double synapses
							g = floor(rnd()*size_pop1);  // don't cross if pop 1, do cross if pop2
							if (!(conn.at(g).at(s))) { break; } // new synapse
						}
					}
					else {
						while(true) { // avoid double synapses
							g = floor(rnd()*size_pop2) + size_pop1;
							if (!(conn.at(g).at(s))) { break; } // new synapse
						}
					}
					conn.at(g).at(s) = 1;
					pout<<g<<" ";
				}
				pout<<std::endl;
			}
    		break;
		}
		case 2: { // competitive model: binomial probability of picking number of GPe/SNr synapses
			if (control_model) {
				size_pop_reg = size;
			} else {
				// size_pop1 = size;
				size_pop1 = round(size/2);
				size_pop2 = round(size/2);
			}
			rand_phases_gpe = 1; // same_delay_modulo_gpe makes them all the same

			int n_conn = 8; // total number of synapses each SNr neuron receives, whether from SNr or GPe (n for binomial distribution)
			// bip_snr (binomial probability) set in input flags
			if (n_delays_gpe) { same_delay_modulo_gpe = n_delays_gpe;} // input flag
    		else { same_delay_modulo_gpe = size; } // all same delay

 			std::default_random_engine bi_generator;
  			std::binomial_distribution<int> bi_dist(n_conn,bip_snr);
			int n_syn[size]={}; // # SNr synapses, n_conn-n_syn is # GPe synapses
			std::vector<bool> cs1(size, 0); // rows = presynaptic GPe, columns postsynaptic SNr

			for (int i=0; i<size; ++i) {
				int number = bi_dist(bi_generator);
				n_syn[i] = number;

				conn_gpe.push_back(cs1); // initialize conn mat of zeros
				conn.push_back(cs1);
			}

			// GPe->SNr connections
			pout<<"SNr neurons receiving GPe input from:"<<std::endl;
			for (int s = 0; s < size; ++s) { // s = SNr neuron index
				for (int i = 0; i < n_conn-n_syn[s]; ++i) {
					int g; // GPe neuron index
					while(true) { // avoid double synapses
						g = floor(rnd()*size);
						if (!(conn_gpe.at(g).at(s))) { break; } // new synapse
					}
					conn_gpe.at(g).at(s) = 1;
					pout<<g<<" ";
				}
				pout<<std::endl;
			}

			// SNr->SNr connections   separated only to print to .prp in separate sections...
			pout<<"SNr neurons receiving SNr input from:"<<std::endl;
			for (int s = 0; s < size; ++s) { // s = post neuron index
				for (int i = 0; i < n_syn[s]; ++i) {
					int g; // pre neuron index
					while(true) { // avoid double synapses
						g = floor(rnd()*size);  // don't cross if pop 1, do cross if pop2
						if (!(conn.at(g).at(s))) { break; } // new synapse
					}
					conn.at(g).at(s) = 1;
					pout<<g<<" ";
				}
				pout<<std::endl;
			}
    		break;
		}
		cout<<"Invalid or unspecified model architecture (arch = " << arch;
	}
	int n_ext_gpe = n_import_gpe + size_pop1 + size_pop2 + size_pop_reg; // total # of GPe neurons
	int n_ext_stn = n_import_stn + n_pois_stn + n_osc_stn;

	if(seed !=0) srand (seed);
	//seed = time(NULL);
	//srand(seed); // line impelmented to test randomness
	int freq=int(.05/dt);
	double g_dend_sum = 0;
	double g_dend_ave = 0;

	gkcc2 = gkcc2*kcc2BLK;
	gkcc2_dend = gkcc2_dend*kcc2BLK;

	if(beta_A !=0){beta_p =rnd()*10000;}

	Population pop(size,gca[0],gca[1],gcal[0],gcal[1],gleak[0],gleak[1],gnap[0],gnap[1],gcan[0],gcan[1],gchr2[0],gchr2[1],garch[0],garch[1],snr_wght[0],snr_wght[1],gpe_wght[0],gpe_wght[1],stn_wght[0],stn_wght[1],stimnum,gsk[0],gsk[1],ghcn[0],ghcn[1],gtrpc3[0],gtrpc3[1],gsd,gkcc2,gkcc2_dend, n_ext_gpe, n_ext_stn, conn, conn_gpe, conn_stn, &pout, gauss_rand, seed );

	ofstream out(fname+".sp");

	int sp=0;
	int h[size];
	double drive = 0;
	double lastSP[size];
	double isi[size];
	double spf[size];

	//pv
	double pv_lastSP=0;
	double pv_isi=0;
	double pv_spf=0;
	double lhx6_lastSP=0;
	double lhx6_isi=0;
	double lhx6_spf=0;
	double stn_lastSP=0;
	double stn_isi=0;
	double stn_spf=0;
	double pop_sp = 0;

	double irr=0;
	for(int i=0;i<size;i++) h[i]=0;
	for(int i=0;i<size;i++) lastSP[i]=0;
	for(int i=0;i<size;i++) isi[i]=0;
	for(int i=0;i<size;i++) spf[i]=0.00001;

	//// Populate external neuron vector
	std::vector<StimulusNeuron*> ext_gpe;
	std::vector<StimulusNeuron*> ext_stn;

	//// Oscillating neurons
	std::default_random_engine generator(seed+42);

	// phases TODO: find a way to make an arbitrary distribution so both don't have to be defined but one unused
	std::normal_distribution<double> delay_dist_gauss_gpe(0,delay_std_gpe); // delay distribution, ms
	std::uniform_real_distribution<double> delay_dist_uniform_gpe(0,delay_uniform_max_gpe);

	std::normal_distribution<double> delay_dist_gauss_stn(0,delay_std_stn);
	std::uniform_real_distribution<double> delay_dist_uniform_stn(0,delay_uniform_max_stn);

	// modulations
	std::uniform_real_distribution<double> mod_dist_gpe(osc_mod_range_gpe[0], osc_mod_range_gpe[1]);
	std::uniform_real_distribution<double> mod_dist_gpe2(osc_mod_range_gpe2[0], osc_mod_range_gpe2[1]);
	std::uniform_real_distribution<double> mod_dist_stn(osc_mod_range_stn[0], osc_mod_range_stn[1]);
	pout<<"# of Oscillating GPe Neurons: " << size_pop1 << std::endl;
	pout<<"GPe Oscillation Shape: " << osc_shape_gpe << std::endl;
	pout<<"GPe Oscillating Neuron Delays, Mod Indices:"<<std::endl;
	double delay;
	double mod;
	for(int i = 0; i < size_pop1; i++){
		if (rand_phases_gpe) {
			if (i % same_delay_modulo_gpe == 0){
				if (phase_dist_type_gpe == "gaussian") {
					delay = delay_dist_gauss_gpe(generator);
				} else if (phase_dist_type_gpe == "uniform") {
					delay = delay_dist_uniform_gpe(generator);
				}
				else {throw std::invalid_argument("phase_dist_type_gpe must be \"gaussian\" or \"uniform\"");}
			}
		} else {
			delay = phases_gpe.at(i);
		}

		mod = mod_dist_gpe(generator);
		pout<<delay<<", "<<mod<<std::endl;

		PoissonFunctionNeuron* n; // TODO: make new lambda func that can just take mod, delay, etc. but has rect, sine, etc. already determined before loop
		//PoissonFunctionNeuronJitter* k;
		if (osc_shape_gpe.compare("sine")) {
			//k = new PoissonFunctionNeuronJitter(0, T+1, 1, osc_freq, osc_cent_gpe,mod,delay,fracup_gpe, i, &gout);
			n = new PoissonFunctionNeuron(0, T+1, 1, [&, osc_cent_gpe, mod, osc_freq, delay](double t) {return sineWave(t, osc_cent_gpe, mod, osc_freq, delay);}, i, &gout);

		}
		else if (osc_shape_gpe.compare("rect")) {
			//double rnd_shift = ((double) rand()) / (double) RAND_MAX;
			//double new_freq = osc_freq * (0.5*rnd_shift + 0.75);
			//k = new PoissonFunctionNeuronJitter(0, T+1, 1, osc_freq, osc_cent_gpe,mod,delay,fracup_gpe, i, &gout);
			//k = new PoissonFunctionNeuronJitter(0, T+1, 1, [&, osc_cent_gpe, mod, osc_freq, delay, fracup_gpe](double t) {return squareWave(t, osc_cent_gpe, mod, osc_freq, delay, fracup_gpe);}, i, &gout);
			n = new PoissonFunctionNeuron(0, T+1, 1, [&, osc_cent_gpe, mod, osc_freq, delay, fracup_gpe](double t) {return squareWave(t, osc_cent_gpe, mod, osc_freq, delay, fracup_gpe);}, i, &gout);

		}
		else {throw std::invalid_argument("osc_shape_gpe must be \"sine\" or \"rect\"");}

		ext_gpe.push_back(n);
		//ext_gpe.push_back(k);
	}
	pout<<std::endl;

	//// Poisson neurons (or a 2nd oscillaitng population)
	if (!pop2_type) { // standard Poisson pop
		pout<<"# of Poisson GPe Neurons: " << size_pop2 << std::endl;
		for(int i = 0; i < size_pop2; i++){
			PoissonNeuron* n = new PoissonNeuron(0, 200000, 25, 1, size_pop1+i, &gout);
			ext_gpe.push_back(n);
		}
	} else { // 2nd osc pop instead of Poisson
		pout<<"# of Oscillating GPe Neurons (in 2nd pop): " << size_pop2 << std::endl;
		pout<<"GPe 2 Oscillation Shape: " << osc_shape_gpe << std::endl;
		pout<<"GPe 2 Oscillating Neuron Delays, Mod Indices:"<<std::endl;
		for(int i = 0; i < size_pop2; i++){
			if (rand_phases_gpe) {
				if ((i % same_delay_modulo_gpe == 0) && (same_delay_modulo_gpe!=size)){ // extra cond to keep same delay as first pop if all same delay
					if (phase_dist_type_gpe == "gaussian") {
						delay = delay_dist_gauss_gpe(generator);
					} else if (phase_dist_type_gpe == "uniform") {
						delay = delay_dist_uniform_gpe(generator);
					}
					else {throw std::invalid_argument("phase_dist_type_gpe must be \"gaussian\" or \"uniform\"");}
				}
			} else {
				delay = phases_gpe.at(i);
			}

			mod = mod_dist_gpe2(generator);
			pout<<delay<<", "<<mod<<std::endl;
			PoissonFunctionNeuron* n;
		//	PoissonFunctionNeuronJitter* k;
			if (osc_shape_gpe.compare("sine")) {
				n = new PoissonFunctionNeuron(0, T+1, 1, [&, osc_cent_gpe2, mod, osc_freq, delay](double t) {return sineWave(t, osc_cent_gpe2, mod, osc_freq, delay);}, i+size_pop1, &gout);
			}
			else if (osc_shape_gpe.compare("rect")) {
				//double rnd_shift = ((double) rand()) / (double) RAND_MAX;
				//double new_freq = osc_freq * (0.5*rnd_shift + 0.75);
				//n = new PoissonFunctionNeuron(0, T+1, 1, [&, osc_cent_gpe2, mod, new_freq, delay, fracup_gpe2](double t) {return squareWave(t, osc_cent_gpe2, mod, new_freq, delay, fracup_gpe2);}, i+size_pop1, &gout);
				n = new PoissonFunctionNeuron(0, T+1, 1, [&, osc_cent_gpe2, mod, osc_freq, delay, fracup_gpe2](double t) {return squareWave(t, osc_cent_gpe2, mod, osc_freq, delay, fracup_gpe2);}, i+size_pop1, &gout);

			}
			else {throw std::invalid_argument("osc_shape_gpe must be \"sine\" or \"rect\"");}
			ext_gpe.push_back(n);
		//	ext_gpe.push_back(k);
		}
	}

	//// Regular neurons
	std::uniform_real_distribution<double> reg_rate_dist_gpe(reg_rate_range_gpe[0], reg_rate_range_gpe[1]);
	double reg_period = 0;

	if (size_pop_reg>0) {
		pout<<"# of Regular GPe Neurons: " << size_pop_reg << std::endl;
		for(int i = 0; i < size_pop_reg; i++){
			reg_period = 1000/reg_rate_dist_gpe(generator);
			RegularNeuron* n = new RegularNeuron(0, 200000, reg_period, size_pop1+size_pop2+i, &gout);
			ext_gpe.push_back(n);
		}
	}

	pout<<"# of Oscillating STN Neurons: " << n_osc_stn << std::endl;
	pout<<"STN Oscillation Shape: " << osc_shape_stn << std::endl;
	pout<<"STN Oscillating Neuron Delays, Mod Indices:"<<std::endl;
	for(int i = 0; i < n_osc_stn; i++){
		if (rand_stn_phases) {
			if (i % same_delay_modulo_stn == 0){
				if (phase_dist_type_stn == "gaussian") {
					delay = delay_dist_gauss_stn(generator);
				} else if (phase_dist_type_stn == "uniform") {
					delay = delay_dist_uniform_stn(generator);
				}
				else {throw std::invalid_argument("phase_dist_type_stn must be \"gaussian\" or \"uniform\"");}
			}
		} else {
			delay = phases_stn.at(i);
		}

		mod = mod_dist_stn(generator);
		pout<<delay<<", "<<mod<<std::endl;
		PoissonFunctionNeuron* n;
		if (osc_shape_stn.compare("sine")) {
			n = new PoissonFunctionNeuron(0, T+1, 1, [&, osc_cent_stn, mod, osc_freq, delay](double t) {return sineWave(t, osc_cent_stn, mod, osc_freq, delay);}, i, &stout);
		}
		else if (osc_shape_stn.compare("rect")) {
			n = new PoissonFunctionNeuron(0, T+1, 1, [&, osc_cent_stn, mod, osc_freq, delay, fracup_stn](double t) {return squareWave(t, osc_cent_stn, mod, osc_freq, delay, fracup_stn);}, i, &stout);
		}

		else {throw std::invalid_argument( "osc_shape_stn must be \"sine\" or \"rect\"");}
		ext_stn.push_back(n);
	}
	pout<<std::endl;

	// Poisson neurons
	pout<<"# of Poisson STN Neurons: " << n_pois_stn << std::endl;
	for(int i = 0; i < n_pois_stn; i++){
		PoissonNeuron* n = new PoissonNeuron(0, 200000, 15, 1, n_osc_stn+i, &stout);
		ext_stn.push_back(n);
	}

	bool ext_gpe_spike[n_ext_gpe] = {0}; // whether external neuron spiked at current t
	bool ext_stn_spike[n_ext_stn] = {0};

	//// Main timestep loop

	// Jitter variables
	int num_steps = 0;
	int res1, res2;
	double applied_shift = 0;
	double new_shift; double new_freq = osc_freq + 0;
	ofstream frequencytracker;
	frequencytracker.open("frequency_tracker.txt");
	frequencytracker << "# " << "implement_jitter:" << "\t"<< implement_jitter << endl;
	frequencytracker << "# "<<"jitt1:" << "\t"<< jitt1 << endl;
	frequencytracker << "# "<< "jitt2:" << "\t"<< jitt2 << endl;
	frequencytracker << "# "<< "shift:" << "\t"<<shift << endl;
	frequencytracker << "# "<< "freq:" << "\t"<< osc_freq << endl;
	frequencytracker << "# "<< "seed:" << "\t"<< seed <<endl;
	// End jitter variables
	for(double t=0;t<=T;t+=dt)
	{
		double fcan= fc[1];
		double fw=fcw[1];
		double fsynCa=fsca[1]-t/T*(fsca[1]-fsca[0]);
		iapplied = 0;
		// int iapp_on = 1500; int iapp_off = 3500;
		// if (t>=iapp_on & t<iapp_off)
		// {
		// 	iapplied = iapp[0] + (t-iapp_on)/(iapp_off-iapp_on) * (iapp[1]-iapp[0]);
		// }

		//reinitialize variables
		for(int i=0;i<size;i++) {isi[i]=0; h[i]=0;}
		sp=0;
		int spk[size];

		/* JITTER CODE WILL GO HERE */
		if (implement_jitter){
			res1 = num_steps % jitt1;
			res2 = num_steps % jitt2;
			if (res2 == 0){
				new_shift = (osc_freq*(1*shift))*rnd();
				if (rand()%2){
					applied_shift = 1*new_shift;
				} else {
					applied_shift = -1*new_shift;
				}
			}
			if (res1 == 0){
				//osc_freq = osc_freq + applied_shift;
				new_freq = new_freq + applied_shift;
			}
		}
		/* END JITTER CODE */

		//step
		for(int i=0;i<n_ext_gpe;i++){
			if (implement_jitter){
				ext_gpe[i]->rate_func = [&, osc_cent_gpe, mod, new_freq, delay, fracup_gpe](double t) {return squareWave(t, osc_cent_gpe, mod, new_freq, delay, fracup_gpe);};
			}

			// the following IF statement records the new osc_freq for the GPe neuron
			if ((num_steps % 1000) && (i == 0)){
				frequencytracker << t << "\t" << new_freq << endl;
			}
			ext_gpe_spike[i] = ext_gpe.at(i)->step(t);
		}
		for(int i=0;i<n_ext_stn;i++){
			ext_stn_spike[i] = ext_stn.at(i)->step(t);
		}

		auto ss =pop.step(dt,drive,spk,fcan,fw,fsynCa,irr,iapplied, ext_gpe_spike, n_ext_gpe, ext_stn_spike, n_ext_stn);

		sp = sp+ss.sp;

		// write synaptic variables
		if (write_syn){
			int tind = round(t/dt);
			if (tind % 40 == 0) {
				extiout<<"\n"<<t/1000<<" ";
				intiout<<"\n"<<t/1000<<" ";
				gsynout<<"\n"<<t/1000<<" ";
				for(int s = 0; s<size; s++){
					extiout<<ss.ext_inh[s]<<" ";
					intiout<<ss.int_inh[s]<<" ";
					gsynout<<pop.net[s].gsyn<<" ";
				}
			}
		}
		// write calcium variables
		if (write_ca){
			int tind = round(t/dt);
			if (tind % 40 == 0) {
				caout<<"\n"<<t/1000<<" ";
				kskout<<"\n"<<t/1000<<" ";
				for(int s = 0; s<size; s++){
					caout<<pop.net[s].ICa<<" ";
					kskout<<pop.net[s].ISK<<" ";
				}
			}
		}

		for(int j=0;j<size;j++){
					mean_v+=(pop.net[j].v/size);
					mean_gaba_som+=(pop.net[j].EGABA/size);
					mean_gaba_dend+=(pop.net[j].EGABA_dend/size);
					}

		//calculate spike isi/freq and write each spike to file
		for(int i=0;i<size;i++)
		{
			h[i]+=spk[i];
			double PV = pop.net[0].pv_sp;
			double LHX6 = pop.net[0].lhx6_sp;
			double STN = pop.net[0].stn_sp;
			pv_spf = 1/0.0;
			lhx6_spf = 1/0.0;
			stn_spf = 1/0.0;
			if(h[i]>=1)
			{
				double cellnum = i;
				if(h[i]<1){cellnum = 1/0.0;}
				if(h[i]>=1){
				isi[i] = (t-lastSP[i])/1000;
				lastSP[i] = t;
				spf[i] =1/isi[i]; //get new spike frequency
				}
				if(PV==1){
				pv_isi = (t-pv_lastSP)/1000;
				pv_lastSP = t;
				pv_spf =1/pv_isi; //get new spike frequency
				}
				if(LHX6==1){
				lhx6_isi = (t-lhx6_lastSP)/1000;
				lhx6_lastSP = t;
				lhx6_spf =1/lhx6_isi; //get new spike frequency
				}
				if(STN==1){
				stn_isi = (t-stn_lastSP)/1000;
				stn_lastSP = t;
				stn_spf =1/stn_isi; //get new spike frequency
				}
				pop_sp +=1;
				last_net_sp = t;
				cout<<(t/1000)<<'\t'<<cellnum<<'\t'<<endl;
			}
		}

		if(int(t/dt)%(freq)==0)
		{
			out<<(t/1000);
			for(int i=0;i<size;i++) {
				out<<'\t'<<pop.net[i].v;
			}
			out<<endl;
		}
		mean_v=0;
		mean_gaba_som=0;
		mean_gaba_dend=0;
		num_steps += 1;
	}
	frequencytracker.close();
	return 0;
}
