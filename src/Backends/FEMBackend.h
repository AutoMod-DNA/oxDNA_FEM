/*
 * FEMBackend.h
 *
 *  Created on: 5 Aug 2026
 *      Author: elina
 */

#ifndef FEMBACKEND_H_
#define FEMBACKEND_H_

#include "MDBackend.h"
#include "../Eigen/Dense"
#include "../Eigen/Core"

class FEMBackend: public MDBackend
{
protected:
	double _n_iter_in;
	double _alpha;
	double _alpha_step;
	double _max_f;
	double _max_ff;
	double _r_max;
	double _mod_B;
	double _mod_C;

public:
	FEMBackend();
	virtual ~FEMBackend();

	// helper functions
	LR_vector get_descriptor_vector(BaseParticle *p, BaseParticle *q, int N);
	double get_descriptor_scalar(BaseParticle *p, BaseParticle *q, int N, int M);

	// 1st and 2nd order derivatives term 1: Backbone potential
	double get_D11_G1(BaseParticle *p, BaseParticle *q);
	double get_D21_G1(BaseParticle *p, BaseParticle *q);

	// helper functions
	double F1(double d, double RLOW, double RHIGH, double EPS, double A, double R0, double RC, double RCLOW, double RCHIGH, double BLOW, double BHIGH);
	double dF1(double d, double RLOW, double RHIGH, double EPS, double A, double R0, double RCLOW, double RCHIGH, double BLOW, double BHIGH);
	double ddF1(double d, double RLOW, double RHIGH, double EPS, double A, double R0, double RCLOW, double RCHIGH, double BLOW, double BHIGH);
	double F4(double T, double A, double B, double T0, double TS, double TC);
	double dF4(double T, double A, double B, double T0, double TS, double TC);
	double ddF4(double T, double A, double B, double T0, double TS, double TC);

	// 1st and 2nd order derivatives term 2: Hydrogen bonding
	double get_D12_G1(BaseParticle *p, BaseParticle *q);
	double get_D12_G2(BaseParticle *p, BaseParticle *q);
	double get_D12_G3(BaseParticle *p, BaseParticle *q);
	double get_D12_G4(BaseParticle *p, BaseParticle *q);
	double get_D12_G5(BaseParticle *p, BaseParticle *q);
	double get_D12_G6(BaseParticle *p, BaseParticle *q);
	double get_D12_G7(BaseParticle *p, BaseParticle *q);

	double get_D22_G1(BaseParticle *p, BaseParticle *q);
	double get_D22_G2(BaseParticle *p, BaseParticle *q);
	double get_D22_G3(BaseParticle *p, BaseParticle *q);
	double get_D22_G4(BaseParticle *p, BaseParticle *q);
	double get_D22_G5(BaseParticle *p, BaseParticle *q);
	double get_D22_G6(BaseParticle *p, BaseParticle *q);
	double get_D22_G7(BaseParticle *p, BaseParticle *q);

	void get_quaternion(BaseParticle *p, Eigen::Vector4d *q);

	// gradients of geometric descriptors term 1: Backbone potential
	void get_T1_G1(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G1, int which);

	// gradients of geometric descriptors term 2: Hydrogen bonding
	void get_T2_G1(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G1, int which);
	void get_T2_G2(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G1, int which);
	void get_T2_G3(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G1, int which);
	void get_T2_G4(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G1, int which);
	void get_T2_G5(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G1, int which);
	void get_T2_G6(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G1, int which);
	void get_T2_G7(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G1, int which);

	// construction of global stiffness matrix and global internal force vector
	void get_KF_local_bonded(BaseParticle *p, BaseParticle *q, Eigen::MatrixXd *K, Eigen::VectorXd *F);
	void get_KF_local_nonbonded(BaseParticle *p, BaseParticle *q, Eigen::MatrixXd *K, Eigen::VectorXd *F);
	void compute_KF(Eigen::MatrixXd *K_global, Eigen::VectorXd *F_global);


	void init();
	void get_settings(input_file &inp);
	void evolve(Eigen::VectorXd *dU);
	void sim_step();
};


#endif /* FEMBACKEND_H_ */
