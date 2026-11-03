/*
 * FEMBackend.cpp
 *
 *  Created on: 5 Aug 2026
 *      Author: elina
 */


#include "FEMBackend.h"

#include "../Utilities/Utils.h"
#include "../Particles/DNANucleotide.h"

#include "../model.h"
#include "../defs.h"

#include <sstream>
#include <cmath>
#include <iostream>


FEMBackend::FEMBackend()
: MDBackend()
{
	_alpha = 0.1;
	_alpha_step = 0.1;
	_n_iter_in = 10.0;
	_max_f = 0.5;
	_max_ff = 0.04;
	_r_max = 0;
	_mod_B = 0;
	_mod_C = 0;
}

FEMBackend::~FEMBackend()
{}

LR_vector FEMBackend::get_descriptor_vector(BaseParticle *p, BaseParticle *q, int N)
{
	LR_vector descriptor;
	switch ( N )
	{
	case 1:		// BACK - BACK
		descriptor = (LR_vector) ((p->pos + p->int_centers[DNANucleotide::BACK]) - (q->pos + q->int_centers[DNANucleotide::BACK]));
		break;
	case 2:		// BASE - BASE
		descriptor = (LR_vector) ((p->pos + p->int_centers[DNANucleotide::BASE]) - (q->pos + q->int_centers[DNANucleotide::BASE]));
		break;
	case 3:		// STACK - STACK
		descriptor = (LR_vector) ((p->pos + p->int_centers[DNANucleotide::STACK]) - (q->pos + q->int_centers[DNANucleotide::STACK]));
		break;
	case 4:		// BACK - BASE
		descriptor = (LR_vector) ((p->pos + p->int_centers[DNANucleotide::BACK]) - (q->pos + q->int_centers[DNANucleotide::BASE]));
		break;
	case 5:		// BASE - BACK
		descriptor = (LR_vector) ((p->pos + p->int_centers[DNANucleotide::BASE]) - (q->pos + q->int_centers[DNANucleotide::BACK]));
		break;
	// default:
		// error message here
	}
	return descriptor;
}

double FEMBackend::get_descriptor_scalar(BaseParticle *p, BaseParticle *q, int N, int M)
{
	double descriptor;
	LR_vector dir_norm;

	LR_vector &a1 = p->orientationT.v1;
	LR_vector &a2 = p->orientationT.v2;
	LR_vector &a3 = p->orientationT.v3;
	LR_vector &b1 = q->orientationT.v1;
	LR_vector &b2 = q->orientationT.v2;
	LR_vector &b3 = q->orientationT.v3;

	switch ( M )
	{
	case 1:
		// BASE - BASE
		dir_norm = get_descriptor_vector(p, q, 2);
		dir_norm = dir_norm / dir_norm.module();
		break;
	// case 2:
		// other interactions here as needed
	}

	switch ( N )
	{
	case 6:		// THETA 1
		descriptor = std::acos(-a1 * b1);
		break;
	case 7:		// THETA 2
		descriptor = std::acos(-b1 * dir_norm);
		break;
	case 8:		// THETA 3
		descriptor = std::acos(a1 * dir_norm);
		break;
	case 9:		// THETA 4
		descriptor = std::acos(a3 * b3);
		break;
	case 10:	// THETA 5
		descriptor = std::acos(a3 * dir_norm);
		break;
	case 11:	// THETA 6
		descriptor = std::acos(-b3 * dir_norm);
		break;
	case 12:	// THETA 7
		descriptor = std::acos(-b3 * dir_norm);
		break;
	case 13:	// THETA 8
		descriptor = std::acos(a3 * dir_norm);
		break;

	case 14:	// THETA 5'
		break;
	case 15:	// THETA 6'
		break;

	case 16:	// COS PHI 1
		break;
	case 17:	// COS PHI 2
		break;
	case 18:	// COS PHI 3
		break;
	case 19:	// COS PHI 4
		break;

	// default:
		// error message here
	}
	return descriptor;
}

double FEMBackend::get_D11_G1(BaseParticle *p, BaseParticle *q)
{
	// D11: DERIVATIVE 1 TERM 1: Backbone potential once differentiated with respect to
	// G1: GEOMETRIC DESCRIPTOR 1: Vector between interaction sites
	LR_vector r_back = get_descriptor_vector(p, q, 1);
	double D11, d_back(r_back.module());
	if (d_back > _r_max && _r_max != 0) {
		D11 = _max_f + _mod_B / d_back;
	}
	else {
		D11 = FENE_EPS*(d_back - FENE_R0_OXDNA2) / (std::pow(FENE_DELTA2, 2) - std::pow(d_back - FENE_R0_OXDNA2, 2));
	}
	return D11;
}

double FEMBackend::get_D21_G1(BaseParticle *p, BaseParticle *q)
{
	// D21: DERIVATIVE 2 TERM 1: Backbone potential twice differentiated with respect to
	// G1: GEOMETRIC DESCRIPTOR 1: Vector between interaction sites
	LR_vector r_back = get_descriptor_vector(p, q, 1);
	double D21, d_back(r_back.module());
	if (d_back > _r_max && _r_max != 0) {
		D21 = (-1)*_mod_B / std::pow(d_back, 2);
	}
	else {
		D21 = FENE_EPS*(std::pow(FENE_DELTA2, 2) + std::pow(d_back - FENE_R0_OXDNA2, 2)) / (std::pow((std::pow(FENE_DELTA2, 2) - std::pow(d_back - FENE_R0_OXDNA2, 2)), 2));
	}
	return D21;
}

double FEMBackend::F1(double d, double RLOW, double RHIGH, double EPS, double A, double R0, double RC, double RCLOW, double RCHIGH, double BLOW, double BHIGH)
{
	double F;
	if (d > RLOW && d < RHIGH)
	{
		F = EPS * std::pow(1 - std::exp((-1)*(d - R0)*A), 2) - EPS * std::pow(1 - std::exp((-1)*(RC - R0)*A), 2);
	}
	else if (d > RCLOW && d < RLOW)
	{
		F = EPS * BLOW * std::pow((RCLOW - d), 2);
	}
	else if (d > RHIGH && d < RCHIGH)
	{
		F = EPS * BHIGH * std::pow((RCHIGH), 2);
	}
	else
	{
		F = 0;
	}
	return F;
}

double FEMBackend::dF1(double d, double RLOW, double RHIGH, double EPS, double A, double R0, double RCLOW, double RCHIGH, double BLOW, double BHIGH)
{
	double D;
	if (d > RLOW && d < RHIGH)
	{
		D = 2 * EPS * A * std::exp((-1)*(d - R0)*A) * (1 - std::exp((-1)*(d - R0)*A));
	}
	else if (d > RCLOW && d < RLOW)
	{
		D = 2 * EPS * BLOW * (d - RCLOW);
	}
	else if (d > RHIGH && d < RCHIGH)
	{
		D = 2 * EPS * BHIGH * (d - RCHIGH);
	}
	else
	{
		D = 0;
	}
	return D;
}

double FEMBackend::ddF1(double d, double RLOW, double RHIGH, double EPS, double A, double R0, double RCLOW, double RCHIGH, double BLOW, double BHIGH)
{
	double D;
	if (d > RLOW && d < RHIGH)
	{
		D = (-1)*2*EPS*A*A * std::exp((-1)*(d - R0)*A) * (1 - 2*std::exp((-1)*(d - R0)*A));
	}
	else if (d > RCLOW && d < RLOW)
	{
		D = 2*EPS*BLOW;
	}
	else if (d > RHIGH && d < RCHIGH)
	{
		D = 2*EPS*BHIGH;
	}
	else
	{
		D = 0;
	}
	return D;
}

double FEMBackend::F4(double T, double A, double B, double T0, double TS, double TC)
{
	double F4;
	if (T > T0 - TS && T < T0 + TS)
	{
		F4 = 1 - A*std::pow(T - T0, 2);
	}
	else if (T > T0 - TC && T < T0 - TS)
	{
		F4 = B*std::pow((T0 - TC) - T, 2);
	}
	else if (T > T0 + TS && T < T0 + TC)
	{
		F4 = B*std::pow((T0 + TC) - T, 2);
	}
	else
	{
		F4 = 0;
	}
	return F4;
}

double FEMBackend::dF4(double T, double A, double B, double T0, double TS, double TC)
{
	double F4;
	if (T > T0 - TS && T < T0 + TS)
	{
		F4 = (-1)*2*A*(T - T0);
	}
	else if (T > T0 - TC && T < T0 - TS)
	{
		F4 = 2*B*(T - (T0 - TC));
	}
	else if (T > T0 + TS && T < T0 + TC)
	{
		F4 = 2*B*(T - (T0 + TC));
	}
	else
	{
		F4 = 0;
	}
	return F4;
}

double FEMBackend::ddF4(double T, double A, double B, double T0, double TS, double TC)
{
	double F4;
	if (T > T0 - TS && T < T0 + TS)
	{
		F4 = (-1)*2*A;
	}
	else if (T > T0 - TC && T < T0 - TS)
	{
		F4 = 2*B;
	}
	else if (T > T0 + TS && T < T0 + TC)
	{
		F4 = 2*B;
	}
	else
	{
		F4 = 0;
	}
	return F4;
}

double FEMBackend::get_D12_G1(BaseParticle *p, BaseParticle *q)
{
	// D22: DERIVATIVE 1 TERM 2: Hydrogen bonding potential once differentiated with respect to
	// G1: GEOMETRIC DESCRIPTOR 1: Vector between interaction sites
	double D12;
	LR_vector r_hb = get_descriptor_vector(p, q, 2);

	D12 = dF1(r_hb.module(), HYDR_RLOW, HYDR_RHIGH, HYDR_EPS_OXDNA2, HYDR_A, HYDR_R0, HYDR_RCLOW, HYDR_RCHIGH, HYDR_BLOW, HYDR_BHIGH);

	D12 *= F4(get_descriptor_scalar(p, q, 6, 1), HYDR_THETA1_A, HYDR_THETA1_B, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 7, 1), HYDR_THETA2_A, HYDR_THETA2_B, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 8, 1), HYDR_THETA3_A, HYDR_THETA3_B, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 9, 1), HYDR_THETA4_A, HYDR_THETA4_B, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 12, 1), HYDR_THETA7_A, HYDR_THETA7_B, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 13, 1), HYDR_THETA8_A, HYDR_THETA8_B, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC);

	return D12;
}

double FEMBackend::get_D12_G2(BaseParticle *p, BaseParticle *q)
{
	// D22: DERIVATIVE 1 TERM 2: Hydrogen bonding potential once differentiated with respect to
	// G2: GEOMETRIC DESCRIPTOR 2: Angle theta 1
	double D12;
	LR_vector r_hb = get_descriptor_vector(p, q, 2);

	D12 = dF4(get_descriptor_scalar(p, q, 6, 1), HYDR_THETA1_A, HYDR_THETA1_B, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC);

	D12 *= F1(r_hb.module(), HYDR_RLOW, HYDR_RHIGH, HYDR_EPS_OXDNA2, HYDR_A, HYDR_R0, HYDR_RC, HYDR_RCLOW, HYDR_RCHIGH, HYDR_BLOW, HYDR_BHIGH);

	D12 *= F4(get_descriptor_scalar(p, q, 7, 1), HYDR_THETA2_A, HYDR_THETA2_B, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 8, 1), HYDR_THETA3_A, HYDR_THETA3_B, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 9, 1), HYDR_THETA4_A, HYDR_THETA4_B, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 12, 1), HYDR_THETA7_A, HYDR_THETA7_B, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 13, 1), HYDR_THETA8_A, HYDR_THETA8_B, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC);

	return D12;
}

double FEMBackend::get_D12_G3(BaseParticle *p, BaseParticle *q)
{
	// D22: DERIVATIVE 1 TERM 2: Hydrogen bonding potential once differentiated with respect to
	// G3: GEOMETRIC DESCRIPTOR 3: Angle theta 2
	double D12;
	LR_vector r_hb = get_descriptor_vector(p, q, 2);

	D12 = dF4(get_descriptor_scalar(p, q, 7, 1), HYDR_THETA2_A, HYDR_THETA2_B, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC);

	D12 *= F1(r_hb.module(), HYDR_RLOW, HYDR_RHIGH, HYDR_EPS_OXDNA2, HYDR_A, HYDR_R0, HYDR_RC, HYDR_RCLOW, HYDR_RCHIGH, HYDR_BLOW, HYDR_BHIGH);

	D12 *= F4(get_descriptor_scalar(p, q, 6, 1), HYDR_THETA1_A, HYDR_THETA1_B, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 8, 1), HYDR_THETA3_A, HYDR_THETA3_B, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 9, 1), HYDR_THETA4_A, HYDR_THETA4_B, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 12, 1), HYDR_THETA7_A, HYDR_THETA7_B, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 13, 1), HYDR_THETA8_A, HYDR_THETA8_B, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC);

	return D12;
}

double FEMBackend::get_D12_G4(BaseParticle *p, BaseParticle *q)
{
	// D22: DERIVATIVE 1 TERM 2: Hydrogen bonding potential once differentiated with respect to
	// G4: GEOMETRIC DESCRIPTOR 4: Angle theta 3
	double D12;
	LR_vector r_hb = get_descriptor_vector(p, q, 2);

	D12 = dF4(get_descriptor_scalar(p, q, 8, 1), HYDR_THETA3_A, HYDR_THETA3_B, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC);

	D12 *= F1(r_hb.module(), HYDR_RLOW, HYDR_RHIGH, HYDR_EPS_OXDNA2, HYDR_A, HYDR_R0, HYDR_RC, HYDR_RCLOW, HYDR_RCHIGH, HYDR_BLOW, HYDR_BHIGH);

	D12 *= F4(get_descriptor_scalar(p, q, 6, 1), HYDR_THETA1_A, HYDR_THETA1_B, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 7, 1), HYDR_THETA2_A, HYDR_THETA2_B, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 9, 1), HYDR_THETA4_A, HYDR_THETA4_B, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 12, 1), HYDR_THETA7_A, HYDR_THETA7_B, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 13, 1), HYDR_THETA8_A, HYDR_THETA8_B, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC);

	return D12;
}

double FEMBackend::get_D12_G5(BaseParticle *p, BaseParticle *q)
{
	// D22: DERIVATIVE 1 TERM 2: Hydrogen bonding potential once differentiated with respect to
	// G5: GEOMETRIC DESCRIPTOR 5: Angle theta 4
	double D12;
	LR_vector r_hb = get_descriptor_vector(p, q, 2);

	D12 = dF4(get_descriptor_scalar(p, q, 9, 1), HYDR_THETA4_A, HYDR_THETA4_B, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC);

	D12 *= F1(r_hb.module(), HYDR_RLOW, HYDR_RHIGH, HYDR_EPS_OXDNA2, HYDR_A, HYDR_R0, HYDR_RC, HYDR_RCLOW, HYDR_RCHIGH, HYDR_BLOW, HYDR_BHIGH);

	D12 *= F4(get_descriptor_scalar(p, q, 6, 1), HYDR_THETA1_A, HYDR_THETA1_B, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 7, 1), HYDR_THETA2_A, HYDR_THETA2_B, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 8, 1), HYDR_THETA3_A, HYDR_THETA3_B, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 12, 1), HYDR_THETA7_A, HYDR_THETA7_B, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 13, 1), HYDR_THETA8_A, HYDR_THETA8_B, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC);

	return D12;
}

double FEMBackend::get_D12_G6(BaseParticle *p, BaseParticle *q)
{
	// D22: DERIVATIVE 1 TERM 2: Hydrogen bonding potential once differentiated with respect to
	// G6: GEOMETRIC DESCRIPTOR 6: Angle theta 7
	double D12;
	LR_vector r_hb = get_descriptor_vector(p, q, 2);

	D12 = dF4(get_descriptor_scalar(p, q, 12, 1), HYDR_THETA7_A, HYDR_THETA7_B, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC);

	D12 *= F1(r_hb.module(), HYDR_RLOW, HYDR_RHIGH, HYDR_EPS_OXDNA2, HYDR_A, HYDR_R0, HYDR_RC, HYDR_RCLOW, HYDR_RCHIGH, HYDR_BLOW, HYDR_BHIGH);

	D12 *= F4(get_descriptor_scalar(p, q, 6, 1), HYDR_THETA1_A, HYDR_THETA1_B, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 7, 1), HYDR_THETA2_A, HYDR_THETA2_B, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 8, 1), HYDR_THETA3_A, HYDR_THETA3_B, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 9, 1), HYDR_THETA4_A, HYDR_THETA4_B, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 13, 1), HYDR_THETA8_A, HYDR_THETA8_B, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC);

	return D12;
}

double FEMBackend::get_D12_G7(BaseParticle *p, BaseParticle *q)
{
	// D22: DERIVATIVE 1 TERM 2: Hydrogen bonding potential once differentiated with respect to
	// G7: GEOMETRIC DESCRIPTOR 7: Angle theta 8
	double D12;
	LR_vector r_hb = get_descriptor_vector(p, q, 2);

	D12 = dF4(get_descriptor_scalar(p, q, 13, 1), HYDR_THETA8_A, HYDR_THETA8_B, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC);

	D12 *= F1(r_hb.module(), HYDR_RLOW, HYDR_RHIGH, HYDR_EPS_OXDNA2, HYDR_A, HYDR_R0, HYDR_RC, HYDR_RCLOW, HYDR_RCHIGH, HYDR_BLOW, HYDR_BHIGH);

	D12 *= F4(get_descriptor_scalar(p, q, 6, 1), HYDR_THETA1_A, HYDR_THETA1_B, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 7, 1), HYDR_THETA2_A, HYDR_THETA2_B, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 8, 1), HYDR_THETA3_A, HYDR_THETA3_B, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 9, 1), HYDR_THETA4_A, HYDR_THETA4_B, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC);
	D12 *= F4(get_descriptor_scalar(p, q, 12, 1), HYDR_THETA7_A, HYDR_THETA7_B, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC);

	return D12;
}

double FEMBackend::get_D22_G1(BaseParticle *p, BaseParticle *q)
{
	// D22: DERIVATIVE 2 TERM 2: Hydrogen bonding potential twice differentiated with respect to
	// G1: GEOMETRIC DESCRIPTOR 1: Vector between interaction sites
	double D22;
	LR_vector r_hb = get_descriptor_vector(p, q, 2);

	D22 = ddF1(r_hb.module(), HYDR_RLOW, HYDR_RHIGH, HYDR_EPS_OXDNA2, HYDR_A, HYDR_R0, HYDR_RCLOW, HYDR_RCHIGH, HYDR_BLOW, HYDR_BHIGH);

	D22 *= F4(get_descriptor_scalar(p, q, 6, 1), HYDR_THETA1_A, HYDR_THETA1_B, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 7, 1), HYDR_THETA2_A, HYDR_THETA2_B, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 8, 1), HYDR_THETA3_A, HYDR_THETA3_B, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 9, 1), HYDR_THETA4_A, HYDR_THETA4_B, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 12, 1), HYDR_THETA7_A, HYDR_THETA7_B, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 13, 1), HYDR_THETA8_A, HYDR_THETA8_B, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC);

	return D22;
}

double FEMBackend::get_D22_G2(BaseParticle *p, BaseParticle *q)
{
	// D22: DERIVATIVE 2 TERM 2: Hydrogen bonding potential twice differentiated with respect to
	// G2: GEOMETRIC DESCRIPTOR 2: Angle theta 1
	double D22;
	LR_vector r_hb = get_descriptor_vector(p, q, 2);

	D22 = ddF4(get_descriptor_scalar(p, q, 6, 1), HYDR_THETA1_A, HYDR_THETA1_B, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC);

	D22 *= F1(r_hb.module(), HYDR_RLOW, HYDR_RHIGH, HYDR_EPS_OXDNA2, HYDR_A, HYDR_R0, HYDR_RC, HYDR_RCLOW, HYDR_RCHIGH, HYDR_BLOW, HYDR_BHIGH);

	D22 *= F4(get_descriptor_scalar(p, q, 7, 1), HYDR_THETA2_A, HYDR_THETA2_B, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 8, 1), HYDR_THETA3_A, HYDR_THETA3_B, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 9, 1), HYDR_THETA4_A, HYDR_THETA4_B, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 12, 1), HYDR_THETA7_A, HYDR_THETA7_B, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 13, 1), HYDR_THETA8_A, HYDR_THETA8_B, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC);

	return D22;
}

double FEMBackend::get_D22_G3(BaseParticle *p, BaseParticle *q)
{
	// D22: DERIVATIVE 2 TERM 2: Hydrogen bonding potential twice differentiated with respect to
	// G3: GEOMETRIC DESCRIPTOR 3: Angle theta 2
	double D22;
	LR_vector r_hb = get_descriptor_vector(p, q, 2);

	D22 = ddF4(get_descriptor_scalar(p, q, 7, 1), HYDR_THETA2_A, HYDR_THETA2_B, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC);

	D22 *= F1(r_hb.module(), HYDR_RLOW, HYDR_RHIGH, HYDR_EPS_OXDNA2, HYDR_A, HYDR_R0, HYDR_RC, HYDR_RCLOW, HYDR_RCHIGH, HYDR_BLOW, HYDR_BHIGH);

	D22 *= F4(get_descriptor_scalar(p, q, 6, 1), HYDR_THETA1_A, HYDR_THETA1_B, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 8, 1), HYDR_THETA3_A, HYDR_THETA3_B, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 9, 1), HYDR_THETA4_A, HYDR_THETA4_B, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 12, 1), HYDR_THETA7_A, HYDR_THETA7_B, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 13, 1), HYDR_THETA8_A, HYDR_THETA8_B, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC);

	return D22;
}

double FEMBackend::get_D22_G4(BaseParticle *p, BaseParticle *q)
{
	// D22: DERIVATIVE 2 TERM 2: Hydrogen bonding potential twice differentiated with respect to
	// G4: GEOMETRIC DESCRIPTOR 4: Angle theta 3
	double D22;
	LR_vector r_hb = get_descriptor_vector(p, q, 2);

	D22 = ddF4(get_descriptor_scalar(p, q, 8, 1), HYDR_THETA3_A, HYDR_THETA3_B, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC);

	D22 *= F1(r_hb.module(), HYDR_RLOW, HYDR_RHIGH, HYDR_EPS_OXDNA2, HYDR_A, HYDR_R0, HYDR_RC, HYDR_RCLOW, HYDR_RCHIGH, HYDR_BLOW, HYDR_BHIGH);

	D22 *= F4(get_descriptor_scalar(p, q, 6, 1), HYDR_THETA1_A, HYDR_THETA1_B, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 7, 1), HYDR_THETA2_A, HYDR_THETA2_B, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 9, 1), HYDR_THETA4_A, HYDR_THETA4_B, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 12, 1), HYDR_THETA7_A, HYDR_THETA7_B, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 13, 1), HYDR_THETA8_A, HYDR_THETA8_B, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC);

	return D22;
}

double FEMBackend::get_D22_G5(BaseParticle *p, BaseParticle *q)
{
	// D22: DERIVATIVE 2 TERM 2: Hydrogen bonding potential twice differentiated with respect to
	// G5: GEOMETRIC DESCRIPTOR 5: Angle theta 4
	double D22;
	LR_vector r_hb = get_descriptor_vector(p, q, 2);

	D22 = ddF4(get_descriptor_scalar(p, q, 9, 1), HYDR_THETA4_A, HYDR_THETA4_B, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC);

	D22 *= F1(r_hb.module(), HYDR_RLOW, HYDR_RHIGH, HYDR_EPS_OXDNA2, HYDR_A, HYDR_R0, HYDR_RC, HYDR_RCLOW, HYDR_RCHIGH, HYDR_BLOW, HYDR_BHIGH);

	D22 *= F4(get_descriptor_scalar(p, q, 6, 1), HYDR_THETA1_A, HYDR_THETA1_B, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 7, 1), HYDR_THETA2_A, HYDR_THETA2_B, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 8, 1), HYDR_THETA3_A, HYDR_THETA3_B, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 12, 1), HYDR_THETA7_A, HYDR_THETA7_B, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 13, 1), HYDR_THETA8_A, HYDR_THETA8_B, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC);

	return D22;
}

double FEMBackend::get_D22_G6(BaseParticle *p, BaseParticle *q)
{
	// D22: DERIVATIVE 2 TERM 2: Hydrogen bonding potential twice differentiated with respect to
	// G6: GEOMETRIC DESCRIPTOR 6: Angle theta 7
	double D22;
	LR_vector r_hb = get_descriptor_vector(p, q, 2);

	D22 = ddF4(get_descriptor_scalar(p, q, 12, 1), HYDR_THETA7_A, HYDR_THETA7_B, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC);

	D22 *= F1(r_hb.module(), HYDR_RLOW, HYDR_RHIGH, HYDR_EPS_OXDNA2, HYDR_A, HYDR_R0, HYDR_RC, HYDR_RCLOW, HYDR_RCHIGH, HYDR_BLOW, HYDR_BHIGH);

	D22 *= F4(get_descriptor_scalar(p, q, 6, 1), HYDR_THETA1_A, HYDR_THETA1_B, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 7, 1), HYDR_THETA2_A, HYDR_THETA2_B, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 8, 1), HYDR_THETA3_A, HYDR_THETA3_B, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 9, 1), HYDR_THETA4_A, HYDR_THETA4_B, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 13, 1), HYDR_THETA8_A, HYDR_THETA8_B, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC);

	return D22;
}

double FEMBackend::get_D22_G7(BaseParticle *p, BaseParticle *q)
{
	// D22: DERIVATIVE 2 TERM 2: Hydrogen bonding potential twice differentiated with respect to
	// G7: GEOMETRIC DESCRIPTOR 7: Angle theta 8
	double D22;
	LR_vector r_hb = get_descriptor_vector(p, q, 2);

	D22 = ddF4(get_descriptor_scalar(p, q, 13, 1), HYDR_THETA8_A, HYDR_THETA8_B, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC);

	D22 *= F1(r_hb.module(), HYDR_RLOW, HYDR_RHIGH, HYDR_EPS_OXDNA2, HYDR_A, HYDR_R0, HYDR_RC, HYDR_RCLOW, HYDR_RCHIGH, HYDR_BLOW, HYDR_BHIGH);

	D22 *= F4(get_descriptor_scalar(p, q, 6, 1), HYDR_THETA1_A, HYDR_THETA1_B, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 7, 1), HYDR_THETA2_A, HYDR_THETA2_B, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 8, 1), HYDR_THETA3_A, HYDR_THETA3_B, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 9, 1), HYDR_THETA4_A, HYDR_THETA4_B, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC);
	D22 *= F4(get_descriptor_scalar(p, q, 12, 1), HYDR_THETA7_A, HYDR_THETA7_B, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC);

	return D22;
}

void FEMBackend::get_quaternion(BaseParticle *p, Eigen::Vector4d *q)
{
	LR_matrix Q = p->orientation;
	double trace = Q.v1.x + Q.v2.y + Q.v3.z;
	if (trace > 0)
	{
		double r = std::sqrt(1 + trace);
		(*q)(0) = r/2;
		(*q)(1) = (Q.v3.y - Q.v2.z)/(2*r);
		(*q)(2) = (Q.v1.z - Q.v3.x)/(2*r);
		(*q)(3) = (Q.v2.x - Q.v1.y)/(2*r);
	}
	else if ((Q.v1.x >= Q.v2.y) && (Q.v1.x >= Q.v3.z))
	{
		double r = std::sqrt(1 + Q.v1.x - Q.v2.y - Q.v3.z);
		(*q)(0) = (Q.v3.y - Q.v2.z)/(2*r);
		(*q)(1) = r/2;
		(*q)(2) = (Q.v1.y + Q.v2.x)/(2*r);
		(*q)(3) = (Q.v3.x + Q.v1.z)/(2*r);
	}
	else if ((Q.v2.y >= Q.v1.x) && (Q.v2.y >= Q.v3.z))
	{
		double r = std::sqrt(1 - Q.v1.x + Q.v2.y - Q.v3.z);
		(*q)(0) = (Q.v3.x - Q.v1.z)/(2*r);
		(*q)(1) = (Q.v2.x + Q.v1.y)/(2*r);
		(*q)(2) = r/2;
		(*q)(3) = (Q.v3.x + Q.v1.z)/(2*r);
	}
	else
	{
		double r = std::sqrt(1 - Q.v1.x - Q.v2.y + Q.v3.z);
		(*q)(0) = (Q.v2.x - Q.v1.y)/(2*r);
		(*q)(1) = (Q.v1.z + Q.v3.x)/(2*r);
		(*q)(2) = (Q.v3.y + Q.v2.z)/(2*r);
		(*q)(3) = r/2;
	}
}

void FEMBackend::get_T1_G1(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G, int which)
{
	// T1: TERM 1: Backbone potential
	// G1: GRADIENT OF GEOMETRIC DESCRIPTOR 1: Vector between interaction sites
	// particle p experiences the force/torque (vectors go to p from q)
	// differentiation is done with respect to the degrees of freedom (position and orientation) of the particle indicated by which
	// which: 1 = p, 2 = q
	LR_vector r_back = get_descriptor_vector(p, q, 1);
	Eigen::Vector4d quat;
	Eigen::template Matrix<double, 3, 4, 0, 3, 4> Q;
	double dx, sign;

	r_back.normalize();
	Eigen::Vector3d r;
	r << r_back.x, r_back.y, r_back.z;

	if (which == 1)
	{
		sign = +1;
		get_quaternion(p, &quat);
		dx = ((LR_vector) (p->int_centers[DNANucleotide::BACK])).module();
	}
	else
	{
		sign = -1;
		get_quaternion(q, &quat);
		dx = ((LR_vector) (q->int_centers[DNANucleotide::BACK])).module();
		sign = -1;
	}
	Q << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
	Q = 2*dx*Q;

	(*G).block(0, 0, 3, 1) = sign * r;
	(*G).block(3, 0, 4, 1) = sign * (r.transpose() * Q).transpose();
}

void FEMBackend::get_T2_G1(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G, int which)
{
	// T2: TERM 2: Hydrogen bonding
	// G1: GRADIENT OF GEOMETRIC DESCRIPTOR 1: Vector between interaction sites
	LR_vector r_hb = get_descriptor_vector(p, q, 2);
	Eigen::Vector4d quat;
	Eigen::template Matrix<double, 3, 4, 0, 3, 4> Q;

	r_hb.normalize();
	Eigen::Vector3d r;
	r << r_hb.x, r_hb.y, r_hb.z;

	double sign, dx;
	if (which == 1)
	{
		sign = +1;
		get_quaternion(p, &quat);
		dx = ((LR_vector) (p->int_centers[DNANucleotide::BASE])).module();
	}
	else
	{
		sign = -1;
		get_quaternion(q, &quat);
		dx = ((LR_vector) (q->int_centers[DNANucleotide::BASE])).module();
	}
	Q << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
	Q = 2*dx*Q;

	(*G).block(0, 0, 3, 1) = sign * r;
	(*G).block(3, 0, 4, 1) = sign * (r.transpose() * Q).transpose();
}

void FEMBackend::get_T2_G2(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G, int which)
{
	// T2: TERM 2: Hydrogen bonding
	// G2: GRADIENT OF GEOMETRIC DESCRIPTOR 2: Angle theta 1
	Eigen::Vector4d quat;
	Eigen::template Matrix<double, 3, 4, 0, 3, 4> Q;

	Eigen::Vector3d b_bar, b_tilde;
	b_bar << p->orientationT.v1.x, p->orientationT.v1.y, p->orientationT.v1.z;
	b_tilde << q->orientationT.v1.x, q->orientationT.v1.y, q->orientationT.v1.z;
	Eigen::Matrix3d p_AT, q_AT;
	p_AT << p->orientationT.v1.x, p->orientationT.v1.y, p->orientationT.v1.z;
	p_AT << p->orientationT.v2.x, p->orientationT.v2.y, p->orientationT.v2.z;
	p_AT << p->orientationT.v3.x, p->orientationT.v3.y, p->orientationT.v3.z;

	q_AT << q->orientationT.v1.x, q->orientationT.v1.y, q->orientationT.v1.z;
	q_AT << q->orientationT.v2.x, q->orientationT.v2.y, q->orientationT.v2.z;
	q_AT << q->orientationT.v3.x, q->orientationT.v3.y, q->orientationT.v3.z;

	Eigen::Vector3d b_bar_s(p_AT * b_bar);
	Eigen::Vector3d b_tilde_s(q_AT * b_tilde);

	double sign;
	if (which == 1)
	{
		sign = +1;
		get_quaternion(p, &quat);
		Q << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
		Q = 2*Q;
		(*G).block(3, 0, 4, 1) =  sign * ((-1) / std::sqrt(1 - std::pow(b_bar_s.transpose() * b_tilde_s, 2))) * (b_tilde.transpose() * Q).transpose();
	}
	else
	{
		sign = -1;
		get_quaternion(q, &quat);
		Q << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
		Q = 2*Q;
		(*G).block(3, 0, 4, 1) =  sign * ((-1) / std::sqrt(1 - std::pow(b_tilde_s.transpose() * b_bar_s, 2))) * (b_bar.transpose() * Q).transpose();
	}

	(*G).block(0, 0, 3, 1) << 0, 0, 0;
}

void FEMBackend::get_T2_G3(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G, int which)
{
	// T2: TERM 2: Hydrogen bonding
	// G3: GRADIENT OF GEOMETRIC DESCRIPTOR 3: Angle theta 2
	LR_vector r_hb_lr = get_descriptor_vector(p, q, 2);

	Eigen::Vector3d r_hb, r_hb_norm, b_bar, b_tilde;

	r_hb << r_hb_lr.x, r_hb_lr.y, r_hb_lr.z;
	r_hb_norm = r_hb / r_hb_lr.module();
	b_bar << p->orientationT.v1.x, p->orientationT.v1.y, p->orientationT.v1.z;
	b_tilde << q->orientationT.v1.x, q->orientationT.v1.y, q->orientationT.v1.z;

	Eigen::Vector4d quat;
	Eigen::template Matrix<double, 3, 4, 0, 3, 4> Q;

	double sign, d(r_hb_lr.module()), dx;
	if (which == 1)
	{
		sign = +1;
		get_quaternion(p, &quat);
		Q << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
		dx = ((LR_vector) (p->int_centers[DNANucleotide::BASE])).module();

		(*G).block(0, 0, 3, 1) = sign * (1/std::sqrt(1 - std::pow(b_bar.transpose()*r_hb_norm, 2))) * (1/d) * (b_bar - (b_bar.transpose()*r_hb_norm)*r_hb_norm);
		(*G).block(3, 0, 4, 1) = 2*(r_hb_norm.transpose() * Q).transpose() + 2*dx*(1/d)*(b_bar.transpose() * Q).transpose();
		(*G).block(3, 0, 4, 1) -= (1/std::pow(d, 3)) * (double)(b_bar.transpose()*r_hb).value() * 2*dx*(r_hb.transpose() * Q).transpose();
		(*G).block(3, 0, 4, 1) *= sign / std::sqrt(1 - std::pow(b_bar.transpose()*r_hb_norm, 2));
	}
	else
	{
		sign = -1;
		get_quaternion(q, &quat);
		Q << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
		dx = ((LR_vector) (q->int_centers[DNANucleotide::BASE])).module();

		(*G).block(0, 0, 3, 1) = sign * (-1/std::sqrt(1 - std::pow(b_tilde.transpose()*r_hb_norm, 2))) * (1/d) * (b_tilde - (b_tilde.transpose()*r_hb_norm)*r_hb_norm);
		(*G).block(3, 0, 4, 1) = 2*dx*(1/d)*(b_tilde.transpose() * Q).transpose();
		(*G).block(3, 0, 4, 1) -= (1/std::pow(d, 3)) * (double)(b_tilde.transpose()*r_hb).value() * 2*dx*(r_hb.transpose() * Q).transpose();
		(*G).block(3, 0, 4, 1) *= sign*(-1) / std::sqrt(1 - std::pow(b_tilde.transpose()*r_hb_norm, 2));
	}
	// which 2 is the same as which 1 for theta 3 but with reversed signs
}

void FEMBackend::get_T2_G4(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G, int which)
{
	// T2: TERM 2: Hydrogen bonding
	// G4: GRADIENT OF GEOMETRIC DESCRIPTOR 4: Angle theta 3
	LR_vector r_hb_lr = get_descriptor_vector(p, q, 2);

	Eigen::Vector3d r_hb, r_hb_norm, b_bar, b_tilde;

	r_hb << r_hb_lr.x, r_hb_lr.y, r_hb_lr.z;
	r_hb_norm = r_hb / r_hb_lr.module();
	b_bar << p->orientationT.v1.x, p->orientationT.v1.y, p->orientationT.v1.z;
	b_tilde << q->orientationT.v1.x, q->orientationT.v1.y, q->orientationT.v1.z;

	Eigen::Vector4d quat;
	Eigen::template Matrix<double, 3, 4, 0, 3, 4> Q;

	double sign, d(r_hb_lr.module()), dx;
	if (which == 1)
	{
		sign = +1;
		get_quaternion(p, &quat);
		Q << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
		dx = ((LR_vector) (p->int_centers[DNANucleotide::BASE])).module();

		(*G).block(0, 0, 3, 1) = sign * ((-1)/std::sqrt(1 - std::pow(b_tilde.transpose()*r_hb_norm, 2))) * (1/d) * (b_tilde - (b_tilde.transpose()*r_hb_norm)*r_hb_norm);
		(*G).block(3, 0, 4, 1) = 2*dx*(1/d)*(b_tilde.transpose() * Q).transpose();
		(*G).block(3, 0, 4, 1) -= (1/std::pow(d, 3)) * (double)(b_tilde.transpose()*r_hb).value() * 2*dx*(r_hb.transpose() * Q).transpose();
		(*G).block(3, 0, 4, 1) *= sign*(-1) / std::sqrt(1 - std::pow(b_tilde.transpose()*r_hb_norm, 2));
	}
	else
	{
		sign = -1;
		get_quaternion(q, &quat);
		Q << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
		dx = ((LR_vector) (q->int_centers[DNANucleotide::BASE])).module();

		(*G).block(0, 0, 3, 1) = sign * (1/std::sqrt(1 - std::pow(b_bar.transpose()*r_hb_norm, 2))) * (1/d) * (b_bar - (b_bar.transpose()*r_hb_norm)*r_hb_norm);
		(*G).block(3, 0, 4, 1) = 2*(r_hb_norm.transpose() * Q).transpose() + 2*dx*(1/d)*(b_bar.transpose() * Q).transpose();
		(*G).block(3, 0, 4, 1) -= (1/std::pow(d, 3)) * (double)(b_bar.transpose()*r_hb).value() * 2*dx*(r_hb.transpose() * Q).transpose();
		(*G).block(3, 0, 4, 1) *= sign / std::sqrt(1 - std::pow(b_bar.transpose()*r_hb_norm, 2));
	}
	// which 2 is the same as which 1 for theta 2 but with reversed signs
}

void FEMBackend::get_T2_G5(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G, int which)
{
	// T2: TERM 2: Hydrogen bonding
	// G5: GRADIENT OF GEOMETRIC DESCRIPTOR 5: Angle theta 4
	Eigen::Vector3d n_bar, n_tilde;

	n_bar << p->orientationT.v3.x, p->orientationT.v3.y, p->orientationT.v3.z;
	n_tilde << p->orientationT.v3.x, q->orientationT.v3.y, q->orientationT.v3.z;

	Eigen::Vector4d quat;
	Eigen::template Matrix<double, 3, 4, 0, 3, 4> Q;

	double sign;
	if (which == 1)
	{
		sign = +1;
		get_quaternion(p, &quat);
		Q << quat(2), quat(3), quat(0), quat(1), -quat(1), -quat(0), quat(3), quat(2), quat(0), -quat(1), -quat(2), quat(3);
		Q = 2*Q;
		(*G).block(3, 0, 4, 1) = sign * (-1) / std::sqrt(1 - std::pow(n_bar.transpose()*n_tilde, 2)) * (n_tilde.transpose() * Q).transpose();
	}
	else
	{
		sign = -1;
		get_quaternion(q, &quat);
		Q << quat(2), quat(3), quat(0), quat(1), -quat(1), -quat(0), quat(3), quat(2), quat(0), -quat(1), -quat(2), quat(3);
		Q = 2*Q;
		(*G).block(3, 0, 4, 1) = sign * (-1) / std::sqrt(1 - std::pow(n_tilde.transpose()*n_bar, 2)) * (n_bar.transpose() * Q).transpose();
	}
	(*G).block(0, 0, 3, 1) << 0, 0, 0;
}

void FEMBackend::get_T2_G6(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G, int which)
{
	// T2: TERM 2: Hydrogen bonding
	// G6: GRADIENT OF GEOMETRIC DESCRIPTOR 6: Angle theta 7
	LR_vector r_hb_lr = get_descriptor_vector(p, q, 2);

	Eigen::Vector3d r_hb, r_hb_norm, n_bar, n_tilde;

	r_hb << r_hb_lr.x, r_hb_lr.y, r_hb_lr.z;
	r_hb_norm = r_hb / r_hb_lr.module();
	n_bar << p->orientationT.v3.x, p->orientationT.v3.y, p->orientationT.v3.z;
	n_tilde << q->orientationT.v3.x, q->orientationT.v3.y, q->orientationT.v3.z;

	Eigen::Vector4d quat;
	Eigen::template Matrix<double, 3, 4, 0, 3, 4> Q_r, Q_n;

	double sign, d(r_hb_lr.module()), dx;
	if (which == 1)
	{
		sign = +1;
		get_quaternion(p, &quat);
		Q_r << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
		Q_n << quat(2), quat(3), quat(0), quat(1), -quat(1), -quat(0), quat(3), quat(2), quat(0), -quat(1), -quat(2), quat(3);
		dx = ((LR_vector) (p->int_centers[DNANucleotide::BASE])).module();

		(*G).block(0, 0, 3, 1) = sign * (1/std::sqrt(1 - std::pow(n_bar.transpose()*r_hb_norm, 2))) * (1/d) * (n_bar - (n_bar.transpose()*r_hb_norm)*r_hb_norm);
		(*G).block(3, 0, 4, 1) = 2*(r_hb_norm.transpose() * Q_n).transpose() + 2*dx*(1/d)*(n_bar.transpose() * Q_r).transpose();
		(*G).block(3, 0, 4, 1) -= (1/std::pow(d, 3)) * (double)(n_bar.transpose()*r_hb).value() * 2*dx*(r_hb.transpose() * Q_r).transpose();
		(*G).block(3, 0, 4, 1) *= sign / std::sqrt(1 - std::pow(n_bar.transpose()*r_hb_norm, 2));
	}
	else
	{
		sign = -1;
		get_quaternion(q, &quat);
		Q_r << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
		Q_n << quat(2), quat(3), quat(0), quat(1), -quat(1), -quat(0), quat(3), quat(2), quat(0), -quat(1), -quat(2), quat(3);
		dx = ((LR_vector) (q->int_centers[DNANucleotide::BASE])).module();

		(*G).block(0, 0, 3, 1) = sign * (-1/std::sqrt(1 - std::pow(n_tilde.transpose()*r_hb_norm, 2))) * (1/d) * (n_tilde - (n_tilde.transpose()*r_hb_norm)*r_hb_norm);
		(*G).block(3, 0, 4, 1) = 2*dx*(1/d)*(n_tilde.transpose() * Q_r).transpose();
		(*G).block(3, 0, 4, 1) -= (1/std::pow(d, 3)) * (double)(n_tilde.transpose()*r_hb).value() * 2*dx*(r_hb.transpose() * Q_r).transpose();
		(*G).block(3, 0, 4, 1) *= sign*(-1) / std::sqrt(1 - std::pow(n_tilde.transpose()*r_hb_norm, 2));
	}
	// which 2 is the same as which 1 for theta 8 but with reversed signs
}

void FEMBackend::get_T2_G7(BaseParticle *p, BaseParticle *q, Eigen::VectorXd *G, int which)
{
        // T2: TERM 2: Hydrogen bonding
	// G7: GRADIENT OF GEOMETRIC DESCRIPTOR 7: Angle theta 8
	LR_vector r_hb_lr = get_descriptor_vector(p, q, 2);

	Eigen::Vector3d r_hb, r_hb_norm, n_bar, n_tilde;

	r_hb << r_hb_lr.x, r_hb_lr.y, r_hb_lr.z;
	r_hb_norm = r_hb / r_hb_lr.module();
	n_bar << p->orientationT.v3.x, p->orientationT.v3.y, p->orientationT.v3.z;
	n_tilde << q->orientationT.v3.x, q->orientationT.v3.y, q->orientationT.v3.z;

	Eigen::Vector4d quat;
	Eigen::template Matrix<double, 3, 4, 0, 3, 4> Q_r, Q_n;

	double sign, d(r_hb_lr.module()), dx;
	if (which == 1)
	{
		sign = +1;
		get_quaternion(p, &quat);
		Q_r << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
		Q_n << quat(2), quat(3), quat(0), quat(1), -quat(1), -quat(0), quat(3), quat(2), quat(0), -quat(1), -quat(2), quat(3);
		dx = ((LR_vector) (p->int_centers[DNANucleotide::BASE])).module();

		(*G).block(0, 0, 3, 1) = sign * ((-1)/std::sqrt(1 - std::pow(n_tilde.transpose()*r_hb_norm, 2))) * (1/d) * (n_tilde - (n_tilde.transpose()*r_hb_norm)*r_hb_norm);
		(*G).block(3, 0, 4, 1) = 2*dx*(1/d)*(n_tilde.transpose() * Q_r).transpose();
		(*G).block(3, 0, 4, 1) -= (1/std::pow(d, 3)) * (double)(n_tilde.transpose()*r_hb).value() * 2*dx*(r_hb.transpose() * Q_r).transpose();
		(*G).block(3, 0, 4, 1) *= sign*(-1) / std::sqrt(1 - std::pow(n_tilde.transpose()*r_hb_norm, 2));
	}
	else
	{
		sign = -1;
		get_quaternion(q, &quat);
		Q_r << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
		Q_n << quat(2), quat(3), quat(0), quat(1), -quat(1), -quat(0), quat(3), quat(2), quat(0), -quat(1), -quat(2), quat(3);
		dx = ((LR_vector) (q->int_centers[DNANucleotide::BASE])).module();

		(*G).block(0, 0, 3, 1) = sign * (1/std::sqrt(1 - std::pow(n_bar.transpose()*r_hb_norm, 2))) * (1/d) * (n_bar - (n_bar.transpose()*r_hb_norm)*r_hb_norm);

		(*G).block(3, 0, 4, 1) = 2*(r_hb_norm.transpose() * Q_n).transpose() + 2*dx*(1/d)*(n_bar.transpose() * Q_r).transpose();
		(*G).block(3, 0, 4, 1) -= (1/std::pow(d, 3)) * (double)(n_bar.transpose()*r_hb).value() * 2*dx*(r_hb.transpose() * Q_r).transpose();
		(*G).block(3, 0, 4, 1) *= sign / std::sqrt(1 - std::pow(n_bar.transpose()*r_hb_norm, 2));
	}
	// which 2 is the same as which 1 for theta 7 but with reversed signs
}

void FEMBackend::get_HG_THETA1(BaseParticle *p, BaseParticle *q, Eigen::MatrixXd *K, int which)
{
        Eigen::Vector3d b_bar, b_tilde, b_bar_s, b_tilde_s;
        Eigen::Vector4d quat;
        Eigen::template Matrix<double, 3, 4, 0, 3, 4> Q_bar, Q_tilde;
        Eigen::Matrix3d A_bar, A_tilde;

        b_bar << p->orientationT.v1.x, p->orientationT.v1.y, p->orientationT.v1.z;
        b_tilde << q->orientationT.v1.x, q->orientationT.v1.y, q->orientationT.v1.z;

        A_bar << p->orientation.v1.x, p->orientation.v1.y, p->orientation.v1.z;
        A_bar << p->orientation.v2.x, p->orientation.v2.y, p->orientation.v2.z;
        A_bar << p->orientation.v3.x, p->orientation.v3.y, p->orientation.v3.z;

        A_tilde << q->orientation.v1.x, q->orientation.v1.y, q->orientation.v1.z;
        A_tilde << q->orientation.v2.x, q->orientation.v2.y, q->orientation.v2.z;
        A_tilde << q->orientation.v3.x, q->orientation.v3.y, q->orientation.v3.y;

        b_bar_s = A_bar.transpose() * b_bar;
        b_tilde_s = A_tilde.transpose() * b_tilde;

        get_quaternion(p, quat);
        Q_bar << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
        Q_bar = 2 * Q_bar;
        get_quaternion(q, quat);
        Q_tilde << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
        Q_tilde = 2 * Q_tilde;

        switch ( which )
        {
        case 1:   // partial derivative p, nabla p
                (*K).block(3, 3, 4, 4) = Q_bar.transpose() * A_bar * A_tilde.transpose() * b_tilde * b_tilde.transpose() * Q_bar;
                (*K).block(3, 3, 4, 4) *= (-1)*b_bar_s.transpose()*b_tilde_s / std::sqrt(std::pow(1 - std::pow(b_bar_s.transpose()*b_tilde_s, 2), 3));
	        break;
        case 2:  // partial derivative p, nabla q
                (*K).block(3, 3, 4, 4) = Q_bar.transpose() * A_bar * A_tilde.transpose() * b_tilde * b_bar.transpose() * Q_tilde;
                (*K).block(3, 3, 4, 4) *= (-1)*b_tilde_s.transpose()*b_bar_s / std::sqrt(std::pow(1 - std::pow(b_tilde_s.transpose()*b_bar_s, 2), 3));
                (*K).block(3, 3, 4, 4) += ((-1) / std::sqrt(1 - std::pow(b_tilde_s.transpose()*b_bar_s, 2))) * Q_bar.transpose()*Q_tilde;
	        break;
        case 3:  // partial derivative q, nabla p
                (*K).block(3, 3, 4, 4) = Q_tilde.transpose() * b_bar.transpose() * A_bar * A_tilde.transpose() * b_tilde.transpose() * Q_bar;
                (*K).block(3, 3, 4, 4) *= (-1)*b_bar_s.transpose()*b_tilde_s / std::sqrt(std::pow(1 - std::pow(b_bar_s.transpose()*b_tilde_s, 2), 3));
                (*K).block(3, 3, 4, 4) += ((-1) / std::sqrt(1 - std::pow(b_bar_s.transpose()*b_tilde_s, 2))) * Q_tilde.transpose()*Q_bar;
	        break;
        case 4:  // partial derivative q, nabla q
                (*K).block(3, 3, 4, 4) = Q_tilde.transpose() * b_bar.transpose() * A_bar * A_tilde.transpose() * b_bar.transpose() * Q_tilde;
                (*K).block(3, 3, 4, 4) *=  (-1)*b_tilde_s.transpose()*b_bar_s / std::sqrt(std::pow(1 - std::pow(b_tilde_s.transpose()*b_bar_s, 2), 3));
	        break;
        }
        
        (*K).block(0, 0, 3, 3) << 0, 0, 0, 0, 0, 0, 0, 0, 0;
        (*K).block(0, 3, 3, 4) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        (*K).block(3, 0, 4, 3) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
}

void FEMBackend::get_HG_THETA2(BaseParticle *p, BaseParticle *q, Eigen::MatrixXd *K, int which)
{
        Eigen::Vector3d b_bar, b_tilde, r_hat;
        Eigen::Vector4d quat;
        Eigen::template Matrix<double, 3, 4, 0, 3, 4> Q_bar, Q_tilde, Q_r;
        Eigen::Matrix3d I;
        LRVector r;
        double d, dx;
        
        r = get_descriptor_vector(p, q, 2);
        r_hat << r.x, r.y, r.z;
        r_hat = r_hat / r.module();
        d = r.module();
        dx = ((LR_vector) (p->int_centers[DNANucleotide::BASE])).module();
        
        I.set_identity();

        b_bar << p->orientationT.v1.x, p->orientationT.v1.y, p->orientationT.v1.z;
        b_tilde << q->orientationT.v1.x, q->orientationT.v1.y, q->orientationT.v1.z;

        get_quaternion(p, quat);
        Q_bar << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
        Q_bar = 2*Q_bar;
        Q_r << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
        Q_r = 2*dx*Q_r;
        get_quaternion(q, quat);
        Q_tilde << quat(0), quat(1), -quat(2), -quat(3), quat(3), quat(2), quat(1), quat(0), -quat(2), quat(3), -quat(0), quat(1);
        Q_tilde = 2*Q_tilde;
        
        double s1, s2, s3, s4 dot1, dot2;
        s1 = (-1) / std::sqrt(1 - std::pow(b_bar.transpose()*r_hat, 2));
        s3 = 1 / std::sqrt(1 - std::pow(b_tilde.transpose()*r_hat, 2));
        s2 = b_bar.transpose()*r_hat / std::sqrt(std::pow(1 - std::pow(b_bar.transpose()*r_hat, 2), 3));
        s4 = b_tilde.transpose()*r_hat / std::sqrt(std::pow(1 - std::pow(b_tilde.transpose()*r_hat, 2), 3));
        dot1 = (double) (b_bar.transpose()*r_hat).value();
        dot2 = (double) (b_tilde.transpose()*r_hat).value();

        switch ( which )
        {
        case 1:   // partial derivative p, nabla p
                (*K).block(0, 0, 3, 3) = s1*(I*dot1 + r_hat*b_bar.transpose()) + s2*(b_bar - r_hat*dot1)*b_bar.transpose();
                (*K).block(0, 0, 3, 3) *= (1/std::pow(d, 2)) * (I - r_hat*r_hat.transpose());
                
                (*K).block(0, 3, 3, 4) = (-1)*s1 * (I*Q_bar - (1/d)*b_bar*r_hat.transpose()*Q_r - (dot1/d)*Q_r);
                (*K).block(0, 3, 3, 4) += s2 * b_bar * (r_hat.transpose()*Q_bar + (1/d)*b_bar.transpose()*Q_r - (dot1/d)*r_hat.transpose()*Q_r) ;
                (*K).block(0, 3, 3, 4) = (1/d) * (I - r_hat*r_hat.transpose()) * (*K).block(0, 3, 3, 4);
                
                (*K).block(3, 0, 4, 3) = (-1)*s1 * (I - r_hat*r_hat.transpose()) + s2 * (b_bar - dot1 * r_hat) * r_hat.transpose();
                (*K).block(3, 0, 4, 3) = (1/d)*Q_bar.transpose() * (*K).block(3, 0, 4, 3);
                
                (*K).block(3, 3, 4, 4) = (-1)*s1*(1/d) * (Q_r - r_hat*r_hat.transpose()*Q_r);
                (*K).block(3, 3, 4, 4) += s2 * r_hat * (r_hat.transpose()*Q_bar + (1/d)*b_bar.transpose()*Q_r - (dot1/d)*r_hat.transpose()*Q_r);
                (*K).block(3, 3, 4, 4) = Q_bar.transpose() * (*K).block(3, 3, 4, 4);
	        break;
        case 2:  // partial derivative p, nabla q
                (*K).block(0, 0, 3, 3) = (-1)*s3*(I*dot2 + r_hat*b_tilde.transpose()) + s4*(b_tilde - r_hat*dot2)*b_tilde.transpose();
                (*K).block(0, 0, 3, 3) = (1/std::pow(d, 2)) * (I - r_hat*r_hat.transpose()) * (*K).block(0, 0, 3, 3);
                
                (*K).block(0, 3, 3, 4) = s3*(Q_tilde - (1/d)*b_tilde*r_hat.transpose()*Q_r - (dot2/d)*Q_r);
                (*K).block(0, 3, 3, 4) += s4*b_tilde*(r_hat.transpose()*Q_tilde + (1/d)*b_tilde.transpose()*Q_r - (dot2/d)*r_hat.transpose()*Q_r);
                (*K).block(0, 3, 3, 4) = (1/d)*(I - r_hat*r_hat.transpose()) * (*K).block(0, 3, 3, 4);
                
                (*K).block(3, 0, 4, 3) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                (*K).block(3, 3, 4, 4) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
	        break;
        case 3:  // partial derivative q, nabla p
                (*K).block(0, 0, 3, 3) = s1*(I*dot1 + r_hat*b_bar.transpose()) + s2*(b_bar - r_hat*dot1)*b_bar.transpose();
                (*K).block(0, 0, 3, 3) = (-1)*(1/std::pow(d, 2)) * (I - r_hat*r_hat.transpose()) * (*K).block(0, 0, 3, 3);
                
                (*K).block(0, 3, 3, 4) = (-1)*s1 * (I*Q_bar - (1/d)*b_bar*r_hat.transpose()*Q_r - (dot1/d)*Q_r);
                (*K).block(0, 3, 3, 4) += s2 * b_bar * (r_hat.transpose()*Q_bar + (1/d)*b_bar.transpose()*Q_r - (dot1/d)*r_hat.transpose()*Q_r) ;
                (*K).block(0, 3, 3, 4) = (-1)*(1/d) * (I - r_hat*r_hat.transpose()) * (*K).block(0, 3, 3, 4);
                
                (*K).block(3, 0, 4, 3) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                (*K).block(3, 3, 4, 4) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
	        break;
        case 4:  // partial derivative q, nabla q
                (*K).block(0, 0, 3, 3) = (-1)*s3*(I*dot2 + r_hat*b_tilde.transpose()) + s4*(b_tilde - r_hat*dot2)*b_tilde.transpose();
                (*K).block(0, 0, 3, 3) = (1/std::pow(d, 2)) * (I - r_hat*r_hat.transpose()) * (*K).block(0, 0, 3, 3);
                
                (*K).block(0, 3, 3, 4) = s3*(Q_tilde - (1/d)*b_tilde*r_hat.transpose()*Q_r - (dot2/d)*Q_r);
                (*K).block(0, 3, 3, 4) += s4*b_tilde*(r_hat.transpose()*Q_tilde + (1/d)*b_tilde.transpose()*Q_r - (dot2/d)*r_hat.transpose()*Q_r);
                (*K).block(0, 3, 3, 4) = (1/d)*(I - r_hat*r_hat.transpose()) * (*K).block(0, 3, 3, 4);
                
                (*K).block(3, 0, 4, 3) = s3*(I - r_hat*r_hat.transpose()) + s4*(b_tilde - r_hat*dot2)*r_hat.transpose();
                (*K).block(3, 0, 4, 3) = (1/d)*Q_tilde.transpose() * (*K).block(3, 0, 4, 3);
                
                (*K).block(3, 3, 4, 4) = s3*(1/d)*(Q_r - r_hat*r_hat.transpose()*Q_r) + s4*(1/d)*r_hat*(b_tilde.transpose()*Q_r - dot2*r_hat.transpose()*Q_r);
                (*K).block(3, 3, 4, 4) = Q_tilde.transpose() * (*K).block(3, 3, 4, 4)
	        break;
        }
}

void FEMBackend::get_HG_THETA3(BaseParticle *p, BaseParticle *q, Eigen::MatrixXd *K, int which)
{
        // WRITE THIS
}

void FEMBackend::get_HG_THETA4(BaseParticle *p, BaseParticle *q, Eigen::MatrixXd *K, int which)
{
        // WRITE THIS
}

void FEMBackend::get_HG_THETA7(BaseParticle *p, BaseParticle *q, Eigen::MatrixXd *K, int which)
{
        // WRITE THIS
}

void FEMBackend::get_HG_THETA8(BaseParticle *p, BaseParticle *q, Eigen::MatrixXd *K, int which)
{
        // WRITE THIS
}

void FEMBackend::get_KF_local_bonded(BaseParticle *p, BaseParticle *q, Eigen::MatrixXd *K, Eigen::VectorXd *F)
{
	//std::cout << "Computing K, F local bonded." << std::endl;
	Eigen::MatrixXd I(7, 7), HG(7, 7);
	Eigen::VectorXd Gp(7), Gq(7);

	I.setIdentity();
	
	// BACKBONE BONDS
	
	// G1: R BACK-BACK

	get_T1_G1(p, q, &Gp, 1);
	get_T1_G1(q, p, &Gq, 1);

	(*K).block(p->index*7, p->index*7, 7, 7) += get_D21_G1(p, q) * (Gp * Gp.transpose()) + get_D11_G1(p, q) * (I - Gp * Gp.transpose());
	(*K).block(q->index*7, p->index*7, 7, 7) += get_D21_G1(q, p) * (Gq * Gp.transpose()) + get_D11_G1(q, p) * (I - Gq * Gp.transpose());
	(*F).block(p->index*7, 0, 7, 1) += get_D11_G1(p, q) * Gp;

	get_T1_G1(p, q, &Gp, 2);
	get_T1_G1(q, p, &Gq, 2);

	(*K).block(p->index*7, q->index*7, 7, 7) += get_D21_G1(p, q) * (Gp * Gq.transpose()) + get_D11_G1(p, q) * (I - Gp * Gq.transpose());
	(*K).block(q->index*7, q->index*7, 7, 7) += get_D21_G1(q, p) * (Gq * Gq.transpose()) + get_D11_G1(q, p) * (I - Gq * Gq.transpose());
	(*F).block(q->index*7, 0, 7, 1) += get_D11_G1(q, p) * Gq;

        // STACKING
        
        // EXCLUDED VOLUME
}

void FEMBackend::get_KF_local_nonbonded(BaseParticle *p, BaseParticle *q, Eigen::MatrixXd *K, Eigen::VectorXd *F)
{
	//std::cout << "Computing K, F local nonbonded." << std::endl;
	Eigen::MatrixXd I(7, 7), HG(7, 7);
	Eigen::VectorXd Gp(7), Gq(7);
	
	I.setIdentity();
	
	// HYDROGEN BONDS
	
	// G1: R BASE-BASE
	
	get_T2_G1(p, q, &Gp, 1);
	get_T2_G1(p, q, &Gq, 1);
	
	(*K).block(p->index*7, p->index*7, 7, 7) += get_D22_G1(p, q) * (Gp * Gp.transpose()) + get_D12_G1(p, q) * (I - Gp * Gp.transpose());
	(*K).block(q->index*7, p->index*7, 7, 7) += get_D22_G1(q, p) * (Gq * Gp.transpose()) + get_D12_G1(q, p) * (I - Gq * Gp.transpose());
	(*F).block(p->index*7, 0, 7, 1) += get_D12_G1(p, q) * Gp;
	
	get_T2_G1(p, q, &Gp, 2);
	get_T2_G1(p, q, &Gq, 2);

	(*K).block(p->index*7, q->index*7, 7, 7) += get_D22_G1(p, q) * (Gp * Gq.transpose()) + get_D12_G1(p, q) * (I - Gp * Gq.transpose());
	(*K).block(q->index*7, q->index*7, 7, 7) += get_D22_G1(q, p) * (Gq * Gq.transpose()) + get_D12_G1(q, p) * (I - Gq * Gq.transpose());
	(*F).block(q->index*7, 0, 7, 1) += get_D12_G1(q, p) * Gq;

	// G2: THETA 1
	
	get_T2_G2(p, q, &Gp, 1);
	get_T2_G2(p, q, &Gq, 1);
	
	get_HG_THETA1(p, q, HG, 1)
	(*K).block(p->index*7, p->index*7, 7, 7) += get_D22_G2(p, q) * (Gp * Gp.transpose()) + get_D12_G2(p, q) * HG;
	get_HG_THETA1(p, q, HG, 3)
	(*K).block(q->index*7, p->index*7, 7, 7) += get_D22_G2(q, p) * (Gq * Gp.transpose()) + get_D12_G2(q, p) * HG;
	(*F).block(p->index*7, 0, 7, 1) += get_D12_G2(p, q) * Gp;
		
	get_T2_G2(p, q, &Gp, 2);
	get_T2_G2(p, q, &Gq, 2);

        get_HG_THETA1(p, q, HG, 2)
	(*K).block(p->index*7, q->index*7, 7, 7) += get_D22_G2(p, q) * (Gp * Gq.transpose()) + get_D12_G2(p, q) * HG;
	get_HG_THETA1(p, q, HG, 4)
	(*K).block(q->index*7, q->index*7, 7, 7) += get_D22_G2(q, p) * (Gq * Gq.transpose()) + get_D12_G2(q, p) * HG;
	(*F).block(q->index*7, 0, 7, 1) += get_D12_G2(q, p) * Gq;
	
	// G3: THETA 2
	
	get_T2_G3(p, q, &Gp, 1);
	get_T2_G3(p, q, &Gq, 1);
	
	get_HG_THETA2(p, q, HG, 1)
	(*K).block(p->index*7, p->index*7, 7, 7) += get_D22_G3(p, q) * (Gp * Gp.transpose()) + get_D12_G3(p, q) * HG;
	get_HG_THETA2(p, q, HG, 3)
	(*K).block(q->index*7, p->index*7, 7, 7) += get_D22_G3(q, p) * (Gq * Gp.transpose()) + get_D12_G3(q, p) * HG;
	(*F).block(p->index*7, 0, 7, 1) += get_D12_G3(p, q) * Gp;
			
	get_T2_G3(p, q, &Gp, 2);
	get_T2_G3(p, q, &Gq, 2);

        get_HG_THETA2(p, q, HG, 2)
	(*K).block(p->index*7, q->index*7, 7, 7) += get_D22_G3(p, q) * (Gp * Gq.transpose()) + get_D12_G3(p, q) * HG;
	get_HG_THETA2(p, q, HG, 4)
	(*K).block(q->index*7, q->index*7, 7, 7) += get_D22_G3(q, p) * (Gq * Gq.transpose()) + get_D12_G3(q, p) * HG;
	(*F).block(q->index*7, 0, 7, 1) += get_D12_G3(q, p) * Gq;
	
	// G4: THETA 3
	
	get_T2_G4(p, q, &Gp, 1);
	get_T2_G4(p, q, &Gq, 1);
			
	(*K).block(p->index*7, p->index*7, 7, 7) += get_D22_G4(p, q) * (Gp * Gp.transpose()) // + get_D12_G4(p, q) * (...);
	(*K).block(q->index*7, p->index*7, 7, 7) += get_D22_G4(q, p) * (Gq * Gp.transpose()) // + get_D12_G4(q, p) * (...);
	(*F).block(p->index*7, 0, 7, 1) += get_D12_G4(p, q) * Gp;
			
	get_T2_G4(p, q, &Gp, 2);
	get_T2_G4(p, q, &Gq, 2);

	(*K).block(p->index*7, q->index*7, 7, 7) += get_D22_G4(p, q) * (Gp * Gq.transpose()) // + get_D12_G4(p, q) * (...);
	(*K).block(q->index*7, q->index*7, 7, 7) += get_D22_G4(q, p) * (Gq * Gq.transpose()) // + get_D12_G4(q, p) * (...);
	(*F).block(q->index*7, 0, 7, 1) += get_D12_G4(q, p) * Gq;
	
	// G5: THETA 4
	
	get_T2_G5(p, q, &Gp, 1);
	get_T2_G5(p, q, &Gq, 1);
			
	(*K).block(p->index*7, p->index*7, 7, 7) += get_D22_G5(p, q) * (Gp * Gp.transpose()) // + get_D12_G5(p, q) * (...);
	(*K).block(q->index*7, p->index*7, 7, 7) += get_D22_G5(q, p) * (Gq * Gp.transpose()) // + get_D12_G5(q, p) * (...);
	(*F).block(p->index*7, 0, 7, 1) += get_D12_G5(p, q) * Gp;
			
	get_T2_G5(p, q, &Gp, 2);
	get_T2_G5(p, q, &Gq, 2);

	(*K).block(p->index*7, q->index*7, 7, 7) += get_D22_G5(p, q) * (Gp * Gq.transpose()) // + get_D12_G5(p, q) * (...);
	(*K).block(q->index*7, q->index*7, 7, 7) += get_D22_G5(q, p) * (Gq * Gq.transpose()) // + get_D12_G5(q, p) * (...);
	(*F).block(q->index*7, 0, 7, 1) += get_D12_G5(q, p) * Gq;
	
	// G6: THETA 7
	
	get_T2_G6(p, q, &Gp, 1);
	get_T2_G6(p, q, &Gq, 1);
			
	(*K).block(p->index*7, p->index*7, 7, 7) += get_D22_G6(p, q) * (Gp * Gp.transpose()) // + get_D12_G6(p, q) * (...);
	(*K).block(q->index*7, p->index*7, 7, 7) += get_D22_G6(q, p) * (Gq * Gp.transpose()) // + get_D12_G6(q, p) * (...);
	(*F).block(p->index*7, 0, 7, 1) += get_D12_G6(p, q) * Gp;
			
	get_T2_G6(p, q, &Gp, 2);
	get_T2_G6(p, q, &Gq, 2);

	(*K).block(p->index*7, q->index*7, 7, 7) += get_D22_G6(p, q) * (Gp * Gq.transpose()) // + get_D12_G6(p, q) * (...);
	(*K).block(q->index*7, q->index*7, 7, 7) += get_D22_G6(q, p) * (Gq * Gq.transpose()) // + get_D12_G6(q, p) * (...);
	(*F).block(q->index*7, 0, 7, 1) += get_D12_G6(q, p) * Gq;
	
	// G7: THETA 8
		
	get_T2_G7(p, q, &Gp, 1);
	get_T2_G7(p, q, &Gq, 1);
				
	(*K).block(p->index*7, p->index*7, 7, 7) += get_D22_G7(p, q) * (Gp * Gp.transpose()) // + get_D12_G7(p, q) * (...);
	(*K).block(q->index*7, p->index*7, 7, 7) += get_D22_G7(q, p) * (Gq * Gp.transpose()) // + get_D12_G7(q, p) * (...);
	(*F).block(p->index*7, 0, 7, 1) += get_D12_G7(p, q) * Gp;
				
	get_T2_G7(p, q, &Gp, 2);
	get_T2_G7(p, q, &Gq, 2);

	(*K).block(p->index*7, q->index*7, 7, 7) += get_D22_G7(p, q) * (Gp * Gq.transpose()) // + get_D12_G7(p, q) * (...);
	(*K).block(q->index*7, q->index*7, 7, 7) += get_D22_G7(q, p) * (Gq * Gq.transpose()) // + get_D12_G7(q, p) * (...);
	(*F).block(q->index*7, 0, 7, 1) += get_D12_G7(q, p) * Gq;

	// CROSS STACKING
	
	// COAXIAL STACKING
	
	// EXCLUDED VOLUME
	
	// DEBYE-HUCKEL
}

void FEMBackend::get_settings(input_file &inp)
{
	MDBackend::get_settings(inp);
	getInputDouble(&inp, "alpha", &_alpha, 0);
	getInputDouble(&inp, "alpha_step", &_alpha_step, 0);
	getInputDouble(&inp, "n_iter_in", &_n_iter_in, 0);
	getInputDouble(&inp, "max_backbone_force", &_max_f, 0);
	getInputDouble(&inp, "max_backbone_force_far", &_max_ff, 0);
}

void FEMBackend::init()
{
	MDBackend::init();
	OX_LOG(Logger::LOG_INFO, "Using start alpha = %g, alpha step = %g, and number of inner iterations = %g", _alpha, _alpha_step, _n_iter_in);

	if (_max_f != 0 && _max_ff != 0) {
		double b = FENE_EPS - 2*_max_f*FENE_R0_OXDNA2;
		double c = (-1)*FENE_EPS*FENE_R0_OXDNA2 + _max_f*FENE_DELTA2*FENE_DELTA2 - _max_f*FENE_R0_OXDNA*FENE_R0_OXDNA;
		_r_max = (1/(2*_max_f)) * ((-1)*b + std::sqrt(b*b - 4*_max_f*c));

		_mod_B = (_r_max*FENE_EPS*(_r_max - FENE_R0_OXDNA2)) / (std::pow(FENE_DELTA2, 2) - std::pow(_r_max - FENE_R0_OXDNA2, 2)) - _max_ff;
		_mod_C = (((-1)*FENE_EPS)/2) * std::log( 1 - std::pow(_r_max - FENE_R0_OXDNA2, 2)/std::pow(FENE_DELTA2, 2) ) - _max_ff * _r_max - _mod_B * std::log(_r_max);
	}
}

void FEMBackend::evolve(Eigen::VectorXd *dU)
{
	for(auto p: _particles) {
		LR_vector dpos = LR_vector((double) (*dU)((int) p->index), (double) (*dU)(p->index + 1), (double) (*dU)(p->index + 2));

		double dR_11, dR_12, dR_13, dR_21, dR_22, dR_23, dR_31, dR_32, dR_33, s;
		double n = std::pow((*dU)(p->index + 3), 2) + std::pow((*dU)(p->index + 4), 2) + std::pow((*dU)(p->index + 5), 2) + std::pow((*dU)(p->index + 6), 2);

		if (n*n <= 1e-24)
		{
			s = 0;
		}
		else
		{
			s = 2/n;
		}

		dR_11 = 1 - s*(std::pow((*dU)(p->index + 5), 2) + std::pow((*dU)(p->index + 6), 2));
		dR_12 = s*((*dU)(p->index + 4)*(*dU)(p->index + 5) - (*dU)(p->index + 3)*(*dU)(p->index + 6));
		dR_13 = s*((*dU)(p->index + 4)*(*dU)(p->index + 6) + (*dU)(p->index + 3)*(*dU)(p->index + 5));

		dR_21 = s*((*dU)(p->index + 4)*(*dU)(p->index + 5) + (*dU)(p->index + 3)*(*dU)(p->index + 6));
		dR_22 = 1 - s*(std::pow((*dU)(p->index + 4), 2) + std::pow((*dU)(p->index + 6), 2));
		dR_23 = s*((*dU)(p->index + 5)*(*dU)(p->index + 6) - (*dU)(p->index + 3)*(*dU)(p->index + 4));

		dR_31 = s*((*dU)(p->index + 4)*(*dU)(p->index + 6) - (*dU)(p->index + 3)*(*dU)(p->index + 5));
		dR_32 = s*((*dU)(p->index + 5)*(*dU)(p->index + 6) + (*dU)(p->index + 3)*(*dU)(p->index + 4));
		dR_33 = 1 - s*(std::pow((*dU)(p->index + 4), 2) + std::pow((*dU)(p->index + 5), 2));


		LR_matrix dR = LR_matrix(dR_11, dR_12, dR_13, dR_21, dR_22, dR_23, dR_31, dR_32, dR_33);

		p->pos = p->pos + dpos;
		p->orientation = p->orientation * dR;

		p->set_positions();
		p->orientationT = p->orientation.get_transpose();

		_lists->single_update(p);
	}
	return;
}

void FEMBackend::compute_KF(Eigen::MatrixXd *K_global, Eigen::VectorXd *F_global)
{
	for(auto p: _particles) {
			for (auto q: _particles) {
				if (p->is_bonded(q) && p->index > q->index) {
					get_KF_local_bonded(p, q, K_global, F_global);
				}
				else if (p->index > q->index) {
					get_KF_local_nonbonded(p, q, K_global, F_global);
				}
			}
		}
}

void FEMBackend::sim_step() {

	int N_dim = _particles.size() * 7;

	Eigen::VectorXd F_global(N_dim), dU(N_dim);
	Eigen::MatrixXd K_global(N_dim, N_dim);

	Eigen::VectorXd *F_p, *dU_p;
	Eigen::MatrixXd *K_p;

	F_p = &F_global;
	dU_p = &dU;
	K_p = &K_global;

	_mytimer->resume();

	for(int i = 0; i < _n_iter_in; i++)
	{
		F_global = 0 * F_global;
		K_global = 0 * K_global;

		std::cout << "Inner loop iteration " << i << " with alpha value " << _alpha << std::endl;
		compute_KF(K_p, F_p);
		std::cout << K_global.block(0, 0, 14, 14) << std::endl;

		F_global = _alpha * F_global;

		if (K_global.determinant() != 0)
		{
			dU = K_global.householderQr().solve(F_global);
			double error = (K_global * dU - F_global).norm() / F_global.norm();
			if (error <= 0.1)
			{
				evolve(dU_p);
			}
			else
			{
				std::cout << "Relative error is " << error << ". No update." << std::endl;
			}
		}
		else
		{
			std::cout << "Global stiffness matrix is singular. No update." << std::endl;
		}

		if(!_lists->is_updated())
		{
			_lists->global_update();
			_N_updates++;
		}
	}
	_alpha = _alpha + _alpha_step;
	_mytimer->pause();
}

