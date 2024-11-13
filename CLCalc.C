#include "CLCalc.h"

registerMooseObject("otterApp", CLCalc);

InputParameters
CLCalc::validParams()
{
    InputParameters params = ADKernelGrad::validParams();
    params.addRequiredCoupledVar("C_L", "Concentration C_L.");
    params.addRequiredCoupledVar("d", "Phase field damage variable.");
    params.addRequiredCoupledVar("hidrostatic_stress", "hidrostatic_stress variable.");
    params.addRequiredCoupledVar("N_T", "N_T variable.");
    params.addRequiredCoupledVar("effective_plastic_strain", "Effective plastic strain.");
    params.addParam<Real>("D_L", 12700.0, "Diffusion coefficient D_L (micrometers squared per second).");
    params.addParam<Real>("V_H", 2e12, "Volume or velocity V_H (micrometers cubed per mole).");
    params.addParam<Real>("R", 4.1240, "R value.");
    params.addParam<Real>("Theta", 554.25, "Theta value.");
    params.addParam<Real>("Beta_L", 6, "Number of Lattice sites per solvent");
    params.addParam<Real>("N_L", 8.46e10, "Lattice site density");
    params.addParam<Real>("W_B", 29820.4, "Binding energy");   
    params.addRequiredCoupledVar("CLCalc_2", "Term used by the second part of the equation");
    params.addRequiredCoupledVar("CLCalc_3", "Term used by the third part of the equation");
    return params;
}

CLCalc::CLCalc(const InputParameters & parameters)
    : ADKernelGrad(parameters),
    _C_L(coupledValue("C_L")),
    _grad_C_L(coupledGradient("C_L")),
    _der_C_L(coupledDot("C_L")),
    _d(coupledValue("d")),
    _hidrostatic_stress(coupledValue("hidrostatic_stress")),
    _grad_hidrostatic_stress(coupledGradient("hidrostatic_stress")),
    _N_T(coupledValue("N_T")),
    _effective_plastic_strain(coupledValue("effective_plastic_strain")),
    _der_effective_plastic_strain(coupledDot("effective_plastic_strain")),
    _D_L(getParam<Real>("D_L")),  // Obtener D_L desde los parámetros
    _V_H(getParam<Real>("V_H")),  // Obtener V_H desde los parámetros
    _R(getParam<Real>("R")),
    _Theta(getParam<Real>("Theta")),
    _Beta_L(getParam<Real>("Beta_L")),
    _N_L(getParam<Real>("N_L")),
    _W_B(getParam<Real>("W_B")),
    _div_CLCalc_2(coupledDiv("CLCalc_2")),
    _div_CLCalc_3(coupledDiv("CLCalc_3"))
{
}

ADRealVectorValue
CLCalc::precomputeQpResidual()
{
    Real term_theta_T = ((((_C_L[_qp])/(_Beta_L*_N_L))/(1-(_C_L[_qp])/(_Beta_L*_N_L)))*exp(_W_B/(_R*_Theta)))/(1+((((_C_L[_qp])/(_Beta_L*_N_L))/(1-(_C_L[_qp])/(_Beta_L*_N_L)))*exp(_W_B/(_R*_Theta))));
    Real term1 = ((_C_L[_qp] + (_Beta_L*term_theta_T*_N_T[_qp]) * (1.0 - term_theta_T)) / _C_L[_qp]) * _der_C_L[_qp];
    Real term2 = _div_CLCalc_2[_qp];
    Real term3 = _div_CLCalc_3[_qp];
    Real derivative_N_T = 2.30258 * _N_T[_qp] * 12.815 * std::exp(-5.5 * _effective_plastic_strain[_qp]);  // Derivada de N_T con respecto a effective plastic strain
    Real term4 = term_theta_T * derivative_N_T * _der_effective_plastic_strain[_qp];
    return (-term1 + term2 - term3 - term4);
//  return (_permeability / _viscosity) * _grad_u[_qp]; 
}

// Real
// KernelEq1::computeQpResidual()
// {
    // Término 1: Derivada temporal de C_L multiplicada por un factor
//    Real term1 = ((_C_L[_qp] + _C_T[_qp] * (1.0 - _theta_T[_qp])) / _C_L[_qp]) * _der_C_L[_qp];

    // Término 2: Divergencia del término de difusión
//    Real term2 = -coupledDiv((1.0 - _d[_qp]) * (1.0 - _d[_qp]) * _D_L * _grad_C_L[_qp]);

    // Término 3: Divergencia del término adicional con concentración, velocidad y gradiente de sigma_H
//    Real term3 = coupledDiv((1.0 - _d[_qp]) * (1.0 - _d[_qp]) * _D_L * _C_L[_qp] * _V_H * _grad_hidrostatic_stress[_qp] / (_R * _Theta));

    // Término 4: Derivada temporal de epsilon multiplicada por la derivada de N_T con respecto a epsilon
//    Real derivative_N_T = 2.30258 * _N_T[_qp] * 12.815 * std::exp(-5.5 * _effective_plastic_strain[_qp]);  // Derivada de N_T con respecto a epsilon
//    Real term4 = _theta_T[_qp] * derivative_N_T * _der_effective_plastic_strain[_qp];

//    return term1 + term2 + term3 + term4;
//}

//Real
//KernelEq1::computeQpJacobian()
//{
    // Aquí puedes calcular la derivada respecto a las variables si es necesario
//    return 1.0; // Para simplificar el Jacobiano
//}