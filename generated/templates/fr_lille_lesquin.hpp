// Auto-generated linear harmonic climate model
// Station name : LILLE-LESQUIN
// Station code : 59343001
#pragma once
#include <cmath>
namespace harmoclimat {
static constexpr double longitude_deg = 3.0975000858306885;
static constexpr double latitude_deg = 50.569999694824219;
static constexpr double delta_utc_solar_h = 0.20650000572204597;
namespace detail {
static constexpr double two_pi = 6.2831853071795864769;
static constexpr double solar_year_days = 365.242189;
static constexpr double omega_annual = two_pi / solar_year_days;
static constexpr double omega_diurnal = two_pi / 24.0;
inline double eval_annual(const double* coeffs, int n_annual, double day){
    double value = coeffs[0];
    for(int k = 1; k <= n_annual; ++k){
        double angle = k * omega_annual * day;
        value += coeffs[2*k - 1] * std::cos(angle);
        value += coeffs[2*k] * std::sin(angle);
    }
    return value;
}
inline double wrap_day(double d){
    while (d >= solar_year_days) d -= solar_year_days;
    while (d < 0.0)   d += solar_year_days;
    return d;
}
inline double wrap_hour(double h){
    while (h >= 24.0) h -= 24.0;
    while (h < 0.0)   h += 24.0;
    return h;
}
} // namespace detail
namespace temperature_model {
static constexpr int n_diurnal = 3;
static constexpr int c0_coeffs_n_annual = 3;
static constexpr double c0_coeffs[] = {
    11.408654816501448, -7.1344330783836538, -2.5371986716383326, -0.091709115110473194, 0.45328033834862852, 0.0026837603685719649, -0.061581750793561228
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -2.7465633712497821, 1.6904548056629876, -0.042141091640045994, 0.40636509901327889, -0.068655941724880487, -0.13884030341669135, 0.040753018318074159
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -0.74141510762252605, -0.035416507439799023, -0.25602417389460025, 0.12777496074074618, -0.14608540338394838, 0.15439520795695535, 0.099691984970186442
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    0.43840019600051477, 0.22895145789819604, -0.094036553837696504, -0.23057119287836345, 0.064568197837564525, -0.022572130549296531, -0.0023603405718186808
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.093737502897558195, 0.20433780493914117, 0.083719318383519523, 0.13721036727640723, 0.041940736329872698, -0.087014716315704579, -0.019546002848383884
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.09221857683898578, -0.27177928119979056, 0.032292508344728278, -0.064499220108517386, 0.037153027742423475, 0.058211839622487847, -0.025492009720612836
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -0.080588410186253137, 0.1466931464120971, -0.027562129820183974, -0.055604489452664135, 0.0051169797498820304, -0.051825507280423407, -0.00092902164614900532
};
inline double evaluate(double day_solar, double hour_solar) {
    double result = detail::eval_annual(c0_coeffs, c0_coeffs_n_annual, day_solar);
    const double ome_d = detail::omega_diurnal;
    result += detail::eval_annual(a1_coeffs, a1_coeffs_n_annual, day_solar) * std::cos(1 * ome_d * hour_solar);
    result += detail::eval_annual(b1_coeffs, b1_coeffs_n_annual, day_solar) * std::sin(1 * ome_d * hour_solar);
    result += detail::eval_annual(a2_coeffs, a2_coeffs_n_annual, day_solar) * std::cos(2 * ome_d * hour_solar);
    result += detail::eval_annual(b2_coeffs, b2_coeffs_n_annual, day_solar) * std::sin(2 * ome_d * hour_solar);
    result += detail::eval_annual(a3_coeffs, a3_coeffs_n_annual, day_solar) * std::cos(3 * ome_d * hour_solar);
    result += detail::eval_annual(b3_coeffs, b3_coeffs_n_annual, day_solar) * std::sin(3 * ome_d * hour_solar);
    return result;
}
} // namespace temperature_model
namespace specific_humidity_model {
static constexpr int n_diurnal = 3;
static constexpr int c0_coeffs_n_annual = 3;
static constexpr double c0_coeffs[] = {
    0.0067742570829238914, -0.0021971583411144985, -0.0013682875519197609, 0.00015978738960744295, 0.00030188801859361303, -5.7087226771165781e-05, 5.9316722474003108e-05
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -0.00010708345425870761, -7.7938468316760234e-05, 3.1879455878288038e-05, 3.8282273951182557e-05, 5.3630405983920669e-05, 3.4598354093296022e-05, -4.7633623115655288e-05
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    4.697654232522565e-06, -0.00011989190116445543, -1.4122753014024712e-05, 2.6162266599793299e-05, 1.7770687014171122e-05, 2.6377655600795191e-05, -1.596669777997595e-05
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -4.3969499806769321e-05, 0.00011732290354991743, 1.6716563212460533e-05, 2.0078747679867898e-05, -2.9659754356780019e-05, -3.1926713877598147e-05, 2.1109441714468897e-06
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -6.470371457607242e-05, 4.3137409809967147e-05, 1.2559303123569934e-05, 4.1492796342237407e-05, -1.1405088450888165e-05, -1.0379510129079166e-05, -2.1974707143099903e-06
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    1.0006086944715888e-05, -2.8869517556712432e-06, -1.960288265903049e-05, -4.250121244020914e-05, 9.8311898459003777e-06, 8.6422390671760096e-06, 1.0023672239332394e-05
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -1.4009818377578179e-05, 5.0472725082196644e-05, 1.4814713366999668e-06, -7.4485594630700862e-06, -3.3129444998515559e-06, -1.7797133002649847e-05, 2.6366938554885292e-06
};
inline double evaluate(double day_solar, double hour_solar) {
    double result = detail::eval_annual(c0_coeffs, c0_coeffs_n_annual, day_solar);
    const double ome_d = detail::omega_diurnal;
    result += detail::eval_annual(a1_coeffs, a1_coeffs_n_annual, day_solar) * std::cos(1 * ome_d * hour_solar);
    result += detail::eval_annual(b1_coeffs, b1_coeffs_n_annual, day_solar) * std::sin(1 * ome_d * hour_solar);
    result += detail::eval_annual(a2_coeffs, a2_coeffs_n_annual, day_solar) * std::cos(2 * ome_d * hour_solar);
    result += detail::eval_annual(b2_coeffs, b2_coeffs_n_annual, day_solar) * std::sin(2 * ome_d * hour_solar);
    result += detail::eval_annual(a3_coeffs, a3_coeffs_n_annual, day_solar) * std::cos(3 * ome_d * hour_solar);
    result += detail::eval_annual(b3_coeffs, b3_coeffs_n_annual, day_solar) * std::sin(3 * ome_d * hour_solar);
    return result;
}
} // namespace specific_humidity_model
namespace pressure_model {
static constexpr int n_diurnal = 3;
static constexpr int c0_coeffs_n_annual = 3;
static constexpr double c0_coeffs[] = {
    1010.3909377848545, -0.080941794241282256, 0.26451957251142677, 0.43364188387017738, 0.62208772516102007, 0.4222564687497235, 0.66818786260877283
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    0.03356089319550138, 0.0046092892334771594, 0.023536026745992516, 0.010307771944655986, 0.0038838648418681562, 0.023768447197713982, 0.0014966159564794104
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    0.058277992322525639, -0.090849973221084138, 0.036816091935996345, -0.025111125828640803, -0.0045283400144537477, 0.009830944601582705, -0.0069284660463496692
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -0.10396053131804685, 0.086785062149291728, 0.037217750142345948, 0.053115516940608135, 0.035399732488098502, -0.042678520877882389, -0.0028089153102113431
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.33784088335922219, -0.018579007244784472, -0.03194908855365803, 0.032759838993314951, -0.018047447827538934, 0.0062401216802394745, 0.017421164798204446
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.010516750669683635, 0.092856282665106621, -0.033661351529480896, -0.015954499082804517, -0.010995917147740958, -0.0046139242894680296, 0.0027265025623667591
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    0.038092108213150282, 0.092701748933566488, -0.0022421757616459346, 0.016564615056069973, 0.0090014558232676643, 0.0035219503902614068, 0.0070958139748153121
};
inline double evaluate(double day_solar, double hour_solar) {
    double result = detail::eval_annual(c0_coeffs, c0_coeffs_n_annual, day_solar);
    const double ome_d = detail::omega_diurnal;
    result += detail::eval_annual(a1_coeffs, a1_coeffs_n_annual, day_solar) * std::cos(1 * ome_d * hour_solar);
    result += detail::eval_annual(b1_coeffs, b1_coeffs_n_annual, day_solar) * std::sin(1 * ome_d * hour_solar);
    result += detail::eval_annual(a2_coeffs, a2_coeffs_n_annual, day_solar) * std::cos(2 * ome_d * hour_solar);
    result += detail::eval_annual(b2_coeffs, b2_coeffs_n_annual, day_solar) * std::sin(2 * ome_d * hour_solar);
    result += detail::eval_annual(a3_coeffs, a3_coeffs_n_annual, day_solar) * std::cos(3 * ome_d * hour_solar);
    result += detail::eval_annual(b3_coeffs, b3_coeffs_n_annual, day_solar) * std::sin(3 * ome_d * hour_solar);
    return result;
}
} // namespace pressure_model
namespace stochastic {
// VAR(1) model for residuals: [P, Q, T]
// x_t = x_{t-1} * A + noise
// noise = u * L^T, where u ~ N(0, I)
static constexpr double A[9] = {
    1.1718524122501768, -0.16250827863926495, -0.042868449379527017,
    0.0007948252085985541, 0.98135054428547319, 0.059779552342641915,
    -0.00052961557671421088, -0.0012003433998291832, 1.0018994497855744,
    -0.19052620100771195, 0.21532152282034112, 0.07046400000903752,
    0.0018443882356720676, -0.027374054227459393, -0.024167072881457509,
    0.0026209284631848453, 0.018379461431658188, -0.088526034726346903,
    0.014980343993425446, -0.058886131275230569, -0.017227203990532798,
    -0.0014438204091440261, 0.0088280488884341091, -0.0056164491654602296,
    -0.0023656498637016719, 0.0024640335889273141, 0.042966609168269929
};
static constexpr double L[9] = {
    0.038787388451621996, 0, 0,
    -0.011011358561646678, 0.20522226763063708, 0,
    -0.042084626125372816, 0.026122804129809177, 0.19287369451177205
};

struct State {
    double p_res = 0.0;
    double q_res = 0.0;
    double t_res = 0.0;
};

// Updates state using VAR(1) process.
// u_p, u_q, u_t: Independent standard normal random numbers provided by the caller.
inline void update_state(State& s, double u_p, double u_q, double u_t) {
    // Map state and input to arrays
    double x[3] = {s.p_res, s.q_res, s.t_res};
    double u[3] = {u_p, u_q, u_t};

    // Compute noise n = u * L^T
    // n[j] = sum_k u[k] * L[j][k] (since L^T at [k][j] is L[j][k])
    double n[3] = {0.0, 0.0, 0.0};
    for(int j=0; j<3; ++j) {
        for(int k=0; k<3; ++k) {
             n[j] += u[k] * L[j*3 + k];
        }
    }

    // Compute transition v = x * A
    // v[j] = sum_k x[k] * A[k][j] (assuming row vector x)
    double v[3] = {0.0, 0.0, 0.0};
    for(int j=0; j<3; ++j) {
        for(int k=0; k<3; ++k) {
             v[j] += x[k] * A[k*3 + j];
        }
    }

    s.p_res = v[0] + n[0];
    s.q_res = v[1] + n[1];
    s.t_res = v[2] + n[2];
}

} // namespace stochastic
inline double predict_temperature(double day_utc, double hour_utc){
    double hour_solar = detail::wrap_hour(hour_utc + delta_utc_solar_h);
    double day_solar  = detail::wrap_day(day_utc + (delta_utc_solar_h / 24.0));
    return temperature_model::evaluate(day_solar, hour_solar);
}
inline double predict_specific_humidity(double day_utc, double hour_utc){
    double hour_solar = detail::wrap_hour(hour_utc + delta_utc_solar_h);
    double day_solar  = detail::wrap_day(day_utc + (delta_utc_solar_h / 24.0));
    return specific_humidity_model::evaluate(day_solar, hour_solar);
}
inline double predict_pressure(double day_utc, double hour_utc){
    double hour_solar = detail::wrap_hour(hour_utc + delta_utc_solar_h);
    double day_solar  = detail::wrap_day(day_utc + (delta_utc_solar_h / 24.0));
    return pressure_model::evaluate(day_solar, hour_solar);
}
inline void predict(double day_utc, double hour_utc, double& temperature_c, double& specific_humidity_kg_kg, double& pressure_hpa){
    double hour_solar = detail::wrap_hour(hour_utc + delta_utc_solar_h);
    double day_solar  = detail::wrap_day(day_utc + (delta_utc_solar_h / 24.0));
    temperature_c = temperature_model::evaluate(day_solar, hour_solar);
    specific_humidity_kg_kg = specific_humidity_model::evaluate(day_solar, hour_solar);
    pressure_hpa = pressure_model::evaluate(day_solar, hour_solar);
}
inline void predict_stochastic(double day_utc, double hour_utc, const stochastic::State& state, double& temperature_c, double& specific_humidity_kg_kg, double& pressure_hpa){
    predict(day_utc, hour_utc, temperature_c, specific_humidity_kg_kg, pressure_hpa);
    temperature_c += state.t_res;
    specific_humidity_kg_kg += state.q_res;
    pressure_hpa += state.p_res;
}
} // namespace harmoclimat