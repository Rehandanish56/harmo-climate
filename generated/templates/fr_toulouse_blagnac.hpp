// Auto-generated linear harmonic climate model
// Station name : TOULOUSE-BLAGNAC
// Station code : 31069001
#pragma once
#include <cmath>
namespace harmoclimat {
static constexpr double longitude_deg = 1.3788330554962158;
static constexpr double latitude_deg = 43.620998382568359;
static constexpr double delta_utc_solar_h = 0.091922203699747726;
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
    14.273282385572992, -8.022285117458571, -3.0484153744995028, -0.19553959821959854, 0.81792583054249623, -0.24978308542679609, -0.081234928876341581
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -3.4077878397263524, 1.4131711019302859, 0.20732053923787275, 0.33545841910355112, -0.2085622990237313, -0.04927361587172286, -0.032368163818283231
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -1.309003076152067, 0.047943241532432276, -0.10710427989114617, -0.095824567095055213, -0.36559014878625662, 0.17608959628728874, 0.12029875482223751
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    0.61807257384978254, 0.3108590111843354, -0.10228640952348921, -0.13143031151409845, 0.047611126722333079, -0.032785403524782292, -0.0038405964231750648
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.10187423826799735, 0.2907448850920395, 0.085799482923615722, 0.14406671958558037, 0.058031250317396081, -0.043004659519286746, -0.020436100980948758
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.069138955816887407, -0.25394083541539264, 0.011690616768016684, -0.090522498360907297, 0.044048763595765296, 0.018753711499789834, -0.013613872664433199
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -0.056227894596304832, 0.15248616033135079, -0.029277886928603173, -0.055744511174944729, 0.01368952048566696, -0.06185916515550223, -0.010014890331307173
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
    0.0075655207998298795, -0.0025350630165942387, -0.0014045771500531771, 0.00017681052766573688, 0.0001824116934380476, -0.00012594465968076887, 9.5307365452126989e-05
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -3.7147225712032773e-05, -0.0002093700824025859, -0.00010218912128056653, 6.1721182878553852e-06, 0.00012170647240886344, 1.7641978171082015e-05, -2.717227044609629e-05
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    5.93840272626606e-05, -0.00021872435925770993, -3.8578601695817434e-05, -5.8050521247540954e-05, 6.0721661571403159e-05, 4.2345340864334393e-05, -1.2127563281995592e-05
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -2.6361709318486475e-05, 0.00010699378229114072, 4.396602122773055e-05, 3.4769713896847181e-05, -1.3927979484497944e-05, -1.5504836700063521e-05, -1.9962935186096961e-05
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -9.7522486152387152e-05, 3.464822929328937e-05, 8.8564561439784421e-06, 5.8672049523849726e-05, 4.1759836448112254e-06, -6.7465833015355766e-07, -1.5638995803006916e-05
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    1.2164294820494788e-05, 2.881289862229336e-06, -9.9489376962722196e-06, -3.6142195396713892e-05, -3.233460475332886e-06, -5.1980381687368717e-06, 1.1194649135058731e-05
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -3.0979290164248423e-06, 5.6996056367915023e-05, -5.2480503609936028e-06, -2.1174371879067239e-06, -1.146867597029321e-06, -1.7259704572760959e-05, -1.4809160551240821e-06
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
    999.77430082324588, 1.3865976758112346, -0.21780913932325699, 1.3471993183087752, 0.79490901446281315, 0.33870965825124316, 0.67632081989985671
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    0.24600669337416067, -0.17950923922575379, -0.0062837035398475927, -0.0059035623220044606, 0.030113703548022043, 0.012237877782063235, -0.0098484462217933878
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    0.25741632109790041, -0.21052563877098604, 0.004905420412367789, -0.0012527466491801598, 0.093709032045727544, -0.014024968925405596, -0.028918814618839495
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -0.21382421937578899, 0.10779432978527498, 0.065603041381895114, 0.066989488775622136, 0.060931978185616249, -0.060238978981043891, -0.01866104534716129
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.5001909271578121, -0.039535681499775266, -0.028586115050679699, 0.029766268025567481, -0.020690388864226698, 0.012992424960215682, 0.012915626793775017
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.020391491554601442, 0.13529614448148325, -0.036427313494580622, -0.016853332518542184, -0.0088061556805097037, -0.0059174800570697293, 0.0058312771506858581
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    0.053867247960696592, 0.11591588737666447, 0.0091308496937692867, 0.021919596642768675, 0.01644141031618835, 0.0055830902160303227, 0.0054361645280903868
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
    1.1245245789305516, -0.03744678743235684, -0.028271505669571255,
    0.00068748660437571718, 0.97279740936746395, 0.038680221302962649,
    -0.0041955870300600558, -0.00047107427405150559, 1.0248977283953671,
    -0.11994667031839903, 0.0095254592856595099, 0.025126777150285086,
    5.9740173803482857e-05, -0.024962105319203537, -0.010575708514773923,
    0.0088890114138552717, 0.017126350734933515, -0.11232068039362586,
    -0.010467231837099307, 0.021011872160682343, 0.011532680773991326,
    -0.0005977082919511427, 0.012428291203965144, -0.0086103056487850556,
    -0.003869974651076376, 0.0093004683898576385, 0.049654482419693863
};
static constexpr double L[9] = {
    0.057539036691476766, 0, 0,
    0.0030170892572706407, 0.21324583362800878, 0,
    -0.0396701972190543, 0.025826742961999858, 0.19502753559770825
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