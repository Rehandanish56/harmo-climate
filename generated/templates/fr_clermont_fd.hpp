// Auto-generated linear harmonic climate model
// Station name : CLERMONT-FD
// Station code : 63113001
#pragma once
#include <cmath>
namespace harmoclimat {
static constexpr double longitude_deg = 3.1493330001831055;
static constexpr double latitude_deg = 45.786834716796875;
static constexpr double delta_utc_solar_h = 0.20995553334554037;
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
    12.300994936143741, -7.9787482143095376, -2.7311254528785249, -0.02711119351093828, 0.58642870903363875, -0.25509974254635065, -0.10589948188994526
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -3.6650372948205754, 1.6918074739831701, 0.062516440752354796, 0.4184213042092062, -0.21459375247073748, -0.028524069098807999, 0.065291832058520294
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -1.0054419885899764, -0.1503681089741426, -0.25213088316883736, 0.11047560533606539, -0.3147953855613061, 0.20764506816141715, 0.14490181688924383
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    0.67615779077397142, 0.43806268019821282, -0.14128656659090547, -0.24292710724153491, 0.059685537910504717, -0.057294246180649773, 0.023918840451927957
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.11513113084501404, 0.34443337745265462, 0.095349466880132946, 0.15677283458658622, 0.10307974296532882, -0.087644479189481095, -0.021572911126609241
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.12651512972259332, -0.32543445091387219, 0.038859093819770089, -0.1001545684429786, 0.066176284546065764, 0.027087463661004712, -0.025407071339738772
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -0.10739034952388736, 0.16996983749732616, 0.0017817386370606308, -0.070770224475284471, 0.0091418494883382627, -0.069438398598772003, -0.017429460877056196
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
    0.0066910382775709279, -0.0026560901579774287, -0.0013121378990567285, 0.0002827851030976144, 0.00016749994071474443, -9.3797400599445377e-05, 0.00012063635442915753
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -0.00023579517628408912, 2.6274801514915866e-05, 1.4474916249335259e-05, 5.4187075147912132e-05, 6.2593378731698958e-05, -5.2222160521406125e-06, -3.4549915186705666e-05
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -9.8274721978643118e-05, -2.6872027124952088e-05, 3.8910547873293792e-06, 7.3578812688016753e-07, 2.1051175811957958e-05, 1.4878457782090662e-05, -2.2590178231792452e-05
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -5.1390909567137402e-05, 0.00015743635833123797, 9.6685619968475656e-06, -8.3516955395762958e-06, -1.8226881039126371e-05, -2.5462006316836192e-05, -9.3382889016175585e-06
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -7.5844256892276468e-05, 6.3602972605551129e-05, 3.9688350015677148e-06, 3.620793694258599e-05, -8.6316927610716499e-06, -1.9135889505125025e-06, 8.9478503078143269e-07
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    2.9646628934302435e-05, -1.036658263633777e-05, -1.8287456122021683e-05, -4.4557763000151844e-05, 6.7303934010877052e-06, 3.7536151735508592e-06, 9.7156347178866896e-06
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -9.5849288090997408e-06, 4.7437553462129785e-05, -3.0829333976505489e-06, -6.9309250609284852e-06, 5.9911808640928691e-06, -1.6676728546661402e-05, -2.3354903985320512e-07
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
    978.65786919620859, 0.69817527989537753, -0.51273485298785371, 1.0381092571241686, 0.76686629079802615, 0.31108527747360576, 0.66327722653092736
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    0.29983349395970688, -0.15583793438238189, 0.035610928130596159, -0.028190655175134791, 0.0079494649551140485, -0.0041162803688548277, 0.00073417214327844212
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    0.21807517170027296, -0.1615665889252467, 0.025541212273602929, -0.029944464457922725, 0.072399617191030399, -0.024850648268732369, -0.041552628562379866
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -0.19882485419701981, 0.10885735376375759, 0.047547207334768825, 0.070472529634608275, 0.041703853769156946, -0.044160269029837898, -0.023531163721803576
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.43339032184661846, -0.038518659496591431, -0.027063617399856486, 0.035111955796705785, -0.025252768969160802, 0.019124492443873122, 0.016163064499322258
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.017191580421563928, 0.12134754948363433, -0.038896432840417487, -0.013305403101415694, -0.012169329668268655, -0.0048877798418896267, 0.0067830834413509251
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    0.045759717933063271, 0.10313259803101606, 0.011414341475021792, 0.022378361020474937, 0.020622823795512238, 0.001554335878206315, 0.0068180343594548337
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
    1.1294492854178189, -0.058617075251460246, -0.030898457093026136,
    -0.0029869059814436466, 0.94677687889995699, 0.044130776914078823,
    -0.0024756659755641912, 0.021250412537648825, 1.0009472928817329,
    -0.12387659988485833, 0.054030832514101013, 0.022884843164417037,
    0.0041410506292484071, -0.0014795212547655534, -0.020706596898363648,
    0.0091979251240281411, 0.0023176671172943503, -0.085993634992074464,
    -0.011066962245391582, 0.0006039394683524657, 0.014498865730222819,
    -0.0018728502477758143, 0.012059664479746976, -0.0032149074408171687,
    -0.0047460428281673279, 0.0031469689552073438, 0.044488365588944874
};
static constexpr double L[9] = {
    0.052937032644652202, 0, 0,
    -0.0014918561261965605, 0.23474301890066959, 0,
    -0.049857410682750215, 0.02656148303764953, 0.21289689037052081
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