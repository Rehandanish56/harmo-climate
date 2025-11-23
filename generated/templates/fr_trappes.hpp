// Auto-generated linear harmonic climate model
// Station name : TRAPPES
// Station code : 78621001
#pragma once
#include <cmath>
namespace harmoclimat {
static constexpr double longitude_deg = 2.0098330974578857;
static constexpr double latitude_deg = 48.774333953857422;
static constexpr double delta_utc_solar_h = 0.13398887316385899;
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
    11.638417272450793, -7.3782932062515254, -2.5086140763819555, -0.0060918821725246498, 0.54466731449633721, -0.060069061607810025, -0.11667871143973588
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -2.8390005693877773, 1.6745020055742577, 0.00066828715375652004, 0.45397783842263006, -0.1357719740962666, -0.11275214184262912, 0.033292077592963584
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -0.98053122654412772, 0.17446539919042628, -0.25302612817521358, 0.040413236946483265, -0.22408947127082038, 0.17784672917738092, 0.14541898630878289
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    0.42981410544625692, 0.2471156692680421, -0.10542613832806358, -0.21435220476302844, 0.057121418180870701, -0.038680270081056366, 0.012640041712347623
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.070512022985157949, 0.21690675712795168, 0.054709955022445973, 0.12588275748154798, 0.057211022722459974, -0.073383129619802645, -0.0048341734381682679
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.099982407267709156, -0.26239183641782621, 0.029431072075564793, -0.070334477120593167, 0.039258469923823605, 0.05014935267433273, -0.024290641936622712
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -0.082142485734420359, 0.12750809457361242, -0.0073638876060988121, -0.025338680058218186, 0.0055901951957137075, -0.060564093839071825, -0.018125414545795516
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
    0.0068275037066256941, -0.002186793246854688, -0.0013199803684118363, 0.00021253771449055843, 0.00020229747988541577, -7.4107498701314132e-05, 0.00012095640611953705
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -8.3973133691718127e-05, -7.6457083925289071e-05, 1.9703292485341472e-05, 4.9089678548356538e-05, 5.4024923022088134e-05, 2.9236534956823457e-05, -4.8191037077659593e-05
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -1.0114457263557731e-05, -8.4935609665241332e-05, -3.532993725130242e-05, -8.8213791358582123e-06, 5.5683038078363855e-05, 2.509832152159941e-05, -3.0776756813410895e-05
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -5.1249688005141719e-05, 0.00011429255467430262, 1.5999278832724847e-05, 1.0123730339673972e-05, -2.0715575190824556e-05, -2.7168566273618418e-05, 8.2005914957169169e-07
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -7.2368157692215875e-05, 5.9997697812010312e-05, 1.0390491223274431e-05, 2.5258663349094558e-05, -1.2313708083956927e-05, 1.1385699427268625e-06, 1.9190987825439142e-06
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    6.6618676898670565e-06, 1.4506848264921625e-06, -2.0907370591588932e-05, -3.2507040003657517e-05, 8.6484328127020218e-06, 5.8417675737644973e-07, 8.5627047571332946e-06
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -8.7853882709634574e-06, 4.2004860795453037e-05, -6.7915643775755307e-07, -3.9388938346363978e-06, -2.2838968603717063e-06, -1.4264461066672783e-05, 4.0246624979124246e-07
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
    996.34802359457456, 0.17924505988554773, -0.07560082035139172, 0.73528077788660839, 0.75737319972123407, 0.39281000192045673, 0.67880053978922394
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    0.081588345894084438, -0.042517210104816738, 0.02742601007386563, 0.023541935483631274, 0.018007378275815875, 0.0046122781846075463, -0.011130120433258436
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    0.1324183554847731, -0.16755232128234898, 0.0256399266640145, -0.0049023446722671208, 0.016157904373887033, -0.0062732242580598356, -0.022190307220323752
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -0.11592167918276605, 0.083859386591716406, 0.040688657327454136, 0.058717831830352514, 0.045551251010201535, -0.042242846895098614, -0.012093326475328004
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.38477787018960452, -0.0027458370714378206, -0.026303142568752783, 0.037582905737921904, -0.017942655253648711, 0.0038618371728954865, 0.019471787866456619
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.010142375884566738, 0.10238663684916191, -0.032024164343931641, -0.015970457538759662, -0.0055999623540091033, -0.002165742175722509, 0.0098446705566142354
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    0.043032667653678636, 0.098356153223307305, 0.0037940104188911258, 0.012512679407018069, 0.010861785477886169, -0.0029120250738200125, 0.007997061976681246
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
    1.1694122875776052, -0.14039967802379366, -0.03410032624287207,
    -0.00046590893976720948, 0.9913834796449198, 0.055008678636095217,
    -0.00059252208634516957, -0.005378592233782023, 1.0118830764093423,
    -0.18578203835912963, 0.18385528382843541, 0.057661281137225406,
    0.0028553951882849694, -0.032808443805464022, -0.020646816872410104,
    0.0022374962169729954, 0.024843385760810653, -0.091148920894285071,
    0.012595131567281717, -0.050287879389941059, -0.012605259102434645,
    -0.0012969621203889788, 0.0061769941087087832, -0.0092772017094874346,
    -0.0018958576897605359, 0.00016112765384336536, 0.040816218485560371
};
static constexpr double L[9] = {
    0.039471415250593077, 0, 0,
    -0.0051241339108887155, 0.19362996979130392, 0,
    -0.038551631795136994, 0.020508208197155816, 0.18346998018394098
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