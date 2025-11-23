// Auto-generated linear harmonic climate model
// Station name : PARIS-MONTSOURIS
// Station code : 75114001
#pragma once
#include <cmath>
namespace harmoclimat {
static constexpr double longitude_deg = 2.3378329277038574;
static constexpr double latitude_deg = 48.821666717529297;
static constexpr double delta_utc_solar_h = 0.15585552851359052;
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
    12.908384689303571, -7.5731027290108042, -2.3958958734299047, -0.0096667700235419235, 0.56032868190893581, -0.036610217876375417, -0.13926178656921026
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -2.5031898373722248, 1.6085994021415113, -0.028920315556659297, 0.34491577618657399, -0.055753098632676246, -0.11804478314643117, 0.013765714468710222
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -1.0345038545492955, 0.1687642055046277, -0.20176269312058456, 0.087407607482594241, -0.20837485735212652, 0.16361862711002392, 0.11452316753244655
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    0.44559100344269686, 0.12843576316994143, -0.06765203698791708, -0.23496541633669951, 0.067445703691818806, -0.021843764795717868, -0.0041879886801589571
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.11972298327603043, 0.26358697352212535, 0.04642600090097794, 0.11252800573268937, 0.054975104142701105, -0.097278799692118237, -0.016473901824592942
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.061617503472662441, -0.26369949519984287, 0.025774891904274653, -0.0075549270201693934, 0.01173831249181787, 0.052335489503180815, -0.023237661482805691
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -0.060047254558118948, 0.12884991691935266, -0.045449113088334825, -0.080049923539305248, 0.0022665804437719764, -0.034054361763979843, 0.0042756085294198866
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
    0.0067789508647979335, -0.0021807762653945546, -0.0013320012736019213, 0.00018310484074214974, 0.00021082092732907225, -7.0066467256123214e-05, 0.00012845764962743229
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -3.7586438735782814e-05, -2.62887376668821e-05, 7.6728087868284532e-06, 9.7374274368640949e-06, 4.1250510474508543e-05, 1.9083750836860351e-05, -3.1663005699014358e-05
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    1.1706165115947021e-05, -0.00010931528952790457, -3.0104227894760167e-05, 3.6240287256764505e-06, 6.5415183417196236e-05, 2.7364184124665256e-05, -3.7422470397072011e-05
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -1.314897843318658e-05, 4.3229933164991232e-05, 4.467894875572187e-06, 1.1351767949439079e-07, -6.2870997835582151e-06, -1.3575388912174891e-05, 4.1274587430846082e-06
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -5.8928077995885015e-05, 5.1822361715790932e-05, 1.1312793226090216e-05, 1.8784053779795973e-05, -1.3287144898019315e-05, -2.8797530148433881e-06, -3.8195039126522181e-06
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    2.8540853048179568e-06, -4.4874366457520791e-06, -1.2357985751851086e-05, -1.2296613174177007e-05, 2.2451374314676853e-06, 1.259415141901682e-06, 5.0422156051825986e-06
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -4.5286245929270255e-06, 3.2079629029040555e-05, -5.8011595170388899e-06, -1.2571278480289337e-05, -1.8673712287818164e-06, -9.256270454218634e-06, 5.9192371990240591e-06
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
    1007.5582885282835, 0.50817048394516251, -0.0065260736983468693, 0.71424900395392121, 0.70514645693153144, 0.36196788722196005, 0.70883633914133737
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    0.12772115775033965, -0.080410220684036401, 0.027628478838537918, 0.014296485301854088, 0.0080423869611655954, 0.0085955019376157162, -0.0060023432962823782
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    0.20351834146319386, -0.18752375294913493, 0.033174909025870045, -0.01414289559146005, 0.03172114376655704, -0.0096788490092857751, -0.027269392907007233
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -0.12360156425535095, 0.081725660513302448, 0.042526026857808652, 0.063691347865278719, 0.044213747683844375, -0.039697572177534203, -0.014895522454673306
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.38673371551115154, -0.0059169182626219599, -0.019507090695451576, 0.035817092767074718, -0.024289754269282877, 0.0053793419176254784, 0.02207430732378881
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.010039447250828375, 0.10673086411782905, -0.034419653483326089, -0.013613681163140091, -0.0075431976087895957, -0.0029430248602720301, 0.0074375458394995457
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    0.040396053267163064, 0.1003837944307587, -0.00082070305419456799, 0.016855212216105843, 0.013157120644934342, -0.00082720824674950187, 0.007819149974521343
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
    1.1673119154514502, -0.13798985705540298, -0.022758430870999487,
    -0.0011096004362053519, 0.9819528839175985, 0.047705205384920485,
    6.8584487173691151e-05, 0.0024031952344317001, 1.0106138996947456,
    -0.18270312069222433, 0.17989170059696991, 0.047529031989378648,
    0.0038937931025137748, -0.023810756677061062, -0.017026588530763602,
    0.0023980929846539392, 0.019512742294686337, -0.080509511533776931,
    0.011424223366979236, -0.050180610815153237, -0.013456553711863883,
    -0.0018088618396124613, 0.0048439716290819557, -0.0095213439692080869,
    -0.0026703616218363988, -0.00072076446008860484, 0.037244609635053001
};
static constexpr double L[9] = {
    0.04040022629708008, 0, 0,
    -0.0069089808399735537, 0.19906674984910003, 0,
    -0.040652761635923927, 0.0097030525944209779, 0.17788675455586228
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